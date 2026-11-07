/******************************************************************************
 *
 * Copyright (C) 2026 Ittiam Systems Pvt Ltd, Bangalore
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 ******************************************************************************/

#include "EncHelper.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "BitsFile.h"
#include "Md5Wrapper.h"
#include "RawFile.h"

// Include internal libavc encoder headers here inside the .cpp file to prevent
// public header pollution
// clang-format off
extern "C" {
#include "ih264_defs.h"
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264e.h"
}
// clang-format on

static void* alignedAlloc(size_t alignment, size_t size) {
#if defined(_WIN32)
  return _aligned_malloc(size, alignment);
#else
  void* buf = nullptr;
  size_t align = std::max<size_t>(alignment, sizeof(void*));
  if (0 != posix_memalign(&buf, align, size)) {
    return nullptr;
  }
  return buf;
#endif
}

static void alignedFree(void* pv_buf) {
#if defined(_WIN32)
  _aligned_free(pv_buf);
#else
  free(pv_buf);
#endif
}

static IV_ARCH_T getDefaultArch() { return ARCH_NA; }

// Map architecture to native enum using local mapper helper
static IV_ARCH_T getNativeArch(libavc::test::Arch arch) {
  switch (arch) {
    case libavc::test::Arch::defaultArch:
      return getDefaultArch();
    case libavc::test::Arch::generic:
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || \
    defined(_M_IX86)
      return ARCH_X86_GENERIC;
#else  // defined(__arm__) || defined(_M_ARM) || defined(__aarch64__) ||
       // defined(_M_ARM64)
      return ARCH_ARM_NONEON;
#endif
    case libavc::test::Arch::armv7:
      return ARCH_ARM_A9Q;
    case libavc::test::Arch::armv8:
      return ARCH_ARM_V8_NEON;
    case libavc::test::Arch::x86_sse42:
    case libavc::test::Arch::x86_avx2:
      return ARCH_X86_SSE42;
  }
  return ARCH_NA;
}

static IV_COLOR_FORMAT_T getNativeColorFormat(libavc::test::Format format) {
  switch (format) {
    case libavc::test::Format::yuv400p:
      return IV_GRAY;
    case libavc::test::Format::yuv420p:
      return IV_YUV_420P;
    case libavc::test::Format::yuv420sp_uv:
      return IV_YUV_420SP_UV;
    case libavc::test::Format::yuv420sp_vu:
      return IV_YUV_420SP_VU;
    case libavc::test::Format::yuv422p:
      return IV_YUV_422P;
    case libavc::test::Format::yuv444p:
      return IV_YUV_444P;
  }
  return IV_YUV_420P;
}

static bool setEncMode(void* codec, IVE_ENC_MODE_T encMode) {
  ive_ctl_set_enc_mode_ip_t encModeIp = {};
  ive_ctl_set_enc_mode_op_t encModeOp = {};
  encModeIp.u4_size = sizeof(ive_ctl_set_enc_mode_ip_t);
  encModeOp.u4_size = sizeof(ive_ctl_set_enc_mode_op_t);
  encModeIp.e_cmd = IVE_CMD_VIDEO_CTL;
  encModeIp.e_sub_cmd = IVE_CMD_CTL_SET_ENC_MODE;
  encModeIp.e_enc_mode = encMode;
  encModeIp.u4_timestamp_high = -1;
  encModeIp.u4_timestamp_low = -1;
  return ih264e_api_function(static_cast<iv_obj_t*>(codec), &encModeIp,
                             &encModeOp) == IV_SUCCESS;
}

std::unique_ptr<EncHelper> EncHelper::Builder::build() {
  std::unique_ptr<EncHelper> helper(new EncHelper());
  helper->mFormat = mFormat;
  helper->mCores = mCores;
  helper->mArch = mArch;
  helper->mWidth = mWidth;
  helper->mHeight = mHeight;
  helper->mBitrate = mBitrate;
  helper->mFrameRateNum = mFrameRateNum;
  helper->mFrameRateDenom = mFrameRateDenom;
  helper->mRcMode = mRcMode;
  helper->mFrameQp = mFrameQp;
  helper->mQualityPreset = mQualityPreset;
  helper->mIInterval = mIInterval;
  helper->mIdrInterval = mIdrInterval;
  helper->mBframes = mBframes;
  helper->mProfile = mProfile;
  helper->mLevel = mLevel;
  helper->mInputFilePath = mInputFilePath;
  helper->mOutFilePath = mOutFilePath;
  helper->mReconMd5FilePath = mReconMd5FilePath;

  // Open input raw YUV file
  if (!mInputFilePath.has_value() ||
      !helper->mRawFile.open(*mInputFilePath, true)) {
    return nullptr;
  }

  // Open output bitstream file
  if (!mOutFilePath.has_value() ||
      !helper->mBitsFile.open(*mOutFilePath, false)) {
    return nullptr;
  }

  // Initialize native encoder
  if (!helper->initEncoder()) {
    return nullptr;
  }
  return helper;
}

bool EncHelper::initEncoder() {
  // 1. Get number of memory records
  iv_num_mem_rec_ip_t numMemRecIp = {};
  iv_num_mem_rec_op_t numMemRecOp = {};
  numMemRecIp.u4_size = sizeof(iv_num_mem_rec_ip_t);
  numMemRecOp.u4_size = sizeof(iv_num_mem_rec_op_t);
  numMemRecIp.e_cmd = IV_CMD_GET_NUM_MEM_REC;

  if (ih264e_api_function(nullptr, &numMemRecIp, &numMemRecOp) != IV_SUCCESS) {
    return false;
  }

  mNumMemRecords = numMemRecOp.u4_num_mem_rec;
  iv_mem_rec_t* memRecords = static_cast<iv_mem_rec_t*>(
      malloc(mNumMemRecords * sizeof(iv_mem_rec_t)));
  if (!memRecords) {
    mNumMemRecords = 0;
    return false;
  }
  mMemRecords = memRecords;

  for (size_t i = 0; i < mNumMemRecords; ++i) {
    memRecords[i].u4_size = sizeof(iv_mem_rec_t);
    memRecords[i].pv_base = nullptr;
    memRecords[i].u4_mem_size = 0;
    memRecords[i].u4_mem_alignment = 0;
    memRecords[i].e_mem_type = IV_NA_MEM_TYPE;
  }

  // 2. Fill memory records
  IV_COLOR_FORMAT_T colorFmt = getNativeColorFormat(mFormat);
  iv_fill_mem_rec_ip_t fillMemRecIp = {};
  iv_fill_mem_rec_op_t fillMemRecOp = {};
  fillMemRecIp.u4_size = sizeof(iv_fill_mem_rec_ip_t);
  fillMemRecOp.u4_size = sizeof(iv_fill_mem_rec_op_t);
  fillMemRecIp.e_cmd = IV_CMD_FILL_NUM_MEM_REC;
  fillMemRecIp.ps_mem_rec = memRecords;
  fillMemRecIp.u4_num_mem_rec = mNumMemRecords;
  fillMemRecIp.u4_max_wd = mWidth;
  fillMemRecIp.u4_max_ht = mHeight;
  fillMemRecIp.u4_max_ref_cnt = 2;
  fillMemRecIp.u4_max_reorder_cnt = 0;
  fillMemRecIp.u4_max_level = mLevel;
  fillMemRecIp.e_color_format = IV_YUV_420P;
  fillMemRecIp.u4_max_srch_rng_x = 256;
  fillMemRecIp.u4_max_srch_rng_y = 256;
  fillMemRecIp.u4_keep_threads_active = 1;

  if (ih264e_api_function(nullptr, &fillMemRecIp, &fillMemRecOp) !=
      IV_SUCCESS) {
    deinitEncoder();
    return false;
  }

  // 3. Allocate memory for each memory record
  for (size_t i = 0; i < mNumMemRecords; ++i) {
    memRecords[i].pv_base =
        alignedAlloc(memRecords[i].u4_mem_alignment, memRecords[i].u4_mem_size);
    if (!memRecords[i].pv_base) {
      deinitEncoder();
      return false;
    }
  }

  // 4. Initialize encoder instance
  mCodec = memRecords[0].pv_base;
  iv_obj_t* codecObj = static_cast<iv_obj_t*>(mCodec);
  codecObj->u4_size = sizeof(iv_obj_t);
  codecObj->pv_fxns = reinterpret_cast<void*>(&ih264e_api_function);

  ive_init_ip_t initIp = {};
  ive_init_op_t initOp = {};
  initIp.u4_size = sizeof(ive_init_ip_t);
  initOp.u4_size = sizeof(ive_init_op_t);
  initIp.e_cmd = IV_CMD_INIT;
  initIp.u4_num_mem_rec = mNumMemRecords;
  initIp.ps_mem_rec = memRecords;
  initIp.u4_max_wd = mWidth;
  initIp.u4_max_ht = mHeight;
  initIp.u4_max_ref_cnt = 2;
  initIp.u4_max_reorder_cnt = 0;
  initIp.u4_max_level = mLevel;
  initIp.e_inp_color_fmt = colorFmt;
  initIp.u4_enable_recon = mReconMd5FilePath.has_value() ? 1 : 0;
  initIp.e_recon_color_fmt = IV_YUV_420P;
  initIp.e_rc_mode = static_cast<IVE_RC_MODE_T>(mRcMode);
  initIp.u4_max_framerate = 120000;
  initIp.u4_max_bitrate = 240000000;
  initIp.u4_num_bframes = mBframes;
  initIp.e_content_type = IV_PROGRESSIVE;
  initIp.u4_max_srch_rng_x = 256;
  initIp.u4_max_srch_rng_y = 256;
  initIp.e_slice_mode = IVE_SLICE_MODE_NONE;
  initIp.u4_slice_param = 256;
  initIp.e_arch = getNativeArch(mArch);
  initIp.e_soc = SOC_GENERIC;
  initIp.u4_keep_threads_active = 1;

  if (ih264e_api_function(codecObj, &initIp, &initOp) != IV_SUCCESS) {
    deinitEncoder();
    return false;
  }

  // 5. Set default parameters
  {
    ive_ctl_setdefault_ip_t defaultIp = {};
    ive_ctl_setdefault_op_t defaultOp = {};
    defaultIp.u4_size = sizeof(ive_ctl_setdefault_ip_t);
    defaultOp.u4_size = sizeof(ive_ctl_setdefault_op_t);
    defaultIp.e_cmd = IVE_CMD_VIDEO_CTL;
    defaultIp.e_sub_cmd = IVE_CMD_CTL_SETDEFAULT;
    defaultIp.u4_timestamp_high = -1;
    defaultIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &defaultIp, &defaultOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 6. Set number of cores
  {
    ive_ctl_set_num_cores_ip_t coresIp = {};
    ive_ctl_set_num_cores_op_t coresOp = {};
    coresIp.u4_size = sizeof(ive_ctl_set_num_cores_ip_t);
    coresOp.u4_size = sizeof(ive_ctl_set_num_cores_op_t);
    coresIp.e_cmd = IVE_CMD_VIDEO_CTL;
    coresIp.e_sub_cmd = IVE_CMD_CTL_SET_NUM_CORES;
    coresIp.u4_num_cores = mCores;
    coresIp.u4_timestamp_high = -1;
    coresIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &coresIp, &coresOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 7. Set dimensions
  {
    ive_ctl_set_dimensions_ip_t dimIp = {};
    ive_ctl_set_dimensions_op_t dimOp = {};
    dimIp.u4_size = sizeof(ive_ctl_set_dimensions_ip_t);
    dimOp.u4_size = sizeof(ive_ctl_set_dimensions_op_t);
    dimIp.e_cmd = IVE_CMD_VIDEO_CTL;
    dimIp.e_sub_cmd = IVE_CMD_CTL_SET_DIMENSIONS;
    dimIp.u4_wd = mWidth;
    dimIp.u4_ht = mHeight;
    dimIp.u4_timestamp_high = -1;
    dimIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &dimIp, &dimOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 8. Set frame rate
  {
    ive_ctl_set_frame_rate_ip_t fpsIp = {};
    ive_ctl_set_frame_rate_op_t fpsOp = {};
    fpsIp.u4_size = sizeof(ive_ctl_set_frame_rate_ip_t);
    fpsOp.u4_size = sizeof(ive_ctl_set_frame_rate_op_t);
    fpsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    fpsIp.e_sub_cmd = IVE_CMD_CTL_SET_FRAMERATE;
    uint32_t fps = static_cast<uint32_t>(
        mFrameRateNum / std::max<size_t>(1, mFrameRateDenom));
    fpsIp.u4_src_frame_rate = fps ? fps : 30;
    fpsIp.u4_tgt_frame_rate = fps ? fps : 30;
    fpsIp.u4_timestamp_high = -1;
    fpsIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &fpsIp, &fpsOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 9. Set IPE params
  {
    ive_ctl_set_ipe_params_ip_t ipeIp = {};
    ive_ctl_set_ipe_params_op_t ipeOp = {};
    ipeIp.u4_size = sizeof(ive_ctl_set_ipe_params_ip_t);
    ipeOp.u4_size = sizeof(ive_ctl_set_ipe_params_op_t);
    ipeIp.e_cmd = IVE_CMD_VIDEO_CTL;
    ipeIp.e_sub_cmd = IVE_CMD_CTL_SET_IPE_PARAMS;
    ipeIp.u4_enable_intra_4x4 = 1;
    ipeIp.u4_enc_speed_preset = static_cast<IVE_SPEED_CONFIG>(mQualityPreset);
    ipeIp.u4_constrained_intra_pred = 0;
    ipeIp.u4_timestamp_high = -1;
    ipeIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &ipeIp, &ipeOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 10. Set bitrate
  {
    ive_ctl_set_bitrate_ip_t brIp = {};
    ive_ctl_set_bitrate_op_t brOp = {};
    brIp.u4_size = sizeof(ive_ctl_set_bitrate_ip_t);
    brOp.u4_size = sizeof(ive_ctl_set_bitrate_op_t);
    brIp.e_cmd = IVE_CMD_VIDEO_CTL;
    brIp.e_sub_cmd = IVE_CMD_CTL_SET_BITRATE;
    brIp.u4_target_bitrate = mBitrate;
    brIp.u4_timestamp_high = -1;
    brIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &brIp, &brOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 11. Set QP
  {
    ive_ctl_set_qp_ip_t qpIp = {};
    ive_ctl_set_qp_op_t qpOp = {};
    qpIp.u4_size = sizeof(ive_ctl_set_qp_ip_t);
    qpOp.u4_size = sizeof(ive_ctl_set_qp_op_t);
    qpIp.e_cmd = IVE_CMD_VIDEO_CTL;
    qpIp.e_sub_cmd = IVE_CMD_CTL_SET_QP;
    qpIp.u4_i_qp = mFrameQp;
    qpIp.u4_i_qp_min = 4;
    qpIp.u4_i_qp_max = MAX_H264_QP;
    qpIp.u4_p_qp = mFrameQp;
    qpIp.u4_p_qp_min = 4;
    qpIp.u4_p_qp_max = MAX_H264_QP;
    qpIp.u4_b_qp = mFrameQp;
    qpIp.u4_b_qp_min = 4;
    qpIp.u4_b_qp_max = MAX_H264_QP;
    qpIp.u4_timestamp_high = -1;
    qpIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &qpIp, &qpOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 12. Set AIR params
  {
    ive_ctl_set_air_params_ip_t airIp = {};
    ive_ctl_set_air_params_op_t airOp = {};
    airIp.u4_size = sizeof(ive_ctl_set_air_params_ip_t);
    airOp.u4_size = sizeof(ive_ctl_set_air_params_op_t);
    airIp.e_cmd = IVE_CMD_VIDEO_CTL;
    airIp.e_sub_cmd = IVE_CMD_CTL_SET_AIR_PARAMS;
    airIp.e_air_mode = IVE_AIR_MODE_NONE;
    airIp.u4_air_refresh_period = 30;
    airIp.u4_timestamp_high = -1;
    airIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &airIp, &airOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 13. Set VBV params
  {
    ive_ctl_set_vbv_params_ip_t vbvIp = {};
    ive_ctl_set_vbv_params_op_t vbvOp = {};
    vbvIp.u4_size = sizeof(ive_ctl_set_vbv_params_ip_t);
    vbvOp.u4_size = sizeof(ive_ctl_set_vbv_params_op_t);
    vbvIp.e_cmd = IVE_CMD_VIDEO_CTL;
    vbvIp.e_sub_cmd = IVE_CMD_CTL_SET_VBV_PARAMS;
    vbvIp.u4_vbv_buf_size = 0;
    vbvIp.u4_vbv_buffer_delay = 1000;
    vbvIp.u4_timestamp_high = -1;
    vbvIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &vbvIp, &vbvOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 14. Set ME params
  {
    ive_ctl_set_me_params_ip_t meIp = {};
    ive_ctl_set_me_params_op_t meOp = {};
    meIp.u4_size = sizeof(ive_ctl_set_me_params_ip_t);
    meOp.u4_size = sizeof(ive_ctl_set_me_params_op_t);
    meIp.e_cmd = IVE_CMD_VIDEO_CTL;
    meIp.e_sub_cmd = IVE_CMD_CTL_SET_ME_PARAMS;
    meIp.u4_enable_fast_sad = 0;
    meIp.u4_enable_alt_ref = 0;
    meIp.u4_enable_hpel = 1;
    meIp.u4_enable_qpel = 1;
    meIp.u4_me_speed_preset = 100;
    meIp.u4_srch_rng_x = 16;
    meIp.u4_srch_rng_y = 16;
    meIp.u4_timestamp_high = -1;
    meIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &meIp, &meOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 15. Set GOP params
  {
    ive_ctl_set_gop_params_ip_t gopIp = {};
    ive_ctl_set_gop_params_op_t gopOp = {};
    gopIp.u4_size = sizeof(ive_ctl_set_gop_params_ip_t);
    gopOp.u4_size = sizeof(ive_ctl_set_gop_params_op_t);
    gopIp.e_cmd = IVE_CMD_VIDEO_CTL;
    gopIp.e_sub_cmd = IVE_CMD_CTL_SET_GOP_PARAMS;
    gopIp.u4_i_frm_interval = mIInterval;
    gopIp.u4_idr_frm_interval = mIdrInterval;
    gopIp.u4_timestamp_high = -1;
    gopIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &gopIp, &gopOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 16. Set deblock params
  {
    ive_ctl_set_deblock_params_ip_t deblkIp = {};
    ive_ctl_set_deblock_params_op_t deblkOp = {};
    deblkIp.u4_size = sizeof(ive_ctl_set_deblock_params_ip_t);
    deblkOp.u4_size = sizeof(ive_ctl_set_deblock_params_op_t);
    deblkIp.e_cmd = IVE_CMD_VIDEO_CTL;
    deblkIp.e_sub_cmd = IVE_CMD_CTL_SET_DEBLOCK_PARAMS;
    deblkIp.u4_disable_deblock_level = 0;
    deblkIp.u4_timestamp_high = -1;
    deblkIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &deblkIp, &deblkOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  // 17. Set profile params
  {
    ive_ctl_set_profile_params_ip_t profIp = {};
    ive_ctl_set_profile_params_op_t profOp = {};
    profIp.u4_size = sizeof(ive_ctl_set_profile_params_ip_t);
    profOp.u4_size = sizeof(ive_ctl_set_profile_params_op_t);
    profIp.e_cmd = IVE_CMD_VIDEO_CTL;
    profIp.e_sub_cmd = IVE_CMD_CTL_SET_PROFILE_PARAMS;
    profIp.e_profile = static_cast<IV_PROFILE_T>(mProfile);
    profIp.u4_entropy_coding_mode =
        (profIp.e_profile == IV_PROFILE_BASE) ? 0 : 1;
    profIp.u4_timestamp_high = -1;
    profIp.u4_timestamp_low = -1;
    if (ih264e_api_function(codecObj, &profIp, &profOp) != IV_SUCCESS) {
      deinitEncoder();
      return false;
    }
  }

  if (!mInputBuf.allocBuffer(mWidth, mHeight, /* bitDepth */ 8, mFormat)) {
    deinitEncoder();
    return false;
  }

  if (mReconMd5FilePath.has_value() &&
      !mReconBuf.allocBuffer(mWidth, mHeight, /* bitDepth */ 8,
                             libavc::test::Format::yuv420p)) {
    deinitEncoder();
    return false;
  }

  size_t outBufSize = std::max<size_t>((mWidth * mHeight * 3) >> 1, 2048);
  if (!mOutputBuf.allocBuffer(outBufSize)) {
    deinitEncoder();
    return false;
  }

  mIsLast = false;
  mTimestampLow = 0;
  mReconMd5s.clear();
  return true;
}

void EncHelper::deinitEncoder() {
  if (mMemRecords) {
    iv_mem_rec_t* memRecords = static_cast<iv_mem_rec_t*>(mMemRecords);
    for (size_t i = 0; i < mNumMemRecords; ++i) {
      if (memRecords[i].pv_base) {
        alignedFree(memRecords[i].pv_base);
        memRecords[i].pv_base = nullptr;
      }
    }
    free(memRecords);
    mMemRecords = nullptr;
  }
  mNumMemRecords = 0;
  mCodec = nullptr;
}

void EncHelper::addReconMd5(const std::string& md5) {
  std::lock_guard<std::mutex> lock(mReconMutex);
  mReconMd5s.push_back(md5);
}

std::string EncHelper::computeReconFrameMd5(void* pv_recon_buf) {
  const iv_raw_buf_t* reconBuf = static_cast<const iv_raw_buf_t*>(pv_recon_buf);
  Md5Wrapper md5;
  md5.init();

  size_t numComp = (reconBuf->e_color_fmt == IV_YUV_420P) ? 3 : 2;
  for (size_t c = 0; c < numComp; ++c) {
    const uint8_t* planePtr = static_cast<const uint8_t*>(reconBuf->apv_bufs[c]);
    if (!planePtr) {
      continue;
    }
    size_t wd = reconBuf->au4_wd[c];
    size_t ht = reconBuf->au4_ht[c];
    size_t strd = reconBuf->au4_strd[c];
    for (size_t r = 0; r < ht; ++r) {
      md5.update(planePtr + r * strd, wd);
    }
  }

  return md5.finalize();
}

bool EncHelper::encodeHeader() {
  if (!mCodec) {
    return false;
  }

  if (!setEncMode(mCodec, IVE_ENC_MODE_HEADER)) {
    return false;
  }

  ih264e_video_encode_ip_t encodeIp = {};
  ih264e_video_encode_op_t encodeOp = {};
  encodeIp.s_ive_ip.u4_size = sizeof(ih264e_video_encode_ip_t);
  encodeOp.s_ive_op.u4_size = sizeof(ih264e_video_encode_op_t);
  encodeIp.s_ive_ip.e_cmd = IVE_CMD_VIDEO_ENCODE;
  encodeIp.s_ive_ip.u4_is_last = 0;
  encodeIp.s_ive_ip.s_out_buf.u4_size = sizeof(iv_bits_buf_t);
  encodeIp.s_ive_ip.s_out_buf.pv_buf = mOutputBuf.data();
  encodeIp.s_ive_ip.s_out_buf.u4_bytes = 0;
  encodeIp.s_ive_ip.s_out_buf.u4_bufsize = mOutputBuf.capacity();

  if (ih264e_api_function(static_cast<iv_obj_t*>(mCodec), &encodeIp,
                          &encodeOp) != IV_SUCCESS) {
    return false;
  }

  if (encodeOp.s_ive_op.output_present == 1 &&
      encodeOp.s_ive_op.s_out_buf.u4_bytes > 0) {
    mOutputBuf.setSize(encodeOp.s_ive_op.s_out_buf.u4_bytes);
    mOutputBuf.setOffset(0);
    mBitsFile.write(mOutputBuf, encodeOp.s_ive_op.s_out_buf.u4_bytes);
  }

  return setEncMode(mCodec, IVE_ENC_MODE_PICTURE);
}

bool EncHelper::encodeFrame(bool& bytesProduced, bool noInput) {
  bytesProduced = false;
  if (!mCodec) {
    return false;
  }

  ih264e_video_encode_ip_t encodeIp = {};
  ih264e_video_encode_op_t encodeOp = {};
  encodeIp.s_ive_ip.u4_size = sizeof(ih264e_video_encode_ip_t);
  encodeOp.s_ive_op.u4_size = sizeof(ih264e_video_encode_op_t);
  encodeIp.s_ive_ip.e_cmd = IVE_CMD_VIDEO_ENCODE;

  iv_raw_buf_t& inpRawBuf = encodeIp.s_ive_ip.s_inp_buf;
  inpRawBuf.u4_size = sizeof(iv_raw_buf_t);
  inpRawBuf.e_color_fmt = getNativeColorFormat(mFormat);

  if (!noInput) {
    size_t planesCount = mInputBuf.numPlanes();
    for (size_t i = 0; i < planesCount; ++i) {
      inpRawBuf.apv_bufs[i] = mInputBuf.planeData(i);
      inpRawBuf.au4_wd[i] = mInputBuf.planeWidth(i);
      inpRawBuf.au4_ht[i] = mInputBuf.planeHeight(i);
      inpRawBuf.au4_strd[i] = mInputBuf.stride(i);
    }
    encodeIp.s_ive_ip.u4_is_last = 0;
  } else {
    for (size_t i = 0; i < IV_MAX_RAW_COMPONENTS; ++i) {
      inpRawBuf.apv_bufs[i] = nullptr;
    }
    encodeIp.s_ive_ip.u4_is_last = 1;
  }

  encodeIp.s_ive_ip.u4_timestamp_high = 0;
  encodeIp.s_ive_ip.u4_timestamp_low = mTimestampLow++;

  encodeIp.s_ive_ip.s_out_buf.u4_size = sizeof(iv_bits_buf_t);
  encodeIp.s_ive_ip.s_out_buf.pv_buf = mOutputBuf.data();
  encodeIp.s_ive_ip.s_out_buf.u4_bytes = 0;
  encodeIp.s_ive_ip.s_out_buf.u4_bufsize = mOutputBuf.capacity();

  if (mReconMd5FilePath.has_value()) {
    iv_raw_buf_t& reconRawBuf = encodeIp.s_ive_ip.s_recon_buf;
    reconRawBuf.u4_size = sizeof(iv_raw_buf_t);
    reconRawBuf.e_color_fmt = IV_YUV_420P;
    for (size_t i = 0; i < mReconBuf.numPlanes(); ++i) {
      reconRawBuf.apv_bufs[i] = mReconBuf.planeData(i);
      reconRawBuf.au4_wd[i] = mReconBuf.planeWidth(i);
      reconRawBuf.au4_ht[i] = mReconBuf.planeHeight(i);
      reconRawBuf.au4_strd[i] = mReconBuf.stride(i);
    }
  }

  if (ih264e_api_function(static_cast<iv_obj_t*>(mCodec), &encodeIp,
                          &encodeOp) != IV_SUCCESS) {
    return false;
  }

  if (encodeOp.s_ive_op.output_present == 1 &&
      encodeOp.s_ive_op.s_out_buf.u4_bytes > 0) {
    mOutputBuf.setSize(encodeOp.s_ive_op.s_out_buf.u4_bytes);
    mOutputBuf.setOffset(0);
    bytesProduced = true;
    mBitsFile.write(mOutputBuf, encodeOp.s_ive_op.s_out_buf.u4_bytes);
  }

  if (mReconMd5FilePath.has_value() && encodeOp.s_ive_op.dump_recon == 1 &&
      encodeOp.s_ive_op.s_recon_buf.au4_wd[0] > 0) {
    std::string md5 = computeReconFrameMd5(&encodeOp.s_ive_op.s_recon_buf);
    addReconMd5(md5);
  }

  mIsLast = (encodeOp.s_ive_op.u4_is_last == 1);
  return true;
}

bool EncHelper::flushEncoder() {
  if (!mCodec) {
    return false;
  }

  while (!mIsLast) {
    bool bytesProduced = false;
    if (!encodeFrame(bytesProduced, true)) {
      return false;
    }
    if (!bytesProduced && mIsLast) {
      break;
    }
  }

  return true;
}

bool EncHelper::encodeFile() {
  mReconMd5s.clear();
  mIsLast = false;
  mTimestampLow = 0;

  if (!encodeHeader()) {
    return false;
  }

  while (mRawFile.read(mInputBuf)) {
    bool bytesProduced = false;
    if (!encodeFrame(bytesProduced)) {
      return false;
    }
  }

  if (!flushEncoder()) {
    return false;
  }

  mBitsFile.close();

  if (mReconMd5FilePath.has_value()) {
    std::ofstream md5File(*mReconMd5FilePath);
    if (!md5File.is_open()) {
      return false;
    }
    std::lock_guard<std::mutex> lock(mReconMutex);
    for (const auto& md5 : mReconMd5s) {
      md5File << md5 << " 0\n";
    }
    mReconMd5s.clear();
  }

  return true;
}
