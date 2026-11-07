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

#include "DecHelper.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "BitsFile.h"
#include "Md5Wrapper.h"
#include "RawFile.h"

// Include internal libavc decoder headers here inside the .cpp file to prevent
// public header pollution
// clang-format off
extern "C" {
#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"
#include "ih264d.h"
}
// clang-format on

using ::testing::Eq;

// Standard aligned memory allocators for libavc callbacks
static void* alignedAlloc(void* pv_ctxt, WORD32 alignment, WORD32 size) {
  (void)pv_ctxt;
#if defined(_WIN32)
  return _aligned_malloc(size, alignment);
#else
  void* buf = nullptr;
  if (0 != posix_memalign(&buf, alignment, size)) {
    return nullptr;
  }
  return buf;
#endif
}

static void alignedFree(void* pv_ctxt, void* pv_buf) {
  (void)pv_ctxt;
#if defined(_WIN32)
  _aligned_free(pv_buf);
#else
  free(pv_buf);
#endif
}

static IVD_ARCH_T getDefaultArch() {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || \
    defined(_M_IX86)
  return ARCH_X86_SSE42;
#elif defined(__arm__) || defined(_M_ARM)
  return ARCH_ARM_A9Q;
#elif defined(__aarch64__) || defined(_M_ARM64)
  return ARCH_ARMV8_GENERIC;
#else
  return ARCH_ARM_NONEON;
#endif
}

// Maps high-level Arch enum to stable native Ittiam API architecture values
static IVD_ARCH_T getNativeArch(libavc::test::Arch arch) {
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
      return ARCH_ARMV8_GENERIC;
    case libavc::test::Arch::x86_sse42:
      return ARCH_X86_SSE42;
    case libavc::test::Arch::x86_avx2:
      return ARCH_X86_AVX2;
  }
  return ARCH_ARM_NONEON;
}

// Helper to set decoder configuration parameters and decode mode
static bool setDecoderConfig(void* codec, IVD_VIDEO_DECODE_MODE_T decMode,
                             size_t dispWd) {
  ih264d_ctl_set_config_ip_t setConfigIp = {};
  ih264d_ctl_set_config_op_t setConfigOp = {};
  setConfigIp.s_ivd_ctl_set_config_ip_t.e_cmd = IVD_CMD_VIDEO_CTL;
  setConfigIp.s_ivd_ctl_set_config_ip_t.e_sub_cmd = IVD_CMD_CTL_SETPARAMS;
  setConfigIp.s_ivd_ctl_set_config_ip_t.u4_disp_wd = dispWd;
  setConfigIp.s_ivd_ctl_set_config_ip_t.e_frm_skip_mode = IVD_SKIP_NONE;
  setConfigIp.s_ivd_ctl_set_config_ip_t.e_frm_out_mode = IVD_DISPLAY_FRAME_OUT;
  setConfigIp.s_ivd_ctl_set_config_ip_t.e_vid_dec_mode = decMode;
  setConfigIp.s_ivd_ctl_set_config_ip_t.u4_size =
      sizeof(ih264d_ctl_set_config_ip_t);
  setConfigOp.s_ivd_ctl_set_config_op_t.u4_size =
      sizeof(ih264d_ctl_set_config_op_t);

  IV_API_CALL_STATUS_T ret = ih264d_api_function(
      static_cast<iv_obj_t*>(codec), &setConfigIp, &setConfigOp);
  return (ret == IV_SUCCESS);
}

std::unique_ptr<DecHelper> DecHelper::Builder::build() {
  if (mNumLoops == 0) {
    return nullptr;
  }

  std::unique_ptr<DecHelper> helper(new DecHelper());
  helper->mFormat = mFormat;
  helper->mCores = mCores;
  helper->mArch = mArch;
  helper->mInputFilePath = mInputFilePath;
  helper->mOutputFilePath = mOutputFilePath;
  helper->mRefMd5Path = mRefMd5Path;
  helper->mDisableRefRetention = mDisableRefRetention;
  helper->mNumLoops = mNumLoops;

  // Load and parse expected reference MD5s if the path is provided
  if (mRefMd5Path.has_value()) {
    std::ifstream md5File(*mRefMd5Path);
    if (!md5File.is_open()) {
      return nullptr;
    }
    std::string line;
    while (std::getline(md5File, line)) {
      if (!line.empty()) {
        std::stringstream ss(line);
        std::string md5;
        ss >> md5;
        if (!md5.empty()) {
          helper->mRefMd5s.push_back(md5);
        }
      }
    }
  }

  // Open input bitstream file
  if (!mInputFilePath.has_value() ||
      !helper->mBitsFile.open(*mInputFilePath, true)) {
    return nullptr;
  }

  // Open output YUV file if path is provided
  if (mOutputFilePath.has_value()) {
    if (!helper->mOutFile.open(*mOutputFilePath, false)) {
      return nullptr;
    }
  }

  // Open output MD5 file if path is provided
  if (mOutputMd5Path.has_value()) {
    helper->mOutputMd5File.open(*mOutputMd5Path);
    if (!helper->mOutputMd5File.is_open()) {
      return nullptr;
    }
  }

  // Open output sync-frame bitstream file if path is provided
  if (mOutputSyncFramesPath.has_value()) {
    helper->mOutputSyncFramesFile.open(*mOutputSyncFramesPath,
                                       std::ios::binary);
    if (!helper->mOutputSyncFramesFile.is_open()) {
      return nullptr;
    }
  }

  // Initialize native decoder (allocates mInputBuf internally)
  if (!helper->initDecoder()) {
    return nullptr;
  }
  return helper;
}

bool DecHelper::initDecoder() {
  ih264d_create_ip_t createIp = {};
  ih264d_create_op_t createOp = {};

  createIp.s_ivd_create_ip_t.e_cmd = IVD_CMD_CREATE;
  createIp.s_ivd_create_ip_t.u4_share_disp_buf =
      0;  // Non-sharing mode is safer for tests

  IV_COLOR_FORMAT_T colorFmt = IV_YUV_420P;
  switch (mFormat) {
    case libavc::test::Format::yuv400p:
      colorFmt = IV_GRAY;
      break;
    case libavc::test::Format::yuv420p:
      colorFmt = IV_YUV_420P;
      break;
    case libavc::test::Format::yuv420sp_uv:
      colorFmt = IV_YUV_420SP_UV;
      break;
    case libavc::test::Format::yuv420sp_vu:
      colorFmt = IV_YUV_420SP_VU;
      break;
    case libavc::test::Format::yuv422p:
      colorFmt = IV_YUV_422P;
      break;
    case libavc::test::Format::yuv444p:
      colorFmt = IV_YUV_444P;
      break;
  }
  createIp.s_ivd_create_ip_t.e_output_format = colorFmt;
  createIp.s_ivd_create_ip_t.pf_aligned_alloc = alignedAlloc;
  createIp.s_ivd_create_ip_t.pf_aligned_free = alignedFree;
  createIp.s_ivd_create_ip_t.pv_mem_ctxt = nullptr;
  createIp.s_ivd_create_ip_t.u4_size = sizeof(ih264d_create_ip_t);
  createOp.s_ivd_create_op_t.u4_size = sizeof(ih264d_create_op_t);

  createIp.u4_enable_frame_info = 0;
  createIp.u4_disable_ref_retention = mDisableRefRetention;
  createIp.u4_keep_threads_active = 1;

  IV_API_CALL_STATUS_T ret =
      ih264d_api_function(nullptr, &createIp, &createOp);
  if (ret != IV_SUCCESS) {
    return false;
  }

  mCodec = createOp.s_ivd_create_op_t.pv_handle;
  if (mCodec) {
    iv_obj_t* ps_handle = static_cast<iv_obj_t*>(mCodec);
    ps_handle->pv_fxns = (void*)&ih264d_api_function;
    ps_handle->u4_size = sizeof(iv_obj_t);
  }

  // Set number of cores
  ih264d_ctl_set_num_cores_ip_t setCoresIp = {};
  ih264d_ctl_set_num_cores_op_t setCoresOp = {};
  setCoresIp.e_cmd = IVD_CMD_VIDEO_CTL;
  setCoresIp.e_sub_cmd = static_cast<IVD_CONTROL_API_COMMAND_TYPE_T>(
      IH264D_CMD_CTL_SET_NUM_CORES);
  setCoresIp.u4_num_cores = mCores;
  setCoresIp.u4_size = sizeof(ih264d_ctl_set_num_cores_ip_t);
  setCoresOp.u4_size = sizeof(ih264d_ctl_set_num_cores_op_t);

  ret = ih264d_api_function(static_cast<iv_obj_t*>(mCodec), &setCoresIp,
                            &setCoresOp);
  if (ret != IV_SUCCESS) {
    deinitDecoder();
    return false;
  }

  // Set processor architecture using local mapper helper
  ih264d_ctl_set_processor_ip_t setProcIp = {};
  ih264d_ctl_set_processor_op_t setProcOp = {};
  setProcIp.e_cmd = IVD_CMD_VIDEO_CTL;
  setProcIp.e_sub_cmd = static_cast<IVD_CONTROL_API_COMMAND_TYPE_T>(
      IH264D_CMD_CTL_SET_PROCESSOR);
  setProcIp.u4_arch = getNativeArch(mArch);
  setProcIp.u4_soc = SOC_GENERIC;
  setProcIp.u4_size = sizeof(ih264d_ctl_set_processor_ip_t);
  setProcOp.u4_size = sizeof(ih264d_ctl_set_processor_op_t);

  ret = ih264d_api_function(static_cast<iv_obj_t*>(mCodec), &setProcIp,
                            &setProcOp);
  if (ret != IV_SUCCESS) {
    deinitDecoder();
    return false;
  }

  // Initialize to header decode mode using config helper
  if (!setDecoderConfig(mCodec, IVD_DECODE_HEADER, 0)) {
    deinitDecoder();
    return false;
  }

  // Allocate chunk-sized input bitstream buffer
  if (!mInputBuf.allocBuffer(1024 * 1024)) {
    deinitDecoder();
    return false;
  }

  mHeaderDecoded = false;
  mInFlushMode = false;
  return true;
}

void DecHelper::deinitDecoder() {
  if (mCodec) {
    ih264d_delete_ip_t deleteIp = {};
    ih264d_delete_op_t deleteOp = {};
    deleteIp.s_ivd_delete_ip_t.e_cmd = IVD_CMD_DELETE;
    deleteIp.s_ivd_delete_ip_t.u4_size = sizeof(ih264d_delete_ip_t);
    deleteOp.s_ivd_delete_op_t.u4_size = sizeof(ih264d_delete_op_t);

    ih264d_api_function(static_cast<iv_obj_t*>(mCodec), &deleteIp, &deleteOp);
    mCodec = nullptr;
  }
}

bool DecHelper::decodeHeader(size_t& bytesConsumed) {
  bytesConsumed = 0;
  if (!mCodec) {
    return false;
  }

  ih264d_video_decode_ip_t decodeIp = {};
  ih264d_video_decode_op_t decodeOp = {};

  decodeIp.s_ivd_video_decode_ip_t.e_cmd = IVD_CMD_VIDEO_DECODE;
  decodeIp.s_ivd_video_decode_ip_t.u4_ts = 0;
  decodeIp.s_ivd_video_decode_ip_t.u4_num_Bytes = mInputBuf.size();
  decodeIp.s_ivd_video_decode_ip_t.pv_stream_buffer = mInputBuf.data();
  decodeIp.s_ivd_video_decode_ip_t.u4_size =
      sizeof(ih264d_video_decode_ip_t);
  decodeOp.s_ivd_video_decode_op_t.u4_size =
      sizeof(ih264d_video_decode_op_t);

  IV_API_CALL_STATUS_T ret = ih264d_api_function(
      static_cast<iv_obj_t*>(mCodec), &decodeIp, &decodeOp);

  bytesConsumed = decodeOp.s_ivd_video_decode_op_t.u4_num_bytes_consumed;

  if (ret != IV_SUCCESS ||
      decodeOp.s_ivd_video_decode_op_t.u4_error_code != 0) {
    return false;
  }

  if (mOutputSyncFramesFile.is_open() && bytesConsumed > 0) {
    if (bytesConsumed > mInputBuf.size()) {
      bytesConsumed = 0;
      return false;
    }
    mOutputSyncFramesFile.write(
        reinterpret_cast<const char*>(mInputBuf.data()),
        static_cast<std::streamsize>(bytesConsumed));
    if (!mOutputSyncFramesFile.good()) {
      bytesConsumed = 0;
      return false;
    }
  }

  if (decodeOp.s_ivd_video_decode_op_t.u4_pic_wd == 0 ||
      decodeOp.s_ivd_video_decode_op_t.u4_pic_ht == 0) {
    return false;
  }

  mWidth = decodeOp.s_ivd_video_decode_op_t.u4_pic_wd;
  mHeight = decodeOp.s_ivd_video_decode_op_t.u4_pic_ht;
  mBitDepth = 8;

  // Allocate the output reconstructed YUV frame buffer
  if (!mOutputBuf.allocBuffer(mWidth, mHeight, mBitDepth, mFormat)) {
    return false;
  }

  // Ensure input buffer is large enough for a frame
  size_t requiredInputSize = std::max<size_t>(
      mWidth * mHeight * 3 * ((mBitDepth + 7) / 8), 1024 * 1024);
  if (mInputBuf.capacity() < requiredInputSize) {
    if (!mInputBuf.allocBuffer(requiredInputSize)) {
      return false;
    }
  }

  // Transition decoder to frame decode mode using config helper
  if (!setDecoderConfig(mCodec, IVD_DECODE_FRAME, 0)) {
    return false;
  }

  mHeaderDecoded = true;
  return true;
}

bool DecHelper::decodeFrame(size_t& bytesConsumed, bool& frameReady,
                            size_t& frameIndex, bool noInput) {
  bytesConsumed = 0;
  frameReady = false;
  if (!mCodec) {
    return false;
  }

  ih264d_video_decode_ip_t decodeIp = {};
  ih264d_video_decode_op_t decodeOp = {};

  decodeIp.s_ivd_video_decode_ip_t.e_cmd = IVD_CMD_VIDEO_DECODE;
  decodeIp.s_ivd_video_decode_ip_t.u4_ts = 0;
  if (noInput) {
    decodeIp.s_ivd_video_decode_ip_t.pv_stream_buffer = nullptr;
    decodeIp.s_ivd_video_decode_ip_t.u4_num_Bytes = 0;
  } else {
    decodeIp.s_ivd_video_decode_ip_t.pv_stream_buffer =
        static_cast<uint8_t*>(mInputBuf.data()) + mInputBuf.offset();
    decodeIp.s_ivd_video_decode_ip_t.u4_num_Bytes =
        mInputBuf.size() - mInputBuf.offset();
  }
  decodeIp.s_ivd_video_decode_ip_t.u4_size =
      sizeof(ih264d_video_decode_ip_t);
  decodeOp.s_ivd_video_decode_op_t.u4_size =
      sizeof(ih264d_video_decode_op_t);

  // Map destination output buffers for reconstruction
  size_t pixelSize = (mOutputBuf.bitDepth() + 7) / 8;
  decodeIp.s_ivd_video_decode_ip_t.s_out_buffer.u4_num_bufs =
      mOutputBuf.numPlanes();
  for (size_t i = 0; i < mOutputBuf.numPlanes(); ++i) {
    decodeIp.s_ivd_video_decode_ip_t.s_out_buffer.pu1_bufs[i] =
        mOutputBuf.planeData(i);
    decodeIp.s_ivd_video_decode_ip_t.s_out_buffer.u4_min_out_buf_size[i] =
        mOutputBuf.stride(i) * mOutputBuf.planeHeight(i) * pixelSize;
  }

  IV_API_CALL_STATUS_T ret = ih264d_api_function(
      static_cast<iv_obj_t*>(mCodec), &decodeIp, &decodeOp);
  bytesConsumed = decodeOp.s_ivd_video_decode_op_t.u4_num_bytes_consumed;
  if (ret != IV_SUCCESS) {
    if (mDisableRefRetention && !noInput &&
        decodeOp.s_ivd_video_decode_op_t.u4_output_present == 0 &&
        decodeOp.s_ivd_video_decode_op_t.u4_error_code == IVD_DEC_FRM_SKIPPED) {
      return true;
    }
    return false;
  }

  bool isSyncFrame =
      !noInput &&
      (decodeOp.s_ivd_video_decode_op_t.e_pic_type == IV_I_FRAME ||
       decodeOp.s_ivd_video_decode_op_t.e_pic_type == IV_IDR_FRAME);

  if (mOutputSyncFramesFile.is_open() && isSyncFrame && bytesConsumed > 0) {
    if (bytesConsumed > mInputBuf.size() - mInputBuf.offset()) {
      return false;
    }
    mOutputSyncFramesFile.write(
        reinterpret_cast<const char*>(mInputBuf.data() + mInputBuf.offset()),
        static_cast<std::streamsize>(bytesConsumed));
    if (!mOutputSyncFramesFile.good()) {
      return false;
    }
  }

  // Check if a frame has been decoded and is ready for display
  if (decodeOp.s_ivd_video_decode_op_t.u4_output_present == 1) {
    frameReady = true;

    if (mDisableRefRetention && !isSyncFrame) {
      return false;
    }

    // Process the decoded frame directly inside decodeFrame
    processDecodedFrame(mOutputBuf, frameIndex);
    mFramesDecoded = frameIndex;
  }

  return true;
}

bool DecHelper::flushDecoder(size_t& frameIndex) {
  if (!mCodec || !mHeaderDecoded) {
    return false;
  }

  // Put decoder in flush mode if not already done
  if (!mInFlushMode) {
    ivd_ctl_flush_ip_t flushIp = {};
    ivd_ctl_flush_op_t flushOp = {};
    flushIp.e_cmd = IVD_CMD_VIDEO_CTL;
    flushIp.e_sub_cmd = IVD_CMD_CTL_FLUSH;
    flushIp.u4_size = sizeof(ivd_ctl_flush_ip_t);
    flushOp.u4_size = sizeof(ivd_ctl_flush_op_t);

    IV_API_CALL_STATUS_T ret = ih264d_api_function(
        static_cast<iv_obj_t*>(mCodec), &flushIp, &flushOp);
    if (ret != IV_SUCCESS || flushOp.u4_error_code != 0) {
      return false;
    }
    mInFlushMode = true;
  }

  // Loop internally and delegate to decodeFrame with noInput = true
  while (true) {
    size_t consumed = 0;
    bool frameReady = false;
    if (!decodeFrame(consumed, frameReady, frameIndex, true)) {
      if (frameReady) {
        return false;
      }
      // Decoder returns an error once all the frames are flushed
      // which is expected. So return true here.
      return true;
    }
    if (!frameReady) {
      break;
    }
  }

  return true;
}

std::string DecHelper::computeFrameMd5(const RawBuf& buf) {
  Md5Wrapper md5;
  md5.init();

  size_t planesCount = buf.numPlanes();
  for (size_t p = 0; p < planesCount; ++p) {
    const uint8_t* planePtr = buf.planeData(p);
    size_t pixelSize = (buf.bitDepth() + 7) / 8;
    size_t planeW = buf.planeWidth(p);
    size_t planeH = buf.planeHeight(p);
    size_t planeStride = buf.stride(p);

    for (size_t r = 0; r < planeH; ++r) {
      md5.update(planePtr + (r * planeStride * pixelSize), planeW * pixelSize);
    }
  }

  return md5.finalize();
}

void DecHelper::processDecodedFrame(RawBuf& buf, size_t& frameIndex) {
  if (mOutputFilePath.has_value()) {
    mOutFile.write(buf);
  }

  bool needRefMd5 = frameIndex < mRefMd5s.size();
  bool needOutMd5 = mOutputMd5File.is_open();
  if (needRefMd5 || needOutMd5) {
    std::string computedMd5 = computeFrameMd5(buf);
    if (needRefMd5) {
      EXPECT_THAT(computedMd5, Eq(mRefMd5s[frameIndex]))
          << "MD5 mismatch at decoded frame " << frameIndex << " of file "
          << *mInputFilePath;
    }
    if (needOutMd5) {
      mOutputMd5File << computedMd5 << " " << (frameIndex + 1) << "\n";
    }
  }
  frameIndex++;
}

bool DecHelper::decodeFile() {
  size_t frameIndex = 0;
  mFramesDecoded = 0;

  size_t inputFrameSize = mInputBuf.capacity();
  for (size_t loop = 0; loop < mNumLoops; ++loop) {
    size_t inputFileOffset = 0;
    mBitsFile.seek(0);
    while (!mBitsFile.isEof()) {
      size_t bytesRead = mBitsFile.read(mInputBuf, inputFrameSize);
      if (bytesRead == 0) break;

      size_t consumed = 0;
      if (!mHeaderDecoded) {
        if (!decodeHeader(consumed)) {
          if (consumed == 0) {
            return false;
          }
        }
        inputFrameSize = mInputBuf.capacity();
      } else {
        bool frameReady = false;

        if (!decodeFrame(consumed, frameReady, frameIndex)) {
          return false;
        }
        if (consumed == 0 && !frameReady) {
          break;
        }
      }
      inputFileOffset += consumed;
      mBitsFile.seek(inputFileOffset);
    }

    if (!mHeaderDecoded) {
      return false;
    }
  }

  // Flush remaining frames (invokes processDecodedFrame internally for all
  // remaining frames)
  if (!flushDecoder(frameIndex)) {
    return false;
  }

  mOutFile.close();

  if (mOutputMd5File.is_open()) {
    mOutputMd5File.flush();
    if (!mOutputMd5File.good()) {
      return false;
    }
  }

  if (mOutputSyncFramesFile.is_open()) {
    mOutputSyncFramesFile.flush();
    if (!mOutputSyncFramesFile.good()) {
      return false;
    }
  }

  return true;
}
