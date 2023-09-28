/*
 * Copyright (C) 2021 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "ih264_defs.h"
#include "ih264_typedefs.h"
#include "ih264e.h"
#include "ih264e_error.h"

#include "TestArgs.h"

#define MAX_FRAME_HEIGHT 1080
#define MAX_FRAME_WIDTH 1920
#define MAX_OUTPUT_BUFFER_SIZE (MAX_FRAME_HEIGHT * MAX_FRAME_WIDTH)

#define ive_api_function ih264e_api_function

constexpr int16_t kCompressionRatio = 1;
constexpr size_t kMinQP = 4;
constexpr uint32_t kHeaderLength = 0x800;

static TestArgs* gArgs = nullptr;

class AvcEncTest
    : public ::testing::TestWithParam<tuple<string, int32_t, int32_t, float, int32_t>> {
  private:
    void setRawBuf(iv_raw_buf_t* psInpRawBuf, const uint8_t* data);
    void setFrameType(IV_PICTURE_CODING_TYPE_T eFrameType);
    void setQp();
    void setEncMode(IVE_ENC_MODE_T eEncMode);
    void setDimensions();
    void setNumCores();
    void setFrameRate();
    void setIpeParams();
    void setBitRate();
    void setAirParams();
    void setMeParams();
    void setGopParams();
    void setProfileParams();
    void setDeblockParams();
    void setVbvParams();
    void setDefault();
    void setVuiParams();
    void getBufInfo();
    void setSeiMdcvParams();
    void setSeiCllParams();
    void setSeiAveParams();
    void setSeiCcvParams();
    void logVersion();
    bool mHalfPelEnable = true;
    bool mQPelEnable = true;
    bool mIntra4x4 = true;
    bool mEnableFastSad = false;
    bool mEnableAltRef = false;
    bool mConstrainedIntraFlag = false;
    bool mSeiCllFlag = true;
    bool mSeiAveFlag = true;
    bool mSeiCcvFlag = true;
    bool mSeiMdcvFlag = true;
    bool mAspectRatioFlag = false;
    bool mNalHrdFlag = false;
    bool mVclHrdFlag = false;
    bool mIsForceIdrEnabled = false;
    bool mIsDynamicBitRateChangeEnabled = true;
    bool mIsDynamicFrameRateChangeEnabled = true;
    uint32_t mAvcEncLevel = 41;
    uint32_t mNumMemRecords = 0;
    uint32_t mNumCores = 4;
    uint32_t mBframes = 0;
    uint32_t mSliceParam = 256;
    uint32_t mMeSpeedPreset = 100;
    uint32_t mIInterval = 1000;
    uint32_t mIDRInterval = 60;
    uint32_t mDisableDeblockLevel = 0;
    uint32_t mIQp = 24;
    uint32_t mPQp = 27;
    uint32_t mBQp = 29;
    uint32_t mIntraRefresh = 30;
    uint32_t mSearchRangeX = 16;
    uint32_t mSearchRangeY = 16;
    uint32_t mForceIdrInterval = 0;          // in number of frames
    uint32_t mDynamicBitRateInterval = 0;    // in number of frames
    uint32_t mDynamicFrameRateInterval = 0;  // in number of frame
    float mFrameRate = 30;
    iv_obj_t* mCodecCtx = nullptr;
    iv_mem_rec_t* mMemRecords = nullptr;
    IVE_AIR_MODE_T mAirMode = IVE_AIR_MODE_NONE;
    IVE_SPEED_CONFIG mEncSpeed = IVE_NORMAL;
    IVE_RC_MODE_T mRCMode = IVE_RC_STORAGE;
    IV_ARCH_T mArch = ARCH_NA;
    IVE_SLICE_MODE_T mSliceMode = IVE_SLICE_MODE_NONE;
    IV_COLOR_FORMAT_T mIvVideoColorFormat = IV_YUV_420P;
    IV_COLOR_FORMAT_T mReconFormat = IV_YUV_420P;
    IV_PROFILE_T mProfile = IV_PROFILE_BASE;

  public:
    AvcEncTest()
        : mInputBuffer(nullptr), mOutputBuffer(nullptr), mFpInput(nullptr), mFpOutput(nullptr) {}

    ~AvcEncTest() {
        iv_mem_rec_t* ps_mem_rec = mMemRecords;
        for (size_t i = 0; i < mNumMemRecords; ++i) {
            if (ps_mem_rec) {
                free(ps_mem_rec->pv_base);
            }
            ++ps_mem_rec;
        }
        if (mMemRecords) {
            free(mMemRecords);
        }
        mCodecCtx = nullptr;
    }

    void SetUp() override {
        tuple<string /* fileName */, int32_t /* frameWidth */, int32_t /* frameHeight */,
              float /* frameRate */, int32_t /* bitRate */>
                params = GetParam();
        mFileName = gArgs->getRes() + get<0>(params);
        mFrameWidth = get<1>(params);
        mFrameHeight = get<2>(params);
        mFrameRate = get<3>(params);
        mBitRate = get<4>(params);
        mOutFileName = gArgs->getRes() + "out.bin";

        ASSERT_LE(mFrameWidth, 1080) << "Frame Width <= 1080";

        ASSERT_LE(mFrameHeight, 1920) << "Frame Height <= 1920";

        mOutputBufferSize = (mFrameWidth * mFrameHeight * 3 / 2) / kCompressionRatio;
        mBitRate = mBitRate * 1024;  // Conversion to bytes per sec

        mInputBuffer = (uint8_t*)malloc((mFrameWidth * mFrameHeight * 3) / 2);
        ASSERT_NE(mInputBuffer, nullptr) << "Failed to allocate the input buffer!";

        mOutputBuffer = (uint8_t*)malloc(mOutputBufferSize);
        ASSERT_NE(mOutputBuffer, nullptr) << "Failed to allocate the output buffer!";

        mFpInput = fopen(mFileName.c_str(), "rb");
        ASSERT_NE(mFpInput, nullptr) << "Failed to open the input file: " << mFileName;

        mFpOutput = fopen(mOutFileName.c_str(), "wb");
        ASSERT_NE(mFpOutput, nullptr) << "Failed to open the output file:" << mOutFileName;

        /* Getting Number of MemRecords */
        iv_num_mem_rec_ip_t sNumMemRecIp = {};
        iv_num_mem_rec_op_t sNumMemRecOp = {};

        sNumMemRecIp.u4_size = sizeof(iv_num_mem_rec_ip_t);
        sNumMemRecOp.u4_size = sizeof(iv_num_mem_rec_op_t);
        sNumMemRecIp.e_cmd = IV_CMD_GET_NUM_MEM_REC;

        status = ive_api_function(nullptr, &sNumMemRecIp, &sNumMemRecOp);
        ASSERT_EQ(status, IV_SUCCESS) << "Error in IV_CMD_GET_NUM_MEM_REC!";

        mNumMemRecords = sNumMemRecOp.u4_num_mem_rec;
        mMemRecords = (iv_mem_rec_t*)malloc(mNumMemRecords * sizeof(iv_mem_rec_t));
        ASSERT_NE(mMemRecords, nullptr) << "Failed to allocate memory to nMemRecords!";

        iv_mem_rec_t* psMemRec;
        psMemRec = mMemRecords;
        for (size_t i = 0; i < mNumMemRecords; ++i) {
            psMemRec->u4_size = sizeof(iv_mem_rec_t);
            psMemRec->pv_base = nullptr;
            psMemRec->u4_mem_size = 0;
            psMemRec->u4_mem_alignment = 0;
            psMemRec->e_mem_type = IV_NA_MEM_TYPE;
            ++psMemRec;
        }

        /* Getting MemRecords Attributes */
        iv_fill_mem_rec_ip_t sFillMemRecIp = {};
        iv_fill_mem_rec_op_t sFillMemRecOp = {};

        sFillMemRecIp.u4_size = sizeof(iv_fill_mem_rec_ip_t);
        sFillMemRecOp.u4_size = sizeof(iv_fill_mem_rec_op_t);

        sFillMemRecIp.e_cmd = IV_CMD_FILL_NUM_MEM_REC;
        sFillMemRecIp.ps_mem_rec = mMemRecords;
        sFillMemRecIp.u4_num_mem_rec = mNumMemRecords;
        sFillMemRecIp.u4_max_wd = mFrameWidth;
        sFillMemRecIp.u4_max_ht = mFrameHeight;
        sFillMemRecIp.u4_max_level = mAvcEncLevel;
        sFillMemRecIp.e_color_format = IV_YUV_420SP_UV;
        sFillMemRecIp.u4_max_ref_cnt = 2;
        sFillMemRecIp.u4_max_reorder_cnt = 0;
        sFillMemRecIp.u4_max_srch_rng_x = 256;
        sFillMemRecIp.u4_max_srch_rng_y = 256;

        status = ive_api_function(nullptr, &sFillMemRecIp, &sFillMemRecOp);
        ASSERT_EQ(status, IV_SUCCESS) << "Failed to fill memory records!";

        /* Allocating Memory for Mem Records */
        psMemRec = mMemRecords;
        for (size_t i = 0; i < mNumMemRecords; ++i) {
            posix_memalign(&psMemRec->pv_base, psMemRec->u4_mem_alignment, psMemRec->u4_mem_size);
            ASSERT_NE(psMemRec->pv_base, nullptr)
                    << "Failed to allocate for size " << psMemRec->u4_mem_size;

            ++psMemRec;
        }

        /* Codec Instance Creation */
        ive_init_ip_t sInitIp = {};
        ive_init_op_t sInitOp = {};

        mCodecCtx = (iv_obj_t*)mMemRecords[0].pv_base;
        mCodecCtx->u4_size = sizeof(iv_obj_t);
        mCodecCtx->pv_fxns = (void*)ive_api_function;

        sInitIp.u4_size = sizeof(ive_init_ip_t);
        sInitOp.u4_size = sizeof(ive_init_op_t);

        sInitIp.e_cmd = IV_CMD_INIT;
        sInitIp.u4_num_mem_rec = mNumMemRecords;
        sInitIp.ps_mem_rec = mMemRecords;
        sInitIp.u4_max_wd = mFrameWidth;
        sInitIp.u4_max_ht = mFrameHeight;
        sInitIp.u4_max_ref_cnt = 2;
        sInitIp.u4_max_reorder_cnt = 0;
        sInitIp.u4_max_level = mAvcEncLevel;
        sInitIp.e_inp_color_fmt = mIvVideoColorFormat;
        sInitIp.u4_enable_recon = 0;
        sInitIp.e_recon_color_fmt = mReconFormat;
        sInitIp.e_rc_mode = mRCMode;
        sInitIp.u4_max_framerate = 120000;
        sInitIp.u4_max_bitrate = 240000000;
        sInitIp.u4_num_bframes = mBframes;
        sInitIp.e_content_type = IV_PROGRESSIVE;
        sInitIp.u4_max_srch_rng_x = 256;
        sInitIp.u4_max_srch_rng_y = 256;
        sInitIp.e_slice_mode = mSliceMode;
        sInitIp.u4_slice_param = mSliceParam;
        sInitIp.e_arch = mArch;
        sInitIp.e_soc = SOC_GENERIC;

        status = ive_api_function(mCodecCtx, &sInitIp, &sInitOp);
        ASSERT_EQ(status, IV_SUCCESS) << "Failed to create Codec Instance!";

        mFrameSize = (mIvVideoColorFormat == IV_YUV_422ILE)
                             ? (mFrameWidth * mFrameHeight * 2)
                             : ((mFrameWidth * mFrameHeight * 3) / 2);

        mTotalFrames = getTotalFrames();

        ASSERT_NO_FATAL_FAILURE(logVersion());

        ASSERT_NO_FATAL_FAILURE(setDefault());

        ASSERT_NO_FATAL_FAILURE(getBufInfo());

        ASSERT_NO_FATAL_FAILURE(setNumCores());

        ASSERT_NO_FATAL_FAILURE(setDimensions());

        ASSERT_NO_FATAL_FAILURE(setFrameRate());

        ASSERT_NO_FATAL_FAILURE(setIpeParams());

        ASSERT_NO_FATAL_FAILURE(setBitRate());

        ASSERT_NO_FATAL_FAILURE(setQp());

        ASSERT_NO_FATAL_FAILURE(setAirParams());

        ASSERT_NO_FATAL_FAILURE(setVbvParams());

        ASSERT_NO_FATAL_FAILURE(setMeParams());

        ASSERT_NO_FATAL_FAILURE(setGopParams());

        ASSERT_NO_FATAL_FAILURE(setDeblockParams());

        ASSERT_NO_FATAL_FAILURE(setVuiParams());

        ASSERT_NO_FATAL_FAILURE(setSeiMdcvParams());

        ASSERT_NO_FATAL_FAILURE(setSeiCllParams());

        ASSERT_NO_FATAL_FAILURE(setSeiAveParams());

        ASSERT_NO_FATAL_FAILURE(setSeiCcvParams());

        ASSERT_NO_FATAL_FAILURE(setProfileParams());

        ASSERT_NO_FATAL_FAILURE(setEncMode(IVE_ENC_MODE_PICTURE));
    }

    void TearDown() override {
        if (mInputBuffer) free(mInputBuffer);
        if (mOutputBuffer) free(mOutputBuffer);
        if (mFpInput) fclose(mFpInput);
        if (mFpOutput) fclose(mFpOutput);
    }

    void encodeFrames(int64_t);
    int64_t getTotalFrames();
    int32_t mFrameWidth = MAX_FRAME_WIDTH;
    int32_t mFrameHeight = MAX_FRAME_HEIGHT;
    int32_t mFrameSize = (mFrameWidth * mFrameHeight * 3) / 2;
    int64_t mTotalFrames = 0;
    int32_t mBitRate = 256000;
    int64_t mOutputBufferSize = MAX_OUTPUT_BUFFER_SIZE;
    string mFileName;
    string mOutFileName;
    uint8_t* mInputBuffer = nullptr;
    uint8_t* mOutputBuffer = nullptr;
    FILE* mFpInput = nullptr;
    FILE* mFpOutput = nullptr;
    IV_STATUS_T status;
};

void AvcEncTest::setDimensions() {
    ive_ctl_set_dimensions_ip_t sDimensionsIp = {};
    ive_ctl_set_dimensions_op_t sDimensionsOp = {};

    sDimensionsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sDimensionsIp.e_sub_cmd = IVE_CMD_CTL_SET_DIMENSIONS;
    sDimensionsIp.u4_ht = mFrameHeight;
    sDimensionsIp.u4_wd = mFrameWidth;

    sDimensionsIp.u4_timestamp_high = -1;
    sDimensionsIp.u4_timestamp_low = -1;

    sDimensionsIp.u4_size = sizeof(ive_ctl_set_dimensions_ip_t);
    sDimensionsOp.u4_size = sizeof(ive_ctl_set_dimensions_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sDimensionsIp, &sDimensionsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set dimensions!\n";

    return;
}

void AvcEncTest::setNumCores() {
    ive_ctl_set_num_cores_ip_t sNumCoresIp = {};
    ive_ctl_set_num_cores_op_t sNumCoresOp = {};

    sNumCoresIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sNumCoresIp.e_sub_cmd = IVE_CMD_CTL_SET_NUM_CORES;
    sNumCoresIp.u4_num_cores = mNumCores;

    sNumCoresIp.u4_timestamp_high = -1;
    sNumCoresIp.u4_timestamp_low = -1;

    sNumCoresIp.u4_size = sizeof(ive_ctl_set_num_cores_ip_t);
    sNumCoresOp.u4_size = sizeof(ive_ctl_set_num_cores_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, (void*)&sNumCoresIp, (void*)&sNumCoresOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set number of cores!\n";

    return;
}

void AvcEncTest::setDefault() {
    ive_ctl_setdefault_ip_t sDefaultIp = {};
    ive_ctl_setdefault_op_t sDefaultOp = {};

    sDefaultIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sDefaultIp.e_sub_cmd = IVE_CMD_CTL_SETDEFAULT;

    sDefaultIp.u4_timestamp_high = -1;
    sDefaultIp.u4_timestamp_low = -1;

    sDefaultIp.u4_size = sizeof(ive_ctl_setdefault_ip_t);
    sDefaultOp.u4_size = sizeof(ive_ctl_setdefault_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sDefaultIp, &sDefaultOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set default encoder parameters!\n";

    return;
}

void AvcEncTest::getBufInfo() {
    ih264e_ctl_getbufinfo_ip_t sGetBufInfoIp = {};
    ih264e_ctl_getbufinfo_op_t sGetBufInfoOp = {};

    sGetBufInfoIp.s_ive_ip.u4_size = sizeof(ih264e_ctl_getbufinfo_ip_t);
    sGetBufInfoOp.s_ive_op.u4_size = sizeof(ih264e_ctl_getbufinfo_op_t);

    sGetBufInfoIp.s_ive_ip.e_cmd = IVE_CMD_VIDEO_CTL;
    sGetBufInfoIp.s_ive_ip.e_sub_cmd = IVE_CMD_CTL_GETBUFINFO;
    sGetBufInfoIp.s_ive_ip.u4_max_ht = mFrameHeight;
    sGetBufInfoIp.s_ive_ip.u4_max_wd = mFrameWidth;
    sGetBufInfoIp.s_ive_ip.e_inp_color_fmt = mIvVideoColorFormat;

    IV_STATUS_T status = ih264e_api_function(mCodecCtx, &sGetBufInfoIp, &sGetBufInfoOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to get buffer info!\n";

    return;
}

void AvcEncTest::setFrameRate() {
    ive_ctl_set_frame_rate_ip_t sFrameRateIp = {};
    ive_ctl_set_frame_rate_op_t sFrameRateOp = {};

    sFrameRateIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sFrameRateIp.e_sub_cmd = IVE_CMD_CTL_SET_FRAMERATE;
    sFrameRateIp.u4_src_frame_rate = mFrameRate;
    sFrameRateIp.u4_tgt_frame_rate = mFrameRate;

    sFrameRateIp.u4_timestamp_high = -1;
    sFrameRateIp.u4_timestamp_low = -1;

    sFrameRateIp.u4_size = sizeof(ive_ctl_set_frame_rate_ip_t);
    sFrameRateOp.u4_size = sizeof(ive_ctl_set_frame_rate_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sFrameRateIp, &sFrameRateOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set frame rate!\n";

    return;
}

void AvcEncTest::setIpeParams() {
    ive_ctl_set_ipe_params_ip_t sIpeParamsIp = {};
    ive_ctl_set_ipe_params_op_t sIpeParamsOp = {};

    sIpeParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sIpeParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_IPE_PARAMS;
    sIpeParamsIp.u4_enable_intra_4x4 = mIntra4x4;
    sIpeParamsIp.u4_enc_speed_preset = mEncSpeed;
    sIpeParamsIp.u4_constrained_intra_pred = mConstrainedIntraFlag;

    sIpeParamsIp.u4_timestamp_high = -1;
    sIpeParamsIp.u4_timestamp_low = -1;

    sIpeParamsIp.u4_size = sizeof(ive_ctl_set_ipe_params_ip_t);
    sIpeParamsOp.u4_size = sizeof(ive_ctl_set_ipe_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sIpeParamsIp, &sIpeParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set IPE params!\n";

    return;
}

void AvcEncTest::setBitRate() {
    ive_ctl_set_bitrate_ip_t sBitrateIp = {};
    ive_ctl_set_bitrate_op_t sBitrateOp = {};

    sBitrateIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sBitrateIp.e_sub_cmd = IVE_CMD_CTL_SET_BITRATE;
    sBitrateIp.u4_target_bitrate = mBitRate;

    sBitrateIp.u4_timestamp_high = -1;
    sBitrateIp.u4_timestamp_low = -1;

    sBitrateIp.u4_size = sizeof(ive_ctl_set_bitrate_ip_t);
    sBitrateOp.u4_size = sizeof(ive_ctl_set_bitrate_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sBitrateIp, &sBitrateOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set bit rate!\n";

    return;
}

void AvcEncTest::setFrameType(IV_PICTURE_CODING_TYPE_T eFrameType) {
    ive_ctl_set_frame_type_ip_t sFrameTypeIp = {};
    ive_ctl_set_frame_type_op_t sFrameTypeOp = {};

    sFrameTypeIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sFrameTypeIp.e_sub_cmd = IVE_CMD_CTL_SET_FRAMETYPE;
    sFrameTypeIp.e_frame_type = eFrameType;

    sFrameTypeIp.u4_timestamp_high = -1;
    sFrameTypeIp.u4_timestamp_low = -1;

    sFrameTypeIp.u4_size = sizeof(ive_ctl_set_frame_type_ip_t);
    sFrameTypeOp.u4_size = sizeof(ive_ctl_set_frame_type_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sFrameTypeIp, &sFrameTypeOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set Frame Type!\n";
    return;
}

void AvcEncTest::setQp() {
    ive_ctl_set_qp_ip_t s_QpIp = {};
    ive_ctl_set_qp_op_t s_QpOp = {};

    s_QpIp.e_cmd = IVE_CMD_VIDEO_CTL;
    s_QpIp.e_sub_cmd = IVE_CMD_CTL_SET_QP;

    s_QpIp.u4_i_qp = mIQp;
    s_QpIp.u4_i_qp_max = MAX_H264_QP;
    s_QpIp.u4_i_qp_min = kMinQP;

    s_QpIp.u4_p_qp = mPQp;
    s_QpIp.u4_p_qp_max = MAX_H264_QP;
    s_QpIp.u4_p_qp_min = kMinQP;

    s_QpIp.u4_b_qp = mBQp;
    s_QpIp.u4_b_qp_max = MAX_H264_QP;
    s_QpIp.u4_b_qp_min = kMinQP;

    s_QpIp.u4_timestamp_high = -1;
    s_QpIp.u4_timestamp_low = -1;

    s_QpIp.u4_size = sizeof(ive_ctl_set_qp_ip_t);
    s_QpOp.u4_size = sizeof(ive_ctl_set_qp_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &s_QpIp, &s_QpOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set QP!\n";

    return;
}

void AvcEncTest::setEncMode(IVE_ENC_MODE_T eEncMode) {
    ive_ctl_set_enc_mode_ip_t sEncModeIp = {};
    ive_ctl_set_enc_mode_op_t sEncModeOp = {};

    sEncModeIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sEncModeIp.e_sub_cmd = IVE_CMD_CTL_SET_ENC_MODE;
    sEncModeIp.e_enc_mode = eEncMode;

    sEncModeIp.u4_timestamp_high = -1;
    sEncModeIp.u4_timestamp_low = -1;

    sEncModeIp.u4_size = sizeof(ive_ctl_set_enc_mode_ip_t);
    sEncModeOp.u4_size = sizeof(ive_ctl_set_enc_mode_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sEncModeIp, &sEncModeOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set encode mode!\n";

    return;
}

void AvcEncTest::setVbvParams() {
    ive_ctl_set_vbv_params_ip_t sVbvIp = {};
    ive_ctl_set_vbv_params_op_t sVbvOp = {};

    sVbvIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sVbvIp.e_sub_cmd = IVE_CMD_CTL_SET_VBV_PARAMS;
    sVbvIp.u4_vbv_buf_size = 0;
    sVbvIp.u4_vbv_buffer_delay = 1000;

    sVbvIp.u4_timestamp_high = -1;
    sVbvIp.u4_timestamp_low = -1;

    sVbvIp.u4_size = sizeof(ive_ctl_set_vbv_params_ip_t);
    sVbvOp.u4_size = sizeof(ive_ctl_set_vbv_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sVbvIp, &sVbvOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set VBV params!\n";

    return;
}

void AvcEncTest::setAirParams() {
    ive_ctl_set_air_params_ip_t sAirIp = {};
    ive_ctl_set_air_params_op_t sAirOp = {};

    sAirIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sAirIp.e_sub_cmd = IVE_CMD_CTL_SET_AIR_PARAMS;
    sAirIp.e_air_mode = mAirMode;
    sAirIp.u4_air_refresh_period = mIntraRefresh;

    sAirIp.u4_timestamp_high = -1;
    sAirIp.u4_timestamp_low = -1;

    sAirIp.u4_size = sizeof(ive_ctl_set_air_params_ip_t);
    sAirOp.u4_size = sizeof(ive_ctl_set_air_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sAirIp, &sAirOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set AIR params!\n";

    return;
}

void AvcEncTest::setMeParams() {
    ive_ctl_set_me_params_ip_t sMeParamsIp = {};
    ive_ctl_set_me_params_op_t sMeParamsOp = {};

    sMeParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sMeParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_ME_PARAMS;
    sMeParamsIp.u4_enable_fast_sad = mEnableFastSad;
    sMeParamsIp.u4_enable_alt_ref = mEnableAltRef;

    sMeParamsIp.u4_enable_hpel = mHalfPelEnable;
    sMeParamsIp.u4_enable_qpel = mQPelEnable;
    sMeParamsIp.u4_me_speed_preset = mMeSpeedPreset;
    sMeParamsIp.u4_srch_rng_x = mSearchRangeX;
    sMeParamsIp.u4_srch_rng_y = mSearchRangeY;

    sMeParamsIp.u4_timestamp_high = -1;
    sMeParamsIp.u4_timestamp_low = -1;

    sMeParamsIp.u4_size = sizeof(ive_ctl_set_me_params_ip_t);
    sMeParamsOp.u4_size = sizeof(ive_ctl_set_me_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sMeParamsIp, &sMeParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set ME params!\n";

    return;
}

void AvcEncTest::setGopParams() {
    ive_ctl_set_gop_params_ip_t sGopParamsIp = {};
    ive_ctl_set_gop_params_op_t sGopParamsOp = {};

    sGopParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sGopParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_GOP_PARAMS;

    sGopParamsIp.u4_i_frm_interval = mIInterval;
    sGopParamsIp.u4_idr_frm_interval = mIDRInterval;

    sGopParamsIp.u4_timestamp_high = -1;
    sGopParamsIp.u4_timestamp_low = -1;

    sGopParamsIp.u4_size = sizeof(ive_ctl_set_gop_params_ip_t);
    sGopParamsOp.u4_size = sizeof(ive_ctl_set_gop_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sGopParamsIp, &sGopParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set GOP params!\n";

    return;
}

void AvcEncTest::setProfileParams() {
    ive_ctl_set_profile_params_ip_t sProfileParamsIp = {};
    ive_ctl_set_profile_params_op_t sProfileParamsOp = {};

    sProfileParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sProfileParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_PROFILE_PARAMS;

    sProfileParamsIp.e_profile = mProfile;
    if (sProfileParamsIp.e_profile == IV_PROFILE_BASE) {
        sProfileParamsIp.u4_entropy_coding_mode = 0;
    } else {
        sProfileParamsIp.u4_entropy_coding_mode = 1;
    }
    sProfileParamsIp.u4_timestamp_high = -1;
    sProfileParamsIp.u4_timestamp_low = -1;

    sProfileParamsIp.u4_size = sizeof(ive_ctl_set_profile_params_ip_t);
    sProfileParamsOp.u4_size = sizeof(ive_ctl_set_profile_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sProfileParamsIp, &sProfileParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set profile params!\n";

    return;
}

void AvcEncTest::setDeblockParams() {
    ive_ctl_set_deblock_params_ip_t sDeblockParamsIp = {};
    ive_ctl_set_deblock_params_op_t sDeblockParamsOp = {};

    sDeblockParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sDeblockParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_DEBLOCK_PARAMS;

    sDeblockParamsIp.u4_disable_deblock_level = mDisableDeblockLevel;

    sDeblockParamsIp.u4_timestamp_high = -1;
    sDeblockParamsIp.u4_timestamp_low = -1;

    sDeblockParamsIp.u4_size = sizeof(ive_ctl_set_deblock_params_ip_t);
    sDeblockParamsOp.u4_size = sizeof(ive_ctl_set_deblock_params_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sDeblockParamsIp, &sDeblockParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set deblock params!\n";

    return;
}

void AvcEncTest::setVuiParams() {
    ih264e_vui_ip_t sVuiParamsIp = {};
    ih264e_vui_op_t sVuiParamsOp = {};

    sVuiParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sVuiParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_VUI_PARAMS;

    sVuiParamsIp.u1_aspect_ratio_info_present_flag = mAspectRatioFlag;
    sVuiParamsIp.u1_video_signal_type_present_flag = 1;
    sVuiParamsIp.u1_colour_description_present_flag = 1;
    sVuiParamsIp.u1_nal_hrd_parameters_present_flag = mNalHrdFlag;
    sVuiParamsIp.u1_vcl_hrd_parameters_present_flag = mVclHrdFlag;

    sVuiParamsIp.u4_size = sizeof(ih264e_vui_ip_t);
    sVuiParamsOp.u4_size = sizeof(ih264e_vui_op_t);

    IV_STATUS_T status = ive_api_function(mCodecCtx, &sVuiParamsIp, &sVuiParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set VUI params!\n";

    return;
}

void AvcEncTest::setSeiMdcvParams() {
    ih264e_ctl_set_sei_mdcv_params_ip_t sSeiMdcvParamsIp = {};
    ih264e_ctl_set_sei_mdcv_params_op_t sSeiMdcvParamsOp = {};

    sSeiMdcvParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sSeiMdcvParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_SEI_MDCV_PARAMS;
    sSeiMdcvParamsIp.u1_sei_mdcv_params_present_flag = mSeiMdcvFlag;
    if (mSeiMdcvFlag) {
        for (int i4_count = 0; i4_count < NUM_SEI_MDCV_PRIMARIES; ++i4_count) {
            sSeiMdcvParamsIp.au2_display_primaries_x[i4_count] = 30000;
            sSeiMdcvParamsIp.au2_display_primaries_y[i4_count] = 35000;
        }
        sSeiMdcvParamsIp.u2_white_point_x = 30000;
        sSeiMdcvParamsIp.u2_white_point_y = 35000;
        sSeiMdcvParamsIp.u4_max_display_mastering_luminance = 100000000;
        sSeiMdcvParamsIp.u4_min_display_mastering_luminance = 50000;
    }

    sSeiMdcvParamsIp.u4_timestamp_high = -1;
    sSeiMdcvParamsIp.u4_timestamp_low = -1;

    sSeiMdcvParamsIp.u4_size = sizeof(ih264e_ctl_set_sei_mdcv_params_ip_t);
    sSeiMdcvParamsOp.u4_size = sizeof(ih264e_ctl_set_sei_mdcv_params_op_t);
    IV_STATUS_T status = ih264e_api_function(mCodecCtx, &sSeiMdcvParamsIp, &sSeiMdcvParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set SEI MDCV params!\n";

    return;
}

void AvcEncTest::setSeiCllParams() {
    ih264e_ctl_set_sei_cll_params_ip_t sSeiCllParamsIp = {};
    ih264e_ctl_set_sei_cll_params_op_t sSeiCllParamsOp = {};

    sSeiCllParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sSeiCllParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_SEI_CLL_PARAMS;
    sSeiCllParamsIp.u1_sei_cll_params_present_flag = mSeiCllFlag;
    if (mSeiCllFlag) {
        sSeiCllParamsIp.u2_max_content_light_level = 0;
        sSeiCllParamsIp.u2_max_pic_average_light_level = 0;
    }

    sSeiCllParamsIp.u4_timestamp_high = -1;
    sSeiCllParamsIp.u4_timestamp_low = -1;

    sSeiCllParamsIp.u4_size = sizeof(ih264e_ctl_set_sei_cll_params_ip_t);
    sSeiCllParamsOp.u4_size = sizeof(ih264e_ctl_set_sei_cll_params_op_t);

    IV_STATUS_T status = ih264e_api_function(mCodecCtx, &sSeiCllParamsIp, &sSeiCllParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set SEI CLL params!\n";

    return;
}

void AvcEncTest::setSeiAveParams() {
    ih264e_ctl_set_sei_ave_params_ip_t sSeiAveParamsIp = {};
    ih264e_ctl_set_sei_ave_params_op_t sSeiAveParamsOp = {};

    sSeiAveParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sSeiAveParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_SEI_AVE_PARAMS;
    sSeiAveParamsIp.u1_sei_ave_params_present_flag = mSeiAveFlag;
    if (mSeiAveFlag) {
        sSeiAveParamsIp.u4_ambient_illuminance = 1;
        sSeiAveParamsIp.u2_ambient_light_x = 0;
        sSeiAveParamsIp.u2_ambient_light_y = 0;
    }

    sSeiAveParamsIp.u4_timestamp_high = -1;
    sSeiAveParamsIp.u4_timestamp_low = -1;

    sSeiAveParamsIp.u4_size = sizeof(ih264e_ctl_set_sei_ave_params_ip_t);
    sSeiAveParamsOp.u4_size = sizeof(ih264e_ctl_set_sei_ave_params_op_t);

    IV_STATUS_T status = ih264e_api_function(mCodecCtx, &sSeiAveParamsIp, &sSeiAveParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set SEI AVE params!\n";

    return;
}

void AvcEncTest::setSeiCcvParams() {
    ih264e_ctl_set_sei_ccv_params_ip_t sSeiCcvParamsIp = {};
    ih264e_ctl_set_sei_ccv_params_op_t sSeiCcvParamsOp = {};

    sSeiCcvParamsIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sSeiCcvParamsIp.e_sub_cmd = IVE_CMD_CTL_SET_SEI_CCV_PARAMS;
    sSeiCcvParamsIp.u1_sei_ccv_params_present_flag = mSeiCcvFlag;
    if (mSeiCcvFlag) {
        sSeiCcvParamsIp.u1_ccv_cancel_flag = 0;
        sSeiCcvParamsIp.u1_ccv_persistence_flag = 1;
        sSeiCcvParamsIp.u1_ccv_primaries_present_flag = 1;
        sSeiCcvParamsIp.u1_ccv_min_luminance_value_present_flag = 1;
        sSeiCcvParamsIp.u1_ccv_max_luminance_value_present_flag = 1;
        sSeiCcvParamsIp.u1_ccv_avg_luminance_value_present_flag = 1;
        sSeiCcvParamsIp.u1_ccv_reserved_zero_2bits = 0;
        for (int i4_count = 0; i4_count < NUM_SEI_CCV_PRIMARIES; ++i4_count) {
            sSeiCcvParamsIp.ai4_ccv_primaries_x[i4_count] = 1;
            sSeiCcvParamsIp.ai4_ccv_primaries_y[i4_count] = 1;
        }
        sSeiCcvParamsIp.u4_ccv_min_luminance_value = 1;
        sSeiCcvParamsIp.u4_ccv_max_luminance_value = 1;
        sSeiCcvParamsIp.u4_ccv_avg_luminance_value = 1;
    }

    sSeiCcvParamsIp.u4_timestamp_high = -1;
    sSeiCcvParamsIp.u4_timestamp_low = -1;

    sSeiCcvParamsIp.u4_size = sizeof(ih264e_ctl_set_sei_ccv_params_ip_t);
    sSeiCcvParamsOp.u4_size = sizeof(ih264e_ctl_set_sei_ccv_params_op_t);

    IV_STATUS_T status = ih264e_api_function(mCodecCtx, &sSeiCcvParamsIp, &sSeiCcvParamsOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to set SEI CCV params!\n";

    return;
}

void AvcEncTest::logVersion() {
    ive_ctl_getversioninfo_ip_t sCtlIp = {};
    ive_ctl_getversioninfo_op_t sCtlOp = {};
    UWORD8 au1Buf[512];

    sCtlIp.e_cmd = IVE_CMD_VIDEO_CTL;
    sCtlIp.e_sub_cmd = IVE_CMD_CTL_GETVERSION;

    sCtlIp.u4_size = sizeof(ive_ctl_getversioninfo_ip_t);
    sCtlOp.u4_size = sizeof(ive_ctl_getversioninfo_op_t);
    sCtlIp.pu1_version = au1Buf;
    sCtlIp.u4_version_bufsize = sizeof(au1Buf);

    IV_STATUS_T status = ive_api_function(mCodecCtx, (void*)&sCtlIp, (void*)&sCtlOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to get encoder version!\n";

    return;
}

void AvcEncTest::encodeFrames(int64_t numFramesToEncode) {
    ih264e_video_encode_ip_t ih264e_video_encode_ip = {};
    ih264e_video_encode_op_t ih264e_video_encode_op = {};

    ive_video_encode_ip_t* sEncodeIp = &ih264e_video_encode_ip.s_ive_ip;
    ive_video_encode_op_t* sEncodeOp = &ih264e_video_encode_op.s_ive_op;

    uint8_t header[kHeaderLength];
    iv_raw_buf_t* psInpRawBuf = &sEncodeIp->s_inp_buf;

    sEncodeIp->s_out_buf.pv_buf = header;
    sEncodeIp->s_out_buf.u4_bytes = 0;
    sEncodeIp->s_out_buf.u4_bufsize = kHeaderLength;
    sEncodeIp->u4_size = sizeof(ih264e_video_encode_ip_t);
    sEncodeOp->u4_size = sizeof(ih264e_video_encode_op_t);

    sEncodeIp->e_cmd = IVE_CMD_VIDEO_ENCODE;
    sEncodeIp->pv_bufs = nullptr;
    sEncodeIp->pv_mb_info = nullptr;
    sEncodeIp->pv_pic_info = nullptr;
    sEncodeIp->u4_mb_info_type = 0;
    sEncodeIp->u4_pic_info_type = 0;
    sEncodeIp->u4_is_last = 0;
    sEncodeOp->s_out_buf.pv_buf = nullptr;

    /* Initialize color formats */
    memset(psInpRawBuf, 0, sizeof(iv_raw_buf_t));
    psInpRawBuf->u4_size = sizeof(iv_raw_buf_t);
    psInpRawBuf->e_color_fmt = mIvVideoColorFormat;

    IV_STATUS_T status = ive_api_function(mCodecCtx, sEncodeIp, sEncodeOp);
    ASSERT_EQ(status, IV_SUCCESS) << "Failed to Initialize Color Formats!\n";

    uint32_t numFrame = 0;

    while (numFramesToEncode > 0) {
        int32_t bytesRead;
        bytesRead = fread(mInputBuffer, 1, mFrameSize, mFpInput);

        if (bytesRead != mFrameSize) {
            break;
        }

        setRawBuf(psInpRawBuf, mInputBuffer);

        sEncodeIp->s_out_buf.pv_buf = mOutputBuffer;
        sEncodeIp->s_out_buf.u4_bufsize = mFrameSize;
        if (mIsForceIdrEnabled) {
            if (numFrame == mForceIdrInterval) {
                ASSERT_NO_FATAL_FAILURE(setFrameType(IV_IDR_FRAME));
            }
        }
        if (mIsDynamicBitRateChangeEnabled) {
            if (numFrame == mDynamicBitRateInterval) {
                mBitRate *= 2;
            }
            ASSERT_NO_FATAL_FAILURE(setBitRate());
        }
        if (mIsDynamicFrameRateChangeEnabled) {
            if (numFrame == mDynamicFrameRateInterval) {
                mFrameRate *= 2;
            }
            ASSERT_NO_FATAL_FAILURE(setFrameRate());
        }

        status = ive_api_function(mCodecCtx, &ih264e_video_encode_ip, &ih264e_video_encode_op);
        ASSERT_EQ(status, IV_SUCCESS) << "Failed to encode frame!\n";

        int32_t numOutputBytes = fwrite((UWORD8*)sEncodeOp->s_out_buf.pv_buf, sizeof(UWORD8),
                                        sEncodeOp->s_out_buf.u4_bytes, mFpOutput);
        ASSERT_NE(numOutputBytes, 0) << "Failed to write the output!" << mOutFileName;

        numFramesToEncode--;
        numFrame++;
    }

    sEncodeIp->u4_is_last = 1;
    psInpRawBuf->apv_bufs[0] = nullptr;
    psInpRawBuf->apv_bufs[1] = nullptr;
    psInpRawBuf->apv_bufs[2] = nullptr;

    status = ive_api_function(mCodecCtx, &ih264e_video_encode_ip, &ih264e_video_encode_op);
    ASSERT_EQ(status, IV_SUCCESS) << "Failure after encoding last frame!\n";

    if (sEncodeOp->output_present) {
        int32_t numOutputBytes = fwrite((UWORD8*)sEncodeOp->s_out_buf.pv_buf, sizeof(UWORD8),
                                        sEncodeOp->s_out_buf.u4_bytes, mFpOutput);
        ASSERT_NE(numOutputBytes, 0) << "Failed to write the output!" << mOutFileName;
    }
}

void AvcEncTest::setRawBuf(iv_raw_buf_t* psInpRawBuf, const uint8_t* data) {
    switch (mIvVideoColorFormat) {
        case IV_YUV_420SP_UV:
            [[fallthrough]];
        case IV_YUV_420SP_VU: {
            uint8_t* yPlane = const_cast<uint8_t*>(data);
            uint8_t* uPlane = const_cast<uint8_t*>(data + (mFrameWidth * mFrameHeight));
            int32_t yStride = mFrameWidth;
            int32_t uStride = mFrameWidth / 2;
            psInpRawBuf->apv_bufs[0] = yPlane;
            psInpRawBuf->apv_bufs[1] = uPlane;

            psInpRawBuf->au4_wd[0] = mFrameWidth;
            psInpRawBuf->au4_wd[1] = mFrameWidth;

            psInpRawBuf->au4_ht[0] = mFrameHeight;
            psInpRawBuf->au4_ht[1] = mFrameHeight / 2;

            psInpRawBuf->au4_strd[0] = yStride;
            psInpRawBuf->au4_strd[1] = uStride;
            break;
        }
        case IV_YUV_422ILE: {
            uint8_t* yPlane = const_cast<uint8_t*>(data);
            psInpRawBuf->apv_bufs[0] = yPlane;

            psInpRawBuf->au4_wd[0] = mFrameWidth * 2;

            psInpRawBuf->au4_ht[0] = mFrameHeight;

            psInpRawBuf->au4_strd[0] = mFrameWidth * 2;
            break;
        }
        case IV_YUV_420P:
            [[fallthrough]];
        default: {
            uint8_t* yPlane = const_cast<uint8_t*>(data);
            uint8_t* uPlane = const_cast<uint8_t*>(data + (mFrameWidth * mFrameHeight));
            uint8_t* vPlane = const_cast<uint8_t*>(data + ((mFrameWidth * mFrameHeight) * 5) / 4);
            int32_t yStride = mFrameWidth;
            int32_t uStride = mFrameWidth / 2;
            int32_t vStride = mFrameWidth / 2;

            psInpRawBuf->apv_bufs[0] = yPlane;
            psInpRawBuf->apv_bufs[1] = uPlane;
            psInpRawBuf->apv_bufs[2] = vPlane;

            psInpRawBuf->au4_wd[0] = mFrameWidth;
            psInpRawBuf->au4_wd[1] = mFrameWidth / 2;
            psInpRawBuf->au4_wd[2] = mFrameWidth / 2;

            psInpRawBuf->au4_ht[0] = mFrameHeight;
            psInpRawBuf->au4_ht[1] = mFrameHeight / 2;
            psInpRawBuf->au4_ht[2] = mFrameHeight / 2;

            psInpRawBuf->au4_strd[0] = yStride;
            psInpRawBuf->au4_strd[1] = uStride;
            psInpRawBuf->au4_strd[2] = vStride;
            break;
        }
    }
    return;
}

int64_t AvcEncTest::getTotalFrames() {
    struct stat buf;
    stat(mFileName.c_str(), &buf);
    size_t fileSize = buf.st_size;
    int64_t totalFrames = (int64_t)(fileSize / mFrameSize);
    return totalFrames;
}

TEST_P(AvcEncTest, EncodeTest) {
    ASSERT_NO_FATAL_FAILURE(encodeFrames(mTotalFrames)) << "Failed to Encode: " << mFileName;
}

INSTANTIATE_TEST_SUITE_P(EncodeTest, AvcEncTest,
                         ::testing::Values(make_tuple("bbb_352x288_420p_30fps_32frames.yuv", 352,
                                                      288, 30, 2048),
                                           make_tuple("football_qvga.yuv", 320, 240, 30, 1024)));

int32_t main(int argc, char** argv) {
    gArgs = new TestArgs();
    ::testing::AddGlobalTestEnvironment(gArgs);
    ::testing::InitGoogleTest(&argc, argv);
    uint8_t status = gArgs->initFromOptions(argc, argv);
    if (status == 0) {
        status = RUN_ALL_TESTS();
    }
    return status;
}
