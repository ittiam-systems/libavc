/******************************************************************************
 *
 * Copyright (C) 2022 The Android Open Source Project
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
 *****************************************************************************
 * Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
 */
/**
*******************************************************************************
* @file
*  isvce_defs.h
*
* @brief
*  Definitions used in the encoder
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_DEFS_H_
#define _ISVCE_DEFS_H_

#include "ih264e_defs.h"

#define SVC_MAX_NUM_BFRAMES 0

#define DEFAULT_INIT_QP 1

#define SVC_MAX_NUM_INP_FRAMES ((SVC_MAX_NUM_BFRAMES) + 2)

#define LOG2_MAX_FRAME_NUM_MINUS4 12

#define ENC_MAX_PU_IN_MB ((MB_SIZE / ENC_MIN_PU_SIZE) * (MB_SIZE / ENC_MIN_PU_SIZE))

#define MAX_REF_FRAMES_PER_PRED_DIR 1

#define SVC_MAX_SLICE_HDR_CNT 1

#define MAX_LAYER_REFERENCE_PICS 1

#define ENABLE_RESIDUAL_PREDICTION 1

#define ENABLE_ILP_MV 1

#define USE_ILP_MV_IN_ME (1 && (ENABLE_ILP_MV))

#define USE_ILP_MV_AS_MVP (1 && (ENABLE_ILP_MV))

#define MAX_MVP_IDX (USE_ILP_MV_AS_MVP ? 1 : 0)

#define ENABLE_IBL_MODE 1

#define ENABLE_INTRA_BASE_DEBLOCK (0 && (ENABLE_IBL_MODE))

#define ENABLE_MODE_STAT_VISUALISER 0

#define FORCE_FAST_INTRA4X4 0

#define FORCE_DISTORTION_BASED_INTRA_4X4_GATING 1

#define ENABLE_INTRA16X16_BASED_INTRA4X4_GATING 0

#define ENABLE_ILP_BASED_INTRA4X4_GATING 0

#define DISABLE_POST_ENC_SKIP 1

#define ENABLE_RE_ENC_AS_SKIP 1

#define MAX_ILP_MV_IN_NBR_RGN 4

/* L, T, TL, TR, Zero, Skip, 'Temporal Skip', ILP */
#define MAX_FPEL_SEARCH_CANDIDATES (7 + MAX_PU_IN_MB + MAX_ILP_MV_IN_NBR_RGN)

#define NUM_SVCE_RC_MEMTABS 45

#define SVCE_MAX_INP_DIM 1920

#define SVCE_MAX_INP_FRAME_SIZE (1920 * 1088)

/**
 ***************************************************************************
 * Enum to hold various mem records being request
 ****************************************************************************
 */
typedef enum ISVCE_MEMREC_TYPES_T
{
    /**
     * Codec Object at API level
     */
    ISVCE_MEM_REC_IV_OBJ,

    /**
     * Codec context
     */
    ISVCE_MEM_REC_CODEC,

    /**
     * Cabac context
     */
    ISVCE_MEM_REC_CABAC,

    /**
     * Cabac context_mb_info
     */
    ISVCE_MEM_REC_CABAC_MB_INFO,

    /**
     * entropy context
     */
    ISVCE_MEM_REC_ENTROPY,

    /**
     * Buffer to hold coeff data
     */
    ISVCE_MEM_REC_MB_COEFF_DATA,

    /**
     * Buffer to hold coeff data
     */
    ISVCE_MEM_REC_MB_HEADER_DATA,

    /**
     * Motion vector bank
     */
    ISVCE_MEM_REC_MVBANK,

    /**
     * Motion vector bits
     */
    ISVCE_MEM_REC_MVBITS,

    /**
     * Holds mem records passed to the codec.
     */
    ISVCE_MEM_REC_BACKUP,

    /**
     * Holds SPS
     */
    ISVCE_MEM_REC_SPS,

    /**
     * Holds PPS
     */
    ISVCE_MEM_REC_PPS,

    /**
     * Holds SVC NALU Extension data
     */
    ISVCE_MEM_REC_SVC_NALU_EXT,

    /**
     * Holds subset SPS data
     */
    ISVCE_MEM_REC_SUBSET_SPS,

    /**
     * Holds Slice Headers
     */
    ISVCE_MEM_REC_SLICE_HDR,

    /**
     * Holds SVC Slice Headers
     */
    ISVCE_MEM_REC_SVC_SLICE_HDR,

    /**
     * Contains map indicating slice index per MB basis
     */
    ISVCE_MEM_REC_SLICE_MAP,

    /**
     * Holds thread handles
     */
    ISVCE_MEM_REC_THREAD_HANDLE,

    /**
     * Holds control call mutex
     */
    ISVCE_MEM_REC_CTL_MUTEX,

    /**
     * Holds entropy call mutex
     */
    ISVCE_MEM_REC_ENTROPY_MUTEX,

    /**
     * Holds memory for Process JOB Queue
     */
    ISVCE_MEM_REC_PROC_JOBQ,

    /**
     * Holds memory for Entropy JOB Queue
     */
    ISVCE_MEM_REC_ENTROPY_JOBQ,

    /**
     * Contains status map indicating processing status per MB basis
     */
    ISVCE_MEM_REC_PROC_MAP,

    /**
     * Contains status map indicating deblocking status per MB basis
     */
    ISVCE_MEM_REC_DBLK_MAP,

    /*
     * Contains AIR map and mask
     */
    ISVCE_MEM_REC_AIR_MAP,

    /**
     * Contains status map indicating ME status per MB basis
     */
    ISVCE_MEM_REC_ME_MAP,

    /**
     * Holds dpb manager context
     */
    ISVCE_MEM_REC_DPB_MGR,

    /**
     * Holds intermediate buffers needed during processing stage
     * Memory for process contexts is allocated in this memtab
     */
    ISVCE_MEM_REC_PROC_SCRATCH,

    /**
     * Holds buffers for vert_bs, horz_bs and QP (all frame level)
     */
    ISVCE_MEM_REC_QUANT_PARAM,

    /**
     * Holds top row syntax information
     */
    ISVCE_MEM_REC_TOP_ROW_SYN_INFO,

    /**
     * Holds buffers for vert_bs, horz_bs and QP (all frame level)
     */
    ISVCE_MEM_REC_BS_QP,

    /**
     * Holds input buffer manager context
     */
    ISVCE_MEM_REC_INP_PIC,

    /**
     * Holds output buffer manager context
     */
    ISVCE_MEM_REC_OUT,

    /**
     * Holds picture buffer manager context and array of pic_buf_ts
     * Also holds reference picture buffers in non-shared mode
     */
    ISVCE_MEM_REC_REF_PIC,

    /*
     * Mem record for color space conversion
     */
    ISVCE_MEM_REC_CSC,

    /**
     * NMB info struct
     */
    ISVCE_MEM_REC_MB_INFO_NMB,

    /**
     * SVC Spatial layer Inputs
     */
    ISVCE_MEM_SVC_SPAT_INP,

    /**
     * Downscaler memory records
     */
    ISVCE_MEM_DOWN_SCALER,

    /**
     * SVC ILP data
     */
    ISVCE_MEM_SVC_ILP_DATA,

    /**
     * SVC ILP MV Context
     */
    ISVCE_MEM_SVC_ILP_MV_CTXT,

    /**
     * SVC ResPred Context
     */
    ISVCE_MEM_SVC_RES_PRED_CTXT,

    /**
     * SVC inter-layer intra pred context
     */
    ISVCE_MEM_SVC_INTRA_PRED_CTXT,

    /**
     * RC Utils Context
     */
    ISVCE_MEM_SVC_RC_UTILS_CTXT,

    /**
     * SubPic RC Context
     */
    ISVCE_MEM_SVC_SUB_PIC_RC_CTXT,

#if ENABLE_MODE_STAT_VISUALISER
    ISVCE_MEM_MODE_STAT_VISUALISER_BUF,
#endif

    /**
     * Rate control of memory records.
     */
    ISVCE_MEM_REC_RC,

    /**
     * Place holder to compute number of memory records.
     */
    ISVCE_MEM_REC_CNT = ISVCE_MEM_REC_RC + NUM_SVCE_RC_MEMTABS,

    /*
     * Do not add anything below
     */
} ISVCE_MEMREC_TYPES_T;

#endif
