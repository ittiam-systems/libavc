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
*  isvce_pred_structs.h
*
* @brief
*  Contains struct definition used for prediction
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_PRED_STRUCTS_H_
#define _ISVCE_PRED_STRUCTS_H_

#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_defs.h"

/**
 * PU information
 */
typedef struct
{
    /**
     *  Motion Vector
     */
    mv_t s_mv;

    /**
     *  Ref index
     */
    WORD8 i1_ref_idx;

} isvce_enc_pu_mv_t;

/*
 * Total Pu info for an MB
 */
typedef struct isvce_enc_pu_t
{
    /* Array with ME info for all lists */
    isvce_enc_pu_mv_t as_me_info[NUM_PRED_DIRS];

    UWORD8 au1_mvp_idx[NUM_PRED_DIRS];

    /**
     *  PU X position in terms of min PU (4x4) units
     */
    UWORD8 u1_pos_x_in_4x4;

    /**
     *  PU Y position in terms of min PU (4x4) units
     */
    UWORD8 u1_pos_y_in_4x4;

    /**
     *  PU width in pixels = (u1_wd_in_4x4_m1 + 1) << 2
     */
    UWORD8 u1_wd_in_4x4_m1;

    /**
     *  PU height in pixels = (u1_ht_in_4x4_m1 + 1) << 2
     */
    UWORD8 u1_ht_in_4x4_m1;

    /**
     *  PRED_L0, PRED_L1, PRED_BI
     */
    UWORD8 u1_pred_mode;

} isvce_enc_pu_t;

typedef struct intra4x4_mode_data_t
{
    UWORD8 u1_predicted_mode;

    UWORD8 u1_mode;

} intra4x4_mode_data_t;

typedef intra4x4_mode_data_t intra8x8_mode_data_t;

typedef struct intra16x16_mode_data_t
{
    UWORD8 u1_mode;

} intra16x16_mode_data_t;

typedef struct enc_intra_pu_t
{
    intra4x4_mode_data_t as_i4x4_mode_data[MAX_TU_IN_MB];

    intra8x8_mode_data_t as_i8x8_mode_data[MIN_TU_IN_MB];

    intra16x16_mode_data_t s_i16x16_mode_data;

    UWORD8 u1_chroma_intra_mode;

} enc_intra_pu_t;

typedef struct isvce_mb_info_t
{
    isvce_enc_pu_t as_pu[ENC_MAX_PU_IN_MB];

    enc_intra_pu_t s_intra_pu;

    UWORD32 u4_cbp;

    UWORD32 u4_csbp;

    UWORD32 u4_res_csbp;

    UWORD16 u2_mb_type;

    WORD32 i4_mb_distortion;

    UWORD8 u1_base_mode_flag;

    UWORD8 u1_residual_prediction_flag;

    UWORD8 u1_tx_size;

    UWORD8 u1_mb_qp;

    UWORD8 u1_is_intra;

} isvce_mb_info_t;

#endif
