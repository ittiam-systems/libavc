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
 *  isvce_cabac_structs.h
 *
 * @brief
 *  This file contains cabac related structure definitions.
 *
 * @author
 *  Doney Alex
 *
 * @remarks
 *  none
 *
 *******************************************************************************
 */

#ifndef _ISVCE_CABAC_STRUCTS_H_
#define _ISVCE_CABAC_STRUCTS_H_

#include "ih264_typedefs.h"
#include "isvc_cabac_tables.h"
#include "ih264e_bitstream.h"
#include "ih264e_cabac_structs.h"

/**
 ******************************************************************************
 *  @brief      MB info for cabac
 ******************************************************************************
 */
typedef struct isvce_mb_info_ctxt_t
{
    /* Neighbour availability Variables needed to get CtxtInc, for CABAC */
    UWORD8 u1_mb_type; /* !< macroblock type: I/P/B/SI/SP */

    UWORD8 u1_cbp; /* !< Coded Block Pattern */
    UWORD8 u1_intrapred_chroma_mode;

    /*************************************************************************/
    /*               Arrangnment of AC CSBP                                  */
    /*        bits:  b7 b6 b5 b4 b3 b2 b1 b0                                 */
    /*        CSBP:  V1 V0 U1 U0 Y3 Y2 Y1 Y0                                 */
    /*************************************************************************/
    UWORD8 u1_yuv_ac_csbp;
    /*************************************************************************/
    /*               Arrangnment of DC CSBP                                  */
    /*        bits:  b7  b6  b5  b4  b3  b2  b1  b0                          */
    /*        CSBP:   x   x   x   x   x  Vdc Udc Ydc                         */
    /*************************************************************************/
    UWORD8 u1_yuv_dc_csbp;

    WORD8 i1_ref_idx[4];
    UWORD8 u1_mv[4][4];

    UWORD8 u1_base_mode_flag;
} isvce_mb_info_ctxt_t;

/**
 ******************************************************************************
 *  @brief      CABAC Context structure : Variables to handle Cabac
 ******************************************************************************
 */
typedef struct isvce_cabac_ctxt_t
{
    /*  Base pointer to all the cabac contexts  */
    bin_ctxt_model au1_cabac_ctxt_table[NUM_SVC_CABAC_CTXTS];

    cab_csbp_t s_lft_csbp;

    /**
     * pointer to Bitstream structure
     */
    bitstrm_t *ps_bitstrm;

    /* Pointer to mb_info_ctxt_t map_base */
    isvce_mb_info_ctxt_t *ps_mb_map_ctxt_inc_base;

    /* Pointer to encoding_envirnoment_t */
    encoding_envirnoment_t s_cab_enc_env;

    /* These things need to be updated at each MbLevel */

    /* Prev ps_mb_qp_delta_ctxt */
    WORD8 i1_prevps_mb_qp_delta_ctxt;

    /* Pointer to mb_info_ctxt_t map */
    isvce_mb_info_ctxt_t *ps_mb_map_ctxt_inc;

    /* Pointer to default mb_info_ctxt_t */
    isvce_mb_info_ctxt_t *ps_def_ctxt_mb_info;

    /* Pointer to current mb_info_ctxt_t */
    isvce_mb_info_ctxt_t *ps_curr_ctxt_mb_info;

    /* Pointer to left mb_info_ctxt_t */
    isvce_mb_info_ctxt_t *ps_left_ctxt_mb_info;

    /* Pointer to top mb_info_ctxt_t  */
    isvce_mb_info_ctxt_t *ps_top_ctxt_mb_info;

    /* Poniter to left csbp structure */
    cab_csbp_t *ps_lft_csbp;
    UWORD8 *pu1_left_y_ac_csbp;
    UWORD8 *pu1_left_uv_ac_csbp;
    UWORD8 *pu1_left_yuv_dc_csbp;

    /***************************************************************************/
    /*       Ref_idx contexts  are stored in the following way                 */
    /*  Array Idx 0,1 for reference indices in Forward direction               */
    /*  Array Idx 2,3 for reference indices in backward direction              */
    /***************************************************************************/
    /* Dimensions for u1_left_ref_ctxt_inc_arr is [2][4] for Mbaff:Top and Bot */
    WORD8 i1_left_ref_idx_ctx_inc_arr[2][4];
    WORD8 *pi1_left_ref_idx_ctxt_inc;

    /* Dimensions for u1_left_mv_ctxt_inc_arr is [2][4][4] for Mbaff case */
    UWORD8 u1_left_mv_ctxt_inc_arr[2][4][4];
    UWORD8 (*pu1_left_mv_ctxt_inc)[4];

} isvce_cabac_ctxt_t;

#endif
