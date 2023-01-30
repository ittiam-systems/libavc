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
*  isvce_cabac_utils.h
*
* @brief
*  Contains function declarations for function declared in
*  isvce_svc_cabac_utils.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_CABAC_UTILS_H_
#define _ISVCE_CABAC_UTILS_H_

#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "isvc_defs.h"
#include "isvc_cabac_tables.h"
#include "isvce_cabac_structs.h"
#include "isvce_cabac.h"

static FORCEINLINE void isvce_cabac_enc_base_mode_flag(isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                                       UWORD8 u1_base_mode_flag)
{
    UWORD8 u1_ctx_inc;
    UWORD8 u1_a, u1_b;

    const UWORD32 u4_ctxidx_offset = BASE_MODE_FLAG;

    u1_a = !ps_cabac_ctxt->ps_left_ctxt_mb_info->u1_base_mode_flag;
    u1_b = !ps_cabac_ctxt->ps_top_ctxt_mb_info->u1_base_mode_flag;

    u1_ctx_inc = u1_a + u1_b;

    isvce_cabac_encode_bin(ps_cabac_ctxt, u1_base_mode_flag,
                           ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctxidx_offset + u1_ctx_inc);
}

static FORCEINLINE void isvce_cabac_enc_residual_prediction_flag(isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                                                 UWORD8 u1_base_mode_flag,
                                                                 UWORD8 u1_residual_prediction_flag)
{
    const UWORD32 u4_ctxidx_offset = RESIDUAL_PREDICTION_FLAG;
    UWORD8 u1_ctx_inc = !u1_base_mode_flag;

    isvce_cabac_encode_bin(ps_cabac_ctxt, u1_residual_prediction_flag,
                           ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctxidx_offset + u1_ctx_inc);
}

static FORCEINLINE void isvce_cabac_enc_motion_prediction_flag(isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                                               UWORD8 u1_motion_prediction_flag,
                                                               UWORD8 u1_is_l0_mvp)
{
    const UWORD32 u4_ctxidx_offset =
        u1_is_l0_mvp ? MOTION_PREDICTION_FLAG_L0 : MOTION_PREDICTION_FLAG_L1;

    isvce_cabac_encode_bin(ps_cabac_ctxt, u1_motion_prediction_flag,
                           ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctxidx_offset);
}

#endif
