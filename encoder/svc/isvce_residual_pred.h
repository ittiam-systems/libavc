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
*  isvce_residual_pred.h
*
* @brief
*  Contains function declarations for function declared in
*isvce_residual_pred.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_RESIDUAL_PRED_H_
#define _ISVCE_RESIDUAL_PRED_H_

#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"
#include "isvce_structs.h"

/* Structs */
typedef struct res_pred_constants_t
{
    void *pv_state;
} res_pred_constants_t;

typedef struct res_pred_outputs_t
{
    yuv_buf_props_t s_res_pred;
} res_pred_outputs_t;

typedef struct res_pred_variables_t
{
    svc_ilp_data_t *ps_svc_ilp_data;

    coordinates_t s_mb_pos;

    UWORD8 u1_spatial_layer_id;
} res_pred_variables_t;

typedef struct svc_res_pred_ctxt_t
{
    res_pred_constants_t s_res_pred_constants;

    res_pred_variables_t s_res_pred_variables;

    res_pred_outputs_t s_res_pred_outputs;

} svc_res_pred_ctxt_t;

extern UWORD32 isvce_get_svc_res_pred_ctxt_size(UWORD8 u1_num_spatial_layers,
                                                DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                                UWORD32 u4_ht);

extern void isvce_svc_res_pred_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_get_mb_residual_pred(svc_res_pred_ctxt_t *ps_res_pred_ctxt);

extern void isvce_get_mb_residual_pred_non_dyadic(svc_res_pred_ctxt_t *ps_res_pred_ctxt);

extern void isvce_residual_pred_eval(svc_res_pred_ctxt_t *ps_res_pred_ctxt, yuv_buf_props_t *ps_src,
                                     yuv_buf_props_t *ps_pred, yuv_buf_props_t *ps_res,
                                     UWORD32 *pu4_res_pred_sad,
                                     UWORD8 *pu1_residual_prediction_flag, UWORD32 u4_winning_sad);

extern void isvce_update_res_pred_info(isvce_process_ctxt_t *ps_proc);

#endif
