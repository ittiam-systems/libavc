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
*  isvce_intra_pred.h
*
* @brief
*  Contains function declarations for function declared in
*isvce_intra_pred.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/
#ifndef _ISVCE_IBL_EVAL_H_
#define _ISVCE_IBL_EVAL_H_

#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvc_intra_resample.h"
#include "isvce_structs.h"
#include "isvce_structs.h"

#define TEMP_BUF_SIZE_LUMA (REF_ARRAY_WIDTH * REF_ARRAY_WIDTH)
#define TEMP_BUF_SIZE_CB (REF_ARRAY_WIDTH * REF_ARRAY_WIDTH)
#define TEMP_BUF_SIZE_CR (DYADIC_REF_W_C * DYADIC_REF_H_C)

#define INTERMEDIATE_BUFF_WIDTH 48
#define INTERMEDIATE_BUFF_HEIGHT (MB_SIZE + 4)
#define TEMP_INTERPOLATION_BUF_SIZE (INTERMEDIATE_BUFF_WIDTH * INTERMEDIATE_BUFF_HEIGHT)

/* Structs */
typedef struct intra_pred_constants_t
{
    void *pv_state;
} intra_pred_constants_t;

typedef struct intra_pred_outputs_t
{
    yuv_buf_props_t s_pred_buf;
} intra_pred_outputs_t;

typedef struct intra_pred_variables_t
{
    svc_ilp_data_t *ps_svc_ilp_data;

    coordinates_t s_mb_pos;

    UWORD8 u1_spatial_layer_id;
} intra_pred_variables_t;

typedef struct svc_intra_pred_ctxt_t
{
    intra_pred_constants_t s_intra_pred_constants;

    intra_pred_variables_t s_intra_pred_variables;

    intra_pred_outputs_t s_intra_pred_outputs;

} svc_intra_pred_ctxt_t;

extern UWORD32 isvce_get_svc_intra_pred_ctxt_size(UWORD8 u1_num_spatial_layers,
                                                  DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                                  UWORD32 u4_ht);

extern void isvce_intra_pred_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_update_ibl_info(svc_intra_pred_ctxt_t *ps_intra_pred_ctxt,
                                  UWORD8 u1_num_spatial_layers, UWORD8 u1_spatial_layer_id,
                                  UWORD16 u2_mb_type, WORD32 i4_mb_x, WORD32 i4_mb_y,
                                  WORD8 u1_base_mode_flag);

extern void isvce_evaluate_IBL_mode(isvce_process_ctxt_t *ps_proc);

extern void isvce_pad_mb_mode_buf(svc_intra_pred_ctxt_t *ps_intra_pred_ctxt,
                                  UWORD8 u1_spatial_layer_id, UWORD8 u1_num_spatial_layers,
                                  DOUBLE d_spatial_res_ratio, UWORD32 u4_wd, UWORD32 u4_ht);

#endif
