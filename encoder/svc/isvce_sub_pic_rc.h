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
*  isvce_sub_pic_rc.h
*
* @brief
*  Contains typdefs and externs used for invoking sub-pic RC
*
*******************************************************************************
*/

#ifndef _ISVCE_SUB_PIC_RC_H_
#define _ISVCE_SUB_PIC_RC_H_

#include "ih264_typedefs.h"
#include "isvce_pred_structs.h"
#include "isvce_defs.h"

/* Structs */
typedef struct svc_sub_pic_rc_constants_t
{
    void *pv_state;

} svc_sub_pic_rc_constants_t;

typedef struct mb_bits_info_t
{
    WORD64 i8_header_bits;

    WORD64 i8_texture_bits;
} mb_bits_info_t;

typedef struct svc_sub_pic_rc_entropy_variables_t
{
    coordinates_t s_mb_pos;

    mb_bits_info_t s_mb_bits;

    UWORD8 u1_spatial_layer_id;
} svc_sub_pic_rc_entropy_variables_t;

typedef struct svc_sub_pic_rc_layer_variables_t
{
    WORD32 i4_max_num_reference_frames;

    WORD32 i4_slice_type;

    WORD32 i4_frame_num;

    UWORD8 u1_frame_qp;

    UWORD8 u1_min_qp;

    UWORD8 u1_max_qp;

    UWORD8 u1_spatial_layer_id;
} svc_sub_pic_rc_layer_variables_t;

typedef struct svc_sub_pic_rc_mb_variables_t
{
    buffer_container_t as_quant_coeffs[NUM_SP_COMPONENTS];

    isvce_enc_pu_mv_t *aps_mvps[MAX_MVP_IDX + 1];

    coordinates_t s_mb_pos;

    isvce_mb_info_t *ps_mb_info;

    UWORD8 *apu1_nnzs[NUM_SP_COMPONENTS];

    UWORD32 u4_cbp;
} svc_sub_pic_rc_mb_variables_t;

typedef struct svc_sub_pic_rc_variables_t
{
    svc_sub_pic_rc_layer_variables_t s_layer_variables;

    svc_sub_pic_rc_mb_variables_t s_mb_variables;

} svc_sub_pic_rc_variables_t;

typedef struct svc_sub_pic_rc_ctxt_t
{
    svc_sub_pic_rc_constants_t s_sub_pic_rc_constants;

    svc_sub_pic_rc_variables_t s_sub_pic_rc_variables;

    svc_sub_pic_rc_entropy_variables_t s_sub_pic_rc_entropy_variables;
} svc_sub_pic_rc_ctxt_t;

/* Function declarations */
extern UWORD32 isvce_get_sub_pic_rc_ctxt_size(UWORD8 u1_num_spatial_layers,
                                              DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                              UWORD32 u4_ht);

extern void isvce_sub_pic_rc_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_sub_pic_rc_ctxt_layer_init(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt);

extern void isvce_sub_pic_rc_ctxt_delete(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt);

extern void isvce_sub_pic_rc_ctxt_update(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt);

extern UWORD8 isvce_sub_pic_rc_get_mb_qp(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt,
                                         UWORD8 u1_cur_mb_qp);

extern void isvce_sub_pic_rc_get_entropy_data(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt);

extern void isvce_sub_pic_rc_dump_data(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt);

#endif
