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
*  isvce_ilp_mv.h
*
* @brief
*  Contains function declarations for function declared in
*  isvce_ilp_mv.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_ILP_MV_H_
#define _ISVCE_ILP_MV_H_

#include "ih264_typedefs.h"
#include "iv2.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_defs.h"
#include "isvce_pred_structs.h"
#include "isvce_structs.h"
#include "isvce_structs.h"
#include "isvce_utils.h"

/* Structs */
typedef struct ilp_mv_constants_t
{
    void *pv_state;
} ilp_mv_constants_t;

typedef struct ilp_mv_outputs_t
{
    ilp_mv_t s_ilp_mv;

    ilp_me_cands_t s_ilp_me_cands;

} ilp_mv_outputs_t;

typedef struct ilp_mv_variables_t
{
    svc_ilp_data_t *ps_svc_ilp_data;

    coordinates_t s_mb_pos;

    UWORD8 u1_spatial_layer_id;
} ilp_mv_variables_t;

typedef struct svc_ilp_mv_ctxt_t
{
    ilp_mv_constants_t s_ilp_mv_constants;

    ilp_mv_variables_t s_ilp_mv_variables;

    ilp_mv_outputs_t s_ilp_mv_outputs;

} svc_ilp_mv_ctxt_t;

/* Function declarations */
extern UWORD32 isvce_get_ilp_mv_ctxt_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                          UWORD32 u4_wd, UWORD32 u4_ht);

extern void isvce_ilp_mv_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_get_mb_ilp_mv(svc_ilp_mv_ctxt_t *ps_ilp_mv_ctxt);

extern void isvce_mvp_idx_eval(isvce_mb_info_t *ps_mb_info, isvce_enc_pu_mv_t *ps_spatial_mvp,
                               isvce_enc_pu_mv_t *ps_ilp_mvp, UWORD8 *pu1_mvd_costs);

static FORCEINLINE UWORD8 isvce_is_ilp_mv_winning_mv(isvce_mb_info_t *ps_mb_info,
                                                     ilp_mv_t *ps_ilp_mv)
{
    if(ENABLE_ILP_MV && ps_ilp_mv && (ps_mb_info->u2_mb_type != PSKIP) &&
       (ps_mb_info->u2_mb_type != BSKIP))
    {
        if((ps_mb_info->u2_mb_type == ps_ilp_mv->e_mb_type) &&
           (((PRED_MODE_T) ps_mb_info->as_pu->u1_pred_mode) == ps_ilp_mv->ae_pred_mode[0]))
        {
            return isvce_check_identical_mv(ps_mb_info->as_pu->as_me_info, ps_ilp_mv->as_mv[0],
                                            ps_ilp_mv->ae_pred_mode[0]);
        }
    }

    return 0;
}

#endif
