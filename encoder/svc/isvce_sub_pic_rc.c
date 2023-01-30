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
*  isvce_sub_pic_rc.c
*
* @brief
*  Contains functions used in sub-pic RC
*
*******************************************************************************
*/
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ih264_typedefs.h"
#include "ih264_cavlc_tables.h"
#include "ih264_platform_macros.h"
#include "ithread.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"
#include "isvce_defs.h"
#include "isvce_sub_pic_rc.h"
#include "isvce_sub_pic_rc_private_defs.h"

/* Dependencies of 'irc_picture_type.h' */
#include "irc_mem_req_and_acq.h"

/* Dependencies of 'irc_rate_control_api_structs' */
#include "irc_picture_type.h"
#include "irc_rd_model.h"
#include "irc_vbr_storage_vbv.h"
#include "irc_est_sad.h"
#include "irc_bit_allocation.h"
#include "irc_mb_model_based.h"
#include "irc_cbr_buffer_control.h"
#include "irc_vbr_str_prms.h"
#include "irc_common.h"

#include "irc_rate_control_api_structs.h"
#include "irc_rate_control_api.h"

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing subPicRC ctxt
*
* @returns  Size of buffers
*
*******************************************************************************
*/
UWORD32 isvce_get_sub_pic_rc_ctxt_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                       UWORD32 u4_wd, UWORD32 u4_ht)
{
    WORD32 i;

    UWORD32 u4_size = MAX_PROCESS_CTXT * sizeof(svc_sub_pic_rc_ctxt_t);

    u4_size += sizeof(sub_pic_rc_state_t);
    u4_size += ithread_get_mutex_struct_size();

    for(i = u1_num_spatial_layers - 1; i >= 0; i--)
    {
        WORD32 i4_layer_wd =
            (WORD32) ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
            0.99;
        WORD32 i4_layer_ht =
            ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
        WORD32 i4_layer_mbs = (i4_layer_wd / MB_SIZE) * (i4_layer_ht / MB_SIZE);

        /* ps_mb_bits_info */
        u4_size += i4_layer_mbs * sizeof(mb_bits_info_t);

#if DUMP_SUB_PIC_RC_DATA
        /* ps_mb_bits_actual */
        u4_size += i4_layer_mbs * sizeof(mb_bits_info_t);
#endif
    }

    return u4_size;
}

void isvce_sub_pic_rc_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    sub_pic_rc_state_t *ps_sub_pic_rc_state;

    WORD32 i, j;

    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size =
        isvce_get_sub_pic_rc_ctxt_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    for(i = 0; i < MAX_PROCESS_CTXT; i++)
    {
        svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt = ps_codec->as_process[i].ps_sub_pic_rc_ctxt =
            (svc_sub_pic_rc_ctxt_t *) pu1_buf;

        pu1_buf += sizeof(ps_sub_pic_rc_ctxt[0]);
        i8_alloc_mem_size -= sizeof(ps_sub_pic_rc_ctxt[0]);

        if(0 == i)
        {
            ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants.pv_state = ps_sub_pic_rc_state =
                (sub_pic_rc_state_t *) pu1_buf;
            pu1_buf += sizeof(ps_sub_pic_rc_state[0]);
            i8_alloc_mem_size -= sizeof(ps_sub_pic_rc_state[0]);

            ASSERT(i8_alloc_mem_size >= 0);
            ASSERT(NULL != ps_codec->s_rate_control.apps_rate_control_api);
            ASSERT(NULL != ps_codec->as_process->s_me_ctxt.pu1_mv_bits);

            ps_sub_pic_rc_state->s_svc_params = ps_codec->s_cfg.s_svc_params;
            ps_sub_pic_rc_state->pu1_uev_codeword_to_bits_map = gau1_uev_codeword_to_bits_map;
            ps_sub_pic_rc_state->pu1_sev_codeword_to_bits_map =
                ps_codec->as_process->s_me_ctxt.pu1_mv_bits;
            ps_sub_pic_rc_state->e_rc_mode = ps_codec->s_cfg.e_rc_mode;

            ps_sub_pic_rc_state->pv_bits_accumulator_mutex = (void *) pu1_buf;
            pu1_buf += ithread_get_mutex_struct_size();
            i8_alloc_mem_size -= ithread_get_mutex_struct_size();
            ithread_mutex_init(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

            for(j = u1_num_spatial_layers - 1; j >= 0; j--)
            {
                sub_pic_rc_layer_state_t *ps_layer_state =
                    &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[j];

                WORD32 i4_layer_wd =
                    (WORD32) ((DOUBLE) u4_wd /
                              pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                    0.99;
                WORD32 i4_layer_ht =
                    ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                    0.99;
                WORD32 i4_layer_mbs = (i4_layer_wd / MB_SIZE) * (i4_layer_ht / MB_SIZE);

                ps_layer_state->i4_wd = i4_layer_wd;
                ps_layer_state->i4_ht = i4_layer_ht;
                ps_layer_state->i4_num_mbs = i4_layer_mbs;
                ps_layer_state->pv_layer_rc_ctxt =
                    ps_codec->s_rate_control.apps_rate_control_api[j];
                ps_layer_state->ps_mb_bits_info = (mb_bits_info_t *) pu1_buf;
                pu1_buf += i4_layer_mbs * sizeof(ps_layer_state->ps_mb_bits_info[0]);
                i8_alloc_mem_size -= i4_layer_mbs * sizeof(ps_layer_state->ps_mb_bits_info[0]);

                ASSERT(i8_alloc_mem_size >= 0);

#if DUMP_SUB_PIC_RC_DATA
                ps_layer_state->ps_mb_bits_actual = (mb_bits_info_t *) pu1_buf;
                pu1_buf += i4_layer_mbs * sizeof(ps_layer_state->ps_mb_bits_actual[0]);
                i8_alloc_mem_size -= i4_layer_mbs * sizeof(ps_layer_state->ps_mb_bits_actual[0]);

                ASSERT(i8_alloc_mem_size >= 0);

                {
                    UWORD8 au1_file_path[MAX_SUB_PIC_RC_DUMP_FILE_PATH_LENGTH + 1];

                    sprintf((WORD8 *) au1_file_path, "%ssubPicRC%1d.txt", SUB_PIC_RC_DUMP_FILE_PATH,
                            j);

                    ps_layer_state->ps_data_dump_file = fopen(au1_file_path, "w");

                    ASSERT(NULL != ps_layer_state->ps_data_dump_file);
                }
#endif
            }
        }
        else
        {
            svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt_src =
                ps_codec->as_process[0].ps_sub_pic_rc_ctxt;
            svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt_dst =
                ps_codec->as_process[i].ps_sub_pic_rc_ctxt;
            sub_pic_rc_state_t *ps_proc0_state =
                (sub_pic_rc_state_t *) ps_sub_pic_rc_ctxt_src->s_sub_pic_rc_constants.pv_state;

            ps_sub_pic_rc_ctxt_dst->s_sub_pic_rc_constants.pv_state = ps_proc0_state;
        }
    }
}

static FORCEINLINE void isvce_sub_pic_rc_qp_params_init(sub_pic_rc_qp_params_t *ps_qp_params,
                                                        UWORD8 u1_min_qp, UWORD8 u1_max_qp)
{
    ps_qp_params->u1_min_qp = u1_min_qp;
    ps_qp_params->u1_max_qp = u1_max_qp;
    ps_qp_params->pu4_qp_to_qscale_map = gau4_qp_to_qscale_map;
    ps_qp_params->pu1_qscale_to_qp_map = gau1_qscale_to_qp_map;
}

void isvce_sub_pic_rc_ctxt_layer_init(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
    sub_pic_rc_layer_state_t *ps_layer_state;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->s_layer_variables.u1_spatial_layer_id;

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];

    memset(&ps_layer_state->s_cumulative_mb_bits, 0, sizeof(ps_layer_state->s_cumulative_mb_bits));
    ps_layer_state->u4_num_mbs_sampled = 0;

    /* Frames with frameNum=0 are usually IDR's. RC model will be reset for IDR's.
     */
    /* Hence, using VBVBufSize as a proxy for estimated bits */
    if(0 == ps_sub_pic_rc_variables->s_layer_variables.i4_frame_num)
    {
        ps_layer_state->u4_allocated_bits =
            irc_get_vbv_buf_size(ps_layer_state->pv_layer_rc_ctxt) / 10.;
    }
    else
    {
        ps_layer_state->u4_allocated_bits =
            irc_get_prev_frm_est_bits(ps_layer_state->pv_layer_rc_ctxt);
    }

    isvce_sub_pic_rc_qp_params_init(&ps_layer_state->s_qp_params,
                                    ps_sub_pic_rc_variables->s_layer_variables.u1_min_qp,
                                    ps_sub_pic_rc_variables->s_layer_variables.u1_max_qp);
}

static FORCEINLINE UWORD32 isvce_sub_pic_rc_get_res_pred_flag_bits(
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables, sub_pic_rc_state_t *ps_sub_pic_rc_state)
{
    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UNUSED(ps_sub_pic_rc_state);

    return (ENABLE_RESIDUAL_PREDICTION && !ps_mb_info->u1_is_intra);
}

static FORCEINLINE UWORD32 isvce_sub_pic_rc_get_cbp_bits(
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables, sub_pic_rc_state_t *ps_sub_pic_rc_state)
{
    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD32 u4_cbp = ps_sub_pic_rc_variables->s_mb_variables.u4_cbp;
    bool b_use_inter_cbp_map = !ps_mb_info->u1_is_intra || ps_mb_info->u1_base_mode_flag;

    return ps_sub_pic_rc_state
        ->pu1_uev_codeword_to_bits_map[gu1_cbp_map_tables[u4_cbp][b_use_inter_cbp_map]];
}

static FORCEINLINE UWORD32 isvce_sub_pic_rc_get_mb_type_bits(
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables, sub_pic_rc_state_t *ps_sub_pic_rc_state)
{
    UWORD32 u4_mb_type;

    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD32 u4_cbp = ps_sub_pic_rc_variables->s_mb_variables.u4_cbp;
    UWORD32 au4_cbps[NUM_SP_COMPONENTS] = {u4_cbp & 15, u4_cbp >> 4};

    switch(ps_mb_info->u2_mb_type)
    {
        case I16x16:
        {
            u4_mb_type = ps_mb_info->s_intra_pu.s_i16x16_mode_data.u1_mode + 1 +
                         (au4_cbps[UV] << 2) + (au4_cbps[Y] == 15) * 12;

            break;
        }
        case I4x4:
        {
            u4_mb_type = 5 * (ps_sub_pic_rc_variables->s_layer_variables.i4_slice_type != ISLICE);

            break;
        }
        case P16x16:
        {
            u4_mb_type = 0;

            break;
        }
        default:
        {
            return 0;
        }
    }

    return ps_sub_pic_rc_state->pu1_uev_codeword_to_bits_map[u4_mb_type];
}

static FORCEINLINE UWORD32 isvce_sub_pic_rc_get_mb_pred_bits(
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables, sub_pic_rc_state_t *ps_sub_pic_rc_state)
{
    WORD32 i;

    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD32 u4_bits = 0;

    switch(ps_mb_info->u2_mb_type)
    {
        case I16x16:
        {
            /* intra_chroma_pred_mode */
            u4_bits +=
                ps_sub_pic_rc_state
                    ->pu1_uev_codeword_to_bits_map[ps_mb_info->s_intra_pu.u1_chroma_intra_mode];

            break;
        }
        case I4x4:
        {
            intra4x4_mode_data_t *ps_i4x4_mode_data = ps_mb_info->s_intra_pu.as_i4x4_mode_data;

            for(i = 0; i < MAX_TU_IN_MB; i++)
            {
                /* prev_intra4x4_pred_mode_flag */
                u4_bits += 1;

                /* rem_intra4x4_pred_mode */
                u4_bits +=
                    3 * (ps_i4x4_mode_data[i].u1_mode != ps_i4x4_mode_data[i].u1_predicted_mode);
            }

            /* intra_chroma_pred_mode */
            u4_bits +=
                ps_sub_pic_rc_state
                    ->pu1_uev_codeword_to_bits_map[ps_mb_info->s_intra_pu.u1_chroma_intra_mode];

            break;
        }
        case P16x16:
        {
            mv_t s_mvd;

            /* motion_prediction_flag_l0 */
            u4_bits += USE_ILP_MV_AS_MVP;

            /* ref_idx_l0 */
            if(2 == ps_sub_pic_rc_variables->s_layer_variables.i4_max_num_reference_frames)
            {
                u4_bits += 1;
            }
            else if(2 < ps_sub_pic_rc_variables->s_layer_variables.i4_max_num_reference_frames)
            {
                u4_bits += ps_sub_pic_rc_state->pu1_uev_codeword_to_bits_map
                               [ps_mb_info->as_pu->as_me_info[L0].i1_ref_idx];
            }

            /* mvd_l0 */
            s_mvd.i2_mvx = ps_mb_info->as_pu->as_me_info[L0].s_mv.i2_mvx -
                           ps_sub_pic_rc_variables->s_mb_variables
                               .aps_mvps[ps_mb_info->as_pu->au1_mvp_idx[L0]]
                               ->s_mv.i2_mvx;
            s_mvd.i2_mvy = ps_mb_info->as_pu->as_me_info[L0].s_mv.i2_mvy -
                           ps_sub_pic_rc_variables->s_mb_variables
                               .aps_mvps[ps_mb_info->as_pu->au1_mvp_idx[L0]]
                               ->s_mv.i2_mvy;
            u4_bits += ps_sub_pic_rc_state->pu1_sev_codeword_to_bits_map[s_mvd.i2_mvx];
            u4_bits += ps_sub_pic_rc_state->pu1_sev_codeword_to_bits_map[s_mvd.i2_mvy];

            break;
        }
        default:
        {
            break;
        }
    }

    return u4_bits;
}

static void ihevce_svc_sub_pic_rc_set_header_bits(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
    sub_pic_rc_layer_state_t *ps_layer_state;
    mb_bits_info_t *ps_mb_bits_info;

    UWORD32 u4_mb_idx;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;
    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->s_layer_variables.u1_spatial_layer_id;

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];
    u4_mb_idx = ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_abscissa +
                ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_ordinate *
                    (ps_layer_state->i4_wd / MB_SIZE);
    ps_mb_bits_info = &ps_layer_state->ps_mb_bits_info[u4_mb_idx];

    /* Hypotheses used for header bits estimation - */
    /* 1. mb_skip_run, base_mode_flag, mb_type, mb_pred, residual_prediction_flag,
     * and cbp */
    /*    are considered as contibuting to header bits. */
    /* 2. mb_skip_run = 1 bit */
    /* 3. base_mode_flag = 1 bit */
    /* 4. mb_type = LUT mapping mbType to corresponding ue(v) */
    /* 5. mb_pred.I4x4 = 1 bit for 16 'prev_intra4x4_pred_mode_flag';  */
    /*                   3 bits for each explicitly signaled
     * 'rem_intra4x4_pred_mode' */
    /* 6. mb_pred.Inter = 1 bit for 'motion_prediction_flag_l0' and
     * 'motion_prediction_flag_l1', when necessary; */
    /*                    mvbits LUT for 'mvd_l0' and 'mvd_l1' */
    /* 7. mb_pred.intra_chroma_pred_mode = LUT mapping intra_chroma_pred_mode to
     * corresponding ue(v) */
    /* 8. residual_prediction_flag = 1 bit */
    /* 9. coded_block_pattern = LUT mapping mbType to corresponding me(v) */

    /* mb_skip_run is assumed to be either 0 or 1 */
    ps_mb_bits_info->i8_header_bits += 1;

    /* 'base_mode_flag' */
    if((ENABLE_ILP_MV || ENABLE_IBL_MODE) && u1_spatial_layer_id)
    {
        ps_mb_bits_info->i8_header_bits += 1;

        if(ps_mb_info->u1_base_mode_flag)
        {
            /* 'residual_prediction_flag' */
            ps_mb_bits_info->i8_header_bits += isvce_sub_pic_rc_get_res_pred_flag_bits(
                ps_sub_pic_rc_variables, ps_sub_pic_rc_state);

            /* 'coded_block_pattern' */
            ps_mb_bits_info->i8_header_bits +=
                isvce_sub_pic_rc_get_cbp_bits(ps_sub_pic_rc_variables, ps_sub_pic_rc_state);

            return;
        }
    }

    /* 'mb_type' */
    ps_mb_bits_info->i8_header_bits +=
        isvce_sub_pic_rc_get_mb_type_bits(ps_sub_pic_rc_variables, ps_sub_pic_rc_state);

    if(PSKIP == ps_mb_info->u2_mb_type)
    {
        return;
    }

    /* 'mb_pred' */
    ps_mb_bits_info->i8_header_bits +=
        isvce_sub_pic_rc_get_mb_pred_bits(ps_sub_pic_rc_variables, ps_sub_pic_rc_state);

    /* 'residual_prediction_flag' */
    ps_mb_bits_info->i8_header_bits +=
        isvce_sub_pic_rc_get_res_pred_flag_bits(ps_sub_pic_rc_variables, ps_sub_pic_rc_state);
}

static FORCEINLINE UWORD32 isvce_sub_pic_rc_get_tu_residual_bits(
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables, WORD32 i4_coeff_start_idx,
    UWORD8 u1_num_coded_coeffs, UWORD8 u1_num_coeffs, bool b_is_chroma)
{
    WORD32 i;
    UWORD32 u4_num_bits;

    UWORD32 u4_bits = 0;
    WORD16 *pi2_coeff =
        ((WORD16 *) ps_sub_pic_rc_variables->s_mb_variables.as_quant_coeffs[b_is_chroma ? UV : Y]
             .pv_data) +
        i4_coeff_start_idx;

    if(0 == u1_num_coded_coeffs)
    {
        return 0;
    }

    GETRANGE(u4_num_bits, u1_num_coded_coeffs);
    u4_bits += u4_num_bits;

    for(i = 0; i < u1_num_coeffs; i++)
    {
        if(pi2_coeff[i])
        {
            GETRANGE(u4_num_bits, pi2_coeff[i]);
            u4_bits += u4_num_bits;
        }
    }
    return u4_bits;
}

static void ihevce_svc_sub_pic_rc_set_texture_bits(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
    sub_pic_rc_layer_state_t *ps_layer_state;
    mb_bits_info_t *ps_mb_bits_info;

    UWORD32 u4_mb_idx;
    WORD32 i, j;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;
    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->s_layer_variables.u1_spatial_layer_id;
    UWORD32 au4_cbps[NUM_SP_COMPONENTS] = {ps_sub_pic_rc_variables->s_mb_variables.u4_cbp & 15,
                                           ps_sub_pic_rc_variables->s_mb_variables.u4_cbp >> 4};

    if(0 == ps_sub_pic_rc_variables->s_mb_variables.u4_cbp)
    {
        return;
    }

    if(MIN_TU_SIZE != ps_mb_info->u1_tx_size)
    {
        return;
    }

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];
    u4_mb_idx = ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_abscissa +
                ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_ordinate *
                    (ps_layer_state->i4_wd / MB_SIZE);
    ps_mb_bits_info = &ps_layer_state->ps_mb_bits_info[u4_mb_idx];

    /* Hypotheses used for texture bits estimation - */
    /* 1. Only level information is considered. */
    /* 2. nnz is used as a proxy for coeff_token. */
    /* 3. Both of the above are assumed coded via i(n). */
    if(au4_cbps[Y])
    {
        /* Y - DC */
        if(I16x16 == ps_mb_info->u2_mb_type)
        {
            ps_mb_bits_info->i8_texture_bits += isvce_sub_pic_rc_get_tu_residual_bits(
                ps_sub_pic_rc_variables, 0, ps_sub_pic_rc_variables->s_mb_variables.apu1_nnzs[Y][0],
                NUM_COEFFS_IN_MIN_TU, false);
        }

        for(i = 0; i < MIN_TU_IN_MB; i++)
        {
            if(au4_cbps[Y] & (1 << i))
            {
                UWORD32 u4_csbp = (ps_mb_info->u4_csbp >> (4 * i)) & 15;

                for(j = 0; j < NUM_4x4_IN_8x8; j++)
                {
                    if(u4_csbp & (1 << j))
                    {
                        /* 1 added to account for DC TU */
                        UWORD8 u1_blk_id = 1 + gau4_tu_zscan_id_to_rasterscan_id_map[i][j];
                        UWORD8 u1_nnz =
                            ps_sub_pic_rc_variables->s_mb_variables.apu1_nnzs[Y][u1_blk_id];

                        if(u1_nnz && (I16x16 == ps_mb_info->u2_mb_type))
                        {
                            u1_nnz -= !!(((WORD16 *) (ps_sub_pic_rc_variables->s_mb_variables
                                                          .as_quant_coeffs[Y]
                                                          .pv_data))[u1_blk_id - 1]);

                            ps_mb_bits_info->i8_texture_bits +=
                                isvce_sub_pic_rc_get_tu_residual_bits(
                                    ps_sub_pic_rc_variables,
                                    u1_blk_id * ps_sub_pic_rc_variables->s_mb_variables
                                                    .as_quant_coeffs[Y]
                                                    .i4_data_stride +
                                        (I16x16 == ps_mb_info->u2_mb_type),
                                    u1_nnz,
                                    NUM_COEFFS_IN_MIN_TU - (I16x16 == ps_mb_info->u2_mb_type),
                                    false);
                        }
                    }
                }
            }
        }
    }

    if(au4_cbps[UV])
    {
        for(i = ((WORD32) U); i <= ((WORD32) V); i++)
        {
            bool b_is_v = (i == ((WORD32) V));

            ps_mb_bits_info->i8_texture_bits += isvce_sub_pic_rc_get_tu_residual_bits(
                ps_sub_pic_rc_variables, b_is_v * NUM_4x4_IN_8x8,
                ps_sub_pic_rc_variables->s_mb_variables
                    .apu1_nnzs[UV][0 + b_is_v * (1 + NUM_4x4_IN_8x8)],
                NUM_4x4_IN_8x8, true);

            for(j = 0; j < NUM_4x4_IN_8x8; j++)
            {
                UWORD8 u1_nnz = ps_sub_pic_rc_variables->s_mb_variables
                                    .apu1_nnzs[UV][j + b_is_v * (1 + NUM_4x4_IN_8x8) + 1];

                if(u1_nnz)
                {
                    u1_nnz -=
                        !!(((WORD16 *) (ps_sub_pic_rc_variables->s_mb_variables.as_quant_coeffs[UV]
                                            .pv_data))[j + b_is_v * NUM_4x4_IN_8x8]);

                    ps_mb_bits_info->i8_texture_bits += isvce_sub_pic_rc_get_tu_residual_bits(
                        ps_sub_pic_rc_variables,
                        (j + b_is_v * NUM_4x4_IN_8x8 + 1) *
                                ps_sub_pic_rc_variables->s_mb_variables.as_quant_coeffs[UV]
                                    .i4_data_stride +
                            1,
                        u1_nnz, NUM_COEFFS_IN_MIN_TU - 1, true);
                }
            }
        }
    }
}

void isvce_sub_pic_rc_ctxt_update(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
    sub_pic_rc_layer_state_t *ps_layer_state;
    mb_bits_info_t *ps_mb_bits_info;

    UWORD32 u4_mb_idx;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;
    isvce_mb_info_t *ps_mb_info = ps_sub_pic_rc_variables->s_mb_variables.ps_mb_info;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->s_layer_variables.u1_spatial_layer_id;
    bool b_is_skip_mb = (PSKIP == ps_mb_info->u2_mb_type) || (BSKIP == ps_mb_info->u2_mb_type);

    if(!ENABLE_IN_FRAME_RC || (IVE_RC_NONE == ps_sub_pic_rc_state->e_rc_mode))
    {
        return;
    }

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];
    u4_mb_idx = ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_abscissa +
                ps_sub_pic_rc_variables->s_mb_variables.s_mb_pos.i4_ordinate *
                    (ps_layer_state->i4_wd / MB_SIZE);
    ps_mb_bits_info = &ps_layer_state->ps_mb_bits_info[u4_mb_idx];

    memset(ps_mb_bits_info, 0, sizeof(ps_mb_bits_info[0]));

    if(!b_is_skip_mb)
    {
        ihevce_svc_sub_pic_rc_set_header_bits(ps_sub_pic_rc_ctxt);

        ihevce_svc_sub_pic_rc_set_texture_bits(ps_sub_pic_rc_ctxt);
    }

    ithread_mutex_lock(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

    ps_layer_state->s_cumulative_mb_bits.i8_header_bits += ps_mb_bits_info->i8_header_bits;
    ps_layer_state->s_cumulative_mb_bits.i8_texture_bits += ps_mb_bits_info->i8_texture_bits;
    ps_layer_state->u4_num_mbs_sampled++;

    ithread_mutex_unlock(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);
}

UWORD8 isvce_sub_pic_rc_get_mb_qp(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt, UWORD8 u1_cur_mb_qp)
{
    sub_pic_rc_layer_state_t *ps_layer_state;

    DOUBLE d_bit_consumption_ratio;
    UWORD32 u4_frame_qscale;
    UWORD8 u1_mb_qp;
    UWORD32 u4_num_mbs_sampled;
    WORD32 i4_cumulative_mb_bits;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->s_layer_variables.u1_spatial_layer_id;
    UWORD8 u1_frame_qp = ps_sub_pic_rc_variables->s_layer_variables.u1_frame_qp;

    if(!ENABLE_IN_FRAME_RC || (IVE_RC_NONE == ps_sub_pic_rc_state->e_rc_mode))
    {
        return u1_cur_mb_qp;
    }

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];

    ithread_mutex_lock(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

    u4_num_mbs_sampled = ps_layer_state->u4_num_mbs_sampled;

    if(u4_num_mbs_sampled < (MIN_SAMPLED_MB_RATIO * ps_layer_state->i4_num_mbs))
    {
        ithread_mutex_unlock(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

        return u1_cur_mb_qp;
    }

    i4_cumulative_mb_bits = (WORD32) (ps_layer_state->s_cumulative_mb_bits.i8_header_bits +
                                      ps_layer_state->s_cumulative_mb_bits.i8_texture_bits);

    d_bit_consumption_ratio =
        (((DOUBLE) i4_cumulative_mb_bits) * ((DOUBLE) ps_layer_state->i4_num_mbs)) /
        (((DOUBLE) ps_layer_state->u4_allocated_bits) * ((DOUBLE) u4_num_mbs_sampled));

    ithread_mutex_unlock(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

    if((d_bit_consumption_ratio > BIT_RATIO_FOR_OVERCONSUMPTION) ||
       (d_bit_consumption_ratio < BIT_RATIO_FOR_UNDERCONSUMPTION))
    {
        u4_frame_qscale = ps_layer_state->s_qp_params.pu4_qp_to_qscale_map[u1_frame_qp] *
                              d_bit_consumption_ratio +
                          0.5;
        u4_frame_qscale = CLIP3(ps_layer_state->s_qp_params.pu4_qp_to_qscale_map[0], MAX_SVC_QSCALE,
                                u4_frame_qscale);
        u1_mb_qp = ps_layer_state->s_qp_params.pu1_qscale_to_qp_map[u4_frame_qscale];
        u1_mb_qp = CLIP3(ps_layer_state->s_qp_params.u1_min_qp,
                         ps_layer_state->s_qp_params.u1_max_qp, u1_mb_qp);
        u1_mb_qp = CLIP3(MAX(MIN_H264_QP, ((WORD16) u1_cur_mb_qp) - MAX_MB_QP_DECREMENT),
                         MIN(MAX_H264_QP, ((WORD16) u1_cur_mb_qp) + MAX_MB_QP_INCREMENT),
                         ((WORD16) u1_mb_qp));
        /* This ensures mb_qp_delta stays within the interval [-26, 25] */
        u1_mb_qp = CLIP3(MAX(MIN_H264_QP, ((WORD16) u1_frame_qp) - MAX_FRAME_QP_DECREMENT),
                         MIN(MAX_H264_QP, ((WORD16) u1_frame_qp) + MAX_FRAME_QP_INCREMENT),
                         ((WORD16) u1_mb_qp));
    }
    else
    {
        u1_mb_qp = u1_cur_mb_qp;
    }

    {
        vbv_buf_status_e e_vbv_buf_status;
        picture_type_e e_rc_pic_type;

        DOUBLE d_est_frame_bits;

        WORD32 i4_num_bits_to_prevent_vbv_underflow;

        d_est_frame_bits = ((DOUBLE) i4_cumulative_mb_bits) * ((DOUBLE) ps_layer_state->i4_num_mbs);
        d_est_frame_bits /= u4_num_mbs_sampled;

        switch(ps_sub_pic_rc_variables->s_layer_variables.i4_slice_type)
        {
            case ISLICE:
            {
                e_rc_pic_type = I_PIC;
                break;
            }
            case PSLICE:
            {
                e_rc_pic_type = P_PIC;
                break;
            }
            default:
            {
                e_rc_pic_type = B_PIC;
                break;
            }
        }

        e_vbv_buf_status =
            irc_get_buffer_status(ps_layer_state->pv_layer_rc_ctxt, (WORD32) d_est_frame_bits,
                                  e_rc_pic_type, &i4_num_bits_to_prevent_vbv_underflow);

        /* This models dec VBV buffer */
        if(VBV_OVERFLOW == e_vbv_buf_status)
        {
            u1_mb_qp--;
        }
        else if(VBV_UNDERFLOW == e_vbv_buf_status)
        {
            u1_mb_qp++;
        }

        /* This ensures mb_qp_delta stays within the interval [-26, 25] */
        u1_mb_qp = CLIP3(ps_layer_state->s_qp_params.u1_min_qp,
                         ps_layer_state->s_qp_params.u1_max_qp, u1_mb_qp);
        u1_mb_qp = CLIP3(MAX(MIN_H264_QP, ((WORD16) u1_frame_qp) - MAX_FRAME_QP_DECREMENT),
                         MIN(MAX_H264_QP, ((WORD16) u1_frame_qp) + MAX_FRAME_QP_INCREMENT),
                         ((WORD16) u1_mb_qp));
    }

    return u1_mb_qp;
}

void isvce_sub_pic_rc_get_entropy_data(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
#if DUMP_SUB_PIC_RC_DATA
    sub_pic_rc_layer_state_t *ps_layer_state;

    UWORD32 u4_mb_idx;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    svc_sub_pic_rc_entropy_variables_t *ps_sub_pic_rc_variables =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_entropy_variables;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;

    UWORD8 u1_spatial_layer_id = ps_sub_pic_rc_variables->u1_spatial_layer_id;

    if(!ENABLE_IN_FRAME_RC || (IVE_RC_NONE == ps_sub_pic_rc_state->e_rc_mode))
    {
        return;
    }

    ps_layer_state = &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[u1_spatial_layer_id];
    u4_mb_idx = ps_sub_pic_rc_variables->s_mb_pos.i4_abscissa +
                ps_sub_pic_rc_variables->s_mb_pos.i4_ordinate * (ps_layer_state->i4_wd / MB_SIZE);

    ps_layer_state->ps_mb_bits_actual[u4_mb_idx] = ps_sub_pic_rc_variables->s_mb_bits;
#else
    UNUSED(ps_sub_pic_rc_ctxt);
#endif
}

void isvce_sub_pic_rc_dump_data(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
#if DUMP_SUB_PIC_RC_DATA
    WORD32 i, j, k;

    svc_sub_pic_rc_constants_t *ps_sub_pic_rc_constants =
        &ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants;
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_constants->pv_state;

    if(!ENABLE_IN_FRAME_RC || (IVE_RC_NONE == ps_sub_pic_rc_state->e_rc_mode))
    {
        return;
    }

    for(i = 0; i < ps_sub_pic_rc_state->s_svc_params.u1_num_spatial_layers; i++)
    {
        sub_pic_rc_layer_state_t *ps_layer_state =
            &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[i];

        for(j = 0; j < (ps_layer_state->i4_ht / MB_SIZE); j++)
        {
            for(k = 0; k < (ps_layer_state->i4_wd / MB_SIZE); k++)
            {
                mb_bits_info_t *ps_mb_bits_est =
                    &ps_layer_state->ps_mb_bits_info[k + j * (ps_layer_state->i4_wd / MB_SIZE)];
                mb_bits_info_t *ps_mb_bits_actual =
                    &ps_layer_state->ps_mb_bits_actual[k + j * (ps_layer_state->i4_wd / MB_SIZE)];

                fprintf(ps_layer_state->ps_data_dump_file, "%ld,%ld,%ld,%ld,\n",
                        ps_mb_bits_est->i8_header_bits, ps_mb_bits_est->i8_texture_bits,
                        ps_mb_bits_actual->i8_header_bits, ps_mb_bits_actual->i8_texture_bits);
            }
        }
    }
#else
    UNUSED(ps_sub_pic_rc_ctxt);
#endif
}

void isvce_sub_pic_rc_ctxt_delete(svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt)
{
    sub_pic_rc_state_t *ps_sub_pic_rc_state =
        (sub_pic_rc_state_t *) ps_sub_pic_rc_ctxt->s_sub_pic_rc_constants.pv_state;

    ithread_mutex_destroy(ps_sub_pic_rc_state->pv_bits_accumulator_mutex);

#if DUMP_SUB_PIC_RC_DATA
    {
        WORD32 i;

        UWORD8 u1_num_spatial_layers = ps_sub_pic_rc_state->s_svc_params.u1_num_spatial_layers;

        for(i = u1_num_spatial_layers - 1; i >= 0; i--)
        {
            sub_pic_rc_layer_state_t *ps_layer_state =
                &ps_sub_pic_rc_state->as_sub_pic_rc_layer_states[i];

            if(ps_layer_state->ps_data_dump_file)
            {
                fclose(ps_layer_state->ps_data_dump_file);
            }

            ps_layer_state->ps_data_dump_file = NULL;
        }
    }
#endif
}
