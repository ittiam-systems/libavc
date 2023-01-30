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
*  isvce_ilp_mv.c
*
* @brief
*  Contains functions used for deriving inter_layer MV's
*
*******************************************************************************
*/
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "isvc_macros.h"
#include "isvc_defs.h"
#include "isvce_defs.h"
#include "isvce_structs.h"
#include "isvce_ilp_mv_private_defs.h"
#include "isvce_ilp_mv.h"
#include "isvce_ilp_mv_utils.h"

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing ILP MV ctxt
*
* @param[in] u1_num_spatial_layers
*  Num Spatial Layers
*
* @param[in] d_spatial_res_ratio
*  Resolution Ratio b/w spatial layers
*
* @param[in] u4_wd
*  Input Width
*
* @param[in] u4_ht
*  Input Height
*
* @returns  Size of buffers
*
*******************************************************************************
*/
UWORD32 isvce_get_ilp_mv_ctxt_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                   UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 u4_size = 0;

    if(u1_num_spatial_layers > 1)
    {
        WORD32 i;

        u4_size += MAX_PROCESS_CTXT * sizeof(svc_ilp_mv_ctxt_t);
        u4_size += MAX_PROCESS_CTXT * sizeof(ilp_mv_state_t);

        u4_size += u1_num_spatial_layers * sizeof(ilp_mv_layer_state_t);

        for(i = u1_num_spatial_layers - 1; i >= 1; i--)
        {
            WORD32 i4_layer_luma_wd =
                (WORD32) ((DOUBLE) u4_wd /
                          pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_mbs = (i4_layer_luma_wd / MB_SIZE) * (i4_layer_luma_ht / MB_SIZE);

            u4_size += i4_layer_luma_mbs * sizeof(ilp_mv_mb_state_t);
        }
    }

    return u4_size;
}

static FORCEINLINE void isvce_ref_layer_pu_and_mb_pos_init(layer_resampler_props_t *ps_layer_props,
                                                           ilp_mv_mb_state_t *ps_mb_state,
                                                           coordinates_t *ps_mb_pos,
                                                           UWORD32 u4_ref_wd, UWORD32 u4_ref_ht,
                                                           UWORD8 u1_field_pic_flag,
                                                           UWORD8 u1_field_mb_flag)
{
    UWORD32 i, j;

    coordinates_t(*aps_pu_positions)[MAX_PU_IN_MB_ROW] = ps_mb_state->as_pu_positions;
    coordinates_t(*aps_mb_positions)[MAX_PU_IN_MB_ROW] = ps_mb_state->as_mb_positions;

    for(i = 0; i < MAX_PU_IN_MB_COL; i++)
    {
        UWORD32 u4_y_ref16;

        UWORD32 u4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht +
                        (4 * i + 1) * (1 + u1_field_mb_flag - u1_field_pic_flag);

        u4_y_ref16 =
            (u4_yc * ps_layer_props->u4_scale_y + (1 << (ps_layer_props->u4_shift_y - 1))) >>
            ps_layer_props->u4_shift_y;
        u4_y_ref16 = MIN(u4_y_ref16, u4_ref_ht - 1);

        for(j = 0; j < MAX_PU_IN_MB_ROW; j++)
        {
            UWORD32 u4_x_ref16;

            UWORD32 u4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + 4 * j + 1;

            u4_x_ref16 =
                (u4_xc * ps_layer_props->u4_scale_x + (1 << (ps_layer_props->u4_shift_x - 1))) >>
                ps_layer_props->u4_shift_x;
            u4_x_ref16 = MIN(u4_x_ref16, u4_ref_wd - 1);

            aps_pu_positions[i][j].i4_abscissa = u4_x_ref16;
            aps_pu_positions[i][j].i4_ordinate = u4_y_ref16;

            aps_mb_positions[i][j].i4_abscissa = (u4_x_ref16 / MB_SIZE);
            aps_mb_positions[i][j].i4_ordinate = (u4_y_ref16 / MB_SIZE);
        }
    }
}

static void isvce_ilp_mv_layer_state_init(ilp_mv_layer_state_t *ps_layer_state,
                                          DOUBLE d_spatial_res_ratio, UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 i, j;

    const UWORD8 u1_ref_layer_field_pic_flag = 0;
    const UWORD8 u1_field_pic_flag = 0;
    const UWORD8 u1_field_mb_flag = 0;

    ilp_mv_mb_state_t *ps_mb_states;
    layer_resampler_props_t *ps_layer_props;

    UWORD32 u4_wd_in_mbs;
    UWORD32 u4_ht_in_mbs;

    UWORD32 u4_ref_wd = (u4_wd / d_spatial_res_ratio);
    UWORD32 u4_ref_ht = (u4_ht / d_spatial_res_ratio) * (1 + u1_ref_layer_field_pic_flag);
    UWORD32 u4_scaled_wd = u4_wd;
    UWORD32 u4_scaled_ht = u4_ht * (1 + u1_field_pic_flag);

    ps_mb_states = ps_layer_state->ps_mb_states;
    ps_layer_props = ps_layer_state->ps_props;

    u4_wd_in_mbs = u4_scaled_wd / ps_layer_props->u4_mb_wd;
    u4_ht_in_mbs = u4_scaled_ht / ps_layer_props->u4_mb_ht;

    ps_layer_state->s_mv_scale.i4_abscissa = ((u4_scaled_wd << 16) + (u4_ref_wd >> 1)) / u4_ref_wd;
    ps_layer_state->s_mv_scale.i4_ordinate = ((u4_scaled_ht << 16) + (u4_ref_ht >> 1)) / u4_ref_ht;

    for(i = 0; i < u4_ht_in_mbs; i++)
    {
        for(j = 0; j < u4_wd_in_mbs; j++)
        {
            coordinates_t s_mb_pos = {j, i};

            isvce_ref_layer_pu_and_mb_pos_init(ps_layer_props, &ps_mb_states[j + i * u4_wd_in_mbs],
                                               &s_mb_pos, u4_ref_wd, u4_ref_ht, u1_field_pic_flag,
                                               u1_field_mb_flag);
        }
    }
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize svc ilp buffers
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_mem_rec
*  Pointer to memory allocated for input buffers
*
*******************************************************************************
*/
void isvce_ilp_mv_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    WORD32 i, j;

    const WORD32 i4_num_proc_ctxts = sizeof(ps_codec->as_process) / sizeof(ps_codec->as_process[0]);
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;

    if(u1_num_spatial_layers > 1)
    {
        ilp_mv_layer_state_t *ps_layer_states;
        ilp_mv_mb_state_t *aps_luma_mb_states[MAX_NUM_SPATIAL_LAYERS];

        DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
        UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
        UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
        UWORD8 *pu1_buf = ps_mem_rec->pv_base;
        WORD64 i8_alloc_mem_size =
            isvce_get_ilp_mv_ctxt_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            ilp_mv_state_t *ps_ilp_mv_state;
            svc_ilp_mv_ctxt_t *ps_ilp_mv_ctxt;

            isvce_process_ctxt_t *ps_proc = ps_codec->as_process + i;

            ps_ilp_mv_ctxt = ps_proc->ps_svc_ilp_mv_ctxt = (svc_ilp_mv_ctxt_t *) pu1_buf;
            pu1_buf += sizeof(svc_ilp_mv_ctxt_t);
            i8_alloc_mem_size -= sizeof(svc_ilp_mv_ctxt_t);

            ps_ilp_mv_ctxt->s_ilp_mv_constants.pv_state = pu1_buf;
            ps_ilp_mv_state = (ilp_mv_state_t *) pu1_buf;
            pu1_buf += sizeof(ilp_mv_state_t);
            i8_alloc_mem_size -= sizeof(ilp_mv_state_t);

            if(0 == i)
            {
                ps_ilp_mv_state->ps_layer_state = (ilp_mv_layer_state_t *) pu1_buf;
                ps_layer_states = ps_ilp_mv_state->ps_layer_state;
                pu1_buf += u1_num_spatial_layers * sizeof(ps_ilp_mv_state->ps_layer_state[0]);
                i8_alloc_mem_size -=
                    u1_num_spatial_layers * sizeof(ps_ilp_mv_state->ps_layer_state[0]);
            }
            else
            {
                ps_ilp_mv_state->ps_layer_state = ps_layer_states;
            }

            ASSERT(i8_alloc_mem_size >= 0);

            if(0 == i)
            {
                for(j = u1_num_spatial_layers - 1; j >= 1; j--)
                {
                    ilp_mv_layer_state_t *ps_layer = &ps_ilp_mv_state->ps_layer_state[j];

                    WORD32 i4_layer_luma_wd =
                        ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_ht =
                        ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_mbs =
                        (i4_layer_luma_wd / MB_SIZE) * (i4_layer_luma_ht / MB_SIZE);

                    ps_layer->ps_mb_states = (ilp_mv_mb_state_t *) pu1_buf;
                    aps_luma_mb_states[j] = ps_layer->ps_mb_states;
                    pu1_buf += i4_layer_luma_mbs * sizeof(ps_layer->ps_mb_states[0]);
                    i8_alloc_mem_size -= u1_num_spatial_layers * sizeof(ps_layer->ps_mb_states[0]);

                    ASSERT(i8_alloc_mem_size >= 0);
                    /* Asserts below verify that
                     * 'ps_codec->s_svc_ilp_data.aps_layer_resampler_props' is initialised
                     */
                    ASSERT(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j].u4_mb_wd ==
                           MB_SIZE);

                    ps_layer->ps_props = &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j];

                    isvce_ilp_mv_layer_state_init(ps_layer, d_spatial_res_ratio, i4_layer_luma_wd,
                                                  i4_layer_luma_ht);
                }
            }
            else
            {
                for(j = u1_num_spatial_layers - 1; j >= 1; j--)
                {
                    ilp_mv_layer_state_t *ps_layer = &ps_ilp_mv_state->ps_layer_state[j];

                    ps_layer->ps_mb_states = aps_luma_mb_states[j];

                    ps_layer->ps_props = &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j];
                }
            }
        }
    }
    else
    {
        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            ps_codec->as_process[i].ps_svc_ilp_mv_ctxt = NULL;
        }
    }
}

static void isvce_get_ilp_mvs_for_me(svc_ilp_mv_ctxt_t *ps_ilp_mv_ctxt)
{
    svc_layer_data_t *ps_ref_layer_data;
    ilp_mv_layer_state_t *ps_layer_state;
    ilp_mv_mb_state_t *ps_mb_state;
    isvce_mb_info_t *ps_ref_mb_info;
    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;
    coordinates_t s_ref_frame_dims;
    coordinates_t s_ref_frame_dims_in_mbs;

    bool b_is_mv_non_identical;
    WORD32 i, j, k;

    ilp_mv_constants_t *ps_ilp_mv_constants = &ps_ilp_mv_ctxt->s_ilp_mv_constants;
    ilp_mv_variables_t *ps_ilp_mv_variables = &ps_ilp_mv_ctxt->s_ilp_mv_variables;
    ilp_mv_outputs_t *ps_ilp_mv_outputs = &ps_ilp_mv_ctxt->s_ilp_mv_outputs;
    ilp_mv_state_t *ps_ilp_mv_state = (ilp_mv_state_t *) ps_ilp_mv_constants->pv_state;
    svc_ilp_data_t *ps_svc_ilp_data = ps_ilp_mv_variables->ps_svc_ilp_data;
    svc_au_data_t *ps_svc_au_data = ps_svc_ilp_data->ps_svc_au_data;
    coordinates_t *ps_mb_pos = &ps_ilp_mv_variables->s_mb_pos;
    const isvce_enc_pu_mv_t s_default_mv = {{0, 0}, -1};

    UWORD8 u1_spatial_layer_id = ps_ilp_mv_variables->u1_spatial_layer_id;
    WORD32 i4_num_ilp_mvs = 0;

    s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;
    s_ref_frame_dims.i4_abscissa =
        ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_width;
    s_ref_frame_dims.i4_ordinate =
        ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_height;
    s_ref_frame_dims_in_mbs.i4_abscissa = s_ref_frame_dims.i4_abscissa / MB_SIZE;
    s_ref_frame_dims_in_mbs.i4_ordinate = s_ref_frame_dims.i4_ordinate / MB_SIZE;

    ps_ref_layer_data = &ps_svc_au_data->ps_svc_layer_data[u1_spatial_layer_id - 1];
    ps_layer_state = &ps_ilp_mv_state->ps_layer_state[u1_spatial_layer_id];
    ps_mb_state =
        &ps_layer_state->ps_mb_states[ps_mb_pos->i4_abscissa +
                                      ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa];

    for(i = 0; i < MAX_PU_IN_MB_COL; i++)
    {
        for(j = 0; j < MAX_PU_IN_MB_ROW; j++)
        {
            b_is_mv_non_identical = true;

            ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0] = s_default_mv;
            ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1] = s_default_mv;

            ps_ref_mb_info =
                &ps_ref_layer_data->ps_mb_info[ps_mb_state->as_mb_positions[i][j].i4_abscissa +
                                               ps_mb_state->as_mb_positions[i][j].i4_ordinate *
                                                   s_ref_frame_dims_in_mbs.i4_abscissa];

            if((ps_ref_mb_info->u2_mb_type == P16x16) || (ps_ref_mb_info->u2_mb_type == B16x16))
            {
                ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs] =
                    ps_ref_mb_info->u2_mb_type;

                ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] =
                    ps_ref_mb_info->as_pu->u1_pred_mode;

                if(ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] != L0)
                {
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1] =
                        ps_ref_mb_info->as_pu->as_me_info[L1];

                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvx =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvx *
                             ps_layer_state->s_mv_scale.i4_abscissa +
                         32768) >>
                        16;
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvy =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvy *
                             ps_layer_state->s_mv_scale.i4_ordinate +
                         32768) >>
                        16;
                }

                if(ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] != L1)
                {
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0] =
                        ps_ref_mb_info->as_pu->as_me_info[L0];

                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvx =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvx *
                             ps_layer_state->s_mv_scale.i4_abscissa +
                         32768) >>
                        16;
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvy =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvy *
                             ps_layer_state->s_mv_scale.i4_ordinate +
                         32768) >>
                        16;
                }

                if(i4_num_ilp_mvs == 0)
                {
                    i4_num_ilp_mvs++;
                }
                else
                {
                    for(k = i4_num_ilp_mvs - 1; k >= 0; k--)
                    {
                        if((ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[k] ==
                            ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs]) &&
                           (ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[k] ==
                            ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs]) &&
                           isvce_check_identical_mv(
                               ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[k],
                               ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs],
                               ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[k]))
                        {
                            b_is_mv_non_identical = false;
                        }
                    }

                    if(b_is_mv_non_identical)
                    {
                        i4_num_ilp_mvs++;
                    }
                }
            }
            else
            {
                ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs] = INVALID_MB_TYPE;
            }
        }
    }

    ps_ilp_mv_outputs->s_ilp_me_cands.u4_num_ilp_mvs = i4_num_ilp_mvs;

    for(i = 0; i < MAX_ILP_MV_IN_NBR_RGN; i++)
    {
        b_is_mv_non_identical = true;

        ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0] = s_default_mv;
        ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1] = s_default_mv;

        if(ps_mb_pos->i4_abscissa + gai1_nbr_ilp_mv_map[i][0] >= 0 &&
           ps_mb_pos->i4_abscissa + gai1_nbr_ilp_mv_map[i][0] < s_frame_dims_in_mbs.i4_abscissa &&
           ps_mb_pos->i4_ordinate + gai1_nbr_ilp_mv_map[i][1] >= 0 &&
           ps_mb_pos->i4_ordinate + gai1_nbr_ilp_mv_map[i][1] < s_frame_dims_in_mbs.i4_ordinate)
        {
            ps_mb_state =
                &ps_layer_state->ps_mb_states[(ps_mb_pos->i4_abscissa + gai1_nbr_ilp_mv_map[i][0]) +
                                              (ps_mb_pos->i4_ordinate + gai1_nbr_ilp_mv_map[i][1]) *
                                                  s_frame_dims_in_mbs.i4_abscissa];

            ps_ref_mb_info =
                &ps_ref_layer_data->ps_mb_info[(ps_mb_state
                                                    ->as_mb_positions[gai1_nbr_ilp_mv_map[i][2]]
                                                                     [gai1_nbr_ilp_mv_map[i][3]]
                                                    .i4_abscissa) +
                                               ps_mb_state
                                                       ->as_mb_positions[gai1_nbr_ilp_mv_map[i][2]]
                                                                        [gai1_nbr_ilp_mv_map[i][3]]
                                                       .i4_ordinate *
                                                   s_ref_frame_dims_in_mbs.i4_abscissa];

            if((ps_ref_mb_info->u2_mb_type == P16x16) || (ps_ref_mb_info->u2_mb_type == B16x16))
            {
                ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs] =
                    ps_ref_mb_info->u2_mb_type;

                ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] =
                    ps_ref_mb_info->as_pu->u1_pred_mode;

                if(ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] != L0)
                {
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1] =
                        ps_ref_mb_info->as_pu->as_me_info[L1];

                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvx =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvx *
                             ps_layer_state->s_mv_scale.i4_abscissa +
                         32768) >>
                        16;
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvy =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L1].s_mv.i2_mvy *
                             ps_layer_state->s_mv_scale.i4_ordinate +
                         32768) >>
                        16;
                }

                if(ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs] != L1)
                {
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0] =
                        ps_ref_mb_info->as_pu->as_me_info[L0];

                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvx =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvx *
                             ps_layer_state->s_mv_scale.i4_abscissa +
                         32768) >>
                        16;
                    ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvy =
                        (ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs][L0].s_mv.i2_mvy *
                             ps_layer_state->s_mv_scale.i4_ordinate +
                         32768) >>
                        16;
                }

                if(i4_num_ilp_mvs == 0)
                {
                    i4_num_ilp_mvs++;
                }
                else
                {
                    for(k = i4_num_ilp_mvs - 1; k >= 0; k--)
                    {
                        if((ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[k] ==
                            ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs]) &&
                           (ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[k] ==
                            ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[i4_num_ilp_mvs]) &&
                           isvce_check_identical_mv(
                               ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[k],
                               ps_ilp_mv_outputs->s_ilp_me_cands.as_mv[i4_num_ilp_mvs],
                               ps_ilp_mv_outputs->s_ilp_me_cands.ae_pred_mode[k]))
                            b_is_mv_non_identical = false;
                    }

                    if(b_is_mv_non_identical)
                    {
                        i4_num_ilp_mvs++;
                    }
                }
            }
            else
            {
                ps_ilp_mv_outputs->s_ilp_me_cands.e_mb_type[i4_num_ilp_mvs] = INVALID_MB_TYPE;
            }
        }
    }

    ps_ilp_mv_outputs->s_ilp_me_cands.u4_num_ilp_mvs_incl_nbrs = i4_num_ilp_mvs;
}

void isvce_get_mb_ilp_mv(svc_ilp_mv_ctxt_t *ps_ilp_mv_ctxt)
{
    svc_layer_data_t *ps_ref_layer_data;
    ilp_mv_layer_state_t *ps_layer_state;
    ilp_mv_mb_state_t *ps_mb_state;
    isvce_mb_info_t *ps_ref_mb_info;
    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;
    coordinates_t s_ref_frame_dims;
    coordinates_t s_ref_frame_dims_in_mbs;

    WORD32 i, j;

    ilp_mv_constants_t *ps_ilp_mv_constants = &ps_ilp_mv_ctxt->s_ilp_mv_constants;
    ilp_mv_variables_t *ps_ilp_mv_variables = &ps_ilp_mv_ctxt->s_ilp_mv_variables;
    ilp_mv_outputs_t *ps_ilp_mv_outputs = &ps_ilp_mv_ctxt->s_ilp_mv_outputs;
    ilp_mv_state_t *ps_ilp_mv_state = (ilp_mv_state_t *) ps_ilp_mv_constants->pv_state;
    svc_ilp_data_t *ps_svc_ilp_data = ps_ilp_mv_variables->ps_svc_ilp_data;
    svc_au_data_t *ps_svc_au_data = ps_svc_ilp_data->ps_svc_au_data;
    coordinates_t *ps_mb_pos = &ps_ilp_mv_variables->s_mb_pos;
    const isvce_enc_pu_mv_t s_default_mv = {{0, 0}, -1};

    UWORD8 u1_spatial_layer_id = ps_ilp_mv_variables->u1_spatial_layer_id;

    s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;
    s_ref_frame_dims.i4_abscissa =
        ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_width;
    s_ref_frame_dims.i4_ordinate =
        ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_height;
    s_ref_frame_dims_in_mbs.i4_abscissa = s_ref_frame_dims.i4_abscissa / MB_SIZE;
    s_ref_frame_dims_in_mbs.i4_ordinate = s_ref_frame_dims.i4_ordinate / MB_SIZE;

    ps_ref_layer_data = &ps_svc_au_data->ps_svc_layer_data[u1_spatial_layer_id - 1];
    ps_layer_state = &ps_ilp_mv_state->ps_layer_state[u1_spatial_layer_id];
    ps_mb_state =
        &ps_layer_state->ps_mb_states[ps_mb_pos->i4_abscissa +
                                      ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa];

    ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0] = s_default_mv;
    ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1] = s_default_mv;

    ps_ref_mb_info = &ps_ref_layer_data->ps_mb_info[ps_mb_state->as_mb_positions[0][0].i4_abscissa +
                                                    ps_mb_state->as_mb_positions[0][0].i4_ordinate *
                                                        s_ref_frame_dims_in_mbs.i4_abscissa];

    if((ps_ref_mb_info->u2_mb_type == P16x16) || (ps_ref_mb_info->u2_mb_type == B16x16))
    {
        ps_ilp_mv_outputs->s_ilp_mv.e_mb_type = ps_ref_mb_info->u2_mb_type;

        ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] = ps_ref_mb_info->as_pu->u1_pred_mode;

        if(ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] != L0)
        {
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1] = ps_ref_mb_info->as_pu->as_me_info[L1];
        }

        if(ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] != L1)
        {
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0] = ps_ref_mb_info->as_pu->as_me_info[L0];
        }
    }
    else
    {
        ps_ilp_mv_outputs->s_ilp_mv.e_mb_type = INVALID_MB_TYPE;
    }

    /* Function call to get non 16x16 ilp mvs for me candidates */
    isvce_get_ilp_mvs_for_me(ps_ilp_mv_ctxt);

    /* Encoder supports only 16x16 partition. */
    /* The code below ensures only 16x16 ILP MV's are used */
    for(i = 0; i < MAX_PU_IN_MB_COL; i++)
    {
        for(j = 0; j < MAX_PU_IN_MB_ROW; j++)
        {
            bool b_unsupported_mv;

            ps_ref_mb_info =
                &ps_ref_layer_data->ps_mb_info[ps_mb_state->as_mb_positions[i][j].i4_abscissa +
                                               ps_mb_state->as_mb_positions[i][j].i4_ordinate *
                                                   s_ref_frame_dims_in_mbs.i4_abscissa];

            b_unsupported_mv =
                (ps_ref_mb_info->u2_mb_type != ps_ilp_mv_outputs->s_ilp_mv.e_mb_type) ||
                (ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] !=
                 ps_ref_mb_info->as_pu->u1_pred_mode) ||
                !isvce_check_identical_mv(ps_ilp_mv_outputs->s_ilp_mv.as_mv[0],
                                          ps_ref_mb_info->as_pu->as_me_info,
                                          ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0]);

            if(b_unsupported_mv)
            {
                ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0] = s_default_mv;
                ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1] = s_default_mv;
                ps_ilp_mv_outputs->s_ilp_mv.e_mb_type = INVALID_MB_TYPE;

                return;
            }
        }
    }

    if(ps_ilp_mv_outputs->s_ilp_mv.e_mb_type != INVALID_MB_TYPE)
    {
        if(ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] != L0)
        {
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1].s_mv.i2_mvx =
                (ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1].s_mv.i2_mvx *
                     ps_layer_state->s_mv_scale.i4_abscissa +
                 32768) >>
                16;
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1].s_mv.i2_mvy =
                (ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L1].s_mv.i2_mvy *
                     ps_layer_state->s_mv_scale.i4_ordinate +
                 32768) >>
                16;
        }

        if(ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] != L1)
        {
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0].s_mv.i2_mvx =
                (ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0].s_mv.i2_mvx *
                     ps_layer_state->s_mv_scale.i4_abscissa +
                 32768) >>
                16;
            ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0].s_mv.i2_mvy =
                (ps_ilp_mv_outputs->s_ilp_mv.as_mv[0][L0].s_mv.i2_mvy *
                     ps_layer_state->s_mv_scale.i4_ordinate +
                 32768) >>
                16;
        }
    }
    else
    {
        ps_ilp_mv_outputs->s_ilp_mv.e_mb_type = INVALID_MB_TYPE;
        ps_ilp_mv_outputs->s_ilp_mv.ae_pred_mode[0] = INVALID_PRED_MODE;
    }
}

void isvce_mvp_idx_eval(isvce_mb_info_t *ps_mb_info, isvce_enc_pu_mv_t *ps_spatial_mvp,
                        isvce_enc_pu_mv_t *ps_ilp_mvp, UWORD8 *pu1_mvd_costs)
{
    if(USE_ILP_MV_AS_MVP && ps_ilp_mvp && !ps_mb_info->u1_is_intra &&
       (ps_mb_info->u2_mb_type != PSKIP) && (ps_mb_info->u2_mb_type != BSKIP) &&
       (ps_mb_info->u2_mb_type != BASE_MODE))
    {
        isvce_enc_pu_mv_t *ps_mv;
        isvce_enc_pu_mv_t *aps_mvps[2];

        WORD32 ai4_mvd_costs[2];
        WORD32 i, j;

        for(i = 0; i < NUM_PRED_DIRS; i++)
        {
            PRED_MODE_T e_pred_mode = (PRED_MODE_T) i;
            PRED_MODE_T e_cmpl_pred_mode = (e_pred_mode == L0) ? L1 : L0;

            if(ps_mb_info->as_pu->u1_pred_mode != e_pred_mode)
            {
                ps_mv = &ps_mb_info->as_pu->as_me_info[e_cmpl_pred_mode];
                aps_mvps[0] = &ps_spatial_mvp[e_cmpl_pred_mode];
                aps_mvps[1] = &ps_ilp_mvp[e_cmpl_pred_mode];

                for(j = 0; j < 2; j++)
                {
                    if((aps_mvps[j]->i1_ref_idx != -1) &&
                       (!j || ((j == 1) && (ps_mv->i1_ref_idx == aps_mvps[j]->i1_ref_idx))))
                    {
                        ai4_mvd_costs[j] =
                            pu1_mvd_costs[ps_mv->s_mv.i2_mvx - aps_mvps[j]->s_mv.i2_mvx] +
                            pu1_mvd_costs[ps_mv->s_mv.i2_mvy - aps_mvps[j]->s_mv.i2_mvy];
                    }
                    else
                    {
                        ai4_mvd_costs[j] = INT32_MAX;
                    }
                }

                ps_mb_info->as_pu->au1_mvp_idx[e_cmpl_pred_mode] =
                    ai4_mvd_costs[0] > ai4_mvd_costs[1];
            }
            else
            {
                ps_mb_info->as_pu->au1_mvp_idx[e_cmpl_pred_mode] = 0;
            }
        }
    }
    else
    {
        ps_mb_info->as_pu->au1_mvp_idx[L0] = 0;
        ps_mb_info->as_pu->au1_mvp_idx[L1] = 0;
    }
}
