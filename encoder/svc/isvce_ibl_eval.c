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
*  isvce_ibl_eval.c
*
* @brief
*  Contains functions used for SVC intra prediction
*
*******************************************************************************
*/
#include <math.h>
#include <limits.h>
#include <stdbool.h>

#include "ih264_typedefs.h"
#include "iv2.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "ih264_padding.h"
#include "isvce_defs.h"
#include "isvce_ibl_private_defs.h"
#include "isvce_ibl_eval.h"
#include "isvce_utils.h"
#include "isvc_intra_resample.h"
#include "isvc_defs.h"

static FORCEINLINE WORD32 isvce_get_num_mb_states(UWORD32 u4_wd, UWORD32 u4_ht)
{
    return (u4_wd / MB_SIZE) * (u4_ht / MB_SIZE);
}

static FORCEINLINE WORD32 isvce_get_phase_array_size(DOUBLE d_spatial_res_ratio, bool b_is_chroma)
{
    return (2 == d_spatial_res_ratio) ? (b_is_chroma ? 3 : 0) : 5;
}

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing residual pred ctxt
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
UWORD32 isvce_get_svc_intra_pred_ctxt_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                           UWORD32 u4_wd, UWORD32 u4_ht)
{
    WORD32 i, j;

    UWORD32 u4_size = 0;

    if(u1_num_spatial_layers > 1)
    {
        u4_size += MAX_PROCESS_CTXT * sizeof(svc_intra_pred_ctxt_t);
        u4_size += MAX_PROCESS_CTXT * sizeof(intra_pred_state_t);
        u4_size += MAX_PROCESS_CTXT * u1_num_spatial_layers * sizeof(intra_pred_layer_state_t);

        for(i = u1_num_spatial_layers - 1; i >= 0; i--)
        {
            WORD32 i4_layer_luma_wd =
                (WORD32) ((DOUBLE) u4_wd /
                          pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_wd_mbs = i4_layer_luma_wd / MB_SIZE;
            WORD32 i4_layer_ht_mbs = i4_layer_luma_ht / MB_SIZE;
            /*Add PAD Mbs */
            WORD32 i4_layer_luma_mbs =
                ((i4_layer_luma_wd / MB_SIZE) + 2) * ((i4_layer_luma_ht / MB_SIZE) + 2);
            WORD32 i4_num_mb_states = isvce_get_num_mb_states(i4_layer_luma_wd, i4_layer_luma_ht);

            for(j = 0; j < NUM_SP_COMPONENTS; j++)
            {
                bool b_is_chroma = ((COMPONENT_TYPE) j) != Y;

                u4_size += i4_num_mb_states * sizeof(intra_pred_mb_state_t);

                /* pi4_ref_array_positions_x */
                u4_size += MAX_REF_ARR_WD_HT * i4_layer_wd_mbs * sizeof(WORD32);

                /* pi4_ref_array_positions_y */
                u4_size += (i4_layer_ht_mbs >> b_is_chroma) * i4_layer_ht_mbs * sizeof(WORD32);

                /* ps_ref_array_phases */
                u4_size += isvce_get_phase_array_size(d_spatial_res_ratio, b_is_chroma) *
                           sizeof(coordinates_t);
            }

            /* pi1_mb_mode */
            u4_size += i4_layer_luma_mbs * sizeof(WORD8);

            /* pu1_refarray_buffer */
            u4_size += MAX_PROCESS_CTXT * TEMP_BUF_SIZE_LUMA * sizeof(UWORD8);

            /* pu1_refarray_cb, pu1_refarray_cr */
            u4_size += MAX_PROCESS_CTXT * (TEMP_BUF_SIZE_CB + TEMP_BUF_SIZE_CR) * sizeof(UWORD8);

            /* pi4_temp_interpolation_buffer */
            u4_size += MAX_PROCESS_CTXT * TEMP_INTERPOLATION_BUF_SIZE * sizeof(WORD32);
        }

        /* intra_pred_outputs_t.s_pred_buf */
        u4_size += MAX_PROCESS_CTXT * MB_SIZE * MB_SIZE * sizeof(UWORD8);

        u4_size += MAX_PROCESS_CTXT * MB_SIZE * MB_SIZE * sizeof(UWORD8);
    }

    return u4_size;
}

static FORCEINLINE WORD32 isvce_get_scaled_pixel_pos(layer_resampler_props_t *ps_layer_props,
                                                     WORD32 i4_pixel_pos, UWORD8 u1_dim_id)
{
    if(1 == u1_dim_id)
    {
        return (((i4_pixel_pos - ps_layer_props->i4_offset_y) *
                     ((WORD64) ps_layer_props->u4_scale_y) +
                 ps_layer_props->i4_add_y) >>
                (ps_layer_props->u4_shift_y - 4)) -
               ps_layer_props->i4_delta_y;
    }
    else
    {
        return (((i4_pixel_pos - ps_layer_props->i4_offset_x) *
                     ((WORD64) ps_layer_props->u4_scale_x) +
                 ps_layer_props->i4_add_x) >>
                (ps_layer_props->u4_shift_x - 4)) -
               ps_layer_props->i4_delta_x;
    }
}

static FORCEINLINE void isvce_ref_array_pos_init(
    layer_resampler_props_t *ps_layer_props, intra_pred_mb_state_t *ps_mb_state,
    coordinates_t *ps_mb_pos, DOUBLE d_spatial_res_ratio, UWORD8 u1_frame_mbs_only_flag,
    UWORD8 u1_field_mb_flag, UWORD8 u1_ref_layer_frame_mbs_only_flag)
{
    if(1.5 == d_spatial_res_ratio)
    {
        UWORD32 i;

        WORD32 *pi4_ref_array_positions_x = ps_mb_state->pi4_ref_array_positions_x;
        WORD32 *pi4_ref_array_positions_y = ps_mb_state->pi4_ref_array_positions_y;
        WORD32 i4_x_offset = ps_mb_state->s_offsets.i4_abscissa;
        WORD32 i4_y_offset = ps_mb_state->s_offsets.i4_ordinate;

        if(0 == ps_mb_pos->i4_abscissa)
        {
            for(i = 0; i < ps_layer_props->u4_mb_ht; i++)
            {
                WORD32 i4_y_ref16;

                WORD32 i4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + i;

                if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
                {
                    i4_yc = i4_yc >> (1 - u1_field_mb_flag);
                }

                i4_y_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_yc, 1);

                pi4_ref_array_positions_y[i] = (i4_y_ref16 >> 4) - i4_y_offset;
            }
        }

        if(0 == ps_mb_pos->i4_ordinate)
        {
            for(i = 0; i < MAX_REF_ARR_WD_HT; i++)
            {
                WORD32 i4_x_ref16;

                WORD32 i4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + i;

                i4_x_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_xc, 0);

                pi4_ref_array_positions_x[i] = (i4_x_ref16 >> 4) - i4_x_offset;
            }
        }
    }
}

static FORCEINLINE void isvce_ref_array_phase_init(
    layer_resampler_props_t *ps_layer_props, intra_pred_mb_state_t *ps_mb_state,
    coordinates_t *ps_mb_pos, DOUBLE d_spatial_res_ratio, UWORD8 u1_frame_mbs_only_flag,
    UWORD8 u1_field_mb_flag, UWORD8 u1_ref_layer_frame_mbs_only_flag)
{
    UWORD32 i, j;

    coordinates_t *ps_ref_array_phases = ps_mb_state->ps_ref_array_phases;

    WORD32 i4_x_offset = ps_mb_state->s_offsets.i4_abscissa;
    WORD32 i4_y_offset = ps_mb_state->s_offsets.i4_ordinate;
    UWORD32 u4_phase_array_idx = 0;

    if(1.5 == d_spatial_res_ratio)
    {
        for(i = 0; i < 3; i++)
        {
            WORD32 i4_y_ref16;

            WORD32 i4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + i;

            if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
            {
                i4_yc = i4_yc >> (1 - u1_field_mb_flag);
            }

            i4_y_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_yc, 1);

            for(j = 0; j < ((0 == i) ? 3 : 1); j++)
            {
                WORD32 i4_x_ref16;

                WORD32 i4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + j;

                i4_x_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_xc, 0);

                ps_ref_array_phases[u4_phase_array_idx].i4_abscissa = i4_x_ref16 & 15;
                ps_ref_array_phases[u4_phase_array_idx].i4_ordinate = i4_y_ref16 & 15;

                u4_phase_array_idx++;
            }
        }
    }
    else
    {
        for(i = 0; i < 2; i++)
        {
            WORD32 i4_y_ref16;

            WORD32 i4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + i;

            if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
            {
                i4_yc = i4_yc >> (1 - u1_field_mb_flag);
            }

            i4_y_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_yc, 1);

            for(j = 0; j < ((0 == i) ? 2 : 1); j++)
            {
                WORD32 i4_x_ref16;

                WORD32 i4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + j;

                i4_x_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_xc, 0);

                ps_ref_array_phases[u4_phase_array_idx].i4_abscissa =
                    (i4_x_ref16 - (16 * i4_x_offset)) & 15;
                ps_ref_array_phases[u4_phase_array_idx].i4_ordinate =
                    (i4_y_ref16 - (16 * i4_y_offset)) & 15;

                u4_phase_array_idx++;
            }
        }
    }
}

static FORCEINLINE void isvce_set_mb_states(layer_resampler_props_t *ps_layer_props,
                                            intra_pred_mb_state_t *ps_mb_states,
                                            coordinates_t *ps_mb_pos, DOUBLE d_spatial_res_ratio,
                                            UWORD32 u4_wd_in_mbs, bool b_is_chroma)
{
    WORD32 i4_x_refmin16;
    WORD32 i4_x_refmax16;
    WORD32 i4_y_refmin16;
    WORD32 i4_y_refmax16;
    WORD32 i4_x_offset, i4_y_offset;

    const UWORD8 u1_frame_mbs_only_flag = 1;
    const UWORD8 u1_ref_layer_frame_mbs_only_flag = 1;
    const UWORD8 u1_field_mb_flag = 0;

    i4_x_refmin16 = isvce_get_scaled_pixel_pos(
        ps_layer_props, ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd, 0);
    i4_x_refmax16 = isvce_get_scaled_pixel_pos(
        ps_layer_props,
        ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + ps_layer_props->u4_mb_wd - 1, 0);

    i4_y_refmin16 = isvce_get_scaled_pixel_pos(
        ps_layer_props, ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht, 1);
    i4_y_refmax16 = isvce_get_scaled_pixel_pos(
        ps_layer_props,
        ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + ps_layer_props->u4_mb_ht - 1, 1);

    i4_x_offset = (i4_x_refmin16 >> 4);
    i4_y_offset = (i4_y_refmin16 >> 4);

    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_offsets.i4_abscissa = i4_x_offset;
    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_offsets.i4_ordinate = i4_y_offset;
    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_ref_array_dims.i4_abscissa = (((i4_x_refmax16 + 15) >> 8) << 4) +
                                        ((WORD32) (ps_layer_props->u4_mb_wd >> 1)) - i4_x_offset +
                                        16;
    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_ref_array_dims.i4_ordinate = (((i4_y_refmax16 + 15) >> 8) << 4) +
                                        ((WORD32) (ps_layer_props->u4_mb_ht >> 1)) - i4_y_offset +
                                        16;

    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_max_pos.i4_abscissa = ((i4_x_refmax16 + 15) >> 4) - i4_x_offset;
    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_max_pos.i4_ordinate = ((i4_y_refmax16 + 15) >> 4) - i4_y_offset;

    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_min_pos.i4_abscissa = (i4_x_refmin16 >> 4) - i4_x_offset;
    ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
        .s_min_pos.i4_ordinate = (i4_y_refmin16 >> 4) - i4_y_offset;

    if((1.5 == d_spatial_res_ratio) &&
       ((0 == ps_mb_pos->i4_abscissa) || (0 == ps_mb_pos->i4_ordinate)))
    {
        WORD32 i4_min, i4_max, i4_xr_index, i4_yr_index, i4_ref_array_wd, i4_ref_array_ht;

        i4_x_offset = i4_x_offset - 2;
        i4_ref_array_wd = ((i4_x_refmax16 + 15) >> 4) - (i4_x_refmin16 >> 4) + 1 + 4;

        i4_min = i4_x_offset;
        i4_xr_index = i4_min - ((i4_min / (WORD32) ps_layer_props->u4_mb_wd) *
                                (WORD32) ps_layer_props->u4_mb_wd);

        if(i4_xr_index < (WORD32) (ps_layer_props->u4_mb_wd >> 1))
        {
            i4_ref_array_wd = i4_ref_array_wd + (ps_layer_props->u4_mb_wd >> 1);
            i4_x_offset = i4_x_offset - ((WORD32) (ps_layer_props->u4_mb_wd >> 1));
        }

        i4_max = ((i4_x_refmax16 + 15) >> 4) + 2;
        i4_xr_index = i4_max - ((i4_max / (WORD32) ps_layer_props->u4_mb_wd) *
                                (WORD32) ps_layer_props->u4_mb_wd);

        if(i4_xr_index >= (WORD32) (ps_layer_props->u4_mb_wd >> 1))
        {
            i4_ref_array_wd = i4_ref_array_wd + (ps_layer_props->u4_mb_wd >> 1);
        }

        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_ref_array_dims.i4_abscissa = i4_ref_array_wd;
        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_offsets.i4_abscissa = i4_x_offset;

        i4_ref_array_ht = ((i4_y_refmax16 + 15) >> 4) - (i4_y_refmin16 >> 4) + 1 + 4;

        i4_y_offset = (i4_y_refmin16 >> 4) - 2;

        i4_min = i4_y_offset;

        i4_yr_index = i4_min - ((i4_min / (WORD32) ps_layer_props->u4_mb_ht) *
                                (WORD32) ps_layer_props->u4_mb_ht);

        if(i4_yr_index < (WORD32) (ps_layer_props->u4_mb_ht >> 1))
        {
            i4_ref_array_ht = i4_ref_array_ht + (ps_layer_props->u4_mb_ht >> 1);
            i4_y_offset = i4_y_offset - ((WORD32) (ps_layer_props->u4_mb_ht >> 1));
        }

        i4_max = ((i4_y_refmax16 + 15) >> 4) + 2;
        i4_yr_index = i4_max - ((i4_max / (WORD32) ps_layer_props->u4_mb_ht) *
                                (WORD32) ps_layer_props->u4_mb_ht);

        if(i4_yr_index >= (WORD32) (ps_layer_props->u4_mb_ht >> 1))
        {
            i4_ref_array_ht = i4_ref_array_ht + (ps_layer_props->u4_mb_ht >> 1);
        }

        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_ref_array_dims.i4_ordinate = i4_ref_array_ht;
        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_offsets.i4_ordinate = i4_y_offset;

        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_max_pos.i4_abscissa = ((i4_x_refmax16 + 15) >> 4) - i4_x_offset;
        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_max_pos.i4_ordinate = ((i4_y_refmax16 + 15) >> 4) - i4_y_offset;

        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_min_pos.i4_abscissa = (i4_x_refmin16 >> 4) - i4_x_offset;
        ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs]
            .s_min_pos.i4_ordinate = (i4_y_refmin16 >> 4) - i4_y_offset;

        isvce_ref_array_pos_init(
            ps_layer_props,
            &ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs],
            ps_mb_pos, d_spatial_res_ratio, u1_frame_mbs_only_flag, u1_field_mb_flag,
            u1_ref_layer_frame_mbs_only_flag);

        isvce_ref_array_phase_init(
            ps_layer_props,
            &ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs],
            ps_mb_pos, d_spatial_res_ratio, u1_frame_mbs_only_flag, u1_field_mb_flag,
            u1_ref_layer_frame_mbs_only_flag);
    }
    else if((2. == d_spatial_res_ratio) &&
            ((0 == ps_mb_pos->i4_abscissa) && (0 == ps_mb_pos->i4_ordinate) && b_is_chroma))
    {
        isvce_ref_array_pos_init(
            ps_layer_props,
            &ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs],
            ps_mb_pos, d_spatial_res_ratio, u1_frame_mbs_only_flag, u1_field_mb_flag,
            u1_ref_layer_frame_mbs_only_flag);

        isvce_ref_array_phase_init(
            ps_layer_props,
            &ps_mb_states[ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * u4_wd_in_mbs],
            ps_mb_pos, d_spatial_res_ratio, u1_frame_mbs_only_flag, u1_field_mb_flag,
            u1_ref_layer_frame_mbs_only_flag);
    }
}

static void isvce_ibl_layer_state_init(intra_pred_layer_state_t *ps_layer_state,
                                       DOUBLE d_spatial_res_ratio, UWORD32 u4_wd, UWORD32 u4_ht,
                                       UWORD8 u1_level_idc, IV_COLOR_FORMAT_T e_color_format)
{
    UWORD32 i, j, k;

    const UWORD8 u1_ref_layer_field_pic_flag = 0;
    const UWORD8 u1_field_pic_flag = 0;
    const UWORD8 u1_frame_mbs_only_flag = 1;
    const UWORD8 u1_ref_layer_frame_mbs_only_flag = 1;
    const UWORD8 u1_bot_field_flag = 0;
    const WORD32 i4_scaled_ref_layer_left_offset = 0;
    const WORD32 i4_scaled_ref_layer_top_offset = 0;
    const WORD32 i4_ref_layer_chroma_phase_x_plus1 = 1;
    const WORD32 i4_ref_layer_chroma_phase_y_plus1 = 1;
    const WORD32 i4_chroma_phase_x_plus1 = 1;
    const WORD32 i4_chroma_phase_y_plus1 = 1;
    const WORD32 i4_sub_wd_chroma = 2;
    const WORD32 i4_sub_ht_chroma = 2;

    ASSERT((IV_YUV_420P == e_color_format) || (IV_YUV_420SP_UV == e_color_format));

    UNUSED(e_color_format);

    for(i = 0; i < NUM_SP_COMPONENTS; i++)
    {
        intra_pred_mb_state_t *ps_mb_states;
        layer_resampler_props_t *ps_layer_props;

        UWORD32 u4_wd_in_mbs;
        UWORD32 u4_ht_in_mbs;

        UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));
        UWORD32 u4_ref_wd = (u4_wd / d_spatial_res_ratio);
        UWORD32 u4_ref_ht = (u4_ht / d_spatial_res_ratio) * (1 + u1_ref_layer_field_pic_flag);
        UWORD32 u4_scaled_wd = u4_wd;
        UWORD32 u4_scaled_ht = u4_ht * (1 + u1_field_pic_flag);

        ps_mb_states =
            u1_is_chroma ? ps_layer_state->ps_chroma_mb_states : ps_layer_state->ps_luma_mb_states;
        ps_layer_props =
            u1_is_chroma ? ps_layer_state->ps_chroma_props : ps_layer_state->ps_luma_props;

        u4_ref_wd = u4_ref_wd >> u1_is_chroma;
        u4_ref_ht = u4_ref_ht >> u1_is_chroma;
        u4_scaled_wd = u4_scaled_wd >> u1_is_chroma;
        u4_scaled_ht = u4_scaled_ht >> u1_is_chroma;

        if(u1_is_chroma)
        {
            ps_layer_props->i4_refphase_x = i4_ref_layer_chroma_phase_x_plus1 - 1;
            ps_layer_props->i4_refphase_y = i4_ref_layer_chroma_phase_y_plus1 - 1;
            ps_layer_props->i4_phase_x = i4_chroma_phase_x_plus1 - 1;
            ps_layer_props->i4_phase_y = i4_chroma_phase_y_plus1 - 1;
            ps_layer_props->u4_sub_wd = i4_sub_wd_chroma;
            ps_layer_props->u4_sub_ht = i4_sub_ht_chroma;
            ps_layer_props->u4_mb_wd = MB_SIZE >> 1;
            ps_layer_props->u4_mb_ht = MB_SIZE >> 1;
        }
        else
        {
            ps_layer_props->i4_refphase_x = 0;
            ps_layer_props->i4_refphase_y = 0;
            ps_layer_props->i4_phase_x = 0;
            ps_layer_props->i4_phase_y = 0;
            ps_layer_props->u4_sub_wd = 1;
            ps_layer_props->u4_sub_ht = 1;
            ps_layer_props->u4_mb_wd = MB_SIZE;
            ps_layer_props->u4_mb_ht = MB_SIZE;
        }

        u4_wd_in_mbs = u4_scaled_wd / ps_layer_props->u4_mb_wd;
        u4_ht_in_mbs = u4_scaled_ht / ps_layer_props->u4_mb_ht;

        if(u1_level_idc <= 30)
        {
            ps_layer_props->u4_shift_x = 16;
            ps_layer_props->u4_shift_y = 16;
        }
        else
        {
            ps_layer_props->u4_shift_x = 31 - isvcd_get_ceil_log2(u4_ref_wd);
            ps_layer_props->u4_shift_y = 31 - isvcd_get_ceil_log2(u4_ref_ht);
        }

        if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
        {
            if(1 == u1_ref_layer_frame_mbs_only_flag)
            {
                ps_layer_props->i4_phase_y = ps_layer_props->i4_phase_y + (4 * u1_bot_field_flag) +
                                             3 - ps_layer_props->u4_sub_ht;
                ps_layer_props->i4_refphase_y = (2 * ps_layer_props->i4_refphase_y) + 2;
            }
            else
            {
                ps_layer_props->i4_phase_y = ps_layer_props->i4_phase_y + 4 * u1_bot_field_flag;
                ps_layer_props->i4_refphase_y =
                    ps_layer_props->i4_refphase_y + (4 * u1_bot_field_flag);
            }
        }

        ps_layer_props->u4_scale_x =
            ((u4_ref_wd << ps_layer_props->u4_shift_x) + (u4_scaled_wd >> 1)) / (u4_scaled_wd);
        ps_layer_props->u4_scale_y =
            ((u4_ref_ht << ps_layer_props->u4_shift_y) + (u4_scaled_ht >> 1)) / (u4_scaled_ht);

        ps_layer_props->i4_offset_x = i4_scaled_ref_layer_left_offset / ps_layer_props->u4_sub_wd;
        ps_layer_props->i4_add_x =
            (((u4_ref_wd * (2 + ps_layer_props->i4_phase_x)) << (ps_layer_props->u4_shift_x - 2)) +
             (u4_scaled_wd >> 1)) /
                u4_scaled_wd +
            (1 << (ps_layer_props->u4_shift_x - 5));
        ps_layer_props->i4_delta_x = 4 * (2 + ps_layer_props->i4_refphase_x);

        if((1 == u1_frame_mbs_only_flag) && (1 == u1_ref_layer_frame_mbs_only_flag))
        {
            ps_layer_props->i4_offset_y =
                i4_scaled_ref_layer_top_offset / ps_layer_props->u4_sub_ht;
            ps_layer_props->i4_add_y = (((u4_ref_ht * (2 + ps_layer_props->i4_phase_y))
                                         << (ps_layer_props->u4_shift_y - 2)) +
                                        (u4_scaled_ht >> 1)) /
                                           u4_scaled_ht +
                                       (1 << (ps_layer_props->u4_shift_y - 5));
            ps_layer_props->i4_delta_y = 4 * (2 + ps_layer_props->i4_refphase_y);
        }
        else
        {
            ps_layer_props->i4_offset_y =
                i4_scaled_ref_layer_top_offset / (2 * ps_layer_props->u4_sub_ht);
            ps_layer_props->i4_add_y = (((u4_ref_ht * (2 + ps_layer_props->i4_phase_y))
                                         << (ps_layer_props->u4_shift_y - 3)) +
                                        (u4_scaled_ht >> 1)) /
                                           u4_scaled_ht +
                                       (1 << (ps_layer_props->u4_shift_y - 5));
            ps_layer_props->i4_delta_y = 2 * (2 + ps_layer_props->i4_refphase_y);
        }

        for(j = 0; j < u4_ht_in_mbs; j++)
        {
            for(k = 0; k < u4_wd_in_mbs; k++)
            {
                coordinates_t s_mb_pos = {k, j};

                isvce_set_mb_states(ps_layer_props, ps_mb_states, &s_mb_pos, d_spatial_res_ratio,
                                    u4_wd_in_mbs, u1_is_chroma);
            }
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
void isvce_intra_pred_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    intra_pred_state_t *ps_intra_pred_state;
    svc_intra_pred_ctxt_t *ps_intra_pred_ctxt;
    intra_pred_mb_state_t *aps_luma_mb_states[MAX_NUM_SPATIAL_LAYERS];
    intra_pred_mb_state_t *aps_chroma_mb_states[MAX_NUM_SPATIAL_LAYERS];

    WORD32 i, j, k, l, m;
    WORD8 *api4_mb_modes[MAX_NUM_SPATIAL_LAYERS];

    isvce_process_ctxt_t *ps_proc = ps_codec->as_process;

    const WORD32 i4_num_proc_ctxts = sizeof(ps_codec->as_process) / sizeof(ps_codec->as_process[0]);
    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size = isvce_get_svc_intra_pred_ctxt_size(
        u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    if(u1_num_spatial_layers > 1)
    {
        for(j = 0; j < i4_num_proc_ctxts; j++)
        {
            ps_proc = &ps_codec->as_process[j];
            ps_intra_pred_ctxt = ps_proc->ps_intra_pred_ctxt = (svc_intra_pred_ctxt_t *) pu1_buf;
            pu1_buf += sizeof(svc_intra_pred_ctxt_t);
            i8_alloc_mem_size -= sizeof(svc_intra_pred_ctxt_t);

            ps_intra_pred_ctxt->s_intra_pred_constants.pv_state = pu1_buf;
            ps_intra_pred_state = (intra_pred_state_t *) pu1_buf;
            pu1_buf += sizeof(intra_pred_state_t);
            i8_alloc_mem_size -= sizeof(intra_pred_state_t);

            ps_intra_pred_state->ps_layer_state = (intra_pred_layer_state_t *) pu1_buf;
            pu1_buf += u1_num_spatial_layers * sizeof(ps_intra_pred_state->ps_layer_state[0]);
            i8_alloc_mem_size -=
                u1_num_spatial_layers * sizeof(ps_intra_pred_state->ps_layer_state[0]);

            ASSERT(i8_alloc_mem_size >= 0);

            for(i = u1_num_spatial_layers - 1; i >= 0; i--)
            {
                intra_pred_layer_state_t *ps_layer_state = &ps_intra_pred_state->ps_layer_state[i];

                WORD32 i4_layer_luma_wd =
                    ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                    0.99;
                WORD32 i4_layer_luma_ht =
                    ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                    0.99;
                WORD32 i4_layer_wd_mbs = i4_layer_luma_wd / MB_SIZE;
                WORD32 i4_layer_ht_mbs = i4_layer_luma_ht / MB_SIZE;
                /* Add PAD MBs on all directions */
                WORD32 i4_layer_luma_mbs =
                    ((i4_layer_luma_wd / MB_SIZE) + 2) * ((i4_layer_luma_ht / MB_SIZE) + 2);
                WORD32 i4_num_mb_states =
                    isvce_get_num_mb_states(i4_layer_luma_wd, i4_layer_luma_ht);

                if(0 == j)
                {
                    UWORD32 au4_ref_xpos_array_size[NUM_SP_COMPONENTS];
                    UWORD32 au4_ref_ypos_array_size[NUM_SP_COMPONENTS];
                    UWORD32 au4_ref_phase_array_size[NUM_SP_COMPONENTS];

                    for(k = 0; k < NUM_SP_COMPONENTS; k++)
                    {
                        bool b_is_chroma = ((COMPONENT_TYPE) k) != Y;

                        au4_ref_xpos_array_size[k] = MAX_REF_ARR_WD_HT;
                        au4_ref_ypos_array_size[k] = (i4_layer_ht_mbs >> b_is_chroma);
                        au4_ref_phase_array_size[k] =
                            isvce_get_phase_array_size(d_spatial_res_ratio, b_is_chroma);
                    }

                    ps_layer_state->ps_luma_mb_states = (intra_pred_mb_state_t *) pu1_buf;
                    aps_luma_mb_states[i] = ps_layer_state->ps_luma_mb_states;
                    pu1_buf += i4_num_mb_states * sizeof(ps_layer_state->ps_luma_mb_states[0]);
                    i8_alloc_mem_size -=
                        i4_num_mb_states * sizeof(ps_layer_state->ps_luma_mb_states[0]);

                    ps_layer_state->ps_chroma_mb_states = (intra_pred_mb_state_t *) pu1_buf;
                    aps_chroma_mb_states[i] = ps_layer_state->ps_chroma_mb_states;
                    pu1_buf += i4_num_mb_states * sizeof(ps_layer_state->ps_chroma_mb_states[0]);
                    i8_alloc_mem_size -=
                        i4_num_mb_states * sizeof(ps_layer_state->ps_chroma_mb_states[0]);

                    if(1.5 == d_spatial_res_ratio)
                    {
                        for(k = 0; k < NUM_SP_COMPONENTS; k++)
                        {
                            bool b_is_chroma = ((COMPONENT_TYPE) k) != Y;

                            WORD32 *pi4_ref_array_positions_x = (WORD32 *) pu1_buf;
                            WORD32 *pi4_ref_array_positions_y =
                                pi4_ref_array_positions_x + MAX_REF_ARR_WD_HT * i4_layer_wd_mbs;
                            coordinates_t *ps_ref_array_phases =
                                (coordinates_t *) (pi4_ref_array_positions_y +
                                                   (i4_layer_ht_mbs >> b_is_chroma) *
                                                       i4_layer_ht_mbs);
                            intra_pred_mb_state_t *ps_mb_state =
                                b_is_chroma ? ps_layer_state->ps_chroma_mb_states
                                            : ps_layer_state->ps_luma_mb_states;

                            for(l = 0; l < i4_layer_ht_mbs; l++)
                            {
                                for(m = 0; m < i4_layer_wd_mbs; m++)
                                {
                                    ps_mb_state[l * i4_layer_wd_mbs + m].pi4_ref_array_positions_x =
                                        pi4_ref_array_positions_x + m * au4_ref_xpos_array_size[k];
                                    ps_mb_state[l * i4_layer_wd_mbs + m].pi4_ref_array_positions_y =
                                        pi4_ref_array_positions_y + l * au4_ref_ypos_array_size[k];

                                    ps_mb_state[l * i4_layer_wd_mbs + m].ps_ref_array_phases =
                                        ps_ref_array_phases;
                                }
                            }

                            pu1_buf += i4_layer_wd_mbs * au4_ref_xpos_array_size[k] *
                                       sizeof(pi4_ref_array_positions_x[0]);
                            pu1_buf += i4_layer_ht_mbs * au4_ref_ypos_array_size[k] *
                                       sizeof(pi4_ref_array_positions_y[0]);
                            pu1_buf += au4_ref_phase_array_size[k] * sizeof(ps_ref_array_phases[0]);
                            i8_alloc_mem_size -= i4_layer_wd_mbs * au4_ref_xpos_array_size[k] *
                                                 sizeof(pi4_ref_array_positions_x[0]);
                            i8_alloc_mem_size -= i4_layer_ht_mbs * au4_ref_ypos_array_size[k] *
                                                 sizeof(pi4_ref_array_positions_y[0]);
                            i8_alloc_mem_size -=
                                au4_ref_phase_array_size[k] * sizeof(ps_ref_array_phases[0]);
                        }
                    }
                    else
                    {
                        intra_pred_mb_state_t *ps_mb_state;
                        coordinates_t *ps_ref_array_phases;

                        for(k = 0; k < NUM_SP_COMPONENTS; k++)
                        {
                            bool b_is_chroma = ((COMPONENT_TYPE) k) != Y;

                            ps_mb_state = b_is_chroma ? ps_layer_state->ps_chroma_mb_states
                                                      : ps_layer_state->ps_luma_mb_states;
                            ps_ref_array_phases = b_is_chroma ? ((coordinates_t *) pu1_buf) : NULL;

                            for(l = 0; l < i4_num_mb_states; l++)
                            {
                                ps_mb_state[l].pi4_ref_array_positions_x = NULL;
                                ps_mb_state[l].pi4_ref_array_positions_y = NULL;
                                ps_mb_state[l].ps_ref_array_phases = ps_ref_array_phases;
                            }
                        }

                        pu1_buf += au4_ref_phase_array_size[U] * sizeof(ps_ref_array_phases[0]);
                        i8_alloc_mem_size -=
                            au4_ref_phase_array_size[U] * sizeof(ps_ref_array_phases[0]);
                    }

                    ps_layer_state->i4_mb_mode_stride = (i4_layer_luma_wd / MB_SIZE) + 2;
                    ps_layer_state->pi1_mb_mode = (WORD8 *) pu1_buf;
                    ps_layer_state->pi1_mb_mode += ps_layer_state->i4_mb_mode_stride + 1;
                    api4_mb_modes[i] = ps_layer_state->pi1_mb_mode;
                    pu1_buf += i4_layer_luma_mbs * sizeof(ps_layer_state->pi1_mb_mode[0]);
                    i8_alloc_mem_size -=
                        u1_num_spatial_layers * sizeof(ps_layer_state->pi1_mb_mode[0]);
                    memset(ps_layer_state->pi1_mb_mode, -1, i4_layer_luma_mbs);

                    if(i > 0)
                    {
                        /* Asserts below verify that
                         * 'ps_codec->s_svc_ilp_data.aps_layer_resampler_props' is initialised
                         */
                        ASSERT(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][i].u4_mb_wd ==
                               MB_SIZE);
                        ASSERT(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][i].u4_mb_wd ==
                               (MB_SIZE / 2));

                        ps_layer_state->ps_luma_props =
                            &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][i];
                        ps_layer_state->ps_chroma_props =
                            &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][i];

                        isvce_ibl_layer_state_init(
                            ps_layer_state, d_spatial_res_ratio, i4_layer_luma_wd, i4_layer_luma_ht,
                            ps_codec->s_cfg.u4_max_level, ps_codec->s_cfg.e_inp_color_fmt);
                    }
                    else
                    {
                        ps_layer_state->ps_luma_props = NULL;
                        ps_layer_state->ps_chroma_props = NULL;
                    }
                }
                else
                {
                    ps_layer_state->ps_luma_mb_states = aps_luma_mb_states[i];
                    ps_layer_state->ps_chroma_mb_states = aps_chroma_mb_states[i];

                    ps_layer_state->i4_mb_mode_stride = (i4_layer_luma_wd / MB_SIZE) + 2;
                    ps_layer_state->pi1_mb_mode = api4_mb_modes[i];

                    if(i > 0)
                    {
                        ps_layer_state->ps_luma_props =
                            &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][i];
                        ps_layer_state->ps_chroma_props =
                            &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][i];
                    }
                    else
                    {
                        ps_layer_state->ps_luma_props = NULL;
                        ps_layer_state->ps_chroma_props = NULL;
                    }
                }

                ps_layer_state->pu1_refarray_buffer = (UWORD8 *) pu1_buf;
                memset(ps_layer_state->pu1_refarray_buffer, 0, TEMP_BUF_SIZE_LUMA * sizeof(UWORD8));
                pu1_buf += TEMP_BUF_SIZE_LUMA * sizeof(UWORD8);
                i8_alloc_mem_size -= TEMP_BUF_SIZE_LUMA * sizeof(UWORD8);

                ps_layer_state->pu1_refarray_cb = (UWORD8 *) pu1_buf;
                memset(ps_layer_state->pu1_refarray_cb, 0, TEMP_BUF_SIZE_CB * sizeof(UWORD8));
                pu1_buf += TEMP_BUF_SIZE_CB * sizeof(UWORD8);
                i8_alloc_mem_size -= TEMP_BUF_SIZE_CB * sizeof(UWORD8);

                ps_layer_state->pu1_refarray_cr = (UWORD8 *) pu1_buf;
                memset(ps_layer_state->pu1_refarray_cr, 0, TEMP_BUF_SIZE_CR * sizeof(UWORD8));
                pu1_buf += TEMP_BUF_SIZE_CR * sizeof(UWORD8);
                i8_alloc_mem_size -= TEMP_BUF_SIZE_CR * sizeof(UWORD8);

                ps_layer_state->pi4_temp_interpolation_buffer = (WORD32 *) pu1_buf;
                pu1_buf += (TEMP_INTERPOLATION_BUF_SIZE * sizeof(WORD32));
                i8_alloc_mem_size -= (TEMP_INTERPOLATION_BUF_SIZE * sizeof(WORD32));

                ASSERT(i8_alloc_mem_size >= 0);
            }
        }

        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            isvce_process_ctxt_t *ps_proc = &ps_codec->as_process[i];
            svc_intra_pred_ctxt_t *ps_intra_pred_ctxt = ps_proc->ps_intra_pred_ctxt;
            yuv_buf_props_t *ps_mb_intra_pred_buf =
                &ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf;

            ps_proc->ps_mb_pred_buf = ps_mb_intra_pred_buf;

            for(j = 0; j < NUM_SP_COMPONENTS; j++)
            {
                buffer_container_t *ps_comp_buf = &ps_mb_intra_pred_buf->as_component_bufs[j];

                ps_comp_buf->pv_data = pu1_buf;
                ps_comp_buf->i4_data_stride = MB_SIZE;
                pu1_buf += MB_SIZE * MB_SIZE * sizeof(UWORD8);
                i8_alloc_mem_size -= MB_SIZE * MB_SIZE * sizeof(WORD8);

                ASSERT(i8_alloc_mem_size >= 0);
            }

            ps_mb_intra_pred_buf->as_component_bufs[V].pv_data = NULL;
            ps_mb_intra_pred_buf->e_color_format = IV_YUV_420SP_UV;
            ps_mb_intra_pred_buf->u1_bit_depth = 16;
            ps_mb_intra_pred_buf->u4_width = MB_SIZE;
            ps_mb_intra_pred_buf->u4_height = MB_SIZE;
        }
    }
    else
    {
        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            isvce_process_ctxt_t *ps_proc = &ps_codec->as_process[i];

            ps_proc->ps_intra_pred_ctxt = NULL;
        }
    }
}

void isvce_intra_sampling_function_selector(intra_sampling_ctxt_t *ps_ctxt,
                                            DOUBLE d_spatial_res_ratio, IV_ARCH_T e_arch)
{
    if(2. == d_spatial_res_ratio)
    {
        switch(e_arch)
        {
#if defined(X86)
            case ARCH_X86_SSE42:
            {
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_horz_interpol_chroma =
                    isvc_horz_interpol_chroma_dyadic_sse42;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_vert_interpol_chroma =
                    isvc_vert_interpol_chroma_dyadic_sse42;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_interpolate_luma =
                    isvc_interpolate_base_luma_dyadic_sse42;

                break;
            }
#elif defined(ARMV8)
            case ARCH_ARM_A53:
            case ARCH_ARM_A57:
            case ARCH_ARM_V8_NEON:
            {
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_horz_interpol_chroma =
                    isvc_horz_interpol_chroma_dyadic_neon;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_vert_interpol_chroma =
                    isvc_vert_interpol_chroma_dyadic_neon;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_interpolate_luma =
                    isvc_interpolate_base_luma_dyadic_neon;

                break;
            }
#elif defined(ARM) && !defined(DISABLE_NEON)
            case ARCH_ARM_A9Q:
            case ARCH_ARM_A9A:
            case ARCH_ARM_A9:
            case ARCH_ARM_A7:
            case ARCH_ARM_A5:
            case ARCH_ARM_A15:
            {
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_horz_interpol_chroma =
                    isvc_horz_interpol_chroma_dyadic_neon;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_vert_interpol_chroma =
                    isvc_vert_interpol_chroma_dyadic_neon;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_interpolate_luma =
                    isvc_interpolate_base_luma_dyadic_neon;

                break;
            }
#endif
            default:
            {
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_horz_interpol_chroma =
                    isvc_horz_interpol_chroma_dyadic;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_vert_interpol_chroma =
                    isvc_vert_interpol_chroma_dyadic;
                ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id].pf_interpolate_luma =
                    isvc_interpolate_base_luma_dyadic;

                break;
            }
        }
    }
}

static void isvce_get_mb_intra_pred(isvce_process_ctxt_t *ps_proc)
{
    mem_element_t s_ref_mb_mode;
    mem_element_t s_inp_luma;
    mem_element_t s_inp_chroma;
    mem_element_t s_out_luma;
    mem_element_t s_out_chroma;

    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;

    WORD32 i4_cur_stride;
    WORD32 i4_ref_stride;
    WORD32 i;

    intra_sampling_ctxt_t s_intra_samp_ctxt[NUM_SP_COMPONENTS];
    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    svc_intra_pred_ctxt_t *ps_intra_pred_ctxt = ps_proc->ps_intra_pred_ctxt;
    intra_pred_state_t *ps_intra_pred_state =
        (intra_pred_state_t *) (ps_intra_pred_ctxt->s_intra_pred_constants.pv_state);
    intra_pred_layer_state_t *ps_layer_state =
        &ps_intra_pred_state->ps_layer_state[ps_proc->u1_spatial_layer_id];
    intra_pred_layer_state_t *ps_ref_layer_state =
        &ps_intra_pred_state->ps_layer_state[ps_proc->u1_spatial_layer_id - 1];

    intra_pred_mb_state_t *ps_luma_mb_state;
    intra_pred_mb_state_t *ps_chroma_mb_state;

    coordinates_t *ps_mb_pos = &ps_intra_pred_ctxt->s_intra_pred_variables.s_mb_pos;
    svc_ilp_data_t *ps_svc_ilp_data = ps_intra_pred_ctxt->s_intra_pred_variables.ps_svc_ilp_data;

    s_frame_dims.i4_abscissa =
        ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate =
        ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;

    ps_luma_mb_state = ps_layer_state->ps_luma_mb_states + ps_mb_pos->i4_abscissa +
                       ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;
    ps_chroma_mb_state = ps_layer_state->ps_chroma_mb_states + ps_mb_pos->i4_abscissa +
                         ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;

    for(i = 0; i < NUM_SP_COMPONENTS; i++)
    {
        UWORD32 u4_ref_wd, u4_ref_ht;

        bool b_is_chroma = (Y != ((COMPONENT_TYPE) i));
        mem_element_t *ps_buf = b_is_chroma ? &s_out_chroma : &s_out_luma;
        intra_pred_mb_state_t *ps_mb_state = b_is_chroma ? ps_chroma_mb_state : ps_luma_mb_state;
        layer_resampler_props_t *ps_layer_props =
            b_is_chroma ? ps_layer_state->ps_chroma_props : ps_layer_state->ps_luma_props;

        s_intra_samp_ctxt[i].i4_res_lyr_id = ps_proc->u1_spatial_layer_id;

        s_intra_samp_ctxt[i].i4_refarray_stride = REF_ARRAY_WIDTH;
        s_intra_samp_ctxt[i].i4_ref_width =
            ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_width;
        s_intra_samp_ctxt[i].i4_ref_height =
            ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_height;

        isvce_intra_sampling_function_selector(&s_intra_samp_ctxt[i],
                                               ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio,
                                               ps_codec->s_cfg.e_arch);

        s_intra_samp_ctxt[i].pu1_refarray_buffer = ps_layer_state->pu1_refarray_buffer;
        s_intra_samp_ctxt[i].pu1_refarray_cb = ps_layer_state->pu1_refarray_cb;
        s_intra_samp_ctxt[i].pu1_refarray_cr = ps_layer_state->pu1_refarray_cr;
        s_intra_samp_ctxt[i].pi4_temp_interpolation_buffer =
            ps_layer_state->pi4_temp_interpolation_buffer;

        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].ps_mb_pos = ps_mb_pos;

        /* Phase is used only by chroma functions */
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_x_phase_0 =
            ps_chroma_mb_state->ps_ref_array_phases[0].i4_abscissa;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_x_phase_1 =
            ps_chroma_mb_state->ps_ref_array_phases[1].i4_abscissa;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_y_phase_0 =
            ps_chroma_mb_state->ps_ref_array_phases[0].i4_ordinate;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_y_phase_1 =
            ps_chroma_mb_state->ps_ref_array_phases[2].i4_ordinate;
        s_intra_samp_ctxt[i]
            .as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id]
            .i1_constrained_intra_rsmpl_flag = 0;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_ref_width =
            ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_width;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i4_ref_height =
            ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_height;

        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i2_x_min_pos =
            ps_mb_state->s_min_pos.i4_abscissa;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i2_x_max_pos =
            ps_mb_state->s_max_pos.i4_abscissa;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i2_y_min_pos =
            ps_mb_state->s_min_pos.i4_ordinate;
        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].i2_y_max_pos =
            ps_mb_state->s_max_pos.i4_ordinate;

        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].ps_phase =
            ps_mb_state->ps_ref_array_phases;

        s_intra_samp_ctxt[i]
            .as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id]
            .pi4_ref_array_positions_x = ps_mb_state->pi4_ref_array_positions_x;
        s_intra_samp_ctxt[i]
            .as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id]
            .pi4_ref_array_positions_y = ps_mb_state->pi4_ref_array_positions_y;

        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].ps_offsets =
            &ps_mb_state->s_offsets;

        s_intra_samp_ctxt[i].as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id].ps_ref_array_dims =
            &ps_mb_state->s_ref_array_dims;

        i4_cur_stride =
            ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.as_component_bufs[i].i4_data_stride;
        ps_buf->pv_buffer =
            (UWORD8 *) (ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.as_component_bufs[i]
                            .pv_data);

        ps_buf->i4_element_size = 1;
        ps_buf->i4_num_element_stride = i4_cur_stride;

        ps_buf = b_is_chroma ? &s_inp_chroma : &s_inp_luma;

        i4_ref_stride = ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1]
                            .as_component_bufs[i]
                            .i4_data_stride;

        u4_ref_wd = ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_width;
        u4_ref_ht =
            ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_height;

        /* For chroma, filteringModeFlag=1 */
        /* If filteringModeFlag=1, interpolation requires samples at an offset of -1
         * along both directions */
        if(ps_proc->s_svc_params.d_spatial_res_ratio == 2.0)
        {
            WORD8 i1_x_odd, i1_y_odd;

            ps_buf->pv_buffer =
                (UWORD8 *) ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1]
                    .as_component_bufs[i]
                    .pv_data +
                (ps_mb_state->s_offsets.i4_abscissa << b_is_chroma) +
                ps_mb_state->s_offsets.i4_ordinate * i4_ref_stride;

            if(!b_is_chroma)
            {
                ps_buf->pv_buffer = ((UWORD8 *) ps_buf->pv_buffer) + -1 + -1 * i4_ref_stride;
            }

            i1_x_odd = (ps_proc->i4_mb_x & 1);
            i1_y_odd = (ps_proc->i4_mb_y & 1);

            if(i1_x_odd)
            {
                ps_buf->pv_buffer = (UWORD8 *) ps_buf->pv_buffer - 8;
            }
            if(i1_y_odd)
            {
                ps_buf->pv_buffer =
                    (UWORD8 *) ps_buf->pv_buffer - ((8 >> b_is_chroma) * i4_ref_stride);
            }
        }
        else
        {
            WORD32 i4_horz_dim = 0;
            WORD32 i4_vert_dim = 0;
            WORD32 i4_dim =
                (WORD32) (ps_mb_state->s_max_pos.i4_abscissa - ps_mb_state->s_min_pos.i4_abscissa) +
                (4 >> b_is_chroma);

            if(i4_dim > i4_horz_dim)
            {
                i4_horz_dim = i4_dim;
            }

            i4_dim =
                (WORD32) (ps_mb_state->s_max_pos.i4_ordinate - ps_mb_state->s_min_pos.i4_ordinate) +
                (4 >> b_is_chroma);

            if(i4_dim > i4_vert_dim)
            {
                i4_vert_dim = i4_dim;
            }

            isvc_intra_resamp_generate_segment_lookup(
                &(s_intra_samp_ctxt[i]
                      .as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id]
                      .as_seg_lookup_horz[0]),
                i4_horz_dim, ps_layer_props->u4_mb_wd, 3);

            isvc_intra_resamp_generate_segment_lookup(
                &(s_intra_samp_ctxt[i]
                      .as_res_lyrs[s_intra_samp_ctxt[i].i4_res_lyr_id]
                      .as_seg_lookup_vert[0]),
                i4_vert_dim, ps_layer_props->u4_mb_ht, 4);

            ps_buf->pv_buffer =
                (UWORD8 *) ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1]
                    .as_component_bufs[i]
                    .pv_data +
                (CLIP3(0, (WORD32) u4_ref_wd - 1, ps_mb_state->s_offsets.i4_abscissa)
                 << b_is_chroma) +
                CLIP3(0, (WORD32) u4_ref_ht - 1, ps_mb_state->s_offsets.i4_ordinate) *
                    i4_ref_stride;
        }

        ps_buf->i4_element_size = 1;
        ps_buf->i4_num_element_stride = i4_ref_stride;
    }

    s_ref_mb_mode.i4_element_size = 1;
    s_ref_mb_mode.i4_num_element_stride =
        (ps_svc_ilp_data->ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id - 1].u4_width >> 4) + 2;
    s_ref_mb_mode.pv_buffer = ps_ref_layer_state->pi1_mb_mode;

    if(ps_proc->s_svc_params.d_spatial_res_ratio == 2.0)
    {
        isvc_intra_samp_mb_dyadic(&s_intra_samp_ctxt[Y], &s_inp_luma, &s_inp_chroma, &s_ref_mb_mode,
                                  &s_out_luma, &s_out_chroma, ps_proc->i4_mb_x, ps_proc->i4_mb_y, 0,
                                  0);
    }
    else
    {
        isvc_intra_samp_mb(&s_intra_samp_ctxt[Y], &s_intra_samp_ctxt[UV], &s_inp_luma,
                           &s_inp_chroma, &s_ref_mb_mode, &s_out_luma, &s_out_chroma);
    }
}

static FORCEINLINE void isvce_get_sad(UWORD8 *pu1_src, UWORD8 *pu1_pred, UWORD32 src_strd,
                                      UWORD32 pred_strd, WORD32 *pi4_distortion, UWORD32 u4_width,
                                      UWORD32 u4_height)
{
    UWORD32 i, j;
    *pi4_distortion = 0;
    for(i = 0; i < u4_width; i++)
    {
        for(j = 0; j < u4_height; j++)
        {
            *pi4_distortion += ABS(pu1_src[j] - pu1_pred[j]);
        }
        pu1_src += src_strd;
        pu1_pred += pred_strd;
    }
}

/**
******************************************************************************
*
* @brief
*  evaluate IBL mode
*
* @par Description
*  This function evaluates IBL mode for the macro-block
*
* @param[in]    ps_proc_ctxt
*  pointer to proc ctxt
*
*  @return      none
*
******************************************************************************
*/
void isvce_evaluate_IBL_mode(isvce_process_ctxt_t *ps_proc)
{
    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    svc_intra_pred_ctxt_t *ps_intra_pred_ctxt = ps_proc->ps_intra_pred_ctxt;

    /* SAD(distortion metric) of a block */
    WORD32 i4_mb_distortion_least = INT_MAX;

    /* cost = distortion + lambda*rate */
    WORD32 i4_mb_cost_least = INT_MAX;

    WORD32 i4_src_strd = ps_proc->s_src_buf_props.as_component_bufs[Y].i4_data_stride;

    UWORD8 *pu1_mb_src = (UWORD8 *) (ps_proc->s_src_buf_props.as_component_bufs[Y].pv_data);

    WORD32 u4_cur_stride =
        ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.as_component_bufs[Y].i4_data_stride;

    UWORD8 *pu1_mb_pred =
        (UWORD8 *) (ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.as_component_bufs[Y]
                        .pv_data);

    ps_intra_pred_ctxt->s_intra_pred_variables.ps_svc_ilp_data = &ps_codec->s_svc_ilp_data;
    ps_intra_pred_ctxt->s_intra_pred_variables.s_mb_pos.i4_abscissa = ps_proc->i4_mb_x;
    ps_intra_pred_ctxt->s_intra_pred_variables.s_mb_pos.i4_ordinate = ps_proc->i4_mb_y;
    ps_intra_pred_ctxt->s_intra_pred_variables.u1_spatial_layer_id = ps_proc->u1_spatial_layer_id;

    isvce_get_mb_intra_pred(ps_proc);

    /* Luma cost */
    isvce_get_sad(pu1_mb_src, pu1_mb_pred, i4_src_strd, u4_cur_stride, &i4_mb_distortion_least,
                  ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.u4_width,
                  ps_intra_pred_ctxt->s_intra_pred_outputs.s_pred_buf.u4_height);

    /* cost = distortion + lambda*rate */
    i4_mb_cost_least = i4_mb_distortion_least;

    /* update the type of the mb if necessary */
    if(i4_mb_cost_least < ps_proc->i4_mb_cost)
    {
        ps_proc->i4_mb_cost = i4_mb_cost_least;
        ps_proc->i4_mb_distortion = i4_mb_distortion_least;
        ps_proc->ps_mb_info->i4_mb_distortion = i4_mb_distortion_least;
        ps_proc->ps_mb_info->u2_mb_type = BASE_MODE;
        ps_proc->ps_mb_info->u1_base_mode_flag = 1;
        ps_proc->ps_mb_info->u1_is_intra = 1;
    }
    else if(ps_proc->ps_mb_info->u2_mb_type != BASE_MODE)
    {
        ps_proc->ps_mb_info->u1_base_mode_flag = 0;
    }
}

void isvce_update_ibl_info(svc_intra_pred_ctxt_t *ps_intra_pred_ctxt, UWORD8 u1_num_spatial_layers,
                           UWORD8 u1_spatial_layer_id, UWORD16 u2_mb_type, WORD32 i4_mb_x,
                           WORD32 i4_mb_y, WORD8 u1_base_mode_flag)
{
    if(u1_num_spatial_layers > 1)
    {
        intra_pred_state_t *ps_intra_pred_state =
            (intra_pred_state_t *) (ps_intra_pred_ctxt->s_intra_pred_constants.pv_state);
        intra_pred_layer_state_t *ps_layer_state =
            &ps_intra_pred_state->ps_layer_state[u1_spatial_layer_id];
        WORD8 i1_is_intra = (u2_mb_type == I4x4 || u2_mb_type == I16x16 || u2_mb_type == I8x8);

        WORD8 *pi1_mb_mode =
            &ps_layer_state->pi1_mb_mode[i4_mb_x + (i4_mb_y * (ps_layer_state->i4_mb_mode_stride))];

        if(u1_base_mode_flag == 1)
        {
            *pi1_mb_mode = SVC_IBL_MB;
        }
        else
        {
            if(i1_is_intra)
            {
                *pi1_mb_mode = SVC_INTRA_MB;
            }
            else
            {
                *pi1_mb_mode = SVC_INTER_MB;
            }
        }
    }
}

void isvce_pad_mb_mode_buf(svc_intra_pred_ctxt_t *ps_intra_pred_ctxt, UWORD8 u1_spatial_layer_id,
                           UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                           UWORD32 u4_ht)
{
    if(u1_num_spatial_layers > 1)
    {
        intra_pred_state_t *ps_intra_pred_state =
            (intra_pred_state_t *) (ps_intra_pred_ctxt->s_intra_pred_constants.pv_state);
        intra_pred_layer_state_t *ps_layer_state =
            &ps_intra_pred_state->ps_layer_state[u1_spatial_layer_id];

        WORD32 i4_layer_luma_wd =
            ((DOUBLE) u4_wd /
             pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - u1_spatial_layer_id)) +
            0.99;
        WORD32 i4_layer_luma_ht =
            ((DOUBLE) u4_ht /
             pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - u1_spatial_layer_id)) +
            0.99;

        WORD32 row, src_strd;
        WORD8 *pu1_src;

        WORD8 *pi1_mb_mode = ps_layer_state->pi1_mb_mode;
        WORD32 i4_mb_mode_stride = ps_layer_state->i4_mb_mode_stride;

        /* Add PAD MBs on all directions */
        i4_layer_luma_wd /= MB_SIZE;
        i4_layer_luma_ht /= MB_SIZE;

        if(d_spatial_res_ratio == 2.0)
        {
            UWORD8 *pu1_mb_mode = (UWORD8 *) pi1_mb_mode;
            /* Pad left */
            ih264_pad_left_luma(pu1_mb_mode, i4_mb_mode_stride, i4_layer_luma_ht, 1);

            /* Pad right */
            ih264_pad_right_luma(pu1_mb_mode + i4_layer_luma_wd, i4_mb_mode_stride,
                                 i4_layer_luma_ht, 1);

            /* Pad top */
            ih264_pad_top(pu1_mb_mode - 1, i4_mb_mode_stride, i4_layer_luma_wd + 2, 1);

            /* Pad bottom */
            ih264_pad_bottom(pu1_mb_mode + (i4_layer_luma_ht * i4_mb_mode_stride) - 1,
                             i4_mb_mode_stride, i4_layer_luma_wd + 2, 1);
        }
        else
        {
            /* Pad left */
            pu1_src = pi1_mb_mode;
            src_strd = i4_mb_mode_stride;
            for(row = 0; row < i4_layer_luma_ht; row++)
            {
                memset(pu1_src - 1, -1, 1);
                pu1_src += src_strd;
            }

            /* Pad right */
            pu1_src = pi1_mb_mode + i4_layer_luma_wd;
            for(row = 0; row < i4_layer_luma_ht; row++)
            {
                memset(pu1_src, -1, 1);
                pu1_src += src_strd;
            }

            /* Pad top */
            pu1_src = pi1_mb_mode - 1;
            memset(pu1_src - src_strd, -1, i4_layer_luma_wd + 2);

            /* Pad bottom */
            pu1_src = pi1_mb_mode + (i4_layer_luma_ht * i4_mb_mode_stride) - 1;
            memset(pu1_src, -1, i4_layer_luma_wd + 2);
        }
    }
}
