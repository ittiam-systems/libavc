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
*  isvce_mode_stat_visualiser.c
*
* @brief
*  Contains functions used for synthesising analysis YUV
*
*******************************************************************************
*/
#include "isvce_defs.h"

#if ENABLE_MODE_STAT_VISUALISER
#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"
#include "isvce_structs.h"
#include "ih264e_fmt_conv.h"
#include "isvce_mode_stat_visualiser.h"

#define MAX_NUM_MB_MODE_VISUALISATIONS 1

static const UWORD8 gau1_output_file_path[] = "out.yuv";

static const double gd_alpha = 0.5;

static const UWORD8 gau1_colors[MAX_NUM_MB_MODE_VISUALISATIONS][NUM_COMPONENTS] = {
    /* Red */
    {81, 90, 240},
};

UWORD32 isvce_get_msv_ctxt_size(UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 u4_size = sizeof(mode_stat_visualiser_t);
    WORD32 i4_num_luma_samples = u4_wd * u4_ht;
    WORD32 i4_num_chroma_samples = i4_num_luma_samples / 4;

    u4_size += (i4_num_luma_samples + i4_num_chroma_samples * 2) * sizeof(UWORD8);

    return u4_size;
}

void isvce_msv_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    mode_stat_visualiser_t *ps_mode_stat_visualiser;
    yuv_buf_props_t *ps_frame_buf;

    WORD32 i;

    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    WORD32 i4_num_luma_samples = u4_wd * u4_ht;
    WORD32 i4_num_chroma_samples = i4_num_luma_samples / 4;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size = isvce_get_msv_ctxt_size(u4_wd, u4_ht);

    ps_mode_stat_visualiser = ps_codec->ps_mode_stat_visualiser =
        (mode_stat_visualiser_t *) pu1_buf;
    pu1_buf += sizeof(ps_mode_stat_visualiser[0]);
    i8_alloc_mem_size -= sizeof(ps_mode_stat_visualiser[0]);

    ps_frame_buf = &ps_mode_stat_visualiser->s_frame_buf;

    ps_mode_stat_visualiser->ps_output_file = fopen((const char *) gau1_output_file_path, "w");

    ps_frame_buf->e_color_format = IV_YUV_420P;
    ps_frame_buf->u1_bit_depth = 8;
    ps_frame_buf->u4_width = u4_wd;
    ps_frame_buf->u4_height = u4_ht;

    for(i = 0; i < NUM_COMPONENTS; i++)
    {
        UWORD8 u1_is_chroma = (((COMPONENT_TYPE) i) != Y);
        UWORD32 u4_buf_size = u1_is_chroma ? i4_num_chroma_samples : i4_num_luma_samples;
        UWORD32 u4_stride = u4_wd >> u1_is_chroma;

        ps_frame_buf->as_component_bufs[i].pv_data = pu1_buf;
        ps_frame_buf->as_component_bufs[i].i4_data_stride = u4_stride;

        pu1_buf += u4_buf_size;
        i8_alloc_mem_size -= u4_buf_size;
    }

    ASSERT(i8_alloc_mem_size >= 0);
}

void isvce_msv_ctxt_delete(mode_stat_visualiser_t *ps_mode_stat_visualiser)
{
    fclose(ps_mode_stat_visualiser->ps_output_file);
}

void isvce_msv_get_input_frame(mode_stat_visualiser_t *ps_mode_stat_visualiser,
                               isvce_inp_buf_t *ps_inp_buf)
{
    svc_params_t *ps_svc_params = &ps_inp_buf->s_svc_params;
    yuv_buf_props_t *ps_target_layer_yuv_buf =
        &ps_inp_buf->as_layer_yuv_buf_props[ps_svc_params->u1_num_spatial_layers - 1];
    yuv_buf_props_t *ps_frame_buf = &ps_mode_stat_visualiser->s_frame_buf;

    ASSERT(ps_target_layer_yuv_buf->u4_width == ps_frame_buf->u4_width);
    ASSERT(ps_target_layer_yuv_buf->u4_height == ps_frame_buf->u4_height);
    ASSERT(ps_target_layer_yuv_buf->u1_bit_depth == ps_frame_buf->u1_bit_depth);
    ASSERT(ps_target_layer_yuv_buf->e_color_format == IV_YUV_420SP_UV);
    ASSERT(ps_frame_buf->u1_bit_depth == IV_YUV_420P);
    ASSERT(ps_target_layer_yuv_buf->as_component_bufs[U].i4_data_stride ==
           ps_target_layer_yuv_buf->as_component_bufs[V].i4_data_stride);

    isvce_fmt_conv_420sp_to_420p(
        ps_target_layer_yuv_buf->as_component_bufs[Y].pv_data,
        ps_target_layer_yuv_buf->as_component_bufs[UV].pv_data,
        ps_frame_buf->as_component_bufs[Y].pv_data, ps_frame_buf->as_component_bufs[U].pv_data,
        ps_frame_buf->as_component_bufs[V].pv_data, ps_frame_buf->u4_width, ps_frame_buf->u4_height,
        ps_target_layer_yuv_buf->as_component_bufs[Y].i4_data_stride,
        ps_target_layer_yuv_buf->as_component_bufs[UV].i4_data_stride,
        ps_frame_buf->as_component_bufs[Y].i4_data_stride,
        ps_frame_buf->as_component_bufs[U].i4_data_stride, 1, 0);
}

void isvce_msv_set_mode(mode_stat_visualiser_t *ps_mode_stat_visualiser,
                        isvce_mb_info_t *ps_mb_info, coordinates_t *ps_mb_pos)
{
    UWORD32 i, j, k;

    for(i = 0; i < NUM_COMPONENTS; i++)
    {
        UWORD8 u1_is_chroma = (((COMPONENT_TYPE) i) != Y);
        UWORD32 u4_wd = MB_SIZE >> u1_is_chroma;
        UWORD32 u4_ht = MB_SIZE >> u1_is_chroma;
        UWORD8 *pu1_buf = ps_mode_stat_visualiser->s_frame_buf.as_component_bufs[i].pv_data;
        WORD32 i4_stride = ps_mode_stat_visualiser->s_frame_buf.as_component_bufs[i].i4_data_stride;

        pu1_buf += ps_mb_pos->i4_abscissa * u4_wd + ps_mb_pos->i4_ordinate * u4_ht * i4_stride;

        for(j = 0; j < u4_ht; j++)
        {
            for(k = 0; k < u4_wd; k++)
            {
                if(ps_mb_info->u1_residual_prediction_flag)
                {
                    pu1_buf[k + j * i4_stride] =
                        (UWORD8) (gd_alpha * gau1_colors[0][i] +
                                  (1. - gd_alpha) * pu1_buf[k + j * i4_stride] + 0.5);
                }
            }
        }
    }
}

void isvce_msv_dump_visualisation(mode_stat_visualiser_t *ps_mode_stat_visualiser)
{
    WORD32 i;

    FILE *ps_output_file = ps_mode_stat_visualiser->ps_output_file;
    yuv_buf_props_t *ps_frame_buf = &ps_mode_stat_visualiser->s_frame_buf;

    for(i = 0; i < NUM_COMPONENTS; i++)
    {
        UWORD8 u1_is_chroma = (((COMPONENT_TYPE) i) != Y);
        UWORD32 u4_wd = ps_frame_buf->u4_width >> u1_is_chroma;
        UWORD32 u4_ht = ps_frame_buf->u4_height >> u1_is_chroma;
        UWORD32 u4_size = u4_wd * u4_ht;

        ASSERT(u4_wd == ps_frame_buf->as_component_bufs[i].i4_data_stride);

        fwrite(ps_frame_buf->as_component_bufs[i].pv_data, sizeof(UWORD8), u4_size, ps_output_file);
    }
}
#endif
