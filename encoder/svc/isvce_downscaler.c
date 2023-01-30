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
*  isvce_downscaler.c
*
* @brief
*  Contains downscaler functions required by the SVC encoder
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_get_downscaler_data_size()
*  - isvce_get_downscaler_padding_dims()
*  - isvce_get_downscaler_normalized_filtered_pixel()
*  - isvce_horizontal_downscale_and_transpose()
*  - isvce_process_downscaler()
*  - isvce_initialize_downscaler()
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* system include files */
#include <stdio.h>
#include <stdlib.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "isvc_macros.h"
#include "ih264_platform_macros.h"
#include "iv2.h"
#include "isvc_defs.h"
#include "isvce_defs.h"
#include "isvc_structs.h"
#include "isvc_structs.h"
#include "isvce_downscaler.h"
#include "isvce_downscaler_private_defs.h"

/**
******************************************************************************
* @brief  lanczos filter coefficients for 2x downscaling
* @remarks Though the length of the filter is 8, the
* same coefficients
* are replicated so that 2 rows can be processed at one
* go in SIMD
******************************************************************************
*/
static WORD8 gai1_lanczos_coefficients_2x[NUM_SCALER_FILTER_PHASES][NUM_SCALER_FILTER_TAPS * 2] = {
    {-7, 0, 39, 64, 39, 0, -7, 0, -7, 0, 39, 64, 39, 0, -7, 0},
    {-6, 0, 33, 62, 41, 4, -6, 0, -6, 0, 33, 62, 41, 4, -6, 0},
    {-5, -1, 29, 57, 45, 9, -5, -1, -5, -1, 29, 57, 45, 9, -5, -1},
    {-4, -2, 23, 55, 48, 14, -4, -2, -4, -2, 23, 55, 48, 14, -4, -2},
    {-3, -3, 18, 52, 52, 18, -3, -3, -3, -3, 18, 52, 52, 18, -3, -3},
    {-2, -4, 13, 49, 54, 24, -2, -4, -2, -4, 13, 49, 54, 24, -2, -4},
    {-1, -5, 9, 44, 58, 29, -1, -5, -1, -5, 9, 44, 58, 29, -1, -5},
    {0, -6, 3, 42, 61, 34, 0, -6, 0, -6, 3, 42, 61, 34, 0, -6}};

/**
******************************************************************************
* @brief  lanczos filter coefficients for 1.5x downscaling
* @remarks Though the length of the filter is 8, the same coefficients
* are replicated so that 2 rows can be processed at one go in SIMD.
******************************************************************************
*/
static WORD8 gai1_lanczos_coefficients_3by2x[NUM_SCALER_FILTER_PHASES][NUM_SCALER_FILTER_TAPS * 2] =
    {{0, -11, 32, 86, 32, -11, 0, 0, 0, -11, 32, 86, 32, -11, 0, 0},
     {0, -10, 26, 79, 39, -5, 0, 0, 0, -10, 26, 79, 39, -5, 0, 0},
     {0, -8, 21, 72, 46, 0, -2, 0, 0, -8, 21, 72, 46, 0, -2, 0},
     {0, -6, 15, 66, 52, 3, -3, 0, 0, -6, 15, 66, 52, 3, -3, 0},
     {0, -6, 10, 60, 60, 10, -6, 0, 0, -6, 10, 60, 60, 10, -6, 0},
     {0, -3, 3, 52, 66, 15, -6, 0, 0, -3, 3, 52, 66, 15, -6, 0},
     {0, -2, 0, 46, 72, 21, -8, 0, 0, -2, 0, 46, 72, 21, -8, 0},
     {0, 0, -5, 39, 79, 26, -10, 0, 0, 0, -5, 39, 79, 26, -10, 0}};

/**
*******************************************************************************
*
* @brief
*   gets the memory size required for downscaler
*
* @par Description:
*   returns the memory required by the downscaler context and state structs
*   for allocation.
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

UWORD32 isvce_get_downscaler_data_size(UWORD8 u1_num_spatial_layers, DOUBLE d_scaling_factor,
                                       UWORD32 u4_width, UWORD32 u4_height)
{
    UWORD32 u4_size = 0;

    if(u1_num_spatial_layers > 1)
    {
        u4_size += sizeof(downscaler_state_t);

        u4_size +=
            (u4_height + NUM_SCALER_FILTER_TAPS * 2) * ((UWORD32) (u4_width / d_scaling_factor));
    }

    return u4_size;
}

/**
*******************************************************************************
*
* @brief
*   gets the padding size required for filtering
*
* @par Description:
*   gets the padding size required for filtering
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

void isvce_get_downscaler_padding_dims(padding_dims_t *ps_pad_dims)
{
    ps_pad_dims->u1_left_pad_size = ALIGN8(NUM_SCALER_FILTER_TAPS / 2);
    ps_pad_dims->u1_right_pad_size = ALIGN8(NUM_SCALER_FILTER_TAPS / 2);
    ps_pad_dims->u1_top_pad_size = NUM_SCALER_FILTER_TAPS / 2;
    ps_pad_dims->u1_bottom_pad_size = NUM_SCALER_FILTER_TAPS / 2;
}

/**
*******************************************************************************
*
* @brief
*   processes downscaler
*
* @par Description:
*   calls the function for padding and scaling
*
* @param[in] ps_scaler
*  pointer to downdownscaler context
*
* @param[in] ps_src_buf_props
*  pointer to source buffer props struct
*
* @param[in] u4_blk_wd
*  width of the block to be processed
*
* @param[in] u4_blk_ht
*  height of the block to be processed
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

void isvce_process_downscaler(downscaler_ctxt_t *ps_scaler, yuv_buf_props_t *ps_src_buf_props,
                              yuv_buf_props_t *ps_dst_buf_props, UWORD32 u4_blk_wd,
                              UWORD32 u4_blk_ht)
{
    buffer_container_t s_src_buf;
    buffer_container_t s_dst_buf;

    UWORD32 u4_scaled_block_size_x, u4_scaled_block_size_y;

    downscaler_state_t *ps_scaler_state = (downscaler_state_t *) ps_scaler->pv_scaler_state;

    ASSERT(ps_src_buf_props->e_color_format == IV_YUV_420SP_UV);

    u4_scaled_block_size_x = (UWORD32) (u4_blk_wd / ps_scaler->d_scaling_factor);
    u4_scaled_block_size_y = (UWORD32) (u4_blk_ht / ps_scaler->d_scaling_factor);

    /* luma */
    s_src_buf = ps_src_buf_props->as_component_bufs[Y];
    s_src_buf.pv_data = ((UWORD8 *) s_src_buf.pv_data) - (NUM_SCALER_FILTER_TAPS / 2) -
                        (NUM_SCALER_FILTER_TAPS / 2) * s_src_buf.i4_data_stride;

    s_dst_buf.pv_data = ps_scaler_state->pv_scratch_buf;
    s_dst_buf.i4_data_stride = u4_blk_ht + NUM_SCALER_FILTER_TAPS;

    ps_scaler_state->pf_downscaler(ps_scaler, &s_src_buf, &s_dst_buf, ps_scaler_state->pai1_filters,
                                   u4_scaled_block_size_x, u4_blk_ht + NUM_SCALER_FILTER_TAPS, 0);

    s_src_buf = s_dst_buf;
    s_dst_buf = ps_dst_buf_props->as_component_bufs[Y];

    ps_scaler_state->pf_downscaler(ps_scaler, &s_src_buf, &s_dst_buf, ps_scaler_state->pai1_filters,
                                   u4_scaled_block_size_y, u4_scaled_block_size_x, 0);

    /* chroma */
    u4_blk_ht /= 2;
    u4_scaled_block_size_y /= 2;

    s_src_buf = ps_src_buf_props->as_component_bufs[U];
    s_src_buf.pv_data = ((UWORD8 *) s_src_buf.pv_data) - NUM_SCALER_FILTER_TAPS -
                        (NUM_SCALER_FILTER_TAPS / 2) * s_src_buf.i4_data_stride;

    s_dst_buf.pv_data = ps_scaler_state->pv_scratch_buf;
    s_dst_buf.i4_data_stride = u4_blk_ht + NUM_SCALER_FILTER_TAPS;

    ps_scaler_state->pf_downscaler(ps_scaler, &s_src_buf, &s_dst_buf, ps_scaler_state->pai1_filters,
                                   u4_scaled_block_size_x, u4_blk_ht + NUM_SCALER_FILTER_TAPS, 1);

    s_src_buf = s_dst_buf;
    s_dst_buf = ps_dst_buf_props->as_component_bufs[U];

    ps_scaler_state->pf_downscaler(ps_scaler, &s_src_buf, &s_dst_buf, ps_scaler_state->pai1_filters,
                                   u4_scaled_block_size_y, u4_scaled_block_size_x, 0);
}

/**
*******************************************************************************
*
* @brief
*   normalized dot product computer for downscaler
*
* @par Description:
*   Given the downscaler filter coefficients, source buffer, the function
*   calculates the dot product between them, adds an offset and normalizes it
*
* @param[in] ps_scaler
*  pointer to src buf
*
* @param[in] pi1_filter
*  pointer to filter coefficients
*
* @returns
*
* @remarks
*
*******************************************************************************
*/

static UWORD8 isvce_get_downscaler_normalized_filtered_pixel(UWORD8 *pu1_src, WORD8 *pi1_filter)
{
    WORD32 i;
    WORD32 i4_norm_dot_product;
    UWORD8 u1_out_pixel;
    WORD32 i4_dot_product_sum = 0;
    WORD32 i4_rounding_offset = 1 << (FILTER_COEFF_Q - 1);
    WORD32 i4_normalizing_factor = 1 << FILTER_COEFF_Q;

    for(i = 0; i < NUM_SCALER_FILTER_TAPS; i++)
    {
        i4_dot_product_sum += (pu1_src[i] * pi1_filter[i]);
    }

    i4_norm_dot_product = ((i4_dot_product_sum + i4_rounding_offset) / i4_normalizing_factor);
    u1_out_pixel = (UWORD8) CLIP_U8(i4_norm_dot_product);

    return u1_out_pixel;
}

/**
*******************************************************************************
*
* @brief
*   horizontal scaler function
*
* @par Description:
*   Does horizontal scaling for the given block
*
* @param[in] ps_scaler
*  pointer to downscaler context
*
* @param[in] ps_src
*  pointer to source buffer container
*
* @param[in] ps_dst
*  pointer to destination buffer container
*
* @param[in] pai1_filters
*  pointer to array of downscaler filters
*
* @param[in] u4_blk_wd
*  width of the block after horizontal scaling (output block width)
*
* @param[in] u4_blk_ht
*  height of the current block (input block height)
*
* @param[in] u1_is_chroma
*  flag suggesting whether the buffer is luma or chroma
*
*
* @returns
*
* @remarks
*  The same function is used for vertical scaling too as
*  the horizontally scaled input in stored in transpose fashion.
*
*******************************************************************************
*/

static void isvce_horizontal_downscale_and_transpose(
    downscaler_ctxt_t *ps_scaler, buffer_container_t *ps_src, buffer_container_t *ps_dst,
    FILTER_COEFF_ARRAY pai1_filters, UWORD32 u4_blk_wd, UWORD32 u4_blk_ht, UWORD8 u1_is_chroma)
{
    WORD32 i, j, k;
    UWORD8 u1_phase;
    UWORD8 u1_filtered_out_pixel;
    UWORD8 *pu1_src_j, *pu1_dst_j;
    UWORD8 u1_filtered_out_u_pixel, u1_filtered_out_v_pixel;
    UWORD8 *pu1_in_pixel;
    UWORD8 *pu1_out_pixel;
    WORD8 *pi1_filter_grid;
    UWORD16 u2_full_pixel_inc;
    UWORD8 au1_temp_u_buff[NUM_SCALER_FILTER_TAPS];
    UWORD8 au1_temp_v_buff[NUM_SCALER_FILTER_TAPS];

    downscaler_state_t *ps_scaler_state = (downscaler_state_t *) ps_scaler->pv_scaler_state;

    UWORD32 u4_center_pixel_pos = ps_scaler_state->i4_init_offset;
    UWORD32 u4_src_horz_increments = ps_scaler_state->u4_horz_increment;
    UWORD8 *pu1_src = ps_src->pv_data;
    UWORD32 u4_in_stride = ps_src->i4_data_stride;
    UWORD8 *pu1_dst = ps_dst->pv_data;
    UWORD32 u4_out_stride = ps_dst->i4_data_stride;
    UWORD32 u4_center_pixel_pos_src = u4_center_pixel_pos;

    /* Offset the input so that the input pixel to be processed
    co-incides with the centre of filter (4th coefficient)*/
    pu1_src += (1 + u1_is_chroma);

    ASSERT((1 << DOWNSCALER_Q) == ps_scaler_state->u4_vert_increment);

    if(!u1_is_chroma)
    {
        for(j = 0; j < (WORD32) u4_blk_ht; j++)
        {
            pu1_src_j = pu1_src + (j * u4_in_stride);
            pu1_dst_j = pu1_dst + j;

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_grid = pai1_filters[u1_phase];

                /* Doing the Calculation for current Loop Count  */
                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;
                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);
                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                u1_filtered_out_pixel =
                    isvce_get_downscaler_normalized_filtered_pixel(pu1_in_pixel, pi1_filter_grid);
                *pu1_out_pixel = u1_filtered_out_pixel;

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }
    }
    else
    {
        for(j = 0; j < (WORD32) u4_blk_ht; j++)
        {
            pu1_src_j = pu1_src + (j * u4_in_stride);
            pu1_dst_j = pu1_dst + j;

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_grid = pai1_filters[u1_phase];

                /*Doing the Calculation for current Loop Count  */
                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;
                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);
                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                for(k = 0; k < NUM_SCALER_FILTER_TAPS; k++)
                {
                    au1_temp_u_buff[k] = *(pu1_in_pixel + (2 * k));
                    au1_temp_v_buff[k] = *(pu1_in_pixel + ((2 * k) + 1));
                }

                u1_filtered_out_u_pixel = isvce_get_downscaler_normalized_filtered_pixel(
                    au1_temp_u_buff, pi1_filter_grid);
                u1_filtered_out_v_pixel = isvce_get_downscaler_normalized_filtered_pixel(
                    au1_temp_v_buff, pi1_filter_grid);
                *pu1_out_pixel = u1_filtered_out_u_pixel;
                *(pu1_out_pixel + u4_out_stride) = u1_filtered_out_v_pixel;

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }
    }
}

void isvce_downscaler_function_selector(downscaler_state_t *ps_scaler_state, IV_ARCH_T e_arch)
{
    switch(e_arch)
    {
#if defined(X86)
        case ARCH_X86_SSE42:
        {
            ps_scaler_state->pf_downscaler = isvce_horizontal_downscale_and_transpose_sse42;

            break;
        }
#elif defined(ARMV8)
        case ARCH_ARM_A53:
        case ARCH_ARM_A57:
        case ARCH_ARM_V8_NEON:
        {
            ps_scaler_state->pf_downscaler = isvce_horizontal_downscale_and_transpose_neon;

            break;
        }
#elif !defined(DISABLE_NEON)
        case ARCH_ARM_A9Q:
        case ARCH_ARM_A9A:
        case ARCH_ARM_A9:
        case ARCH_ARM_A7:
        case ARCH_ARM_A5:
        case ARCH_ARM_A15:
        {
            ps_scaler_state->pf_downscaler = isvce_horizontal_downscale_and_transpose_neon;

            break;
        }
#endif
        default:
        {
            ps_scaler_state->pf_downscaler = isvce_horizontal_downscale_and_transpose;

            break;
        }
    }
}

/**
*******************************************************************************
*
* @brief
*   initializes the downscaler context
*
* @par Description:
*   initializes the downscaler context for the given scaling factor
*   with padding size, filter size, etc.
*
* @param[in] ps_scaler
*   pointer downscaler context
*
* @param[in] ps_mem_rec
*   pointer to memory allocated to downscaler process
*
* @param[in] d_scaling_factor
*   scaling reatio of width/ height between two consecutive SVC layers
*
* @param[in] u1_num_spatial_layers
*   scaling reatio of width/ height between two consecutive SVC layers
*
* @param[in] u4_wd
*   width of the input
*
* @param[in] u4_ht
*   height of the input
*
* @param[in] e_arch
*   architecure type
*
* @returns
*
* @remarks
*  when ARM intrinsics are added, update should be done here
*
*******************************************************************************
*/

void isvce_initialize_downscaler(downscaler_ctxt_t *ps_scaler, iv_mem_rec_t *ps_mem_rec,
                                 DOUBLE d_scaling_factor, UWORD8 u1_num_spatial_layers,
                                 UWORD32 u4_in_width, UWORD32 u4_in_height, IV_ARCH_T e_arch)
{
    if(u1_num_spatial_layers > 1)
    {
        downscaler_state_t *ps_scaler_state;

        UWORD8 *pu1_buf = (UWORD8 *) ps_mem_rec->pv_base;

        ps_scaler_state = (downscaler_state_t *) pu1_buf;
        pu1_buf += sizeof(ps_scaler_state[0]);

        ps_scaler_state->pv_scratch_buf = pu1_buf;
        ps_scaler_state->u4_in_wd = u4_in_width;
        ps_scaler_state->u4_in_ht = u4_in_height;

        ps_scaler->pv_scaler_state = ps_scaler_state;
        ps_scaler->d_scaling_factor = d_scaling_factor;
        ps_scaler->u1_num_spatial_layers = u1_num_spatial_layers;

        isvce_downscaler_function_selector(ps_scaler_state, e_arch);

        ps_scaler_state->u4_horz_increment = (UWORD32) (d_scaling_factor * (1 << DOWNSCALER_Q));

        ps_scaler_state->u4_vert_increment = (1 << DOWNSCALER_Q);
        ps_scaler_state->i4_init_offset = 0;
        ps_scaler_state->pai1_filters = (d_scaling_factor == 2.0) ? gai1_lanczos_coefficients_2x
                                                                  : gai1_lanczos_coefficients_3by2x;
    }
}
