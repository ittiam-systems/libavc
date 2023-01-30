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
******************************************************************************
* @file isvce_downscaler_sse42.c
*
* @brief
*  This file contains the x86 SIMD version of the function which does
*  horizontal scaling and transpose
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - isvce_horizontal_downscale_and_transpose_sse42()
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

/* User include files */
#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_platform_macros.h"
#include "isvc_defs.h"
#include "isvce_defs.h"
#include "isvc_structs.h"
#include "isvce_downscaler_private_defs.h"

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

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

void isvce_horizontal_downscale_and_transpose_sse42(
    downscaler_ctxt_t *ps_scaler, buffer_container_t *ps_src, buffer_container_t *ps_dst,
    FILTER_COEFF_ARRAY pai1_filters, UWORD32 u4_blk_wd, UWORD32 u4_blk_ht, UWORD8 u1_is_chroma)
{
    WORD32 i, j;
    UWORD8 u1_phase;
    UWORD8 *pu1_src_j, *pu1_dst_j;
    WORD32 i4_temp_pixel_holder;
    UWORD32 u4_num_iterations_vertical_by_16;
    UWORD32 u4_rem_vert_loop;
    UWORD8 *pu1_in_pixel;
    UWORD8 *pu1_out_pixel;
    WORD8 *pi1_filter_for_grid;
    UWORD16 u2_full_pixel_inc;

    __m128i src_temp_0, src_temp_1, src_temp_2, src_temp_3, src_temp_4, src_temp_5, src_temp_6,
        src_temp_7;

    __m128i reg_all_1s, reg_64val_32bit, reg_all_0s, filt_coeff_grid, reg_shuffle;

    __m128i reg_01_16x8b, reg_02_16x8b, reg_03_16x8b, reg_04_16x8b, reg_05_16x8b;

    downscaler_state_t *ps_scaler_state = (downscaler_state_t *) ps_scaler->pv_scaler_state;

    UWORD32 u4_center_pixel_pos = ps_scaler_state->i4_init_offset;
    UWORD32 u4_src_vert_increments = ps_scaler_state->u4_vert_increment;
    UWORD32 u4_src_horz_increments = ps_scaler_state->u4_horz_increment;

    UWORD8 *pu1_src = ps_src->pv_data;
    UWORD32 u4_in_stride = ps_src->i4_data_stride;
    UWORD8 *pu1_dst = ps_dst->pv_data;
    UWORD32 u4_out_stride = ps_dst->i4_data_stride;
    UWORD32 u4_center_pixel_pos_src = u4_center_pixel_pos;

    ASSERT((1 << DOWNSCALER_Q) == u4_src_vert_increments);

    reg_all_1s = _mm_set1_epi16((short) 1);
    reg_64val_32bit = _mm_set1_epi32((int) 64);
    reg_all_0s = _mm_setzero_si128();
    reg_shuffle = _mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0);

    u4_num_iterations_vertical_by_16 = u4_blk_ht >> 4;
    u4_rem_vert_loop = u4_blk_ht % 16;

    /* Offset the input so that the input pixel to be processed
       co-incides with the centre of filter (4th coefficient)*/
    pu1_src += (1 + u1_is_chroma);

    if(!u1_is_chroma)
    {
        for(j = 0; j < (WORD32) u4_num_iterations_vertical_by_16; j++)
        {
            pu1_src_j = pu1_src + ((j << 4) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 4);

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_for_grid = pai1_filters[u1_phase];

                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                filt_coeff_grid = _mm_loadu_si128((__m128i *) pi1_filter_for_grid);
                /******************************************************/
                /* This loop is going vertically in bottom direction */
                /* but the output pixels are stored in horizontal    */
                /* direction in transpose manner                     */
                /******************************************************/

                /*For row 0,1*/
                src_temp_0 = _mm_loadl_epi64((__m128i *) pu1_in_pixel);
                src_temp_1 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride));
                /*next transfer the 8 pixels from temp_2 to temp_1 higher bits 64-127*/
                src_temp_0 = _mm_unpacklo_epi64(src_temp_0, src_temp_1);

                /*For row 2,3*/
                src_temp_2 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 2));

                src_temp_3 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 3));

                src_temp_2 = _mm_unpacklo_epi64(src_temp_2, src_temp_3);

                reg_01_16x8b = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);

                /*multiply with filter coeffs to get 16 bit results*/
                reg_02_16x8b = _mm_maddubs_epi16(src_temp_2, filt_coeff_grid);

                reg_01_16x8b = _mm_hadd_epi16(reg_01_16x8b, reg_02_16x8b);
                /*add adjacent 16 bit values to get 32 bit values*/
                reg_01_16x8b = _mm_madd_epi16(reg_01_16x8b, reg_all_1s);

                /*Add offset of 64 for rounding each out pixel value*/
                reg_01_16x8b = _mm_add_epi32(reg_01_16x8b, reg_64val_32bit);
                /*Divide by 128 each out pixel value*/
                reg_01_16x8b = _mm_srli_epi32(reg_01_16x8b, 7);

                /*For row 4,5*/
                src_temp_4 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 4));

                src_temp_5 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 5));

                src_temp_4 = _mm_unpacklo_epi64(src_temp_4, src_temp_5);

                /*For row 6,7*/
                src_temp_6 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 6));

                src_temp_7 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 7));

                src_temp_6 = _mm_unpacklo_epi64(src_temp_6, src_temp_7);

                reg_03_16x8b = _mm_maddubs_epi16(src_temp_4, filt_coeff_grid);

                reg_04_16x8b = _mm_maddubs_epi16(src_temp_6, filt_coeff_grid);

                reg_03_16x8b = _mm_hadd_epi16(reg_03_16x8b, reg_04_16x8b);

                reg_03_16x8b = _mm_madd_epi16(reg_03_16x8b, reg_all_1s);

                /*next add 2 adjacent 32 bit values to get a single 32 bit
                **value in each row
                */

                /*Add offset of 64 for rounding each out pixel value*/
                reg_03_16x8b = _mm_add_epi32(reg_03_16x8b, reg_64val_32bit);
                /*Divide by 128 each out pixel value*/
                reg_03_16x8b = _mm_srli_epi32(reg_03_16x8b, 7);

                /*pack the lower 16 bit values corresponding to the 8 output
                pixels from reg1 and reg 2*/
                reg_01_16x8b = _mm_packus_epi32(reg_01_16x8b, reg_03_16x8b);

                /*For row 8,9*/
                src_temp_0 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + 8 * u4_in_stride));

                src_temp_1 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + 9 * u4_in_stride));

                /*next transfer the 8 pixels from temp_2 to temp_1 higher bits 64-127*/
                src_temp_0 = _mm_unpacklo_epi64(src_temp_0, src_temp_1);

                /*For row 10,11*/
                src_temp_2 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 10));

                src_temp_3 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 11));

                src_temp_2 = _mm_unpacklo_epi64(src_temp_2, src_temp_3);

                reg_02_16x8b = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);

                /*multiply with filter coeffs to get 16 bit results*/
                reg_03_16x8b = _mm_maddubs_epi16(src_temp_2, filt_coeff_grid);

                reg_02_16x8b = _mm_hadd_epi16(reg_02_16x8b, reg_03_16x8b);
                /*add adjacent 16 bit values to get 32 bit values*/
                reg_02_16x8b = _mm_madd_epi16(reg_02_16x8b, reg_all_1s);

                /*next add 2 adjacent 32 bit values to get a single
                32 bit value in each row*/

                /*Add offset of 64 for rounding each out pixel value*/
                reg_02_16x8b = _mm_add_epi32(reg_02_16x8b, reg_64val_32bit);
                /*Divide by 128 each out pixel value*/
                reg_02_16x8b = _mm_srli_epi32(reg_02_16x8b, 7);

                /*For row 12,13*/
                src_temp_4 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 12));

                src_temp_5 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 13));

                src_temp_4 = _mm_unpacklo_epi64(src_temp_4, src_temp_5);

                /*For row 14,15*/
                src_temp_6 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 14));

                src_temp_7 = _mm_loadl_epi64((__m128i *) (pu1_in_pixel + u4_in_stride * 15));

                src_temp_6 = _mm_unpacklo_epi64(src_temp_6, src_temp_7);

                reg_04_16x8b = _mm_maddubs_epi16(src_temp_4, filt_coeff_grid);

                reg_05_16x8b = _mm_maddubs_epi16(src_temp_6, filt_coeff_grid);

                reg_04_16x8b = _mm_hadd_epi16(reg_04_16x8b, reg_05_16x8b);
                /*add adjacent 16 bit values to get 32 bit values*/
                reg_04_16x8b = _mm_madd_epi16(reg_04_16x8b, reg_all_1s);

                /*next add 2 adjacent 32 bit values to get a single
                32 bit value in each row*/

                /*Add offset of 64 for rounding each out pixel value*/
                reg_04_16x8b = _mm_add_epi32(reg_04_16x8b, reg_64val_32bit);
                /*Divide by 128 each out pixel value*/
                reg_04_16x8b = _mm_srli_epi32(reg_04_16x8b, 7);

                /*pack the lower 16 bit values corresponding to the 8 output
                pixels from reg1 and reg 2*/
                reg_02_16x8b = _mm_packus_epi32(reg_02_16x8b, reg_04_16x8b);

                /*next get saturated 8 bit output pixel values for row 0-15*/
                reg_01_16x8b = _mm_packus_epi16(reg_01_16x8b, reg_02_16x8b);

                /*Store the 16 output values*/
                _mm_storeu_si128((__m128i *) pu1_out_pixel, reg_01_16x8b);

                pu1_out_pixel += 16;

                pu1_in_pixel += ((u4_src_vert_increments * (u4_in_stride << 4)) >> DOWNSCALER_Q);

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }

        /*if height is not a multiple of 8 process 2 rows at a
        time for the remaining rows*/
        if(u4_rem_vert_loop)
        {
            pu1_src_j = pu1_src + ((j << 4) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 4);

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_for_grid = pai1_filters[u1_phase];

                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                filt_coeff_grid = _mm_loadu_si128((__m128i *) pi1_filter_for_grid);

                for(j = u4_rem_vert_loop; j > 0; j--)
                {
                    src_temp_0 = _mm_loadl_epi64((__m128i const *) pu1_in_pixel);

                    src_temp_0 = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);

                    src_temp_0 = _mm_madd_epi16(src_temp_0, reg_all_1s);

                    reg_01_16x8b = _mm_hadd_epi32(src_temp_0, reg_all_0s);

                    /*Add offset of 64 for rounding each out pixel value*/
                    reg_01_16x8b = _mm_add_epi32(reg_01_16x8b, reg_64val_32bit);
                    /*Divide by 128 each out pixel value*/
                    reg_01_16x8b = _mm_srli_epi32(reg_01_16x8b, (int) 7);

                    reg_01_16x8b = _mm_packus_epi32(reg_01_16x8b, reg_all_0s);

                    /*next get saturated 8 bit output pixel values*/
                    reg_01_16x8b = _mm_packus_epi16(reg_01_16x8b, reg_all_0s);

                    /*Store the 1 output value*/
                    *pu1_out_pixel = (UWORD8) _mm_cvtsi128_si32(reg_01_16x8b);

                    pu1_in_pixel += (u4_src_vert_increments * u4_in_stride) >> DOWNSCALER_Q;

                    pu1_out_pixel++;
                }
                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }
    }

    else /* for chroma */
    {
        for(j = 0; j < (WORD32) u4_num_iterations_vertical_by_16; j++)
        {
            pu1_src_j = pu1_src + ((j << 4) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 4);

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_for_grid = pai1_filters[u1_phase];

                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                filt_coeff_grid = _mm_loadu_si128((__m128i *) pi1_filter_for_grid);
                /******************************************************/
                /* This loop is going vertically in bottom direction */
                /* but the output pixels are stored in horizontal    */
                /* direction in transpose manner                     */
                /******************************************************/

                /*Load 16 values shuffle to separate Cb and Cr and process*/

                src_temp_0 = _mm_loadu_si128((__m128i *) pu1_in_pixel);
                src_temp_1 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride));

                src_temp_2 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 2));

                src_temp_3 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 3));

                src_temp_0 = _mm_shuffle_epi8(src_temp_0, reg_shuffle);
                src_temp_1 = _mm_shuffle_epi8(src_temp_1, reg_shuffle);
                src_temp_2 = _mm_shuffle_epi8(src_temp_2, reg_shuffle);
                src_temp_3 = _mm_shuffle_epi8(src_temp_3, reg_shuffle);

                reg_01_16x8b = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);
                reg_02_16x8b = _mm_maddubs_epi16(src_temp_1, filt_coeff_grid);

                reg_01_16x8b = _mm_hadd_epi16(reg_01_16x8b, reg_02_16x8b);

                reg_01_16x8b = _mm_madd_epi16(reg_01_16x8b, reg_all_1s);

                reg_01_16x8b = _mm_add_epi32(reg_01_16x8b, reg_64val_32bit);

                reg_01_16x8b = _mm_srli_epi32(reg_01_16x8b, (int) 7);

                reg_03_16x8b = _mm_maddubs_epi16(src_temp_2, filt_coeff_grid);
                reg_04_16x8b = _mm_maddubs_epi16(src_temp_3, filt_coeff_grid);

                src_temp_4 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 4));

                src_temp_5 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 5));

                src_temp_6 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 6));

                src_temp_7 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 7));

                src_temp_4 = _mm_shuffle_epi8(src_temp_4, reg_shuffle);
                src_temp_5 = _mm_shuffle_epi8(src_temp_5, reg_shuffle);
                src_temp_6 = _mm_shuffle_epi8(src_temp_6, reg_shuffle);
                src_temp_7 = _mm_shuffle_epi8(src_temp_7, reg_shuffle);

                reg_03_16x8b = _mm_hadd_epi16(reg_03_16x8b, reg_04_16x8b);

                reg_03_16x8b = _mm_madd_epi16(reg_03_16x8b, reg_all_1s);

                reg_03_16x8b = _mm_add_epi32(reg_03_16x8b, reg_64val_32bit);

                reg_03_16x8b = _mm_srli_epi32(reg_03_16x8b, (int) 7);

                reg_01_16x8b = _mm_packus_epi32(reg_01_16x8b, reg_03_16x8b);

                reg_02_16x8b = _mm_maddubs_epi16(src_temp_4, filt_coeff_grid);
                reg_04_16x8b = _mm_maddubs_epi16(src_temp_5, filt_coeff_grid);

                reg_02_16x8b = _mm_hadd_epi16(reg_02_16x8b, reg_04_16x8b);

                reg_02_16x8b = _mm_madd_epi16(reg_02_16x8b, reg_all_1s);

                reg_02_16x8b = _mm_add_epi32(reg_02_16x8b, reg_64val_32bit);

                reg_02_16x8b = _mm_srli_epi32(reg_02_16x8b, (int) 7);

                reg_03_16x8b = _mm_maddubs_epi16(src_temp_6, filt_coeff_grid);
                reg_04_16x8b = _mm_maddubs_epi16(src_temp_7, filt_coeff_grid);

                reg_03_16x8b = _mm_hadd_epi16(reg_03_16x8b, reg_04_16x8b);

                reg_03_16x8b = _mm_madd_epi16(reg_03_16x8b, reg_all_1s);

                reg_03_16x8b = _mm_add_epi32(reg_03_16x8b, reg_64val_32bit);

                reg_03_16x8b = _mm_srli_epi32(reg_03_16x8b, (int) 7);

                reg_02_16x8b = _mm_packus_epi32(reg_02_16x8b, reg_03_16x8b);

                reg_01_16x8b = _mm_packus_epi16(reg_01_16x8b, reg_02_16x8b);

                reg_01_16x8b = _mm_shuffle_epi8(reg_01_16x8b, reg_shuffle);

                src_temp_0 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + 8 * u4_in_stride));

                src_temp_1 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + 9 * u4_in_stride));

                src_temp_2 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 10));

                src_temp_3 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 11));

                src_temp_0 = _mm_shuffle_epi8(src_temp_0, reg_shuffle);
                src_temp_1 = _mm_shuffle_epi8(src_temp_1, reg_shuffle);
                src_temp_2 = _mm_shuffle_epi8(src_temp_2, reg_shuffle);
                src_temp_3 = _mm_shuffle_epi8(src_temp_3, reg_shuffle);

                reg_02_16x8b = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);
                reg_03_16x8b = _mm_maddubs_epi16(src_temp_1, filt_coeff_grid);

                reg_02_16x8b = _mm_hadd_epi16(reg_02_16x8b, reg_03_16x8b);

                reg_02_16x8b = _mm_madd_epi16(reg_02_16x8b, reg_all_1s);

                reg_02_16x8b = _mm_add_epi32(reg_02_16x8b, reg_64val_32bit);

                reg_02_16x8b = _mm_srli_epi32(reg_02_16x8b, (int) 7);

                reg_04_16x8b = _mm_maddubs_epi16(src_temp_2, filt_coeff_grid);
                reg_05_16x8b = _mm_maddubs_epi16(src_temp_3, filt_coeff_grid);

                src_temp_4 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 12));

                src_temp_5 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 13));

                src_temp_6 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 14));

                src_temp_7 = _mm_loadu_si128((__m128i *) (pu1_in_pixel + u4_in_stride * 15));

                src_temp_4 = _mm_shuffle_epi8(src_temp_4, reg_shuffle);
                src_temp_5 = _mm_shuffle_epi8(src_temp_5, reg_shuffle);
                src_temp_6 = _mm_shuffle_epi8(src_temp_6, reg_shuffle);
                src_temp_7 = _mm_shuffle_epi8(src_temp_7, reg_shuffle);

                reg_04_16x8b = _mm_hadd_epi16(reg_04_16x8b, reg_05_16x8b);

                reg_04_16x8b = _mm_madd_epi16(reg_04_16x8b, reg_all_1s);

                reg_04_16x8b = _mm_add_epi32(reg_04_16x8b, reg_64val_32bit);

                reg_04_16x8b = _mm_srli_epi32(reg_04_16x8b, (int) 7);

                reg_02_16x8b = _mm_packus_epi32(reg_02_16x8b, reg_04_16x8b);

                reg_03_16x8b = _mm_maddubs_epi16(src_temp_4, filt_coeff_grid);
                reg_05_16x8b = _mm_maddubs_epi16(src_temp_5, filt_coeff_grid);

                reg_03_16x8b = _mm_hadd_epi16(reg_03_16x8b, reg_05_16x8b);

                reg_03_16x8b = _mm_madd_epi16(reg_03_16x8b, reg_all_1s);

                reg_03_16x8b = _mm_add_epi32(reg_03_16x8b, reg_64val_32bit);

                reg_03_16x8b = _mm_srli_epi32(reg_03_16x8b, (int) 7);

                reg_04_16x8b = _mm_maddubs_epi16(src_temp_6, filt_coeff_grid);
                reg_05_16x8b = _mm_maddubs_epi16(src_temp_7, filt_coeff_grid);

                reg_04_16x8b = _mm_hadd_epi16(reg_04_16x8b, reg_05_16x8b);

                reg_04_16x8b = _mm_madd_epi16(reg_04_16x8b, reg_all_1s);

                reg_04_16x8b = _mm_add_epi32(reg_04_16x8b, reg_64val_32bit);

                reg_04_16x8b = _mm_srli_epi32(reg_04_16x8b, (int) 7);

                reg_03_16x8b = _mm_packus_epi32(reg_03_16x8b, reg_04_16x8b);

                reg_02_16x8b = _mm_packus_epi16(reg_02_16x8b, reg_03_16x8b);

                reg_02_16x8b = _mm_shuffle_epi8(reg_02_16x8b, reg_shuffle);

                reg_03_16x8b = _mm_unpacklo_epi64(reg_01_16x8b, reg_02_16x8b);

                reg_04_16x8b = _mm_unpackhi_epi64(reg_01_16x8b, reg_02_16x8b);

                /*Storing after shuffling again*/

                _mm_storeu_si128((__m128i *) pu1_out_pixel, reg_03_16x8b);
                _mm_storeu_si128((__m128i *) (pu1_out_pixel + u4_out_stride), reg_04_16x8b);

                pu1_out_pixel += 16;

                pu1_in_pixel += (u4_src_vert_increments * (u4_in_stride << 4)) >> DOWNSCALER_Q;

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }

        /*if height is not a multiple of 8 process 2 rows at a
        time for the remaining rows*/
        if(u4_rem_vert_loop)
        {
            pu1_src_j = pu1_src + ((j << 4) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 4);

            u4_center_pixel_pos = u4_center_pixel_pos_src;
            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                UWORD8 u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_for_grid = pai1_filters[u1_phase];

                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                filt_coeff_grid = _mm_loadu_si128((__m128i *) pi1_filter_for_grid);

                for(j = u4_rem_vert_loop; j > 0; j = j - 2)
                {
                    src_temp_0 = _mm_loadu_si128((__m128i const *) pu1_in_pixel);
                    src_temp_0 = _mm_shuffle_epi8(src_temp_0, reg_shuffle);

                    src_temp_1 = _mm_loadu_si128((__m128i const *) (pu1_in_pixel + u4_in_stride));

                    src_temp_1 = _mm_shuffle_epi8(src_temp_1, reg_shuffle);

                    src_temp_0 = _mm_maddubs_epi16(src_temp_0, filt_coeff_grid);
                    src_temp_1 = _mm_maddubs_epi16(src_temp_1, filt_coeff_grid);

                    reg_01_16x8b = _mm_hadd_epi16(src_temp_0, src_temp_1);

                    reg_01_16x8b = _mm_madd_epi16(reg_01_16x8b, reg_all_1s);

                    /*Add offset of 64 for rounding each out pixel value*/
                    reg_01_16x8b = _mm_add_epi32(reg_01_16x8b, reg_64val_32bit);
                    /*Divide by 128 each out pixel value*/
                    reg_01_16x8b = _mm_srli_epi32(reg_01_16x8b, (int) 7);

                    reg_01_16x8b = _mm_packus_epi32(reg_01_16x8b, reg_all_0s);

                    /*next get saturated 8 bit output pixel values*/
                    reg_01_16x8b = _mm_packus_epi16(reg_01_16x8b, reg_all_0s);

                    reg_01_16x8b = _mm_shuffle_epi8(reg_01_16x8b, reg_shuffle);

                    reg_02_16x8b = _mm_srli_si128(reg_01_16x8b, (int) 8);

                    /*Store the 2 output values*/
                    i4_temp_pixel_holder = _mm_cvtsi128_si32(reg_01_16x8b);

                    *pu1_out_pixel = (UWORD8) i4_temp_pixel_holder;
                    i4_temp_pixel_holder >>= 8;

                    *(pu1_out_pixel + 1) = (UWORD8) i4_temp_pixel_holder;

                    i4_temp_pixel_holder = _mm_cvtsi128_si32(reg_02_16x8b);

                    *(pu1_out_pixel + u4_out_stride) = (UWORD8) i4_temp_pixel_holder;
                    i4_temp_pixel_holder >>= 8;

                    *(pu1_out_pixel + u4_out_stride + 1) = (UWORD8) i4_temp_pixel_holder;

                    pu1_in_pixel += (u4_src_vert_increments * (u4_in_stride << 1)) >> DOWNSCALER_Q;
                    pu1_out_pixel += 2;
                }
                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }
    }
}
