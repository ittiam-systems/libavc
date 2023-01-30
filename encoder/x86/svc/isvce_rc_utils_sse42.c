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
* @file isvce_rc_utils_sse42.c
*
* @brief
*  This file contains the x86 SIMD version of the function which computes
*  gradient per pixel value being used in Init Qp
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - isvce_get_gpp_sse42()
*
* @remarks
*  None
*
*******************************************************************************
*/

#include <immintrin.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "isvc_structs.h"
#include "isvce_rc_utils_private_defs.h"

/**
*******************************************************************************
*
* @brief
*   get gpp function
*
* @par Description:
*   computes gradient per pixel value for a given frame
*
* @param[in] ps_input_buf
*  pointer to yuv buffer properties
*
* @returns
*  calculated gpp value
*
* @remarks
*  none
*
*******************************************************************************
*/

DOUBLE isvce_get_gpp_sse42(yuv_buf_props_t *ps_input_buf)
{
    UWORD8 *pu1_input_buf;
    UWORD16 mask_ffff, mask_00ff;
    UWORD32 i, j, k;
    UWORD32 u4_width, u4_height, i4_input_stride;
    DOUBLE d_gpp_y, d_gpp_u, d_gpp_v, d_gpp;

    __m128i u1_src_r0, u1_src_r1, u1_src_r2, u1_src_r3, u1_src_r4;
    __m128i u1_src_right_r0, u1_src_right_r1, u1_src_right_r2, u1_src_right_r3;
    __m128i u2_sad_cur_bot_r01, u2_sad_cur_bot_r12, u2_sad_cur_bot_r23, u2_sad_cur_bot_r34;
    __m128i u2_sad_cur_right_r0, u2_sad_cur_right_r1, u2_sad_cur_right_r2, u2_sad_cur_right_r3;
    __m128i u2_sad_hadd, u1_shuffle_chroma, u2_mask_and_pixY, u2_mask_and_pixUV;

    d_gpp_y = 0;
    d_gpp_u = 0;
    d_gpp_v = 0;
    d_gpp = 0;
    mask_ffff = 0xffff;
    mask_00ff = 0x00ff;
    pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[0].pv_data;
    i4_input_stride = ps_input_buf->as_component_bufs[0].i4_data_stride;
    u4_width = ps_input_buf->u4_width;
    u4_height = ps_input_buf->u4_height;

    u1_shuffle_chroma = _mm_setr_epi8(0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x01, 0x03,
                                      0x05, 0x07, 0x09, 0x0b, 0x0d, 0x0f);
    u2_mask_and_pixY = _mm_setr_epi16(mask_ffff, mask_ffff, mask_ffff, mask_ffff, mask_ffff,
                                      mask_ffff, mask_ffff, mask_00ff);
    u2_mask_and_pixUV = _mm_setr_epi16(mask_ffff, mask_ffff, mask_ffff, mask_00ff, mask_ffff,
                                       mask_ffff, mask_ffff, mask_00ff);

    ASSERT((u4_width % 16) == 0);

    /***********************************************************/
    /* For Luma -                                              */
    /* This code block calculates gpp value for luma by adding */
    /* the absolute difference between the current pixel and   */
    /* it's immediate right pixel with the absolute difference */
    /* between the current pixel and it's immediate bottom     */
    /* pixel and accumulating for every pixel in the frame.    */
    /***********************************************************/
    for(i = 0; i < u4_height - 4; i += 4)
    {
        for(j = 0; j < u4_width - 16; j += 16)
        {
            u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
            u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
            u1_src_r2 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j));
            u1_src_r3 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j));
            u1_src_r4 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 4) + j));
            u1_src_right_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j + 1));
            u1_src_right_r1 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j + 1));
            u1_src_right_r2 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j + 1));
            u1_src_right_r3 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j + 1));

            u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
            u2_sad_cur_bot_r12 = _mm_sad_epu8(u1_src_r1, u1_src_r2);
            u2_sad_cur_bot_r23 = _mm_sad_epu8(u1_src_r2, u1_src_r3);
            u2_sad_cur_bot_r34 = _mm_sad_epu8(u1_src_r3, u1_src_r4);
            u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);
            u2_sad_cur_right_r1 = _mm_sad_epu8(u1_src_r1, u1_src_right_r1);
            u2_sad_cur_right_r2 = _mm_sad_epu8(u1_src_r2, u1_src_right_r2);
            u2_sad_cur_right_r3 = _mm_sad_epu8(u1_src_r3, u1_src_right_r3);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r12);
            u2_sad_cur_bot_r23 = _mm_adds_epu16(u2_sad_cur_bot_r23, u2_sad_cur_bot_r34);
            u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r1);
            u2_sad_cur_right_r2 = _mm_adds_epu16(u2_sad_cur_right_r2, u2_sad_cur_right_r3);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r23);
            u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r2);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

            u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

            d_gpp_y += _mm_extract_epi16(u2_sad_hadd, 0);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 15  */
        /* pixels are getting processed separately by performing    */
        /* and operations with u2_mask_and_pixY mask                */
        /************************************************************/
        u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
        u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
        u1_src_r2 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j));
        u1_src_r3 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j));
        u1_src_r4 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 4) + j));
        u1_src_right_r0 = _mm_srli_si128(u1_src_r0, 1);
        u1_src_right_r1 = _mm_srli_si128(u1_src_r1, 1);
        u1_src_right_r2 = _mm_srli_si128(u1_src_r2, 1);
        u1_src_right_r3 = _mm_srli_si128(u1_src_r3, 1);

        u1_src_r0 = _mm_and_si128(u1_src_r0, u2_mask_and_pixY);
        u1_src_r1 = _mm_and_si128(u1_src_r1, u2_mask_and_pixY);
        u1_src_r2 = _mm_and_si128(u1_src_r2, u2_mask_and_pixY);
        u1_src_r3 = _mm_and_si128(u1_src_r3, u2_mask_and_pixY);
        u1_src_r4 = _mm_and_si128(u1_src_r4, u2_mask_and_pixY);

        u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
        u2_sad_cur_bot_r12 = _mm_sad_epu8(u1_src_r1, u1_src_r2);
        u2_sad_cur_bot_r23 = _mm_sad_epu8(u1_src_r2, u1_src_r3);
        u2_sad_cur_bot_r34 = _mm_sad_epu8(u1_src_r3, u1_src_r4);
        u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);
        u2_sad_cur_right_r1 = _mm_sad_epu8(u1_src_r1, u1_src_right_r1);
        u2_sad_cur_right_r2 = _mm_sad_epu8(u1_src_r2, u1_src_right_r2);
        u2_sad_cur_right_r3 = _mm_sad_epu8(u1_src_r3, u1_src_right_r3);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r12);
        u2_sad_cur_bot_r23 = _mm_adds_epu16(u2_sad_cur_bot_r23, u2_sad_cur_bot_r34);
        u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r1);
        u2_sad_cur_right_r2 = _mm_adds_epu16(u2_sad_cur_right_r2, u2_sad_cur_right_r3);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r23);
        u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r2);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

        u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

        d_gpp_y += _mm_extract_epi16(u2_sad_hadd, 0);

        pu1_input_buf += (i4_input_stride << 2);
    }

    /* Loop for the remaining height */
    for(k = i; k < u4_height - 1; k++)
    {
        for(j = 0; j < u4_width - 16; j += 16)
        {
            u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
            u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
            u1_src_right_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j + 1));

            u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
            u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

            u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

            d_gpp_y += _mm_extract_epi16(u2_sad_hadd, 0);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 15  */
        /* pixels are getting processed separately by performing    */
        /* and operations with u2_mask_and_pixY mask                */
        /************************************************************/
        u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
        u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
        u1_src_right_r0 = _mm_srli_si128(u1_src_r0, 1);

        u1_src_r0 = _mm_and_si128(u1_src_r0, u2_mask_and_pixY);
        u1_src_r1 = _mm_and_si128(u1_src_r1, u2_mask_and_pixY);

        u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
        u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

        u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

        d_gpp_y += _mm_extract_epi16(u2_sad_hadd, 0);

        pu1_input_buf += (i4_input_stride);
    }

    pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[1].pv_data;
    i4_input_stride = ps_input_buf->as_component_bufs[1].i4_data_stride;

    /**************************************************************/
    /* For Chroma -                                               */
    /* This code block first deinterleaves the Cb and Cr values   */
    /* from the loaded registers, calculates gpp value for both   */
    /* Cb and Cr separately by adding the absolute difference     */
    /* between the current pixel and it's immediate right pixel   */
    /* with the absolute difference between the current pixel and */
    /* it's immediate bottom pixel and accumulating for every     */
    /* pixel in the frame.                                        */
    /**************************************************************/
    for(i = 0; i < (u4_height / 2) - 4; i += 4)
    {
        for(j = 0; j < u4_width - 16; j += 16)
        {
            u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
            u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
            u1_src_r2 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j));
            u1_src_r3 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j));
            u1_src_r4 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 4) + j));
            u1_src_right_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j + 2));
            u1_src_right_r1 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j + 2));
            u1_src_right_r2 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j + 2));
            u1_src_right_r3 =
                _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j + 2));

            /* separating u and v */
            u1_src_r0 = _mm_shuffle_epi8(u1_src_r0, u1_shuffle_chroma);
            u1_src_r1 = _mm_shuffle_epi8(u1_src_r1, u1_shuffle_chroma);
            u1_src_r2 = _mm_shuffle_epi8(u1_src_r2, u1_shuffle_chroma);
            u1_src_r3 = _mm_shuffle_epi8(u1_src_r3, u1_shuffle_chroma);
            u1_src_r4 = _mm_shuffle_epi8(u1_src_r4, u1_shuffle_chroma);
            u1_src_right_r0 = _mm_shuffle_epi8(u1_src_right_r0, u1_shuffle_chroma);
            u1_src_right_r1 = _mm_shuffle_epi8(u1_src_right_r1, u1_shuffle_chroma);
            u1_src_right_r2 = _mm_shuffle_epi8(u1_src_right_r2, u1_shuffle_chroma);
            u1_src_right_r3 = _mm_shuffle_epi8(u1_src_right_r3, u1_shuffle_chroma);

            u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
            u2_sad_cur_bot_r12 = _mm_sad_epu8(u1_src_r1, u1_src_r2);
            u2_sad_cur_bot_r23 = _mm_sad_epu8(u1_src_r2, u1_src_r3);
            u2_sad_cur_bot_r34 = _mm_sad_epu8(u1_src_r3, u1_src_r4);
            u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);
            u2_sad_cur_right_r1 = _mm_sad_epu8(u1_src_r1, u1_src_right_r1);
            u2_sad_cur_right_r2 = _mm_sad_epu8(u1_src_r2, u1_src_right_r2);
            u2_sad_cur_right_r3 = _mm_sad_epu8(u1_src_r3, u1_src_right_r3);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r12);
            u2_sad_cur_bot_r23 = _mm_adds_epu16(u2_sad_cur_bot_r23, u2_sad_cur_bot_r34);
            u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r1);
            u2_sad_cur_right_r2 = _mm_adds_epu16(u2_sad_cur_right_r2, u2_sad_cur_right_r3);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r23);
            u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r2);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

            u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

            d_gpp_u += _mm_extract_epi16(u2_sad_hadd, 0);
            d_gpp_v += _mm_extract_epi16(u2_sad_hadd, 1);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 15  */
        /* pixels are getting processed separately by performing    */
        /* and operations with u2_mask_and_pixUV mask               */
        /************************************************************/
        u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
        u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
        u1_src_r2 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 2) + j));
        u1_src_r3 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 3) + j));
        u1_src_r4 = _mm_loadu_si128((__m128i *) (pu1_input_buf + (i4_input_stride * 4) + j));
        u1_src_right_r0 = _mm_srli_si128(u1_src_r0, 2);
        u1_src_right_r1 = _mm_srli_si128(u1_src_r1, 2);
        u1_src_right_r2 = _mm_srli_si128(u1_src_r2, 2);
        u1_src_right_r3 = _mm_srli_si128(u1_src_r3, 2);

        /* separating u and v */
        u1_src_r0 = _mm_shuffle_epi8(u1_src_r0, u1_shuffle_chroma);
        u1_src_r1 = _mm_shuffle_epi8(u1_src_r1, u1_shuffle_chroma);
        u1_src_r2 = _mm_shuffle_epi8(u1_src_r2, u1_shuffle_chroma);
        u1_src_r3 = _mm_shuffle_epi8(u1_src_r3, u1_shuffle_chroma);
        u1_src_r4 = _mm_shuffle_epi8(u1_src_r4, u1_shuffle_chroma);
        u1_src_right_r0 = _mm_shuffle_epi8(u1_src_right_r0, u1_shuffle_chroma);
        u1_src_right_r1 = _mm_shuffle_epi8(u1_src_right_r1, u1_shuffle_chroma);
        u1_src_right_r2 = _mm_shuffle_epi8(u1_src_right_r2, u1_shuffle_chroma);
        u1_src_right_r3 = _mm_shuffle_epi8(u1_src_right_r3, u1_shuffle_chroma);

        u1_src_r0 = _mm_and_si128(u1_src_r0, u2_mask_and_pixUV);
        u1_src_r1 = _mm_and_si128(u1_src_r1, u2_mask_and_pixUV);
        u1_src_r2 = _mm_and_si128(u1_src_r2, u2_mask_and_pixUV);
        u1_src_r3 = _mm_and_si128(u1_src_r3, u2_mask_and_pixUV);
        u1_src_r4 = _mm_and_si128(u1_src_r4, u2_mask_and_pixUV);
        u1_src_right_r0 = _mm_and_si128(u1_src_right_r0, u2_mask_and_pixUV);
        u1_src_right_r1 = _mm_and_si128(u1_src_right_r1, u2_mask_and_pixUV);
        u1_src_right_r2 = _mm_and_si128(u1_src_right_r2, u2_mask_and_pixUV);
        u1_src_right_r3 = _mm_and_si128(u1_src_right_r3, u2_mask_and_pixUV);

        u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
        u2_sad_cur_bot_r12 = _mm_sad_epu8(u1_src_r1, u1_src_r2);
        u2_sad_cur_bot_r23 = _mm_sad_epu8(u1_src_r2, u1_src_r3);
        u2_sad_cur_bot_r34 = _mm_sad_epu8(u1_src_r3, u1_src_r4);
        u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);
        u2_sad_cur_right_r1 = _mm_sad_epu8(u1_src_r1, u1_src_right_r1);
        u2_sad_cur_right_r2 = _mm_sad_epu8(u1_src_r2, u1_src_right_r2);
        u2_sad_cur_right_r3 = _mm_sad_epu8(u1_src_r3, u1_src_right_r3);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r12);
        u2_sad_cur_bot_r23 = _mm_adds_epu16(u2_sad_cur_bot_r23, u2_sad_cur_bot_r34);
        u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r1);
        u2_sad_cur_right_r2 = _mm_adds_epu16(u2_sad_cur_right_r2, u2_sad_cur_right_r3);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r23);
        u2_sad_cur_right_r0 = _mm_adds_epu16(u2_sad_cur_right_r0, u2_sad_cur_right_r2);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

        u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

        d_gpp_u += _mm_extract_epi16(u2_sad_hadd, 0);
        d_gpp_v += _mm_extract_epi16(u2_sad_hadd, 1);

        pu1_input_buf += (i4_input_stride * 4);
    }

    /* Loop for the remaining height */
    for(k = i; k < (u4_height / 2) - 1; k++)
    {
        for(j = 0; j < u4_width - 16; j += 16)
        {
            u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
            u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
            u1_src_right_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j + 2));

            /* separating u and v */
            u1_src_r0 = _mm_shuffle_epi8(u1_src_r0, u1_shuffle_chroma);
            u1_src_r1 = _mm_shuffle_epi8(u1_src_r1, u1_shuffle_chroma);
            u1_src_right_r0 = _mm_shuffle_epi8(u1_src_right_r0, u1_shuffle_chroma);

            u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
            u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);

            u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

            u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
            u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

            d_gpp_u += _mm_extract_epi16(u2_sad_hadd, 0);
            d_gpp_v += _mm_extract_epi16(u2_sad_hadd, 1);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 15  */
        /* pixels are getting processed separately by performing    */
        /* and operations with u2_mask_and_pixUV mask               */
        /************************************************************/
        u1_src_r0 = _mm_loadu_si128((__m128i *) (pu1_input_buf + j));
        u1_src_r1 = _mm_loadu_si128((__m128i *) (pu1_input_buf + i4_input_stride + j));
        u1_src_right_r0 = _mm_srli_si128(u1_src_r0, 2);

        /* separating u and v */
        u1_src_r0 = _mm_shuffle_epi8(u1_src_r0, u1_shuffle_chroma);
        u1_src_r1 = _mm_shuffle_epi8(u1_src_r1, u1_shuffle_chroma);
        u1_src_right_r0 = _mm_shuffle_epi8(u1_src_right_r0, u1_shuffle_chroma);

        u1_src_r0 = _mm_and_si128(u1_src_r0, u2_mask_and_pixUV);
        u1_src_r1 = _mm_and_si128(u1_src_r1, u2_mask_and_pixUV);
        u1_src_right_r0 = _mm_and_si128(u1_src_right_r0, u2_mask_and_pixUV);

        u2_sad_cur_bot_r01 = _mm_sad_epu8(u1_src_r0, u1_src_r1);
        u2_sad_cur_right_r0 = _mm_sad_epu8(u1_src_r0, u1_src_right_r0);

        u2_sad_cur_bot_r01 = _mm_adds_epu16(u2_sad_cur_bot_r01, u2_sad_cur_right_r0);

        u2_sad_hadd = _mm_hadd_epi16(u2_sad_cur_bot_r01, u2_sad_cur_bot_r01);
        u2_sad_hadd = _mm_hadd_epi16(u2_sad_hadd, u2_sad_hadd);

        d_gpp_u += _mm_extract_epi16(u2_sad_hadd, 0);
        d_gpp_v += _mm_extract_epi16(u2_sad_hadd, 1);

        pu1_input_buf += i4_input_stride;
    }

    d_gpp_y /= (u4_width * u4_height);
    d_gpp_u /= ((u4_width / 2) * (u4_height / 2));
    d_gpp_v /= ((u4_width / 2) * (u4_height / 2));

    d_gpp = (DOUBLE) ((WT_LUMA_GPP * d_gpp_y) + d_gpp_u + d_gpp_v) / WT_TOTAL_GPP;

    return d_gpp;
}
