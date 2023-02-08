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
 *  isvcd_residual_resamp_sse42.c
 *
 * @brief
 *  Contains function definitions for intra resampling functions
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_interpolate_residual_sse42
 *  - isvcd_residual_luma_dyadic_sse42
 *  - isvcd_residual_reflayer_const_non_boundary_mb_sse42
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>
/* User include files */
#include "ih264_typedefs.h"
#include "isvcd_structs.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_dyadic_sse42                          */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore         creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_residual_luma_dyadic_sse42(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                      WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                      WORD32 i4_out_res_stride, mem_element_t *ps_ref_mb_mode,
                                      UWORD16 u2_mb_x, UWORD16 u2_mb_y, WORD32 i4_ref_nnz,
                                      WORD32 i4_ref_tx_size)

{
    WORD16 *pi2_refarray_buffer;
    WORD32 i4_blk_ctr;
    residual_sampling_ctxt_t *ps_ctxt;

    UNUSED(ps_ref_mb_mode);
    UNUSED(u2_mb_x);
    UNUSED(u2_mb_y);

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    pi2_refarray_buffer = ps_ctxt->pi2_refarray_buffer;

    /* based on transform size the counter and interpolation width and */
    /* height are intialised as follows                                */

    if((i4_ref_tx_size) && (0 != i4_ref_nnz))
    {
        WORD16 *pi2_ref_data_byte;
        WORD32 i4_i, i4_j;
        WORD16 *pi2_refarray_buffer_tmp = pi2_refarray_buffer;

        __m128i i2_coeff_8x16b_r1_0, i2_coeff_8x16b_r1_1;
        __m128i res_8x16b_r1_0, res_8x16b_r1_1;
        __m128i final_res_8x16b_r1_0, final_res_8x16b_r1_1;

        __m128i coeff_add_8x16b_r1;

        __m128i coeff_add_8x16b_r2;
        __m128i i2_coeff_8x16b_r2_0, i2_coeff_8x16b_r2_1;
        __m128i res_8x16b_r2_0, res_8x16b_r2_1;
        __m128i final_res_8x16b_r2_0, final_res_8x16b_r2_1;

        pi2_ref_data_byte = pi2_inp_data;

        /* ----------- Horizontal Interpolation ---------------- */
        for(i4_i = 0; i4_i < BLOCK_HEIGHT; i4_i += 2)
        {
            i2_coeff_8x16b_r1_0 =
                _mm_loadu_si128((__m128i *) pi2_ref_data_byte);         // a0 a1 a2 a3 a4 a5 a6 a7
            i2_coeff_8x16b_r2_0 = _mm_loadu_si128(
                (__m128i *) (pi2_ref_data_byte + i4_inp_data_stride));  // b0 b1 b2 b3 b4 b5 b6 b7

            i2_coeff_8x16b_r1_1 = _mm_srli_si128(i2_coeff_8x16b_r1_0, 2);  // a1 a2 a3 a4 a5 a6 a7 0
            i2_coeff_8x16b_r2_1 = _mm_srli_si128(i2_coeff_8x16b_r2_0, 2);  // b1 b2 b3 b4 b5 b6 b7 0

            coeff_add_8x16b_r1 = _mm_add_epi16(i2_coeff_8x16b_r1_0, i2_coeff_8x16b_r1_1);
            coeff_add_8x16b_r2 = _mm_add_epi16(i2_coeff_8x16b_r2_0, i2_coeff_8x16b_r2_1);

            i2_coeff_8x16b_r1_0 = _mm_slli_epi16(i2_coeff_8x16b_r1_0, 1);
            i2_coeff_8x16b_r2_0 = _mm_slli_epi16(i2_coeff_8x16b_r2_0, 1);

            i2_coeff_8x16b_r1_1 = _mm_slli_epi16(i2_coeff_8x16b_r1_1, 1);
            i2_coeff_8x16b_r2_1 = _mm_slli_epi16(i2_coeff_8x16b_r2_1, 1);

            res_8x16b_r1_0 = _mm_add_epi16(i2_coeff_8x16b_r1_0, coeff_add_8x16b_r1);
            res_8x16b_r2_0 = _mm_add_epi16(i2_coeff_8x16b_r2_0, coeff_add_8x16b_r2);

            res_8x16b_r1_1 = _mm_add_epi16(i2_coeff_8x16b_r1_1, coeff_add_8x16b_r1);
            res_8x16b_r2_1 = _mm_add_epi16(i2_coeff_8x16b_r2_1, coeff_add_8x16b_r2);

            final_res_8x16b_r1_0 = _mm_unpacklo_epi16(res_8x16b_r1_0, res_8x16b_r1_1);
            final_res_8x16b_r2_0 = _mm_unpacklo_epi16(res_8x16b_r2_0, res_8x16b_r2_1);

            final_res_8x16b_r1_1 = _mm_unpackhi_epi16(res_8x16b_r1_0, res_8x16b_r1_1);
            final_res_8x16b_r2_1 = _mm_unpackhi_epi16(res_8x16b_r2_0, res_8x16b_r2_1);

            _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 1), final_res_8x16b_r1_0);
            _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 9), final_res_8x16b_r1_1);

            _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 17), final_res_8x16b_r2_0);
            _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 25), final_res_8x16b_r2_1);

            pi2_refarray_buffer[0] = (pi2_ref_data_byte[0] << 2);
            pi2_refarray_buffer[15] = (pi2_ref_data_byte[7] << 2);
            pi2_ref_data_byte += i4_inp_data_stride;
            pi2_refarray_buffer[16] = (pi2_ref_data_byte[0] << 2);
            pi2_refarray_buffer[31] = (pi2_ref_data_byte[7] << 2);

            /* vertical loop uopdates */
            pi2_ref_data_byte = pi2_inp_data + ((i4_i + 2) * i4_inp_data_stride);
            pi2_refarray_buffer += 32;
        }

        /* ----------- Vertical Interpolation ---------------- */
        pi2_refarray_buffer = pi2_refarray_buffer_tmp;

        {
            __m128i i4_horz_samp_4x32b_r1_1, i4_horz_samp_4x32b_r1_2, i4_horz_samp_4x32b_r1_3,
                i4_horz_samp_4x32b_r1_4;
            __m128i i4_horz_samp_4x32b_r2_1, i4_horz_samp_4x32b_r2_2, i4_horz_samp_4x32b_r2_3,
                i4_horz_samp_4x32b_r2_4;
            __m128i i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2, i4_res_samp_4x32b_r1_3,
                i4_res_samp_4x32b_r1_4;
            __m128i i4_res_samp_4x32b_r2_1, i4_res_samp_4x32b_r2_2, i4_res_samp_4x32b_r2_3,
                i4_res_samp_4x32b_r2_4;
            __m128i horz_add_4x32b_r2_1, horz_add_4x32b_r2_2, horz_add_4x32b_r2_3,
                horz_add_4x32b_r2_4;

            __m128i i4_horz_samp_8x16b_r1_1, i4_horz_samp_8x16b_r2_1;
            __m128i i4_horz_samp_8x16b_r1_2, i4_horz_samp_8x16b_r2_2;
            __m128i i4_horz_samp_8x16b_r1_3, i4_horz_samp_8x16b_r2_3;
            __m128i i4_horz_samp_8x16b_r1_4, i4_horz_samp_8x16b_r2_4;

            __m128i twos = _mm_set1_epi32(2);
            __m128i eights = _mm_set1_epi32(8);

            WORD16 *pi2_out;

            pi2_out = pi2_out_res;

            i4_horz_samp_8x16b_r1_1 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer));
            i4_horz_samp_8x16b_r1_2 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 4));
            i4_horz_samp_8x16b_r1_3 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 8));
            i4_horz_samp_8x16b_r1_4 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 12));

            i4_horz_samp_4x32b_r1_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_1);
            i4_horz_samp_4x32b_r1_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_2);
            i4_horz_samp_4x32b_r1_3 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_3);
            i4_horz_samp_4x32b_r1_4 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_4);

            /* populate the first inter sample */
            i4_res_samp_4x32b_r1_1 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_1, twos), 2);
            i4_res_samp_4x32b_r1_2 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_2, twos), 2);
            i4_res_samp_4x32b_r1_3 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_3, twos), 2);
            i4_res_samp_4x32b_r1_4 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_4, twos), 2);

            _mm_storeu_si128((__m128i *) pi2_out,
                             _mm_packs_epi32(i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2));
            _mm_storeu_si128((__m128i *) (pi2_out + 8),
                             _mm_packs_epi32(i4_res_samp_4x32b_r1_3, i4_res_samp_4x32b_r1_4));
            pi2_out += i4_out_res_stride;

            for(i4_j = 0; i4_j < 14; i4_j += 2)
            {
                pi2_refarray_buffer += MB_WIDTH;

                i4_horz_samp_8x16b_r2_1 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer));
                i4_horz_samp_8x16b_r2_2 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 4));
                i4_horz_samp_8x16b_r2_3 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 8));
                i4_horz_samp_8x16b_r2_4 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 12));

                i4_horz_samp_4x32b_r2_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_1);
                i4_horz_samp_4x32b_r2_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_2);
                i4_horz_samp_4x32b_r2_3 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_3);
                i4_horz_samp_4x32b_r2_4 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_4);

                horz_add_4x32b_r2_1 =
                    _mm_add_epi32(i4_horz_samp_4x32b_r1_1, i4_horz_samp_4x32b_r2_1);
                horz_add_4x32b_r2_2 =
                    _mm_add_epi32(i4_horz_samp_4x32b_r1_2, i4_horz_samp_4x32b_r2_2);
                horz_add_4x32b_r2_3 =
                    _mm_add_epi32(i4_horz_samp_4x32b_r1_3, i4_horz_samp_4x32b_r2_3);
                horz_add_4x32b_r2_4 =
                    _mm_add_epi32(i4_horz_samp_4x32b_r1_4, i4_horz_samp_4x32b_r2_4);

                i4_res_samp_4x32b_r1_1 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r1_1, 1), horz_add_4x32b_r2_1);
                i4_res_samp_4x32b_r1_2 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r1_2, 1), horz_add_4x32b_r2_2);
                i4_res_samp_4x32b_r1_3 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r1_3, 1), horz_add_4x32b_r2_3);
                i4_res_samp_4x32b_r1_4 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r1_4, 1), horz_add_4x32b_r2_4);

                i4_res_samp_4x32b_r2_1 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r2_1, 1), horz_add_4x32b_r2_1);
                i4_res_samp_4x32b_r2_2 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r2_2, 1), horz_add_4x32b_r2_2);
                i4_res_samp_4x32b_r2_3 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r2_3, 1), horz_add_4x32b_r2_3);
                i4_res_samp_4x32b_r2_4 =
                    _mm_add_epi32(_mm_slli_epi32(i4_horz_samp_4x32b_r2_4, 1), horz_add_4x32b_r2_4);

                i4_res_samp_4x32b_r1_1 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_1, eights), 4);
                i4_res_samp_4x32b_r1_2 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_2, eights), 4);
                i4_res_samp_4x32b_r1_3 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_3, eights), 4);
                i4_res_samp_4x32b_r1_4 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_4, eights), 4);

                i4_res_samp_4x32b_r2_1 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_1, eights), 4);
                i4_res_samp_4x32b_r2_2 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_2, eights), 4);
                i4_res_samp_4x32b_r2_3 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_3, eights), 4);
                i4_res_samp_4x32b_r2_4 =
                    _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_4, eights), 4);

                /* populate 2 samples based on current coeffs */
                _mm_storeu_si128((__m128i *) pi2_out,
                                 _mm_packs_epi32(i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2));
                _mm_storeu_si128((__m128i *) (pi2_out + 8),
                                 _mm_packs_epi32(i4_res_samp_4x32b_r1_3, i4_res_samp_4x32b_r1_4));
                pi2_out += i4_out_res_stride;

                _mm_storeu_si128((__m128i *) pi2_out,
                                 _mm_packs_epi32(i4_res_samp_4x32b_r2_1, i4_res_samp_4x32b_r2_2));
                _mm_storeu_si128((__m128i *) (pi2_out + 8),
                                 _mm_packs_epi32(i4_res_samp_4x32b_r2_3, i4_res_samp_4x32b_r2_4));
                pi2_out += i4_out_res_stride;

                /* store the coeff 2 to coeff 1 */
                /* (used in next iteration)     */
                i4_horz_samp_4x32b_r1_1 = i4_horz_samp_4x32b_r2_1;
                i4_horz_samp_4x32b_r1_2 = i4_horz_samp_4x32b_r2_2;
                i4_horz_samp_4x32b_r1_3 = i4_horz_samp_4x32b_r2_3;
                i4_horz_samp_4x32b_r1_4 = i4_horz_samp_4x32b_r2_4;
            }

            i4_res_samp_4x32b_r1_1 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_1, twos), 2);
            i4_res_samp_4x32b_r1_2 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_2, twos), 2);
            i4_res_samp_4x32b_r1_3 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_3, twos), 2);
            i4_res_samp_4x32b_r1_4 =
                _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r1_4, twos), 2);

            _mm_storeu_si128((__m128i *) pi2_out,
                             _mm_packs_epi32(i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2));
            _mm_storeu_si128((__m128i *) (pi2_out + 8),
                             _mm_packs_epi32(i4_res_samp_4x32b_r1_3, i4_res_samp_4x32b_r1_4));
        }
    }
    else
    {
        /* ----------------------------------------------------------------- */
        /* LOOP over number of blocks                                        */
        /* ----------------------------------------------------------------- */
        for(i4_blk_ctr = 0; i4_blk_ctr < 4; i4_blk_ctr++)
        {
            /* if reference layer is not coded then no processing */
            if(0 != (i4_ref_nnz & 0x1))
            {
                __m128i i2_coeff_8x16b_r1_0, i2_coeff_8x16b_r1_1;
                __m128i i2_coeff_8x16b_r2_0, i2_coeff_8x16b_r2_1;
                __m128i i2_coeff_8x16b_r3_0, i2_coeff_8x16b_r3_1;
                __m128i i2_coeff_8x16b_r4_0, i2_coeff_8x16b_r4_1;

                __m128i res_8x16b_r1_0, res_8x16b_r1_1;
                __m128i res_8x16b_r2_0, res_8x16b_r2_1;
                __m128i res_8x16b_r3_0, res_8x16b_r3_1;
                __m128i res_8x16b_r4_0, res_8x16b_r4_1;
                __m128i final_res_8x16b_r1_0;
                __m128i final_res_8x16b_r2_0;
                __m128i final_res_8x16b_r3_0;
                __m128i final_res_8x16b_r4_0;

                __m128i coeff_add_8x16b_r1;
                __m128i coeff_add_8x16b_r2;
                __m128i coeff_add_8x16b_r3;
                __m128i coeff_add_8x16b_r4;

                /* ----------- Horizontal Interpolation ---------------- */

                i2_coeff_8x16b_r1_0 =
                    _mm_loadu_si128((__m128i *) pi2_inp_data);         // a0 a1 a2 a3 a4 a5 a6 a7
                i2_coeff_8x16b_r2_0 = _mm_loadu_si128(
                    (__m128i *) (pi2_inp_data + i4_inp_data_stride));  // b0 b1 b2 b3 b4 b5 b6 b7
                i2_coeff_8x16b_r3_0 =
                    _mm_loadu_si128((__m128i *) (pi2_inp_data + (i4_inp_data_stride << 1)));
                i2_coeff_8x16b_r4_0 =
                    _mm_loadu_si128((__m128i *) (pi2_inp_data + (i4_inp_data_stride * 3)));

                i2_coeff_8x16b_r1_1 = _mm_srli_si128(i2_coeff_8x16b_r1_0,
                                                     2);  // a1 a2 a3 a4 a5 a6 a7 0
                i2_coeff_8x16b_r2_1 = _mm_srli_si128(i2_coeff_8x16b_r2_0,
                                                     2);  // b1 b2 b3 b4 b5 b6 b7 0
                i2_coeff_8x16b_r3_1 = _mm_srli_si128(i2_coeff_8x16b_r3_0, 2);
                i2_coeff_8x16b_r4_1 = _mm_srli_si128(i2_coeff_8x16b_r4_0, 2);

                coeff_add_8x16b_r1 = _mm_add_epi16(i2_coeff_8x16b_r1_0, i2_coeff_8x16b_r1_1);
                coeff_add_8x16b_r2 = _mm_add_epi16(i2_coeff_8x16b_r2_0, i2_coeff_8x16b_r2_1);
                coeff_add_8x16b_r3 = _mm_add_epi16(i2_coeff_8x16b_r3_0, i2_coeff_8x16b_r3_1);
                coeff_add_8x16b_r4 = _mm_add_epi16(i2_coeff_8x16b_r4_0, i2_coeff_8x16b_r4_1);

                i2_coeff_8x16b_r1_0 = _mm_slli_epi16(i2_coeff_8x16b_r1_0, 1);
                i2_coeff_8x16b_r2_0 = _mm_slli_epi16(i2_coeff_8x16b_r2_0, 1);
                i2_coeff_8x16b_r3_0 = _mm_slli_epi16(i2_coeff_8x16b_r3_0, 1);
                i2_coeff_8x16b_r4_0 = _mm_slli_epi16(i2_coeff_8x16b_r4_0, 1);

                i2_coeff_8x16b_r1_1 = _mm_slli_epi16(i2_coeff_8x16b_r1_1, 1);
                i2_coeff_8x16b_r2_1 = _mm_slli_epi16(i2_coeff_8x16b_r2_1, 1);
                i2_coeff_8x16b_r3_1 = _mm_slli_epi16(i2_coeff_8x16b_r3_1, 1);
                i2_coeff_8x16b_r4_1 = _mm_slli_epi16(i2_coeff_8x16b_r4_1, 1);

                res_8x16b_r1_0 = _mm_add_epi16(i2_coeff_8x16b_r1_0, coeff_add_8x16b_r1);
                res_8x16b_r2_0 = _mm_add_epi16(i2_coeff_8x16b_r2_0, coeff_add_8x16b_r2);
                res_8x16b_r3_0 = _mm_add_epi16(i2_coeff_8x16b_r3_0, coeff_add_8x16b_r3);
                res_8x16b_r4_0 = _mm_add_epi16(i2_coeff_8x16b_r4_0, coeff_add_8x16b_r4);

                res_8x16b_r1_1 = _mm_add_epi16(i2_coeff_8x16b_r1_1, coeff_add_8x16b_r1);
                res_8x16b_r2_1 = _mm_add_epi16(i2_coeff_8x16b_r2_1, coeff_add_8x16b_r2);
                res_8x16b_r3_1 = _mm_add_epi16(i2_coeff_8x16b_r3_1, coeff_add_8x16b_r3);
                res_8x16b_r4_1 = _mm_add_epi16(i2_coeff_8x16b_r4_1, coeff_add_8x16b_r4);

                final_res_8x16b_r1_0 = _mm_unpacklo_epi16(res_8x16b_r1_0, res_8x16b_r1_1);
                final_res_8x16b_r2_0 = _mm_unpacklo_epi16(res_8x16b_r2_0, res_8x16b_r2_1);
                final_res_8x16b_r3_0 = _mm_unpacklo_epi16(res_8x16b_r3_0, res_8x16b_r3_1);
                final_res_8x16b_r4_0 = _mm_unpacklo_epi16(res_8x16b_r4_0, res_8x16b_r4_1);

                _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 1), final_res_8x16b_r1_0);
                _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 9), final_res_8x16b_r2_0);
                _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 17), final_res_8x16b_r3_0);
                _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 25), final_res_8x16b_r4_0);

                pi2_refarray_buffer[0] = (pi2_inp_data[0] << 2);
                pi2_refarray_buffer[7] = (pi2_inp_data[3] << 2);
                pi2_refarray_buffer[8] = (pi2_inp_data[i4_inp_data_stride] << 2);
                pi2_refarray_buffer[15] = (pi2_inp_data[i4_inp_data_stride + 3] << 2);
                pi2_refarray_buffer[16] = (pi2_inp_data[(i4_inp_data_stride << 1)] << 2);
                pi2_refarray_buffer[23] = (pi2_inp_data[(i4_inp_data_stride << 1) + 3] << 2);
                pi2_refarray_buffer[24] = (pi2_inp_data[(i4_inp_data_stride * 3)] << 2);
                pi2_refarray_buffer[31] = (pi2_inp_data[(i4_inp_data_stride * 3) + 3] << 2);

                /* ----------- Vertical Interpolation ---------------- */
                {
                    __m128i i4_horz_samp_8x16b_r0_1, i4_horz_samp_8x16b_r0_2;
                    __m128i i4_horz_samp_8x16b_r1_1, i4_horz_samp_8x16b_r1_2;
                    __m128i i4_horz_samp_8x16b_r2_1, i4_horz_samp_8x16b_r2_2;
                    __m128i i4_horz_samp_8x16b_r3_1, i4_horz_samp_8x16b_r3_2;

                    __m128i i4_horz_samp_4x32b_r0_1, i4_horz_samp_4x32b_r0_2;
                    __m128i i4_horz_samp_4x32b_r1_1, i4_horz_samp_4x32b_r1_2;
                    __m128i i4_horz_samp_4x32b_r2_1, i4_horz_samp_4x32b_r2_2;
                    __m128i i4_horz_samp_4x32b_r3_1, i4_horz_samp_4x32b_r3_2;

                    __m128i i4_res_samp_4x32b_r0_1, i4_res_samp_4x32b_r0_2;
                    __m128i i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2;
                    __m128i i4_res_samp_4x32b_r2_1, i4_res_samp_4x32b_r2_2;
                    __m128i i4_res_samp_4x32b_r3_1, i4_res_samp_4x32b_r3_2;
                    __m128i i4_res_samp_4x32b_r4_1, i4_res_samp_4x32b_r4_2;
                    __m128i i4_res_samp_4x32b_r5_1, i4_res_samp_4x32b_r5_2;
                    __m128i i4_res_samp_4x32b_r6_1, i4_res_samp_4x32b_r6_2;
                    __m128i i4_res_samp_4x32b_r7_1, i4_res_samp_4x32b_r7_2;

                    __m128i horz_add_4x32b_r1_1, horz_add_4x32b_r1_2;
                    __m128i horz_add_4x32b_r2_1, horz_add_4x32b_r2_2;
                    __m128i horz_add_4x32b_r3_1, horz_add_4x32b_r3_2;

                    __m128i twos = _mm_set1_epi32(2);
                    __m128i eights = _mm_set1_epi32(8);

                    i4_horz_samp_8x16b_r0_1 = _mm_loadu_si128((__m128i *) (pi2_refarray_buffer));
                    i4_horz_samp_8x16b_r0_2 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 4));
                    i4_horz_samp_8x16b_r1_1 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + BLOCK_WIDTH));
                    i4_horz_samp_8x16b_r1_2 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + BLOCK_WIDTH + 4));
                    i4_horz_samp_8x16b_r2_1 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLOCK_WIDTH << 1)));
                    i4_horz_samp_8x16b_r2_2 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLOCK_WIDTH << 1) + 4));
                    i4_horz_samp_8x16b_r3_1 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLOCK_WIDTH * 3)));
                    i4_horz_samp_8x16b_r3_2 =
                        _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLOCK_WIDTH * 3) + 4));

                    i4_horz_samp_4x32b_r0_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r0_1);
                    i4_horz_samp_4x32b_r0_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r0_2);
                    i4_horz_samp_4x32b_r1_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_1);
                    i4_horz_samp_4x32b_r1_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r1_2);
                    i4_horz_samp_4x32b_r2_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_1);
                    i4_horz_samp_4x32b_r2_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r2_2);
                    i4_horz_samp_4x32b_r3_1 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r3_1);
                    i4_horz_samp_4x32b_r3_2 = _mm_cvtepi16_epi32(i4_horz_samp_8x16b_r3_2);

                    horz_add_4x32b_r1_1 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r0_1, i4_horz_samp_4x32b_r1_1);
                    horz_add_4x32b_r2_1 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r1_1, i4_horz_samp_4x32b_r2_1);
                    horz_add_4x32b_r3_1 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r2_1, i4_horz_samp_4x32b_r3_1);

                    horz_add_4x32b_r1_2 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r0_2, i4_horz_samp_4x32b_r1_2);
                    horz_add_4x32b_r2_2 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r1_2, i4_horz_samp_4x32b_r2_2);
                    horz_add_4x32b_r3_2 =
                        _mm_add_epi32(i4_horz_samp_4x32b_r2_2, i4_horz_samp_4x32b_r3_2);

                    i4_res_samp_4x32b_r1_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r0_1, 1), horz_add_4x32b_r1_1);
                    i4_res_samp_4x32b_r2_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r1_1, 1), horz_add_4x32b_r1_1);
                    i4_res_samp_4x32b_r3_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r1_1, 1), horz_add_4x32b_r2_1);
                    i4_res_samp_4x32b_r4_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r2_1, 1), horz_add_4x32b_r2_1);
                    i4_res_samp_4x32b_r5_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r2_1, 1), horz_add_4x32b_r3_1);
                    i4_res_samp_4x32b_r6_1 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r3_1, 1), horz_add_4x32b_r3_1);

                    i4_res_samp_4x32b_r1_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r0_2, 1), horz_add_4x32b_r1_2);
                    i4_res_samp_4x32b_r2_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r1_2, 1), horz_add_4x32b_r1_2);
                    i4_res_samp_4x32b_r3_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r1_2, 1), horz_add_4x32b_r2_2);
                    i4_res_samp_4x32b_r4_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r2_2, 1), horz_add_4x32b_r2_2);
                    i4_res_samp_4x32b_r5_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r2_2, 1), horz_add_4x32b_r3_2);
                    i4_res_samp_4x32b_r6_2 = _mm_add_epi32(
                        _mm_slli_epi32(i4_horz_samp_4x32b_r3_2, 1), horz_add_4x32b_r3_2);

                    i4_res_samp_4x32b_r0_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r0_1, twos), 2);
                    i4_res_samp_4x32b_r1_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_1, eights), 4);
                    i4_res_samp_4x32b_r2_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_1, eights), 4);
                    i4_res_samp_4x32b_r3_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r3_1, eights), 4);
                    i4_res_samp_4x32b_r4_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r4_1, eights), 4);
                    i4_res_samp_4x32b_r5_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r5_1, eights), 4);
                    i4_res_samp_4x32b_r6_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r6_1, eights), 4);
                    i4_res_samp_4x32b_r7_1 =
                        _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r3_1, twos), 2);

                    i4_res_samp_4x32b_r0_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r0_2, twos), 2);
                    i4_res_samp_4x32b_r1_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r1_2, eights), 4);
                    i4_res_samp_4x32b_r2_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r2_2, eights), 4);
                    i4_res_samp_4x32b_r3_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r3_2, eights), 4);
                    i4_res_samp_4x32b_r4_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r4_2, eights), 4);
                    i4_res_samp_4x32b_r5_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r5_2, eights), 4);
                    i4_res_samp_4x32b_r6_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_res_samp_4x32b_r6_2, eights), 4);
                    i4_res_samp_4x32b_r7_2 =
                        _mm_srai_epi32(_mm_add_epi32(i4_horz_samp_4x32b_r3_2, twos), 2);

                    /* populate 2 samples based on current coeffs */
                    _mm_storeu_si128(
                        (__m128i *) pi2_out_res,
                        _mm_packs_epi32(i4_res_samp_4x32b_r0_1, i4_res_samp_4x32b_r0_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + i4_out_res_stride),
                        _mm_packs_epi32(i4_res_samp_4x32b_r1_1, i4_res_samp_4x32b_r1_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride << 1)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r2_1, i4_res_samp_4x32b_r2_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride * 3)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r3_1, i4_res_samp_4x32b_r3_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride << 2)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r4_1, i4_res_samp_4x32b_r4_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride * 5)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r5_1, i4_res_samp_4x32b_r5_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride * 6)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r6_1, i4_res_samp_4x32b_r6_2));
                    _mm_storeu_si128(
                        (__m128i *) (pi2_out_res + (i4_out_res_stride * 7)),
                        _mm_packs_epi32(i4_res_samp_4x32b_r7_1, i4_res_samp_4x32b_r7_2));

                    pi2_out_res += BLOCK_WIDTH;
                }
            }
            else
            {
                pi2_out_res += BLOCK_WIDTH;
            }

            /* Block level loop updates */
            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= SUB_BLOCK_WIDTH;
                pi2_inp_data += (i4_inp_data_stride * SUB_BLOCK_HEIGHT);
                pi2_out_res -= MB_WIDTH;
                pi2_out_res += (i4_out_res_stride * BLOCK_HEIGHT);
                i4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += SUB_BLOCK_WIDTH;
            }

            i4_ref_nnz >>= 1;
        } /* end of loop over all the blocks */
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_residual_sse42                          */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore         creation                             */
/*                                                                           */
/*****************************************************************************/

void isvcd_interpolate_residual_sse42(void *pv_residual_samp_ctxt, WORD16 *pi2_out,
                                      WORD32 i4_out_stride, WORD32 i4_refarray_wd, UWORD16 u2_mb_x,
                                      UWORD16 u2_mb_y, WORD32 i4_chroma_flag)
{
    residual_sampling_ctxt_t *ps_ctxt;
    residual_samp_map_ctxt_t *ps_map_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;

    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    WORD32 i4_temp_array_ht;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD16 *pi2_ref_array;
    UWORD8 *pu1_ref_x_ptr_incr, *pu1_ref_y_ptr_incr;

    WORD8 arr_y_ref_pos[16] = {0};
    WORD8 arr_x_ref_pos[16] = {0};
    WORD8 arr_x_phase[32] = {0};
    WORD8 arr_y_phase[32] = {0};
    WORD8 *pi1_y_ref_pos;
    WORD8 *pi1_x_ref_pos;
    WORD8 *pi1_y_phase;
    WORD8 *pi1_x_phase;

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi2_ref_array = ps_ctxt->pi2_refarray_buffer;
    pu1_ref_x_ptr_incr = ps_ctxt->pu1_ref_x_ptr_incr;
    pu1_ref_y_ptr_incr = ps_ctxt->pu1_ref_y_ptr_incr;

    /* --------------------------------------------------------------------- */
    /* Extracting information from the mapping context                       */
    /* --------------------------------------------------------------------- */
    if(1 == i4_chroma_flag)
        ps_map_ctxt = &ps_lyr_ctxt->s_chroma_map_ctxt;
    else
        ps_map_ctxt = &ps_lyr_ctxt->s_luma_map_ctxt;

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

    i4_temp_array_ht = i4_mb_ht;
    i4_frm_mb_y = u2_mb_y * i4_mb_ht;
    i4_frm_mb_x = u2_mb_x * i4_mb_wd;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation                                                */
    /* --------------------------------------------------------------------- */

    if(i4_chroma_flag == 0)
    {
        __m128i const_16_8x16b, const_128, const_ones, const_ones_8x16b, mid_indx_16x8b;
        __m128i ref_arr_8x16b_r0_0;
        __m128i ref_arr_8x16b_r1_0;
        __m128i phs_mask_8x16b_0, phs_mask_16min_8x16b_0, phs_mask_16x8b_0;
        __m128i x_ref_pos_mask_r0, x_ref_rnd_mask_r0_0;
        __m128i x_ref_pos_mask_temp_r0_0;
        __m128i x_ref_pos_mask_temp_r1_0;
        __m128i phs_mask_div8_8x16b_0;
        __m128i u1_incr_8x16b_r0_0, ref_arr_temp0_8x16b_r0_0, res0_8x16b_r0_0,
            u1_incr_not_8x16b_r0_0;
        __m128i u1_incr_8x16b_r1_0, ref_arr_temp1_8x16b_r0_0, res1_8x16b_r0_0;

        __m128i u1_incr_not_8x16b_r0_even, u1_incr_not_8x16b_r1_even, x_ref_pos_mask_temp_r0_even,
            x_ref_pos_mask_temp_r1_even;
        __m128i u1_incr_not_8x16b_r0_odd, u1_incr_not_8x16b_r1_odd, x_ref_pos_mask_temp_r0_odd,
            x_ref_pos_mask_temp_r1_odd;

        __m128i ref_arr_temp0_8x16b_r1_0, res_8x16b_r0_0, res0_8x16b_r1_0, u1_incr_not_8x16b_r1_0;
        __m128i ref_arr_temp1_8x16b_r1_0, res_8x16b_r1_0, res1_8x16b_r1_0;
        __m128i u1_y_incr_8x16b_r0_0, u1_y_incr_8x16b_r0_1, u1_y_incr_8x16b_r0_low,
            u1_y_incr_8x16b_r0_high;

        __m128i prev_res_8x16b_r0_0;
        __m128i prev_res_8x16b_r1_0;
        __m128i prev_res_8x16b_r0_1;
        __m128i prev_res_8x16b_r1_1;

        __m128i u1_prev_y_incr_8x16b_r0_0;
        __m128i u1_prev_y_incr_8x16b_r0_1;

        __m128i ref_arr_8x16b_r0_1;
        __m128i ref_arr_8x16b_r1_1;
        __m128i phs_mask_8x16b_1, phs_mask_div8_8x16b_1, phs_mask_16min_8x16b_1;
        __m128i x_ref_pos_mask_temp_r0_1;
        __m128i x_ref_pos_mask_temp_r1_1;
        __m128i ref_arr_temp0_8x16b_r0_1, res0_8x16b_r0_1, u1_incr_not_8x16b_r0_1;
        __m128i ref_arr_temp1_8x16b_r0_1, res1_8x16b_r0_1;

        __m128i ref_arr_temp0_8x16b_r1_1, res_8x16b_r0_1, res0_8x16b_r1_1, u1_incr_not_8x16b_r1_1;
        __m128i ref_arr_temp1_8x16b_r1_1, res_8x16b_r1_1, res1_8x16b_r1_1;

        __m128i vert_res0_8x16b_r0_0, vert_res0_8x16b_r0_1, res_4x32b_l_0, res_4x32b_h_0;
        __m128i vert_res1_8x16b_r0_0, vert_res1_8x16b_r0_1, res_4x32b_l_1, res_4x32b_h_1;
        __m128i res_8x16b_l, res_8x16b_h;
        __m128i phs_y_mask_16min_8x16b, phs_y_mask_8x16b, phs_y_mask_mix_8x16b;
        __m128i zero_8x16b;
        WORD32 zero_r0_0, zero_r1_0, zero_r0_1, zero_r1_1, zero_r0_r1 = 0;
        WORD32 strt_indx_h;
        WORD16 *pi2_ref_array_temp;
        UWORD8 *pu1_ref_x_ptr_incr_temp, *pu1_ref_y_ptr_incr_temp;
        WORD32 i4_y_phase;
        WORD32 out_stride_temp;
        const_128 = _mm_set1_epi32(128);
        zero_8x16b = _mm_set1_epi16(0);
        const_ones = _mm_set1_epi8(1);
        const_ones_8x16b = _mm_set1_epi16(1);

        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            arr_y_phase[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            arr_y_ref_pos[i4_y] = (WORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
        }
        pi1_y_ref_pos = arr_y_ref_pos;
        pi1_y_phase = arr_y_phase;

        strt_indx_h = 0;
        strt_indx_h = (ps_x_pos_phase[8 + i4_frm_mb_x].i2_ref_pos);
        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            arr_x_ref_pos[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            arr_x_phase[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
        }

        pi1_x_ref_pos = arr_x_ref_pos;
        pi1_x_phase = arr_x_phase;

        x_ref_pos_mask_r0 = _mm_loadu_si128((__m128i *) (pi1_x_ref_pos));
        phs_mask_16x8b_0 = _mm_loadu_si128((__m128i *) (pi1_x_phase));
        phs_mask_8x16b_0 = _mm_cvtepi8_epi16(phs_mask_16x8b_0);
        phs_mask_8x16b_1 = _mm_cvtepi8_epi16(_mm_loadu_si128((__m128i *) (pi1_x_phase + 8)));

        phs_mask_div8_8x16b_0 = _mm_srli_epi16(phs_mask_8x16b_0, 3);
        phs_mask_div8_8x16b_1 = _mm_srli_epi16(phs_mask_8x16b_1, 3);
        phs_mask_div8_8x16b_0 = _mm_packs_epi16(phs_mask_div8_8x16b_0, phs_mask_div8_8x16b_1);
        const_16_8x16b = _mm_set1_epi16(16);

        phs_mask_16min_8x16b_0 = _mm_sub_epi16(const_16_8x16b, phs_mask_8x16b_0);
        phs_mask_16min_8x16b_1 = _mm_sub_epi16(const_16_8x16b, phs_mask_8x16b_1);

        x_ref_rnd_mask_r0_0 = _mm_add_epi8(x_ref_pos_mask_r0, phs_mask_div8_8x16b_0);
        mid_indx_16x8b = _mm_set1_epi8((strt_indx_h << 1));
        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            if((i4_y > 0) && (pi1_y_ref_pos[i4_y] == pi1_y_ref_pos[i4_y - 1]))
            {
                if(zero_r0_r1)
                {
                    res_8x16b_l = _mm_set1_epi16(0);
                    res_8x16b_h = _mm_set1_epi16(0);
                    out_stride_temp = (i4_y * i4_out_stride);
                    _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp), res_8x16b_l);
                    _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp + 8), res_8x16b_h);
                    continue;
                }

                res_8x16b_r0_0 = prev_res_8x16b_r0_0;
                res_8x16b_r1_0 = prev_res_8x16b_r1_0;
                res_8x16b_r0_1 = prev_res_8x16b_r0_1;
                res_8x16b_r1_1 = prev_res_8x16b_r1_1;

                u1_y_incr_8x16b_r0_0 = u1_prev_y_incr_8x16b_r0_0;
                u1_y_incr_8x16b_r0_1 = u1_prev_y_incr_8x16b_r0_1;
            }
            else
            {
                pi2_ref_array_temp = pi2_ref_array + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                pu1_ref_x_ptr_incr_temp =
                    pu1_ref_x_ptr_incr + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                ref_arr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pi2_ref_array_temp));
                ref_arr_8x16b_r1_0 =
                    _mm_loadu_si128((__m128i *) (pi2_ref_array_temp + i4_refarray_wd));
                ref_arr_8x16b_r0_1 =
                    _mm_loadu_si128((__m128i *) (pi2_ref_array_temp + strt_indx_h));
                ref_arr_8x16b_r1_1 = _mm_loadu_si128(
                    (__m128i *) (pi2_ref_array_temp + i4_refarray_wd + strt_indx_h));

                zero_r0_0 = _mm_test_all_ones(_mm_cmpeq_epi16(
                    ref_arr_8x16b_r0_0, zero_8x16b));  // return 1 if all zeros, else 0
                zero_r1_0 = _mm_test_all_ones(_mm_cmpeq_epi16(ref_arr_8x16b_r1_0, zero_8x16b));
                zero_r0_1 = _mm_test_all_ones(_mm_cmpeq_epi16(ref_arr_8x16b_r0_1, zero_8x16b));
                zero_r1_1 = _mm_test_all_ones(_mm_cmpeq_epi16(ref_arr_8x16b_r1_1, zero_8x16b));

                zero_r0_r1 = zero_r0_0 && zero_r1_0 && zero_r0_1 && zero_r1_1;

                if(!zero_r0_r1)
                {
                    u1_incr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pu1_ref_x_ptr_incr_temp));
                    u1_incr_8x16b_r1_0 =
                        _mm_loadu_si128((__m128i *) (pu1_ref_x_ptr_incr_temp + i4_refarray_wd));

                    u1_incr_8x16b_r0_0 = _mm_shuffle_epi8(u1_incr_8x16b_r0_0, x_ref_pos_mask_r0);
                    u1_incr_8x16b_r1_0 = _mm_shuffle_epi8(u1_incr_8x16b_r1_0, x_ref_pos_mask_r0);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_andnot_si128(u1_incr_8x16b_r0_0, phs_mask_div8_8x16b_0);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_andnot_si128(u1_incr_8x16b_r1_0, phs_mask_div8_8x16b_0);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, x_ref_pos_mask_r0);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, x_ref_pos_mask_r0);

                    x_ref_pos_mask_temp_r0_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, u1_incr_8x16b_r0_0);
                    x_ref_pos_mask_temp_r1_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, u1_incr_8x16b_r1_0);

                    /* _mm_slli_epi8(u1_incr_not_8x16b_r0_0, 1)*/
                    u1_incr_not_8x16b_r0_even =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, u1_incr_not_8x16b_r0_0);
                    u1_incr_not_8x16b_r1_even =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, u1_incr_not_8x16b_r1_0);
                    x_ref_pos_mask_temp_r0_even =
                        _mm_add_epi8(x_ref_pos_mask_temp_r0_0, x_ref_pos_mask_temp_r0_0);
                    x_ref_pos_mask_temp_r1_even =
                        _mm_add_epi8(x_ref_pos_mask_temp_r1_0, x_ref_pos_mask_temp_r1_0);

                    u1_incr_not_8x16b_r0_odd = _mm_add_epi8(u1_incr_not_8x16b_r0_even, const_ones);
                    u1_incr_not_8x16b_r1_odd = _mm_add_epi8(u1_incr_not_8x16b_r1_even, const_ones);
                    x_ref_pos_mask_temp_r0_odd =
                        _mm_add_epi8(x_ref_pos_mask_temp_r0_even, const_ones);
                    x_ref_pos_mask_temp_r1_odd =
                        _mm_add_epi8(x_ref_pos_mask_temp_r1_even, const_ones);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_unpacklo_epi8(u1_incr_not_8x16b_r0_even, u1_incr_not_8x16b_r0_odd);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_unpacklo_epi8(u1_incr_not_8x16b_r1_even, u1_incr_not_8x16b_r1_odd);
                    x_ref_pos_mask_temp_r0_0 =
                        _mm_unpacklo_epi8(x_ref_pos_mask_temp_r0_even, x_ref_pos_mask_temp_r0_odd);
                    x_ref_pos_mask_temp_r1_0 =
                        _mm_unpacklo_epi8(x_ref_pos_mask_temp_r1_even, x_ref_pos_mask_temp_r1_odd);

                    u1_incr_not_8x16b_r0_1 =
                        _mm_unpackhi_epi8(u1_incr_not_8x16b_r0_even, u1_incr_not_8x16b_r0_odd);
                    u1_incr_not_8x16b_r1_1 =
                        _mm_unpackhi_epi8(u1_incr_not_8x16b_r1_even, u1_incr_not_8x16b_r1_odd);
                    x_ref_pos_mask_temp_r0_1 =
                        _mm_unpackhi_epi8(x_ref_pos_mask_temp_r0_even, x_ref_pos_mask_temp_r0_odd);
                    x_ref_pos_mask_temp_r1_1 =
                        _mm_unpackhi_epi8(x_ref_pos_mask_temp_r1_even, x_ref_pos_mask_temp_r1_odd);

                    u1_incr_not_8x16b_r0_1 = _mm_sub_epi8(u1_incr_not_8x16b_r0_1, mid_indx_16x8b);
                    u1_incr_not_8x16b_r1_1 = _mm_sub_epi8(u1_incr_not_8x16b_r1_1, mid_indx_16x8b);
                    x_ref_pos_mask_temp_r0_1 =
                        _mm_sub_epi8(x_ref_pos_mask_temp_r0_1, mid_indx_16x8b);
                    x_ref_pos_mask_temp_r1_1 =
                        _mm_sub_epi8(x_ref_pos_mask_temp_r1_1, mid_indx_16x8b);

                    ref_arr_temp0_8x16b_r0_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_0, u1_incr_not_8x16b_r0_0);
                    ref_arr_temp0_8x16b_r1_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_0, u1_incr_not_8x16b_r1_0);
                    ref_arr_temp1_8x16b_r0_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_0, x_ref_pos_mask_temp_r0_0);
                    ref_arr_temp1_8x16b_r1_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_0, x_ref_pos_mask_temp_r1_0);
                    ref_arr_temp0_8x16b_r0_1 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_1, u1_incr_not_8x16b_r0_1);
                    ref_arr_temp0_8x16b_r1_1 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_1, u1_incr_not_8x16b_r1_1);
                    ref_arr_temp1_8x16b_r0_1 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_1, x_ref_pos_mask_temp_r0_1);
                    ref_arr_temp1_8x16b_r1_1 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_1, x_ref_pos_mask_temp_r1_1);

                    res0_8x16b_r0_0 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r0_0, phs_mask_16min_8x16b_0);
                    res0_8x16b_r1_0 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r1_0, phs_mask_16min_8x16b_0);
                    res1_8x16b_r0_0 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r0_0, phs_mask_8x16b_0);
                    res1_8x16b_r1_0 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r1_0, phs_mask_8x16b_0);
                    res0_8x16b_r0_1 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r0_1, phs_mask_16min_8x16b_1);
                    res0_8x16b_r1_1 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r1_1, phs_mask_16min_8x16b_1);
                    res1_8x16b_r0_1 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r0_1, phs_mask_8x16b_1);
                    res1_8x16b_r1_1 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r1_1, phs_mask_8x16b_1);

                    res_8x16b_r0_0 = _mm_add_epi16(res0_8x16b_r0_0, res1_8x16b_r0_0);
                    res_8x16b_r1_0 = _mm_add_epi16(res0_8x16b_r1_0, res1_8x16b_r1_0);
                    res_8x16b_r0_1 = _mm_add_epi16(res0_8x16b_r0_1, res1_8x16b_r0_1);
                    res_8x16b_r1_1 = _mm_add_epi16(res0_8x16b_r1_1, res1_8x16b_r1_1);

                    prev_res_8x16b_r0_0 = res_8x16b_r0_0;
                    prev_res_8x16b_r1_0 = res_8x16b_r1_0;
                    prev_res_8x16b_r0_1 = res_8x16b_r0_1;
                    prev_res_8x16b_r1_1 = res_8x16b_r1_1;

                    pu1_ref_y_ptr_incr_temp =
                        pu1_ref_y_ptr_incr + (pi1_y_ref_pos[i4_y] * i4_refarray_wd);
                    u1_y_incr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pu1_ref_y_ptr_incr_temp));

                    u1_y_incr_8x16b_r0_0 =
                        _mm_shuffle_epi8(u1_y_incr_8x16b_r0_0, x_ref_rnd_mask_r0_0);

                    u1_y_incr_8x16b_r0_low = _mm_cvtepi8_epi16(u1_y_incr_8x16b_r0_0);
                    u1_y_incr_8x16b_r0_high =
                        _mm_cvtepi8_epi16(_mm_unpackhi_epi64(u1_y_incr_8x16b_r0_0, const_ones));

                    u1_y_incr_8x16b_r0_0 =
                        _mm_cmpeq_epi16(u1_y_incr_8x16b_r0_low, const_ones_8x16b);
                    u1_y_incr_8x16b_r0_1 =
                        _mm_cmpeq_epi16(u1_y_incr_8x16b_r0_high, const_ones_8x16b);

                    u1_prev_y_incr_8x16b_r0_0 = u1_y_incr_8x16b_r0_0;
                    u1_prev_y_incr_8x16b_r0_1 = u1_y_incr_8x16b_r0_1;
                }
            }

            if(zero_r0_r1)
            {
                res_8x16b_l = _mm_set1_epi16(0);
                res_8x16b_h = _mm_set1_epi16(0);
            }
            else
            {
                i4_y_phase = pi1_y_phase[i4_y];

                if((i4_y_phase) >> 3)
                {
                    vert_res0_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r1_0, res_8x16b_r0_0, u1_y_incr_8x16b_r0_0);
                    vert_res1_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r1_0, res_8x16b_r1_0, u1_y_incr_8x16b_r0_0);
                    vert_res0_8x16b_r0_1 =
                        _mm_blendv_epi8(res_8x16b_r1_1, res_8x16b_r0_1, u1_y_incr_8x16b_r0_1);
                    vert_res1_8x16b_r0_1 =
                        _mm_blendv_epi8(res_8x16b_r1_1, res_8x16b_r1_1, u1_y_incr_8x16b_r0_1);
                }
                else
                {
                    vert_res0_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r0_0, res_8x16b_r0_0, u1_y_incr_8x16b_r0_0);
                    vert_res1_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r0_0, res_8x16b_r1_0, u1_y_incr_8x16b_r0_0);
                    vert_res0_8x16b_r0_1 =
                        _mm_blendv_epi8(res_8x16b_r0_1, res_8x16b_r0_1, u1_y_incr_8x16b_r0_1);
                    vert_res1_8x16b_r0_1 =
                        _mm_blendv_epi8(res_8x16b_r0_1, res_8x16b_r1_1, u1_y_incr_8x16b_r0_1);
                }
                res0_8x16b_r0_0 = _mm_unpacklo_epi16(vert_res0_8x16b_r0_0, vert_res1_8x16b_r0_0);
                res1_8x16b_r0_0 = _mm_unpackhi_epi16(vert_res0_8x16b_r0_0, vert_res1_8x16b_r0_0);
                res0_8x16b_r0_1 = _mm_unpacklo_epi16(vert_res0_8x16b_r0_1, vert_res1_8x16b_r0_1);
                res1_8x16b_r0_1 = _mm_unpackhi_epi16(vert_res0_8x16b_r0_1, vert_res1_8x16b_r0_1);

                phs_y_mask_16min_8x16b = _mm_set1_epi16(16 - i4_y_phase);
                phs_y_mask_8x16b = _mm_set1_epi16(i4_y_phase);
                phs_y_mask_mix_8x16b = _mm_unpacklo_epi16(phs_y_mask_16min_8x16b, phs_y_mask_8x16b);

                res_4x32b_l_0 = _mm_madd_epi16(res0_8x16b_r0_0, phs_y_mask_mix_8x16b);
                res_4x32b_l_1 = _mm_madd_epi16(res1_8x16b_r0_0, phs_y_mask_mix_8x16b);
                res_4x32b_h_0 = _mm_madd_epi16(res0_8x16b_r0_1, phs_y_mask_mix_8x16b);
                res_4x32b_h_1 = _mm_madd_epi16(res1_8x16b_r0_1, phs_y_mask_mix_8x16b);

                res_4x32b_l_0 = _mm_add_epi32(res_4x32b_l_0, const_128);
                res_4x32b_l_1 = _mm_add_epi32(res_4x32b_l_1, const_128);
                res_4x32b_h_0 = _mm_add_epi32(res_4x32b_h_0, const_128);
                res_4x32b_h_1 = _mm_add_epi32(res_4x32b_h_1, const_128);

                res_4x32b_l_0 = _mm_srai_epi32(res_4x32b_l_0, 8);
                res_4x32b_l_1 = _mm_srai_epi32(res_4x32b_l_1, 8);
                res_4x32b_h_0 = _mm_srai_epi32(res_4x32b_h_0, 8);
                res_4x32b_h_1 = _mm_srai_epi32(res_4x32b_h_1, 8);
                res_8x16b_l = _mm_packs_epi32(res_4x32b_l_0, res_4x32b_l_1);
                res_8x16b_h = _mm_packs_epi32(res_4x32b_h_0, res_4x32b_h_1);
            }

            out_stride_temp = (i4_y * i4_out_stride);
            _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp), res_8x16b_l);
            _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp + 8), res_8x16b_h);
        }
    }
    else
    {
        __m128i const_16_8x16b, const_128, const_ones, const_ones_8x16b;
        __m128i ref_arr_8x16b_r0_0;
        __m128i ref_arr_8x16b_r1_0;
        __m128i phs_mask_8x16b_0, phs_mask_div8_8x16b_0, phs_mask_16min_8x16b_0;
        __m128i x_ref_pos_mask_r0, x_ref_rnd_mask_r0_0;
        __m128i x_ref_pos_mask_temp_r0_0;
        __m128i x_ref_pos_mask_temp_r1_0;

        __m128i u1_incr_8x16b_r0_0, ref_arr_temp0_8x16b_r0_0, res0_8x16b_r0_0,
            u1_incr_not_8x16b_r0_0;
        __m128i u1_incr_8x16b_r1_0, ref_arr_temp1_8x16b_r0_0, res1_8x16b_r0_0;
        __m128i u1_y_incr_8x16b_r0_0;

        __m128i u1_incr_not_8x16b_r0_odd, u1_incr_not_8x16b_r1_odd, x_ref_pos_mask_temp_r0_odd,
            x_ref_pos_mask_temp_r1_odd;
        __m128i u1_incr_not_8x16b_r0_even, u1_incr_not_8x16b_r1_even, x_ref_pos_mask_temp_r0_even,
            x_ref_pos_mask_temp_r1_even;

        __m128i ref_arr_temp0_8x16b_r1_0, res_8x16b_r0_0, res0_8x16b_r1_0, u1_incr_not_8x16b_r1_0;
        __m128i ref_arr_temp1_8x16b_r1_0, res_8x16b_r1_0, res1_8x16b_r1_0;
        __m128i u1_prev_y_incr_8x16b_r0_0;
        __m128i prev_res_8x16b_r0_0;
        __m128i prev_res_8x16b_r1_0;

        __m128i vert_res0_8x16b_r0_0, res_4x32b_l_0, out_4x32b_l;
        __m128i vert_res1_8x16b_r0_0, res_4x32b_l_1, out_4x32b_h;
        __m128i phs_y_mask_16min_8x16b, phs_y_mask_8x16b, phs_y_mask_mix_8x16b;
        __m128i chroma_mask, chroma_mask2;
        __m128i zero_8x16b = _mm_set1_epi16(0);
        WORD32 zero_r0_0, zero_r1_0, zero_r0_r1 = 0;
        WORD16 *pi2_ref_array_temp;
        UWORD8 *pu1_ref_x_ptr_incr_temp, *pu1_ref_y_ptr_incr_temp;
        WORD32 i4_y_phase;
        WORD32 out_stride_temp;
        const_ones = _mm_set1_epi8(1);
        const_ones_8x16b = _mm_set1_epi16(1);
        const_128 = _mm_set1_epi32(128);

        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            arr_y_phase[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            arr_y_ref_pos[i4_y] = (WORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
        }
        pi1_y_ref_pos = arr_y_ref_pos;
        pi1_y_phase = arr_y_phase;

        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            arr_x_ref_pos[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            arr_x_phase[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
        }

        pi1_x_ref_pos = arr_x_ref_pos;
        pi1_x_phase = arr_x_phase;

        phs_mask_8x16b_0 = _mm_cvtepi8_epi16(_mm_loadu_si128((__m128i *) (pi1_x_phase)));
        x_ref_pos_mask_r0 = _mm_loadu_si128((__m128i *) (pi1_x_ref_pos));

        const_16_8x16b = _mm_set1_epi16(16);
        chroma_mask = _mm_set1_epi32(0xFFFF0000);
        chroma_mask2 = _mm_set1_epi32(0x0000FFFF);
        phs_mask_div8_8x16b_0 = _mm_srli_epi16(phs_mask_8x16b_0, 3);
        phs_mask_div8_8x16b_0 = _mm_packs_epi16(phs_mask_div8_8x16b_0, const_ones);

        phs_mask_16min_8x16b_0 = _mm_sub_epi16(const_16_8x16b, phs_mask_8x16b_0);
        x_ref_rnd_mask_r0_0 = _mm_add_epi8(x_ref_pos_mask_r0, phs_mask_div8_8x16b_0);

        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            if((i4_y > 0) && (pi1_y_ref_pos[i4_y] == pi1_y_ref_pos[i4_y - 1]))
            {
                if(zero_r0_r1)
                {
                    res_4x32b_l_0 = _mm_set1_epi32(0);
                    res_4x32b_l_1 = _mm_set1_epi32(0);
                    out_stride_temp = (i4_y * i4_out_stride);

                    out_4x32b_l = _mm_loadu_si128((__m128i *) (pi2_out + out_stride_temp));
                    out_4x32b_h = _mm_loadu_si128((__m128i *) (pi2_out + out_stride_temp + 8));

                    out_4x32b_l = _mm_and_si128(out_4x32b_l, chroma_mask);
                    out_4x32b_h = _mm_and_si128(out_4x32b_h, chroma_mask);

                    res_4x32b_l_0 = _mm_and_si128(res_4x32b_l_0, chroma_mask2);
                    res_4x32b_l_1 = _mm_and_si128(res_4x32b_l_1, chroma_mask2);

                    out_4x32b_l = _mm_add_epi8(res_4x32b_l_0, out_4x32b_l);
                    out_4x32b_h = _mm_add_epi8(res_4x32b_l_1, out_4x32b_h);

                    _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp), out_4x32b_l);
                    _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp + 8), out_4x32b_h);
                    continue;
                }

                res_8x16b_r0_0 = prev_res_8x16b_r0_0;
                res_8x16b_r1_0 = prev_res_8x16b_r1_0;

                u1_y_incr_8x16b_r0_0 = u1_prev_y_incr_8x16b_r0_0;
            }
            else
            {
                pi2_ref_array_temp = pi2_ref_array + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                pu1_ref_x_ptr_incr_temp =
                    pu1_ref_x_ptr_incr + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                ref_arr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pi2_ref_array_temp));
                ref_arr_8x16b_r1_0 =
                    _mm_loadu_si128((__m128i *) (pi2_ref_array_temp + i4_refarray_wd));

                zero_r0_0 = _mm_test_all_ones(_mm_cmpeq_epi16(
                    ref_arr_8x16b_r0_0, zero_8x16b));  // return 1 if all zeros, else 0
                zero_r1_0 = _mm_test_all_ones(_mm_cmpeq_epi16(ref_arr_8x16b_r1_0, zero_8x16b));

                zero_r0_r1 = zero_r0_0 && zero_r1_0;

                if(!zero_r0_r1)
                {
                    u1_incr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pu1_ref_x_ptr_incr_temp));
                    u1_incr_8x16b_r1_0 =
                        _mm_loadu_si128((__m128i *) (pu1_ref_x_ptr_incr_temp + i4_refarray_wd));

                    u1_incr_8x16b_r0_0 = _mm_shuffle_epi8(u1_incr_8x16b_r0_0, x_ref_pos_mask_r0);
                    u1_incr_8x16b_r1_0 = _mm_shuffle_epi8(u1_incr_8x16b_r1_0, x_ref_pos_mask_r0);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_andnot_si128(u1_incr_8x16b_r0_0, phs_mask_div8_8x16b_0);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_andnot_si128(u1_incr_8x16b_r1_0, phs_mask_div8_8x16b_0);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, x_ref_pos_mask_r0);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, x_ref_pos_mask_r0);

                    x_ref_pos_mask_temp_r0_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, u1_incr_8x16b_r0_0);
                    x_ref_pos_mask_temp_r1_0 =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, u1_incr_8x16b_r1_0);

                    u1_incr_not_8x16b_r0_even =
                        _mm_add_epi8(u1_incr_not_8x16b_r0_0, u1_incr_not_8x16b_r0_0);
                    u1_incr_not_8x16b_r1_even =
                        _mm_add_epi8(u1_incr_not_8x16b_r1_0, u1_incr_not_8x16b_r1_0);
                    x_ref_pos_mask_temp_r0_even =
                        _mm_add_epi8(x_ref_pos_mask_temp_r0_0, x_ref_pos_mask_temp_r0_0);
                    x_ref_pos_mask_temp_r1_even =
                        _mm_add_epi8(x_ref_pos_mask_temp_r1_0, x_ref_pos_mask_temp_r1_0);

                    u1_incr_not_8x16b_r0_odd = _mm_add_epi8(u1_incr_not_8x16b_r0_even, const_ones);
                    u1_incr_not_8x16b_r1_odd = _mm_add_epi8(u1_incr_not_8x16b_r1_even, const_ones);
                    x_ref_pos_mask_temp_r0_odd =
                        _mm_add_epi8(x_ref_pos_mask_temp_r0_even, const_ones);
                    x_ref_pos_mask_temp_r1_odd =
                        _mm_add_epi8(x_ref_pos_mask_temp_r1_even, const_ones);

                    u1_incr_not_8x16b_r0_0 =
                        _mm_unpacklo_epi8(u1_incr_not_8x16b_r0_even, u1_incr_not_8x16b_r0_odd);
                    u1_incr_not_8x16b_r1_0 =
                        _mm_unpacklo_epi8(u1_incr_not_8x16b_r1_even, u1_incr_not_8x16b_r1_odd);
                    x_ref_pos_mask_temp_r0_0 =
                        _mm_unpacklo_epi8(x_ref_pos_mask_temp_r0_even, x_ref_pos_mask_temp_r0_odd);
                    x_ref_pos_mask_temp_r1_0 =
                        _mm_unpacklo_epi8(x_ref_pos_mask_temp_r1_even, x_ref_pos_mask_temp_r1_odd);

                    ref_arr_temp0_8x16b_r0_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_0, u1_incr_not_8x16b_r0_0);
                    ref_arr_temp0_8x16b_r1_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_0, u1_incr_not_8x16b_r1_0);
                    ref_arr_temp1_8x16b_r0_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r0_0, x_ref_pos_mask_temp_r0_0);
                    ref_arr_temp1_8x16b_r1_0 =
                        _mm_shuffle_epi8(ref_arr_8x16b_r1_0, x_ref_pos_mask_temp_r1_0);

                    res0_8x16b_r0_0 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r0_0, phs_mask_16min_8x16b_0);
                    res0_8x16b_r1_0 =
                        _mm_mullo_epi16(ref_arr_temp0_8x16b_r1_0, phs_mask_16min_8x16b_0);
                    res1_8x16b_r0_0 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r0_0, phs_mask_8x16b_0);
                    res1_8x16b_r1_0 = _mm_mullo_epi16(ref_arr_temp1_8x16b_r1_0, phs_mask_8x16b_0);

                    res_8x16b_r0_0 = _mm_add_epi16(res0_8x16b_r0_0, res1_8x16b_r0_0);
                    res_8x16b_r1_0 = _mm_add_epi16(res0_8x16b_r1_0, res1_8x16b_r1_0);

                    pu1_ref_y_ptr_incr_temp =
                        pu1_ref_y_ptr_incr + (pi1_y_ref_pos[i4_y] * i4_refarray_wd);
                    u1_y_incr_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pu1_ref_y_ptr_incr_temp));

                    u1_y_incr_8x16b_r0_0 =
                        _mm_shuffle_epi8(u1_y_incr_8x16b_r0_0, x_ref_rnd_mask_r0_0);

                    u1_y_incr_8x16b_r0_0 = _mm_cvtepi8_epi16(u1_y_incr_8x16b_r0_0);
                    u1_y_incr_8x16b_r0_0 = _mm_cmpeq_epi16(u1_y_incr_8x16b_r0_0, const_ones_8x16b);
                    u1_prev_y_incr_8x16b_r0_0 = u1_y_incr_8x16b_r0_0;

                    prev_res_8x16b_r0_0 = res_8x16b_r0_0;
                    prev_res_8x16b_r1_0 = res_8x16b_r1_0;
                }
            }

            if(zero_r0_r1)
            {
                res_4x32b_l_0 = _mm_set1_epi32(0);
                res_4x32b_l_1 = _mm_set1_epi32(0);
            }
            else
            {
                i4_y_phase = pi1_y_phase[i4_y];

                if((i4_y_phase) >> 3)
                {
                    vert_res0_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r1_0, res_8x16b_r0_0, u1_y_incr_8x16b_r0_0);
                    vert_res1_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r1_0, res_8x16b_r1_0, u1_y_incr_8x16b_r0_0);
                }
                else
                {
                    vert_res0_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r0_0, res_8x16b_r0_0, u1_y_incr_8x16b_r0_0);
                    vert_res1_8x16b_r0_0 =
                        _mm_blendv_epi8(res_8x16b_r0_0, res_8x16b_r1_0, u1_y_incr_8x16b_r0_0);
                }

                res0_8x16b_r0_0 = _mm_unpacklo_epi16(vert_res0_8x16b_r0_0, vert_res1_8x16b_r0_0);
                res1_8x16b_r0_0 = _mm_unpackhi_epi16(vert_res0_8x16b_r0_0, vert_res1_8x16b_r0_0);

                phs_y_mask_16min_8x16b = _mm_set1_epi16(16 - i4_y_phase);
                phs_y_mask_8x16b = _mm_set1_epi16(i4_y_phase);
                phs_y_mask_mix_8x16b = _mm_unpacklo_epi16(phs_y_mask_16min_8x16b, phs_y_mask_8x16b);

                res_4x32b_l_0 = _mm_madd_epi16(res0_8x16b_r0_0, phs_y_mask_mix_8x16b);
                res_4x32b_l_1 = _mm_madd_epi16(res1_8x16b_r0_0, phs_y_mask_mix_8x16b);
                res_4x32b_l_0 = _mm_add_epi32(res_4x32b_l_0, const_128);
                res_4x32b_l_1 = _mm_add_epi32(res_4x32b_l_1, const_128);

                res_4x32b_l_0 = _mm_srai_epi32(res_4x32b_l_0, 8);
                res_4x32b_l_1 = _mm_srai_epi32(res_4x32b_l_1, 8);
            }
            out_stride_temp = (i4_y * i4_out_stride);

            out_4x32b_l = _mm_loadu_si128((__m128i *) (pi2_out + out_stride_temp));
            out_4x32b_h = _mm_loadu_si128((__m128i *) (pi2_out + out_stride_temp + 8));

            out_4x32b_l = _mm_and_si128(out_4x32b_l, chroma_mask);
            out_4x32b_h = _mm_and_si128(out_4x32b_h, chroma_mask);

            res_4x32b_l_0 = _mm_and_si128(res_4x32b_l_0, chroma_mask2);
            res_4x32b_l_1 = _mm_and_si128(res_4x32b_l_1, chroma_mask2);

            out_4x32b_l = _mm_add_epi8(res_4x32b_l_0, out_4x32b_l);
            out_4x32b_h = _mm_add_epi8(res_4x32b_l_1, out_4x32b_h);

            _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp), out_4x32b_l);
            _mm_storeu_si128((__m128i *) (pi2_out + out_stride_temp + 8), out_4x32b_h);
        }
    }
    return;
} /* End of Interpolation Function */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_reflayer_const_non_boundary_mb_sse42       */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore         creation                             */
/*                                                                           */
/*****************************************************************************/

void isvcd_residual_reflayer_const_non_boundary_mb_sse42(
    WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride, WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
    WORD32 i4_refarray_ht, WORD32 i4_ref_mb_type_q0, WORD32 i4_ref_mb_type_q1,
    WORD32 i4_ref_mb_type_q2, WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
    WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag)
{
    WORD32 i4_y;

    WORD16 *pi2_ref_data_byte;
    WORD16 *pi2_ref_array_temp;
    if(i4_chroma_flag == 0)
    {
        WORD8 index_0[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        __m128i ref_mb_type_8x16_q0, ref_mb_type_8x16_q1, ref_mb_type_8x16_q2, ref_mb_type_8x16_q3,
            mb_quard1_part_x_8x16;
        __m128i ref_mb_type_8x16_0, ref_mb_type_8x16_1;
        __m128i ref_mb_type_8x16_low_0, ref_mb_type_8x16_low_1;
        __m128i mb_type_mask_8x16_0 = _mm_set1_epi8(-1);
        __m128i mb_type_mask_8x16_1 = _mm_set1_epi8(-1);
        __m128i mb_type_mask_8x16_low_0, mb_type_mask_8x16_low_1;
        __m128i mask_8x16_0;
        __m128i index_arr_0;
        __m128i inp_data_16x8_0, inp_data_16x8_1;
        __m128i res_16x8_0, res_16x8_1;
        __m128i one_8x16 = _mm_set1_epi8(1);
        __m128i zero_8x16 = _mm_set1_epi8(0);

        index_arr_0 = _mm_loadu_si128((__m128i *) index_0);
        ref_mb_type_8x16_q0 = _mm_set1_epi8(i4_ref_mb_type_q0);
        ref_mb_type_8x16_q1 = _mm_set1_epi8(i4_ref_mb_type_q1);
        ref_mb_type_8x16_q2 = _mm_set1_epi8(i4_ref_mb_type_q2);
        ref_mb_type_8x16_q3 = _mm_set1_epi8(i4_ref_mb_type_q3);
        if((i4_mb_quard1_part_x >= i4_refarray_wd) && (i4_mb_quard1_part_y >= i4_refarray_ht))
        {
            // Quard 0
            ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
            ref_mb_type_8x16_1 = ref_mb_type_8x16_q0;
            mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
            mb_type_mask_8x16_1 = mb_type_mask_8x16_0;
        }
        else if((i4_mb_quard1_part_y >= (i4_refarray_ht - 1)) &&
                (i4_mb_quard1_part_x < i4_refarray_wd))
        {
            // Quard 0 & 1
            if(i4_mb_quard1_part_x == 8)
            {
                ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                ref_mb_type_8x16_1 = ref_mb_type_8x16_q1;
            }
            else if(i4_mb_quard1_part_x < 8)
            {
                mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x << 1));
                mask_8x16_0 =
                    _mm_cmplt_epi8(index_arr_0, mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

                ref_mb_type_8x16_0 =
                    _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);
                ref_mb_type_8x16_1 = ref_mb_type_8x16_q1;
            }
            else
            {
                mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x - 8) << 1);
                mask_8x16_0 =
                    _mm_cmplt_epi8(index_arr_0, mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

                ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                ref_mb_type_8x16_1 =
                    _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);
            }

            mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
            mb_type_mask_8x16_1 = _mm_cmpeq_epi8(ref_mb_type_8x16_1, one_8x16);
        }
        else
        {
            if(i4_mb_quard1_part_x >= i4_refarray_wd)
            {
                ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                ref_mb_type_8x16_1 = ref_mb_type_8x16_q0;

                ref_mb_type_8x16_low_0 = ref_mb_type_8x16_q2;
                ref_mb_type_8x16_low_1 = ref_mb_type_8x16_q2;
            }
            else
            {
                // Quard 0, 1, 2, 3
                if(i4_mb_quard1_part_x == 8)
                {
                    ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                    ref_mb_type_8x16_1 = ref_mb_type_8x16_q1;

                    ref_mb_type_8x16_low_0 = ref_mb_type_8x16_q2;
                    ref_mb_type_8x16_low_1 = ref_mb_type_8x16_q3;
                }
                else if(i4_mb_quard1_part_x < 8)
                {
                    mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x << 1));
                    mask_8x16_0 = _mm_cmplt_epi8(index_arr_0,
                                                 mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

                    ref_mb_type_8x16_0 =
                        _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);
                    ref_mb_type_8x16_1 = ref_mb_type_8x16_q1;

                    ref_mb_type_8x16_low_0 =
                        _mm_blendv_epi8(ref_mb_type_8x16_q3, ref_mb_type_8x16_q2, mask_8x16_0);
                    ref_mb_type_8x16_low_1 = ref_mb_type_8x16_q3;
                }
                else
                {
                    mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x - 8) << 1);
                    mask_8x16_0 = _mm_cmplt_epi8(index_arr_0,
                                                 mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

                    ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                    ref_mb_type_8x16_1 =
                        _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);

                    ref_mb_type_8x16_low_0 = ref_mb_type_8x16_q2;
                    ref_mb_type_8x16_low_1 =
                        _mm_blendv_epi8(ref_mb_type_8x16_q3, ref_mb_type_8x16_q2, mask_8x16_0);
                }
                mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
                mb_type_mask_8x16_1 = _mm_cmpeq_epi8(ref_mb_type_8x16_1, one_8x16);

                mb_type_mask_8x16_low_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_low_0, one_8x16);
                mb_type_mask_8x16_low_1 = _mm_cmpeq_epi8(ref_mb_type_8x16_low_1, one_8x16);
            }
        }

        if(i4_mb_quard1_part_y < i4_refarray_ht - 1)
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte));
                inp_data_16x8_1 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte + 8));

                if(i4_y < i4_mb_quard1_part_y)
                {
                    res_16x8_0 = _mm_blendv_epi8(zero_8x16, inp_data_16x8_0, mb_type_mask_8x16_0);
                    res_16x8_1 = _mm_blendv_epi8(zero_8x16, inp_data_16x8_1, mb_type_mask_8x16_1);
                }
                else
                {
                    res_16x8_0 =
                        _mm_blendv_epi8(zero_8x16, inp_data_16x8_0, mb_type_mask_8x16_low_0);
                    res_16x8_1 =
                        _mm_blendv_epi8(zero_8x16, inp_data_16x8_1, mb_type_mask_8x16_low_1);
                }

                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp), res_16x8_0);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp + 8), res_16x8_1);
            }
        }
        else
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte));
                inp_data_16x8_1 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte + 8));

                res_16x8_0 = _mm_blendv_epi8(zero_8x16, inp_data_16x8_0, mb_type_mask_8x16_0);
                res_16x8_1 = _mm_blendv_epi8(zero_8x16, inp_data_16x8_1, mb_type_mask_8x16_1);

                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp), res_16x8_0);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp + 8), res_16x8_1);
            }
        }
    }
    else
    {
        WORD8 index_0[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        WORD8 even_mask_arr[16] = {0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15};
        __m128i ref_mb_type_8x16_q0, ref_mb_type_8x16_q1, ref_mb_type_8x16_q2, ref_mb_type_8x16_q3,
            mb_quard1_part_x_8x16;
        __m128i ref_mb_type_8x16_0;
        __m128i ref_mb_type_8x16_low_0;
        __m128i mb_type_mask_8x16_0 = _mm_set1_epi8(-1);
        __m128i mb_type_mask_8x16_low_0;
        __m128i mask_8x16_0;
        __m128i index_arr_0, even_mask;
        __m128i inp_data_16x8_0, inp_data_16x8_1, inp_data_16x8;
        __m128i res_16x8_0;
        __m128i one_8x16 = _mm_set1_epi8(1);
        __m128i zero_8x16 = _mm_set1_epi8(0);

        index_arr_0 = _mm_loadu_si128((__m128i *) index_0);
        even_mask = _mm_loadu_si128((__m128i *) even_mask_arr);

        ref_mb_type_8x16_q0 = _mm_set1_epi8(i4_ref_mb_type_q0);
        ref_mb_type_8x16_q1 = _mm_set1_epi8(i4_ref_mb_type_q1);
        ref_mb_type_8x16_q2 = _mm_set1_epi8(i4_ref_mb_type_q2);
        ref_mb_type_8x16_q3 = _mm_set1_epi8(i4_ref_mb_type_q3);
        if((i4_mb_quard1_part_x >= i4_refarray_wd) && (i4_mb_quard1_part_y >= i4_refarray_ht))
        {
            // Quard 0
            ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
            mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
        }
        else if((i4_mb_quard1_part_y >= (i4_refarray_ht - 1)) &&
                (i4_mb_quard1_part_x < i4_refarray_wd))
        {
            // Quard 0 & 1
            mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x << 1));
            mask_8x16_0 =
                _mm_cmplt_epi8(index_arr_0, mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

            ref_mb_type_8x16_0 =
                _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);
            mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
        }
        else
        {
            if(i4_mb_quard1_part_x >= i4_refarray_wd)
            {
                // Quard 0 & 2
                ref_mb_type_8x16_0 = ref_mb_type_8x16_q0;
                ref_mb_type_8x16_low_0 = ref_mb_type_8x16_q2;
            }
            else
            {
                // Quard 0, 1, 2, 3
                mb_quard1_part_x_8x16 = _mm_set1_epi8((i4_mb_quard1_part_x << 1));
                mask_8x16_0 =
                    _mm_cmplt_epi8(index_arr_0, mb_quard1_part_x_8x16);  // return 1 if a<b, else 0

                ref_mb_type_8x16_0 =
                    _mm_blendv_epi8(ref_mb_type_8x16_q1, ref_mb_type_8x16_q0, mask_8x16_0);
                ref_mb_type_8x16_low_0 =
                    _mm_blendv_epi8(ref_mb_type_8x16_q3, ref_mb_type_8x16_q2, mask_8x16_0);

                mb_type_mask_8x16_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_0, one_8x16);
                mb_type_mask_8x16_low_0 = _mm_cmpeq_epi8(ref_mb_type_8x16_low_0, one_8x16);
            }
        }

        if(i4_mb_quard1_part_y < i4_refarray_ht - 1)
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte));
                inp_data_16x8_1 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte + 8));

                inp_data_16x8_0 = _mm_shuffle_epi8(inp_data_16x8_0, even_mask);
                inp_data_16x8_1 = _mm_shuffle_epi8(inp_data_16x8_1, even_mask);

                inp_data_16x8 = _mm_unpacklo_epi64(inp_data_16x8_0, inp_data_16x8_1);
                if(i4_y < i4_mb_quard1_part_y)
                {
                    res_16x8_0 = _mm_blendv_epi8(zero_8x16, inp_data_16x8, mb_type_mask_8x16_0);
                }
                else
                {
                    res_16x8_0 = _mm_blendv_epi8(zero_8x16, inp_data_16x8, mb_type_mask_8x16_low_0);
                }

                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp), res_16x8_0);
            }
        }
        else
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte));
                inp_data_16x8_1 = _mm_loadu_si128((__m128i *) (pi2_ref_data_byte + 8));

                inp_data_16x8_0 = _mm_shuffle_epi8(inp_data_16x8_0, even_mask);
                inp_data_16x8_1 = _mm_shuffle_epi8(inp_data_16x8_1, even_mask);
                inp_data_16x8 = _mm_unpacklo_epi64(inp_data_16x8_0, inp_data_16x8_1);

                res_16x8_0 = _mm_blendv_epi8(zero_8x16, inp_data_16x8, mb_type_mask_8x16_0);
                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                _mm_storeu_si128((__m128i *) (pi2_ref_array_temp), res_16x8_0);
            }
        }
    }
}
