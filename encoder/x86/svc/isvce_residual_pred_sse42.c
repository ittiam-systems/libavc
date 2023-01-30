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
*
* @file
*  isvce_residual_pred_sse42.c
*
* @brief
*  Contains functions
* used for SVC residual
* prediction
*
*******************************************************************************
*/
#include <immintrin.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "isvc_structs.h"

void isvce_luma_residual_sampler_2x_sse42(coordinates_t *ps_ref_array_positions,
                                          coordinates_t *ps_ref_array_phases,
                                          buffer_container_t *ps_inp, buffer_container_t *ps_out,
                                          buffer_container_t *ps_scratch, UWORD32 u4_ref_nnz,
                                          UWORD8 u1_ref_tx_size)
{
    WORD16 *pi2_inp_data = (WORD16 *) ps_inp->pv_data;
    WORD16 *pi2_out_res = (WORD16 *) ps_out->pv_data;
    WORD32 i4_inp_data_stride = ps_inp->i4_data_stride;
    WORD32 i4_out_res_stride = ps_out->i4_data_stride;
    WORD16 *pi2_refarray_buffer = (WORD16 *) ps_scratch->pv_data;
    WORD32 i4_blk_ctr;

    UNUSED(ps_ref_array_positions);
    UNUSED(ps_ref_array_phases);

    /* For 2x scaling, offsets always point to TL pixel outside MB */
    /* Hence, refTransBlkIdc will be different and since phase */
    /* for first refArray pos for horiz filtering samples > 8, */
    /* first row and first column from the refArray is never used */
    pi2_inp_data += 1 + i4_inp_data_stride;

    if((u1_ref_tx_size) && (0 != u4_ref_nnz))
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
        for(i4_i = 0; i4_i < BLK8x8SIZE; i4_i += 2)
        {
            /* a0 a1 a2 a3 a4 a5 a6 a7 */
            i2_coeff_8x16b_r1_0 = _mm_loadu_si128((__m128i *) pi2_ref_data_byte);
            /* b0 b1 b2 b3 b4 b5 b6 b7 */
            i2_coeff_8x16b_r2_0 =
                _mm_loadu_si128((__m128i *) (pi2_ref_data_byte + i4_inp_data_stride));

            /* a1 a2 a3 a4 a5 a6 a7 0 */
            i2_coeff_8x16b_r1_1 = _mm_srli_si128(i2_coeff_8x16b_r1_0, 2);
            /* b1 b2 b3 b4 b5 b6 b7 0 */
            i2_coeff_8x16b_r2_1 = _mm_srli_si128(i2_coeff_8x16b_r2_0, 2);

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

            /* vertical loop updates */
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
                pi2_refarray_buffer += MB_SIZE;

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
            if(0 != (u4_ref_nnz & 0x1))
            {
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
                    {
                        /* a0 a1 a2 a3 a4 a5 a6 a7 */
                        i2_coeff_8x16b_r1_0 = _mm_loadu_si128((__m128i *) pi2_inp_data);
                        /* b0 b1 b2 b3 b4 b5 b6 b7 */
                        i2_coeff_8x16b_r2_0 =
                            _mm_loadu_si128((__m128i *) (pi2_inp_data + i4_inp_data_stride));
                        i2_coeff_8x16b_r3_0 =
                            _mm_loadu_si128((__m128i *) (pi2_inp_data + (i4_inp_data_stride << 1)));
                        i2_coeff_8x16b_r4_0 =
                            _mm_loadu_si128((__m128i *) (pi2_inp_data + (i4_inp_data_stride * 3)));

                        /* a1 a2 a3 a4 a5 a6 a7 0 */
                        i2_coeff_8x16b_r1_1 = _mm_srli_si128(i2_coeff_8x16b_r1_0, 2);
                        /* b1 b2 b3 b4 b5 b6 b7 0 */
                        i2_coeff_8x16b_r2_1 = _mm_srli_si128(i2_coeff_8x16b_r2_0, 2);
                        i2_coeff_8x16b_r3_1 = _mm_srli_si128(i2_coeff_8x16b_r3_0, 2);
                        i2_coeff_8x16b_r4_1 = _mm_srli_si128(i2_coeff_8x16b_r4_0, 2);

                        coeff_add_8x16b_r1 =
                            _mm_add_epi16(i2_coeff_8x16b_r1_0, i2_coeff_8x16b_r1_1);
                        coeff_add_8x16b_r2 =
                            _mm_add_epi16(i2_coeff_8x16b_r2_0, i2_coeff_8x16b_r2_1);
                        coeff_add_8x16b_r3 =
                            _mm_add_epi16(i2_coeff_8x16b_r3_0, i2_coeff_8x16b_r3_1);
                        coeff_add_8x16b_r4 =
                            _mm_add_epi16(i2_coeff_8x16b_r4_0, i2_coeff_8x16b_r4_1);

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

                        _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 1),
                                         final_res_8x16b_r1_0);
                        _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 9),
                                         final_res_8x16b_r2_0);
                        _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 17),
                                         final_res_8x16b_r3_0);
                        _mm_storeu_si128((__m128i *) (pi2_refarray_buffer + 25),
                                         final_res_8x16b_r4_0);

                        pi2_refarray_buffer[0] = (pi2_inp_data[0] << 2);
                        pi2_refarray_buffer[7] = (pi2_inp_data[3] << 2);
                        pi2_refarray_buffer[8] = (pi2_inp_data[i4_inp_data_stride] << 2);
                        pi2_refarray_buffer[15] = (pi2_inp_data[i4_inp_data_stride + 3] << 2);
                        pi2_refarray_buffer[16] = (pi2_inp_data[(i4_inp_data_stride << 1)] << 2);
                        pi2_refarray_buffer[23] =
                            (pi2_inp_data[(i4_inp_data_stride << 1) + 3] << 2);
                        pi2_refarray_buffer[24] = (pi2_inp_data[(i4_inp_data_stride * 3)] << 2);
                        pi2_refarray_buffer[31] = (pi2_inp_data[(i4_inp_data_stride * 3) + 3] << 2);
                    }

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

                        i4_horz_samp_8x16b_r0_1 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer));
                        i4_horz_samp_8x16b_r0_2 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + 4));
                        i4_horz_samp_8x16b_r1_1 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + BLK8x8SIZE));
                        i4_horz_samp_8x16b_r1_2 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + BLK8x8SIZE + 4));
                        i4_horz_samp_8x16b_r2_1 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLK8x8SIZE << 1)));
                        i4_horz_samp_8x16b_r2_2 = _mm_loadu_si128(
                            (__m128i *) (pi2_refarray_buffer + (BLK8x8SIZE << 1) + 4));
                        i4_horz_samp_8x16b_r3_1 =
                            _mm_loadu_si128((__m128i *) (pi2_refarray_buffer + (BLK8x8SIZE * 3)));
                        i4_horz_samp_8x16b_r3_2 = _mm_loadu_si128(
                            (__m128i *) (pi2_refarray_buffer + (BLK8x8SIZE * 3) + 4));

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

                        pi2_out_res += BLK8x8SIZE;
                    }
                }
            }
            else
            {
                pi2_out_res += BLK8x8SIZE;
            }

            /* Block level loop updates */
            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= 4;
                pi2_inp_data += (i4_inp_data_stride * 4);
                pi2_out_res -= MB_SIZE;
                pi2_out_res += (i4_out_res_stride * BLK8x8SIZE);
                u4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += 4;
            }

            u4_ref_nnz >>= 1;

        } /* end of loop over all the blocks */
    }
}

UWORD32 isvce_get_sad_with_residual_pred_sse42(buffer_container_t *ps_src,
                                               buffer_container_t *ps_pred,
                                               buffer_container_t *ps_res, UWORD32 u4_mb_wd,
                                               UWORD32 u4_mb_ht)
{
    UWORD32 i, j, u4_sad = 0;
    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    UWORD32 u4_num_rows_per_loop = 8;
    UWORD32 u4_ht_by_8 = u4_mb_ht / u4_num_rows_per_loop;

    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i src_r4, src_r5, src_r6, src_r7;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i pred_r4, pred_r5, pred_r6, pred_r7;
    __m128i res_r0, res_r1, res_r2, res_r3;
    __m128i res_r4, res_r5, res_r6, res_r7;
    __m128i zero_4x32 = _mm_set1_epi32((WORD32) 0);

    if((u4_mb_wd == 16) && (u4_mb_ht % 8 == 0))
    {
        for(i = 0; i < u4_ht_by_8; i++)
        {
            for(j = 0; j < 2; j++)
            {
                src_r0 = _mm_loadl_epi64((__m128i *) (pu1_src));
                src_r1 = _mm_loadl_epi64((__m128i *) (pu1_src + 8));

                pu1_src += i4_src_stride;

                src_r2 = _mm_loadl_epi64((__m128i *) (pu1_src));
                src_r3 = _mm_loadl_epi64((__m128i *) (pu1_src + 8));

                pu1_src += i4_src_stride;

                src_r4 = _mm_loadl_epi64((__m128i *) (pu1_src));
                src_r5 = _mm_loadl_epi64((__m128i *) (pu1_src + 8));

                pu1_src += i4_src_stride;

                src_r6 = _mm_loadl_epi64((__m128i *) (pu1_src));
                src_r7 = _mm_loadl_epi64((__m128i *) (pu1_src + 8));

                pu1_src += i4_src_stride;

                pred_r0 = _mm_loadl_epi64((__m128i *) (pu1_pred));
                pred_r1 = _mm_loadl_epi64((__m128i *) (pu1_pred + 8));

                pu1_pred += i4_pred_stride;

                pred_r2 = _mm_loadl_epi64((__m128i *) (pu1_pred));
                pred_r3 = _mm_loadl_epi64((__m128i *) (pu1_pred + 8));

                pu1_pred += i4_pred_stride;

                pred_r4 = _mm_loadl_epi64((__m128i *) (pu1_pred));
                pred_r5 = _mm_loadl_epi64((__m128i *) (pu1_pred + 8));

                pu1_pred += i4_pred_stride;

                pred_r6 = _mm_loadl_epi64((__m128i *) (pu1_pred));
                pred_r7 = _mm_loadl_epi64((__m128i *) (pu1_pred + 8));

                pu1_pred += i4_pred_stride;

                src_r0 = _mm_cvtepu8_epi16(src_r0);
                src_r1 = _mm_cvtepu8_epi16(src_r1);
                src_r2 = _mm_cvtepu8_epi16(src_r2);
                src_r3 = _mm_cvtepu8_epi16(src_r3);
                src_r4 = _mm_cvtepu8_epi16(src_r4);
                src_r5 = _mm_cvtepu8_epi16(src_r5);
                src_r6 = _mm_cvtepu8_epi16(src_r6);
                src_r7 = _mm_cvtepu8_epi16(src_r7);

                pred_r0 = _mm_cvtepu8_epi16(pred_r0);
                pred_r1 = _mm_cvtepu8_epi16(pred_r1);
                pred_r2 = _mm_cvtepu8_epi16(pred_r2);
                pred_r3 = _mm_cvtepu8_epi16(pred_r3);
                pred_r4 = _mm_cvtepu8_epi16(pred_r4);
                pred_r5 = _mm_cvtepu8_epi16(pred_r5);
                pred_r6 = _mm_cvtepu8_epi16(pred_r6);
                pred_r7 = _mm_cvtepu8_epi16(pred_r7);

                res_r0 = _mm_loadu_si128((__m128i *) (pi2_res));
                res_r1 = _mm_loadu_si128((__m128i *) (pi2_res + 8));

                pi2_res += i4_res_stride;

                res_r2 = _mm_loadu_si128((__m128i *) (pi2_res));
                res_r3 = _mm_loadu_si128((__m128i *) (pi2_res + 8));

                pi2_res += i4_res_stride;

                res_r4 = _mm_loadu_si128((__m128i *) (pi2_res));
                res_r5 = _mm_loadu_si128((__m128i *) (pi2_res + 8));

                pi2_res += i4_res_stride;

                res_r6 = _mm_loadu_si128((__m128i *) (pi2_res));
                res_r7 = _mm_loadu_si128((__m128i *) (pi2_res + 8));

                pi2_res += i4_res_stride;

                src_r0 = _mm_sub_epi16(src_r0, pred_r0);
                src_r1 = _mm_sub_epi16(src_r1, pred_r1);
                src_r2 = _mm_sub_epi16(src_r2, pred_r2);
                src_r3 = _mm_sub_epi16(src_r3, pred_r3);
                src_r4 = _mm_sub_epi16(src_r4, pred_r4);
                src_r5 = _mm_sub_epi16(src_r5, pred_r5);
                src_r6 = _mm_sub_epi16(src_r6, pred_r6);
                src_r7 = _mm_sub_epi16(src_r7, pred_r7);

                src_r0 = _mm_sub_epi16(src_r0, res_r0);
                src_r1 = _mm_sub_epi16(src_r1, res_r1);
                src_r2 = _mm_sub_epi16(src_r2, res_r2);
                src_r3 = _mm_sub_epi16(src_r3, res_r3);
                src_r4 = _mm_sub_epi16(src_r4, res_r4);
                src_r5 = _mm_sub_epi16(src_r5, res_r5);
                src_r6 = _mm_sub_epi16(src_r6, res_r6);
                src_r7 = _mm_sub_epi16(src_r7, res_r7);

                src_r0 = _mm_abs_epi16(src_r0);
                src_r1 = _mm_abs_epi16(src_r1);
                src_r2 = _mm_abs_epi16(src_r2);
                src_r3 = _mm_abs_epi16(src_r3);
                src_r4 = _mm_abs_epi16(src_r4);
                src_r5 = _mm_abs_epi16(src_r5);
                src_r6 = _mm_abs_epi16(src_r6);
                src_r7 = _mm_abs_epi16(src_r7);

                src_r0 = _mm_adds_epu16(src_r0, src_r1);
                src_r1 = _mm_adds_epu16(src_r2, src_r3);
                src_r2 = _mm_adds_epu16(src_r4, src_r5);
                src_r3 = _mm_adds_epu16(src_r6, src_r7);

                src_r0 = _mm_adds_epu16(src_r0, src_r1);
                src_r1 = _mm_adds_epu16(src_r2, src_r3);

                src_r0 = _mm_adds_epu16(src_r0, src_r1);

                src_r1 = _mm_cvtepu16_epi32(src_r0);
                src_r2 = _mm_srli_si128(src_r0, 8);
                src_r2 = _mm_cvtepu16_epi32(src_r2);

                src_r0 = _mm_hadd_epi32(src_r1, src_r2);
                src_r0 = _mm_hadd_epi32(src_r0, zero_4x32);
                src_r0 = _mm_hadd_epi32(src_r0, zero_4x32);

                u4_sad += _mm_extract_epi32(src_r0, 0);
            }
        }
    }
    else
    {
        for(i = 0; i < u4_mb_ht; i++)
        {
            for(j = 0; j < u4_mb_wd; j++)
            {
                WORD16 i2_src = pu1_src[j + i * i4_src_stride];
                WORD16 i2_pred = pu1_pred[j + i * i4_pred_stride];
                WORD16 i2_res = pi2_res[j + i * i4_res_stride];
                u4_sad += ABS(i2_src - i2_pred - i2_res);
            }
        }
    }

    return u4_sad;
}
