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
/*!
 **************************************************************************

 * * \file ih264d_resamp_svc.c
 *
 * \brief
 *    Contains routines that
 * resample for SVC resampling
 *
 * Detailed_description
 *
 * \date
 *
 *
 *
 * \author

 * **************************************************************************

 */
#include <immintrin.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "isvc_intra_resample.h"

void isvc_interpolate_base_luma_dyadic_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                             UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_y;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp, *pu1_out;
    WORD16 *pi2_tmp;

    __m128i i4_samp_16x8b_0, i4_samp_16x8b_1, i4_samp_16x8b_2, i4_samp_16x8b_3;
    __m128i i4_samp_8x16b_0, i4_samp_8x16b_1, i4_samp_8x16b_2, i4_samp_8x16b_3;
    __m128i i4_res_8x16b_r1_1, i4_res_8x16b_r1_2, i4_res_8x16b_r1_3;
    __m128i i4_res_8x16b_r2_1, i4_res_8x16b_r2_2, i4_res_8x16b_r2_3;

    /* Filter coefficient values for phase 4 */
    __m128i i4_coeff_8x16b_0 = _mm_set1_epi16(-3);
    __m128i i4_coeff_8x16b_1 = _mm_set1_epi16(28);
    i4_filt_stride = 12;
    i4_src_stride = DYADIC_REF_W_Y;

    /* Initializing pointers */
    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    pu1_out = pu1_out_buf;

    /* Vertical interpolation */
    /*First 64 bit */
    /* y = 0, y_phase = 12 */
    i4_samp_16x8b_0 = _mm_loadl_epi64((__m128i *) (pu1_inp));
    i4_samp_16x8b_1 = _mm_loadl_epi64((__m128i *) (pu1_inp + i4_src_stride));
    i4_samp_16x8b_2 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1)));
    i4_samp_16x8b_3 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
    pu1_inp += (i4_src_stride << 2);
    i4_samp_8x16b_0 = _mm_cvtepu8_epi16(i4_samp_16x8b_0);
    i4_samp_8x16b_1 = _mm_cvtepu8_epi16(i4_samp_16x8b_1);
    i4_samp_8x16b_2 = _mm_cvtepu8_epi16(i4_samp_16x8b_2);
    i4_samp_8x16b_3 = _mm_cvtepu8_epi16(i4_samp_16x8b_3);

    /* since y_phase 12 for y = 0 */
    /*Multiply by 8 =>  left shift by 3*/
    i4_res_8x16b_r1_1 = _mm_slli_epi16(i4_samp_8x16b_1, 3);
    i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_2, i4_coeff_8x16b_1);
    i4_res_8x16b_r1_3 = _mm_mullo_epi16(i4_samp_8x16b_3, i4_coeff_8x16b_0);

    i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
    i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_0);
    i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);

    _mm_storeu_si128((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);
    pi2_tmp += i4_filt_stride;

    for(i4_y = 1; i4_y < 15; i4_y += 2)
    {
        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y		*/
        //*Multiply by 8 =>  left shift by 3*/
        i4_res_8x16b_r1_1 = _mm_mullo_epi16(i4_samp_8x16b_0, i4_coeff_8x16b_0);
        i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_1, i4_coeff_8x16b_1);
        i4_res_8x16b_r1_3 = _mm_slli_epi16(i4_samp_8x16b_2, 3);

        i4_res_8x16b_r2_1 = _mm_slli_epi16(i4_samp_8x16b_1, 3);
        i4_res_8x16b_r2_2 = _mm_mullo_epi16(i4_samp_8x16b_2, i4_coeff_8x16b_1);
        i4_res_8x16b_r2_3 = _mm_mullo_epi16(i4_samp_8x16b_3, i4_coeff_8x16b_0);

        i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_3);
        i4_res_8x16b_r2_3 = _mm_subs_epi16(i4_res_8x16b_r2_3, i4_samp_8x16b_0);

        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
        i4_res_8x16b_r2_1 = _mm_adds_epi16(i4_res_8x16b_r2_1, i4_res_8x16b_r2_2);

        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);
        i4_res_8x16b_r2_1 = _mm_adds_epi16(i4_res_8x16b_r2_1, i4_res_8x16b_r2_3);

        /* Storing the results */
        _mm_storeu_si128((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);
        _mm_storeu_si128((__m128i *) (pi2_tmp + i4_filt_stride), i4_res_8x16b_r2_1);
        pi2_tmp += (i4_filt_stride << 1);
        pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* y = 15, y_phase = 4 */
        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        i4_res_8x16b_r1_1 = _mm_mullo_epi16(i4_samp_8x16b_0, i4_coeff_8x16b_0);
        i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_1, i4_coeff_8x16b_1);
        i4_res_8x16b_r1_3 = _mm_slli_epi16(i4_samp_8x16b_2, 3);
        i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_3);

        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);

        /* Store the output */
        _mm_storeu_si128((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;

    /*Remaining 32 bit */
    pu1_inp += 8;
    pi2_tmp += 8;

        /* y = 0, y_phase = 12 */
        i4_samp_16x8b_0 = _mm_loadl_epi64((__m128i *) (pu1_inp));
        i4_samp_16x8b_1 = _mm_loadl_epi64((__m128i *) (pu1_inp + i4_src_stride));
        i4_samp_16x8b_2 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1)));
        i4_samp_16x8b_3 =
            _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
        pu1_inp += (i4_src_stride << 2);
        i4_samp_8x16b_0 = _mm_cvtepu8_epi16(i4_samp_16x8b_0);
        i4_samp_8x16b_1 = _mm_cvtepu8_epi16(i4_samp_16x8b_1);
        i4_samp_8x16b_2 = _mm_cvtepu8_epi16(i4_samp_16x8b_2);
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(i4_samp_16x8b_3);

        /* since y_phase 12 for y = 0 */
        /*Multiply by 8 =>  left shift by 3*/
        i4_res_8x16b_r1_1 = _mm_slli_epi16(i4_samp_8x16b_1, 3);
        i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_2, i4_coeff_8x16b_1);
        i4_res_8x16b_r1_3 = _mm_mullo_epi16(i4_samp_8x16b_3, i4_coeff_8x16b_0);

        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
        i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_0);
        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);

        _mm_storel_epi64((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 15; i4_y += 2)
        {
            i4_samp_8x16b_0 = i4_samp_8x16b_1;
            i4_samp_8x16b_1 = i4_samp_8x16b_2;
            i4_samp_8x16b_2 = i4_samp_8x16b_3;
            i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

            /* y_phase is 4 for odd values of y */
            /* and 12 for even values of y		*/
            //*Multiply by 8 =>  left shift by 3*/
            i4_res_8x16b_r1_1 = _mm_mullo_epi16(i4_samp_8x16b_0, i4_coeff_8x16b_0);
            i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_1, i4_coeff_8x16b_1);
            i4_res_8x16b_r1_3 = _mm_slli_epi16(i4_samp_8x16b_2, 3);

            i4_res_8x16b_r2_1 = _mm_slli_epi16(i4_samp_8x16b_1, 3);
            i4_res_8x16b_r2_2 = _mm_mullo_epi16(i4_samp_8x16b_2, i4_coeff_8x16b_1);
            i4_res_8x16b_r2_3 = _mm_mullo_epi16(i4_samp_8x16b_3, i4_coeff_8x16b_0);

            i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_3);
            i4_res_8x16b_r2_3 = _mm_subs_epi16(i4_res_8x16b_r2_3, i4_samp_8x16b_0);

            i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
            i4_res_8x16b_r2_1 = _mm_adds_epi16(i4_res_8x16b_r2_1, i4_res_8x16b_r2_2);

            i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);
            i4_res_8x16b_r2_1 = _mm_adds_epi16(i4_res_8x16b_r2_1, i4_res_8x16b_r2_3);

            /* Storing the results */
            _mm_storel_epi64((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);
            _mm_storel_epi64((__m128i *) (pi2_tmp + i4_filt_stride), i4_res_8x16b_r2_1);
            pi2_tmp += (i4_filt_stride << 1);
            pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* y = 15, y_phase = 4 */
        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        i4_res_8x16b_r1_1 = _mm_mullo_epi16(i4_samp_8x16b_0, i4_coeff_8x16b_0);
        i4_res_8x16b_r1_2 = _mm_mullo_epi16(i4_samp_8x16b_1, i4_coeff_8x16b_1);
        i4_res_8x16b_r1_3 = _mm_slli_epi16(i4_samp_8x16b_2, 3);
        i4_res_8x16b_r1_3 = _mm_subs_epi16(i4_res_8x16b_r1_3, i4_samp_8x16b_3);

        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_2);
        i4_res_8x16b_r1_1 = _mm_adds_epi16(i4_res_8x16b_r1_1, i4_res_8x16b_r1_3);

        /* Store the output */
        _mm_storel_epi64((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;

    {
        __m128i coeff_c0_c1_8x16b = _mm_set_epi16(28, -3, 28, -3, 28, -3, 28, -3);
        __m128i coeff_c2_c3_8x16b = _mm_set_epi16(-1, 8, -1, 8, -1, 8, -1, 8);
        __m128i coeff_c3_c2_8x16b = _mm_set_epi16(8, -1, 8, -1, 8, -1, 8, -1);
        __m128i coeff_c1_c0_8x16b = _mm_set_epi16(-3, 28, -3, 28, -3, 28, -3, 28);

        __m128i i4_samp_8x16b_rpart1_0, i4_samp_8x16b_rpart2_0;
        __m128i i4_samp_8x16b_rpart1_1, i4_samp_8x16b_rpart2_1;
        __m128i i4_samp_8x16b_rpart1_2, i4_samp_8x16b_rpart2_2;
        __m128i i4_samp_8x16b_rpart1_3, i4_samp_8x16b_rpart2_3;
        __m128i i4_samp_8x16b_rpart1_4, i4_samp_8x16b_rpart2_4;

        __m128i i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart2_0;
        __m128i i4_res_4x32b_rpart1_1, i4_res_4x32b_rpart2_1;
        __m128i i4_res_4x32b_rpart1_2, i4_res_4x32b_rpart2_2;
        __m128i i4_res_4x32b_rpart1_3, i4_res_4x32b_rpart2_3;

        __m128i res_512 = _mm_set1_epi32(512);
        /* Horizontal interpolation */
        for(i4_y = 0; i4_y < 16; i4_y++)
        {
            i4_samp_8x16b_rpart1_0 = _mm_loadu_si128((__m128i *) pi2_tmp);
            i4_samp_8x16b_rpart2_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 4));

            i4_samp_8x16b_rpart1_1 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 2);
            i4_samp_8x16b_rpart1_2 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 4);
            i4_samp_8x16b_rpart1_3 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 6);
            i4_samp_8x16b_rpart1_4 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 8);

            i4_samp_8x16b_rpart2_1 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 2);
            i4_samp_8x16b_rpart2_2 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 4);
            i4_samp_8x16b_rpart2_3 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 6);
            i4_samp_8x16b_rpart2_4 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 8);

            i4_samp_8x16b_rpart1_0 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_0, i4_samp_8x16b_rpart1_1);
            i4_samp_8x16b_rpart1_1 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_1, i4_samp_8x16b_rpart1_2);
            i4_samp_8x16b_rpart1_2 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_2, i4_samp_8x16b_rpart1_3);
            i4_samp_8x16b_rpart1_3 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_3, i4_samp_8x16b_rpart1_4);

            i4_samp_8x16b_rpart2_0 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_0, i4_samp_8x16b_rpart2_1);
            i4_samp_8x16b_rpart2_1 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_1, i4_samp_8x16b_rpart2_2);
            i4_samp_8x16b_rpart2_2 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_2, i4_samp_8x16b_rpart2_3);
            i4_samp_8x16b_rpart2_3 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_3, i4_samp_8x16b_rpart2_4);

            i4_res_4x32b_rpart1_0 = _mm_madd_epi16(i4_samp_8x16b_rpart1_0, coeff_c3_c2_8x16b);
            i4_res_4x32b_rpart1_2 = _mm_madd_epi16(i4_samp_8x16b_rpart1_2, coeff_c1_c0_8x16b);

            i4_res_4x32b_rpart1_1 = _mm_madd_epi16(i4_samp_8x16b_rpart1_1, coeff_c0_c1_8x16b);
            i4_res_4x32b_rpart1_3 = _mm_madd_epi16(i4_samp_8x16b_rpart1_3, coeff_c2_c3_8x16b);

            i4_res_4x32b_rpart2_0 = _mm_madd_epi16(i4_samp_8x16b_rpart2_0, coeff_c3_c2_8x16b);
            i4_res_4x32b_rpart2_2 = _mm_madd_epi16(i4_samp_8x16b_rpart2_2, coeff_c1_c0_8x16b);

            i4_res_4x32b_rpart2_1 = _mm_madd_epi16(i4_samp_8x16b_rpart2_1, coeff_c0_c1_8x16b);
            i4_res_4x32b_rpart2_3 = _mm_madd_epi16(i4_samp_8x16b_rpart2_3, coeff_c2_c3_8x16b);

            i4_res_4x32b_rpart1_0 = _mm_add_epi32(i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart1_2);
            i4_res_4x32b_rpart1_1 = _mm_add_epi32(i4_res_4x32b_rpart1_1, i4_res_4x32b_rpart1_3);

            i4_res_4x32b_rpart2_0 = _mm_add_epi32(i4_res_4x32b_rpart2_0, i4_res_4x32b_rpart2_2);
            i4_res_4x32b_rpart2_1 = _mm_add_epi32(i4_res_4x32b_rpart2_1, i4_res_4x32b_rpart2_3);

            i4_res_4x32b_rpart1_2 =
                _mm_unpacklo_epi32(i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart1_1);
            i4_res_4x32b_rpart1_3 =
                _mm_unpackhi_epi32(i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart1_1);

            i4_res_4x32b_rpart2_2 =
                _mm_unpacklo_epi32(i4_res_4x32b_rpart2_0, i4_res_4x32b_rpart2_1);
            i4_res_4x32b_rpart2_3 =
                _mm_unpackhi_epi32(i4_res_4x32b_rpart2_0, i4_res_4x32b_rpart2_1);

            i4_res_4x32b_rpart1_0 = _mm_add_epi32(i4_res_4x32b_rpart1_2, res_512);
            i4_res_4x32b_rpart1_1 = _mm_add_epi32(i4_res_4x32b_rpart1_3, res_512);

            i4_res_4x32b_rpart1_0 = _mm_srai_epi32(i4_res_4x32b_rpart1_0, 10);
            i4_res_4x32b_rpart1_1 = _mm_srai_epi32(i4_res_4x32b_rpart1_1, 10);

            i4_res_4x32b_rpart2_0 = _mm_add_epi32(i4_res_4x32b_rpart2_2, res_512);
            i4_res_4x32b_rpart2_1 = _mm_add_epi32(i4_res_4x32b_rpart2_3, res_512);

            i4_res_4x32b_rpart2_0 = _mm_srai_epi32(i4_res_4x32b_rpart2_0, 10);
            i4_res_4x32b_rpart2_1 = _mm_srai_epi32(i4_res_4x32b_rpart2_1, 10);

            _mm_storeu_si128(
                (__m128i *) pu1_out,
                _mm_packus_epi16(_mm_packus_epi32(i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart1_1),
                                 _mm_packus_epi32(i4_res_4x32b_rpart2_0, i4_res_4x32b_rpart2_1)));

            pi2_tmp += i4_filt_stride;
            pu1_out += i4_out_stride;

        } /* End of loop over y */
    }
}

void isvc_vert_interpol_chroma_dyadic_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                            WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD8 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;
    __m128i i4_samp_16x8b_0, i4_samp_16x8b_1, i4_samp_16x8b_2, i4_samp_16x8b_3, i4_samp_16x8b_4,
        i4_samp_16x8b_5;
    __m128i i4_res_8x16b_r0, i4_res_8x16b_r1, i4_res_8x16b_r2, i4_res_8x16b_r3, i4_res_8x16b_r4,
        i4_res_8x16b_r5, i4_res_8x16b_r6, i4_res_8x16b_r7;
    __m128i i4_res_8x16b_r7_temp;
    __m128i i4_c0_c1_16x8b, i4_c2_c3_16x8b;

    i4_coeff_0 = (WORD8) (16 - i4_phase_0);
    i4_coeff_1 = (WORD8) (i4_phase_0);
    i4_coeff_2 = (WORD8) (16 - i4_phase_1);
    i4_coeff_3 = (WORD8) (i4_phase_1);

    i4_c0_c1_16x8b =
        _mm_set_epi8(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0);
    i4_c2_c3_16x8b =
        _mm_set_epi8(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2);

    /* Initializing pointers */
    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;

    i4_samp_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_inp));
    i4_samp_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_inp + i4_src_stride));
    i4_samp_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1)));
    i4_samp_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
    i4_samp_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 2)));
    i4_samp_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 2) + i4_src_stride));

    i4_samp_16x8b_0 = _mm_unpacklo_epi8(i4_samp_16x8b_0, i4_samp_16x8b_1);
    i4_res_8x16b_r0 = _mm_maddubs_epi16(i4_samp_16x8b_0, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp), i4_res_8x16b_r0);

    i4_samp_16x8b_1 = _mm_unpacklo_epi8(i4_samp_16x8b_1, i4_samp_16x8b_2);
    i4_res_8x16b_r1 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + i4_filt_stride), i4_res_8x16b_r1);

    i4_res_8x16b_r2 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1)), i4_res_8x16b_r2);

    i4_samp_16x8b_2 = _mm_unpacklo_epi8(i4_samp_16x8b_2, i4_samp_16x8b_3);
    i4_res_8x16b_r3 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1) + i4_filt_stride),
                     i4_res_8x16b_r3);

    i4_res_8x16b_r4 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2)), i4_res_8x16b_r4);

    i4_samp_16x8b_3 = _mm_unpacklo_epi8(i4_samp_16x8b_3, i4_samp_16x8b_4);
    i4_res_8x16b_r5 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + i4_filt_stride),
                     i4_res_8x16b_r5);

    i4_res_8x16b_r6 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c0_c1_16x8b);
    _mm_storel_epi64((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1)),
                     i4_res_8x16b_r6);

    i4_res_8x16b_r6 = _mm_shuffle_epi32(i4_res_8x16b_r6, 78);

    i4_samp_16x8b_4 = _mm_unpacklo_epi8(i4_samp_16x8b_4, i4_samp_16x8b_5);

    i4_res_8x16b_r7 = _mm_maddubs_epi16(i4_samp_16x8b_4, i4_c2_c3_16x8b);

    i4_res_8x16b_r7 = _mm_shuffle_epi32(i4_res_8x16b_r7, 147);

    i4_res_8x16b_r7_temp = _mm_blend_epi16(i4_res_8x16b_r6, i4_res_8x16b_r7, 252);

    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1) + 4),
                     i4_res_8x16b_r7_temp);
}

void isvc_horz_interpol_chroma_dyadic_sse42(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                            WORD32 i4_out_stride, WORD32 i4_phase_0,
                                            WORD32 i4_phase_1)
{
    WORD32 i4_dst_stride, i4_dst_stride2, i4_dst_stride4;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    __m128i i4_samp_8x16b_r1_0, i4_samp_8x16b_r1_1, i4_samp_8x16b_r1_2;
    __m128i i4_samp_8x16b_r2_0, i4_samp_8x16b_r2_1, i4_samp_8x16b_r2_2;
    __m128i i4_samp_8x16b_r3_0, i4_samp_8x16b_r3_1, i4_samp_8x16b_r3_2;
    __m128i i4_samp_8x16b_r4_0, i4_samp_8x16b_r4_1, i4_samp_8x16b_r4_2;
    __m128i i4_samp_8x16b_r5_0, i4_samp_8x16b_r5_1, i4_samp_8x16b_r5_2;
    __m128i i4_samp_8x16b_r6_0, i4_samp_8x16b_r6_1, i4_samp_8x16b_r6_2;
    __m128i i4_samp_8x16b_r7_0, i4_samp_8x16b_r7_1, i4_samp_8x16b_r7_2;
    __m128i i4_samp_8x16b_r8_0, i4_samp_8x16b_r8_1, i4_samp_8x16b_r8_2;

    __m128i i4_res_4x32b_r1_0, i4_res_4x32b_r1_1;
    __m128i i4_res_4x32b_r2_0, i4_res_4x32b_r2_1;
    __m128i i4_res_4x32b_r3_0, i4_res_4x32b_r3_1;
    __m128i i4_res_4x32b_r4_0, i4_res_4x32b_r4_1;
    __m128i i4_res_4x32b_r5_0, i4_res_4x32b_r5_1;
    __m128i i4_res_4x32b_r6_0, i4_res_4x32b_r6_1;
    __m128i i4_res_4x32b_r7_0, i4_res_4x32b_r7_1;
    __m128i i4_res_4x32b_r8_0, i4_res_4x32b_r8_1;

    __m128i i4_res_final_8x16b_r1, i4_res_final_8x16b_r2, i4_res_final_8x16b_r3,
        i4_res_final_8x16b_r4, i4_res_final_8x16b_r5, i4_res_final_8x16b_r6, i4_res_final_8x16b_r7,
        i4_res_final_8x16b_r8;

    __m128i out_16x8b_r1, out_16x8b_r2, out_16x8b_r3, out_16x8b_r4, out_16x8b_r5, out_16x8b_r6,
        out_16x8b_r7, out_16x8b_r8;

    __m128i i4_res_final_8x16b_r12_0, i4_res_final_8x16b_r12_1;
    __m128i i4_res_final_8x16b_r34_0, i4_res_final_8x16b_r34_1;
    __m128i i4_res_final_8x16b_r56_0, i4_res_final_8x16b_r56_1;
    __m128i i4_res_final_8x16b_r67_0, i4_res_final_8x16b_r67_1;
    __m128i chroma_mask, chroma_mask2;

    WORD32 i4_coeff_0 = 16 - i4_phase_0;
    WORD32 i4_coeff_1 = i4_phase_0;
    WORD32 i4_coeff_2 = 16 - i4_phase_1;
    WORD32 i4_coeff_3 = i4_phase_1;
    __m128i coeff_c0_c1_8x16b = _mm_set_epi16(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                                              i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0);
    __m128i coeff_c2_c3_8x16b = _mm_set_epi16(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                                              i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2);
    __m128i res_128 = _mm_set1_epi32(128);
    UWORD32 u4_norm_factor = 8;

    /* Initializing pointers */
    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_dst_stride = i4_out_stride;

    i4_dst_stride2 = i4_dst_stride << 1;
    i4_dst_stride4 = i4_dst_stride << 2;

    /* Horizontal interpolation */
    i4_samp_8x16b_r1_0 = _mm_loadu_si128((__m128i *) pi2_tmp);
    i4_samp_8x16b_r2_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 6));
    i4_samp_8x16b_r3_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 12));
    i4_samp_8x16b_r4_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 18));
    i4_samp_8x16b_r5_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 24));
    i4_samp_8x16b_r6_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 30));
    i4_samp_8x16b_r7_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 36));
    i4_samp_8x16b_r8_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 42));

    i4_samp_8x16b_r1_1 = _mm_srli_si128(i4_samp_8x16b_r1_0, 2);
    i4_samp_8x16b_r1_2 = _mm_srli_si128(i4_samp_8x16b_r1_0, 4);

    i4_samp_8x16b_r2_1 = _mm_srli_si128(i4_samp_8x16b_r2_0, 2);
    i4_samp_8x16b_r2_2 = _mm_srli_si128(i4_samp_8x16b_r2_0, 4);

    i4_samp_8x16b_r3_1 = _mm_srli_si128(i4_samp_8x16b_r3_0, 2);
    i4_samp_8x16b_r3_2 = _mm_srli_si128(i4_samp_8x16b_r3_0, 4);

    i4_samp_8x16b_r4_1 = _mm_srli_si128(i4_samp_8x16b_r4_0, 2);
    i4_samp_8x16b_r4_2 = _mm_srli_si128(i4_samp_8x16b_r4_0, 4);

    i4_samp_8x16b_r5_1 = _mm_srli_si128(i4_samp_8x16b_r5_0, 2);
    i4_samp_8x16b_r5_2 = _mm_srli_si128(i4_samp_8x16b_r5_0, 4);

    i4_samp_8x16b_r6_1 = _mm_srli_si128(i4_samp_8x16b_r6_0, 2);
    i4_samp_8x16b_r6_2 = _mm_srli_si128(i4_samp_8x16b_r6_0, 4);

    i4_samp_8x16b_r7_1 = _mm_srli_si128(i4_samp_8x16b_r7_0, 2);
    i4_samp_8x16b_r7_2 = _mm_srli_si128(i4_samp_8x16b_r7_0, 4);

    i4_samp_8x16b_r8_1 = _mm_srli_si128(i4_samp_8x16b_r8_0, 2);
    i4_samp_8x16b_r8_2 = _mm_srli_si128(i4_samp_8x16b_r8_0, 4);

    i4_samp_8x16b_r1_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r1_0, i4_samp_8x16b_r1_1);
    i4_samp_8x16b_r2_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r2_0, i4_samp_8x16b_r2_1);
    i4_samp_8x16b_r3_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r3_0, i4_samp_8x16b_r3_1);
    i4_samp_8x16b_r4_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r4_0, i4_samp_8x16b_r4_1);
    i4_samp_8x16b_r5_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r5_0, i4_samp_8x16b_r5_1);
    i4_samp_8x16b_r6_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r6_0, i4_samp_8x16b_r6_1);
    i4_samp_8x16b_r7_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r7_0, i4_samp_8x16b_r7_1);
    i4_samp_8x16b_r8_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r8_0, i4_samp_8x16b_r8_1);

    i4_samp_8x16b_r1_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r1_1, i4_samp_8x16b_r1_2);
    i4_samp_8x16b_r2_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r2_1, i4_samp_8x16b_r2_2);
    i4_samp_8x16b_r3_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r3_1, i4_samp_8x16b_r3_2);
    i4_samp_8x16b_r4_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r4_1, i4_samp_8x16b_r4_2);
    i4_samp_8x16b_r5_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r5_1, i4_samp_8x16b_r5_2);
    i4_samp_8x16b_r6_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r6_1, i4_samp_8x16b_r6_2);
    i4_samp_8x16b_r7_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r7_1, i4_samp_8x16b_r7_2);
    i4_samp_8x16b_r8_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r8_1, i4_samp_8x16b_r8_2);

    // a0c0 + a1c1  a1c0 + a2c1  a2c0 + a3c1  a3c0 + a4c1
    i4_res_4x32b_r1_0 = _mm_madd_epi16(i4_samp_8x16b_r1_0, coeff_c0_c1_8x16b);
    // b0c0+b1c1  b1c0+b2c1  b2c0+b3c1  b3c0+b4c1
    i4_res_4x32b_r2_0 = _mm_madd_epi16(i4_samp_8x16b_r2_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r3_0 = _mm_madd_epi16(i4_samp_8x16b_r3_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r4_0 = _mm_madd_epi16(i4_samp_8x16b_r4_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r5_0 = _mm_madd_epi16(i4_samp_8x16b_r5_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r6_0 = _mm_madd_epi16(i4_samp_8x16b_r6_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r7_0 = _mm_madd_epi16(i4_samp_8x16b_r7_0, coeff_c0_c1_8x16b);
    i4_res_4x32b_r8_0 = _mm_madd_epi16(i4_samp_8x16b_r8_0, coeff_c0_c1_8x16b);

    // a1c2+a2c3  a2c2+a3c3  a3c2+a4c3  a4c2+a5c3
    i4_res_4x32b_r1_1 = _mm_madd_epi16(i4_samp_8x16b_r1_1, coeff_c2_c3_8x16b);
    // b1c2+b2c3  b2c2+b3c3  b3c2+b4c3  b4c2+b5c3
    i4_res_4x32b_r2_1 = _mm_madd_epi16(i4_samp_8x16b_r2_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r3_1 = _mm_madd_epi16(i4_samp_8x16b_r3_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r4_1 = _mm_madd_epi16(i4_samp_8x16b_r4_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r5_1 = _mm_madd_epi16(i4_samp_8x16b_r5_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r6_1 = _mm_madd_epi16(i4_samp_8x16b_r6_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r7_1 = _mm_madd_epi16(i4_samp_8x16b_r7_1, coeff_c2_c3_8x16b);
    i4_res_4x32b_r8_1 = _mm_madd_epi16(i4_samp_8x16b_r8_1, coeff_c2_c3_8x16b);

    i4_res_4x32b_r1_0 = _mm_add_epi32(i4_res_4x32b_r1_0, res_128);
    i4_res_4x32b_r2_0 = _mm_add_epi32(i4_res_4x32b_r2_0, res_128);
    i4_res_4x32b_r3_0 = _mm_add_epi32(i4_res_4x32b_r3_0, res_128);
    i4_res_4x32b_r4_0 = _mm_add_epi32(i4_res_4x32b_r4_0, res_128);
    i4_res_4x32b_r5_0 = _mm_add_epi32(i4_res_4x32b_r5_0, res_128);
    i4_res_4x32b_r6_0 = _mm_add_epi32(i4_res_4x32b_r6_0, res_128);
    i4_res_4x32b_r7_0 = _mm_add_epi32(i4_res_4x32b_r7_0, res_128);
    i4_res_4x32b_r8_0 = _mm_add_epi32(i4_res_4x32b_r8_0, res_128);

    i4_res_4x32b_r1_1 = _mm_add_epi32(i4_res_4x32b_r1_1, res_128);
    i4_res_4x32b_r2_1 = _mm_add_epi32(i4_res_4x32b_r2_1, res_128);
    i4_res_4x32b_r3_1 = _mm_add_epi32(i4_res_4x32b_r3_1, res_128);
    i4_res_4x32b_r4_1 = _mm_add_epi32(i4_res_4x32b_r4_1, res_128);
    i4_res_4x32b_r5_1 = _mm_add_epi32(i4_res_4x32b_r5_1, res_128);
    i4_res_4x32b_r6_1 = _mm_add_epi32(i4_res_4x32b_r6_1, res_128);
    i4_res_4x32b_r7_1 = _mm_add_epi32(i4_res_4x32b_r7_1, res_128);
    i4_res_4x32b_r8_1 = _mm_add_epi32(i4_res_4x32b_r8_1, res_128);

    i4_res_4x32b_r1_0 = _mm_srai_epi32(i4_res_4x32b_r1_0, u4_norm_factor);
    i4_res_4x32b_r2_0 = _mm_srai_epi32(i4_res_4x32b_r2_0, u4_norm_factor);
    i4_res_4x32b_r3_0 = _mm_srai_epi32(i4_res_4x32b_r3_0, u4_norm_factor);
    i4_res_4x32b_r4_0 = _mm_srai_epi32(i4_res_4x32b_r4_0, u4_norm_factor);
    i4_res_4x32b_r5_0 = _mm_srai_epi32(i4_res_4x32b_r5_0, u4_norm_factor);
    i4_res_4x32b_r6_0 = _mm_srai_epi32(i4_res_4x32b_r6_0, u4_norm_factor);
    i4_res_4x32b_r7_0 = _mm_srai_epi32(i4_res_4x32b_r7_0, u4_norm_factor);
    i4_res_4x32b_r8_0 = _mm_srai_epi32(i4_res_4x32b_r8_0, u4_norm_factor);

    i4_res_4x32b_r1_1 = _mm_srai_epi32(i4_res_4x32b_r1_1, u4_norm_factor);
    i4_res_4x32b_r2_1 = _mm_srai_epi32(i4_res_4x32b_r2_1, u4_norm_factor);
    i4_res_4x32b_r3_1 = _mm_srai_epi32(i4_res_4x32b_r3_1, u4_norm_factor);
    i4_res_4x32b_r4_1 = _mm_srai_epi32(i4_res_4x32b_r4_1, u4_norm_factor);
    i4_res_4x32b_r5_1 = _mm_srai_epi32(i4_res_4x32b_r5_1, u4_norm_factor);
    i4_res_4x32b_r6_1 = _mm_srai_epi32(i4_res_4x32b_r6_1, u4_norm_factor);
    i4_res_4x32b_r7_1 = _mm_srai_epi32(i4_res_4x32b_r7_1, u4_norm_factor);
    i4_res_4x32b_r8_1 = _mm_srai_epi32(i4_res_4x32b_r8_1, u4_norm_factor);

    i4_res_final_8x16b_r12_0 = _mm_packs_epi32(i4_res_4x32b_r1_0, i4_res_4x32b_r2_0);
    i4_res_final_8x16b_r34_0 = _mm_packs_epi32(i4_res_4x32b_r3_0, i4_res_4x32b_r4_0);
    i4_res_final_8x16b_r56_0 = _mm_packs_epi32(i4_res_4x32b_r5_0, i4_res_4x32b_r6_0);
    i4_res_final_8x16b_r67_0 = _mm_packs_epi32(i4_res_4x32b_r7_0, i4_res_4x32b_r8_0);

    i4_res_final_8x16b_r12_1 = _mm_packs_epi32(i4_res_4x32b_r1_1, i4_res_4x32b_r2_1);
    i4_res_final_8x16b_r34_1 = _mm_packs_epi32(i4_res_4x32b_r3_1, i4_res_4x32b_r4_1);
    i4_res_final_8x16b_r56_1 = _mm_packs_epi32(i4_res_4x32b_r5_1, i4_res_4x32b_r6_1);
    i4_res_final_8x16b_r67_1 = _mm_packs_epi32(i4_res_4x32b_r7_1, i4_res_4x32b_r8_1);

    i4_res_final_8x16b_r1 = _mm_unpacklo_epi16(i4_res_final_8x16b_r12_0, i4_res_final_8x16b_r12_1);
    i4_res_final_8x16b_r2 = _mm_unpackhi_epi16(i4_res_final_8x16b_r12_0, i4_res_final_8x16b_r12_1);
    i4_res_final_8x16b_r3 = _mm_unpacklo_epi16(i4_res_final_8x16b_r34_0, i4_res_final_8x16b_r34_1);
    i4_res_final_8x16b_r4 = _mm_unpackhi_epi16(i4_res_final_8x16b_r34_0, i4_res_final_8x16b_r34_1);
    i4_res_final_8x16b_r5 = _mm_unpacklo_epi16(i4_res_final_8x16b_r56_0, i4_res_final_8x16b_r56_1);
    i4_res_final_8x16b_r6 = _mm_unpackhi_epi16(i4_res_final_8x16b_r56_0, i4_res_final_8x16b_r56_1);
    i4_res_final_8x16b_r7 = _mm_unpacklo_epi16(i4_res_final_8x16b_r67_0, i4_res_final_8x16b_r67_1);
    i4_res_final_8x16b_r8 = _mm_unpackhi_epi16(i4_res_final_8x16b_r67_0, i4_res_final_8x16b_r67_1);

    chroma_mask = _mm_set1_epi16(0xFF00);
    chroma_mask2 = _mm_set1_epi16(0x00FF);
    out_16x8b_r1 = _mm_loadu_si128((__m128i *) (&pu1_out[0]));
    out_16x8b_r2 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride]));
    out_16x8b_r3 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride2]));
    out_16x8b_r4 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride2 + i4_dst_stride]));
    out_16x8b_r5 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride4]));
    out_16x8b_r6 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride4 + i4_dst_stride]));
    out_16x8b_r7 = _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride4 + i4_dst_stride2]));
    out_16x8b_r8 =
        _mm_loadu_si128((__m128i *) (&pu1_out[i4_dst_stride4 + i4_dst_stride2 + i4_dst_stride]));

    out_16x8b_r1 = _mm_and_si128(out_16x8b_r1, chroma_mask);
    out_16x8b_r2 = _mm_and_si128(out_16x8b_r2, chroma_mask);
    out_16x8b_r3 = _mm_and_si128(out_16x8b_r3, chroma_mask);
    out_16x8b_r4 = _mm_and_si128(out_16x8b_r4, chroma_mask);
    out_16x8b_r5 = _mm_and_si128(out_16x8b_r5, chroma_mask);
    out_16x8b_r6 = _mm_and_si128(out_16x8b_r6, chroma_mask);
    out_16x8b_r7 = _mm_and_si128(out_16x8b_r7, chroma_mask);
    out_16x8b_r8 = _mm_and_si128(out_16x8b_r8, chroma_mask);

    i4_res_final_8x16b_r1 = _mm_and_si128(i4_res_final_8x16b_r1, chroma_mask2);
    i4_res_final_8x16b_r2 = _mm_and_si128(i4_res_final_8x16b_r2, chroma_mask2);
    i4_res_final_8x16b_r3 = _mm_and_si128(i4_res_final_8x16b_r3, chroma_mask2);
    i4_res_final_8x16b_r4 = _mm_and_si128(i4_res_final_8x16b_r4, chroma_mask2);
    i4_res_final_8x16b_r5 = _mm_and_si128(i4_res_final_8x16b_r5, chroma_mask2);
    i4_res_final_8x16b_r6 = _mm_and_si128(i4_res_final_8x16b_r6, chroma_mask2);
    i4_res_final_8x16b_r7 = _mm_and_si128(i4_res_final_8x16b_r7, chroma_mask2);
    i4_res_final_8x16b_r8 = _mm_and_si128(i4_res_final_8x16b_r8, chroma_mask2);

    out_16x8b_r1 = _mm_add_epi8(i4_res_final_8x16b_r1, out_16x8b_r1);
    out_16x8b_r2 = _mm_add_epi8(i4_res_final_8x16b_r2, out_16x8b_r2);
    out_16x8b_r3 = _mm_add_epi8(i4_res_final_8x16b_r3, out_16x8b_r3);
    out_16x8b_r4 = _mm_add_epi8(i4_res_final_8x16b_r4, out_16x8b_r4);
    out_16x8b_r5 = _mm_add_epi8(i4_res_final_8x16b_r5, out_16x8b_r5);
    out_16x8b_r6 = _mm_add_epi8(i4_res_final_8x16b_r6, out_16x8b_r6);
    out_16x8b_r7 = _mm_add_epi8(i4_res_final_8x16b_r7, out_16x8b_r7);
    out_16x8b_r8 = _mm_add_epi8(i4_res_final_8x16b_r8, out_16x8b_r8);

    _mm_storeu_si128((__m128i *) pu1_out, out_16x8b_r1);
    _mm_storeu_si128((__m128i *) (pu1_out + i4_dst_stride), out_16x8b_r2);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 2)), out_16x8b_r3);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 3)), out_16x8b_r4);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 4)), out_16x8b_r5);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 5)), out_16x8b_r6);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 6)), out_16x8b_r7);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 7)), out_16x8b_r8);
}
