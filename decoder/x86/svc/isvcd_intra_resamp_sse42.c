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
 *  isvcd_intra_resamp_sse42.c
 *
 * @brief
 *  Contains function definitions for intra resampling functions
 *
 * @author
 * Kishore
 *
 * @par List of Functions:
 *  - isvcd_interpolate_base_luma_dyadic_sse42
 *  - isvcd_vert_interpol_chroma_dyadic_1_sse42
 *  - isvcd_vert_interpol_chroma_dyadic_2_sse42
 *  - isvcd_vert_interpol_chroma_dyadic_3_sse42
 *  - isvcd_horz_interpol_chroma_dyadic_1_sse42
 *  - isvcd_horz_interpol_chroma_dyadic_2_sse42
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
/*  Function Name : isvcd_interpolate_base_luma_dyadic_sse42                  */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  intra resampling for dyadic scaling ratios               */
/*  Inputs        : pu1_inp_buf : ptr to the 12x12 reference sample buffer   */
/*                    pi2_tmp_filt_buf : ptr to the 12x16 buffer to hold the */
/*                        vertically interpolated data                       */
/*                  pu1_out_buf : output buffer pointer                      */
/*                  i4_out_stride : output buffer stride                     */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         05 21 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/

void isvcd_interpolate_base_luma_dyadic_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                              UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_x, i4_y;
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

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    pu1_out = pu1_out_buf;

    /* Vertical interpolation */
    /*First 64 bit */
    for(i4_x = 0; i4_x < 1; i4_x++)
    {
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

        _mm_storeu_si128((__m128i *) pi2_tmp, i4_res_8x16b_r1_1);
        pi2_tmp += i4_filt_stride;

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/

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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/

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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/

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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));

        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/
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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));
        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/
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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));
        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/

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

        i4_samp_8x16b_0 = i4_samp_8x16b_1;
        i4_samp_8x16b_1 = i4_samp_8x16b_2;
        i4_samp_8x16b_2 = i4_samp_8x16b_3;
        i4_samp_8x16b_3 = _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *) (pu1_inp)));
        /* y_phase is 4 for odd values of y */
        /* and 12 for even values of y        */
        /*Multiply by 8 =>  left shift by 3*/

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
    } /* End of loop over x */

    /*Remaining 32 bit */
    pu1_inp += 8;
    pi2_tmp += 8;
    for(i4_x = 0; i4_x < 1; i4_x++)
    {
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
            /* and 12 for even values of y        */
            /*Multiply by 8 =>  left shift by 3*/

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
    }

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
            // a0 a1 a2 a3 a4 a5 a6 a7
            i4_samp_8x16b_rpart1_0 = _mm_loadu_si128((__m128i *) pi2_tmp);
            // a4 a5 a6 a7 a8 a9 a10 a11
            i4_samp_8x16b_rpart2_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 4));
            // a1 a2 a3 a4 a5 a6 a7 0
            i4_samp_8x16b_rpart1_1 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 2);
            // a2 a3 a4 a5 a6 a7 0 0
            i4_samp_8x16b_rpart1_2 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 4);
            // a3 a4 a5 a6 a7 0 0 0
            i4_samp_8x16b_rpart1_3 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 6);
            // a4 a5 a6 a7 0 0 0 0
            i4_samp_8x16b_rpart1_4 = _mm_srli_si128(i4_samp_8x16b_rpart1_0, 8);

            // a5 a6 a7 a8 a9 a10 a11 0
            i4_samp_8x16b_rpart2_1 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 2);
            // a6 a7 a8 a9 a10 a11 0 0
            i4_samp_8x16b_rpart2_2 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 4);
            // a7 a8 a9 a10 a11 0 0 0
            i4_samp_8x16b_rpart2_3 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 6);
            // a8 a9 a10 a11 0 0 0 0
            i4_samp_8x16b_rpart2_4 = _mm_srli_si128(i4_samp_8x16b_rpart2_0, 8);
            // a0 a1  a1 a2  a2 a3  a3 a4
            i4_samp_8x16b_rpart1_0 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_0, i4_samp_8x16b_rpart1_1);
            // a1 a2  a2 a3  a3 a4  a4 a5
            i4_samp_8x16b_rpart1_1 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_1, i4_samp_8x16b_rpart1_2);
            // a2 a3  a3 a4  a4 a5  a5 a6
            i4_samp_8x16b_rpart1_2 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_2, i4_samp_8x16b_rpart1_3);
            // a3 a4  a4 a5  a5 a6  a6 a7
            i4_samp_8x16b_rpart1_3 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart1_3, i4_samp_8x16b_rpart1_4);
            // a4 a5  a5 a6  a6 a7  a7 a8
            i4_samp_8x16b_rpart2_0 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_0, i4_samp_8x16b_rpart2_1);
            // a5 a6  a6 a7  a7 a8  a8 a9
            i4_samp_8x16b_rpart2_1 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_1, i4_samp_8x16b_rpart2_2);
            // a6 a7  a7 a8  a8 a9  a9 a10
            i4_samp_8x16b_rpart2_2 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_2, i4_samp_8x16b_rpart2_3);
            // a7 a8  a8 a9  a9 a10 a10 a11
            i4_samp_8x16b_rpart2_3 =
                _mm_unpacklo_epi16(i4_samp_8x16b_rpart2_3, i4_samp_8x16b_rpart2_4);
            // a0c3+a1c2  a1c3+a2c2  a2c3+a3c2  a3c3+a4c2
            i4_res_4x32b_rpart1_0 = _mm_madd_epi16(i4_samp_8x16b_rpart1_0, coeff_c3_c2_8x16b);
            // a2c1+a3c0  a3c1+a4c0  a4c1+a5c0  a5c1+a6c0
            i4_res_4x32b_rpart1_2 = _mm_madd_epi16(i4_samp_8x16b_rpart1_2, coeff_c1_c0_8x16b);
            // a1c0+a2c1  a2c0+a3c1  a3c0+a4c1  a4c0+a5c1
            i4_res_4x32b_rpart1_1 = _mm_madd_epi16(i4_samp_8x16b_rpart1_1, coeff_c0_c1_8x16b);
            // a3c2+a4c3  a5c2+a5c3  a5c2+a6c3  a6c2+a7c3
            i4_res_4x32b_rpart1_3 = _mm_madd_epi16(i4_samp_8x16b_rpart1_3, coeff_c2_c3_8x16b);
            // a4c3+a5c2  a5a3+a6c2  a6c3+a7c2  a7c3+a8c2
            i4_res_4x32b_rpart2_0 = _mm_madd_epi16(i4_samp_8x16b_rpart2_0, coeff_c3_c2_8x16b);
            // a6c1+a7c0  a7c1+a8c0  a8c1+a9c0  a9c1+a10c0
            i4_res_4x32b_rpart2_2 = _mm_madd_epi16(i4_samp_8x16b_rpart2_2, coeff_c1_c0_8x16b);
            // a5c0+a6c1  a6c0+a7c1  a7c0+a8c1  a8c0+a9c1
            i4_res_4x32b_rpart2_1 = _mm_madd_epi16(i4_samp_8x16b_rpart2_1, coeff_c0_c1_8x16b);
            // a7c2+a8c3  a8c2+a9c3  a9c2+a10c3  a10c2+a11c3
            i4_res_4x32b_rpart2_3 = _mm_madd_epi16(i4_samp_8x16b_rpart2_3, coeff_c2_c3_8x16b);
            // a0c3+a1c2 + a2c1+a3c0  a1c3+a2c2 + a3c1+a4c0 a2c3+a3c2 + a4c1+a5c0
            // a3c3+a4c2 +a5c1+a6c0
            i4_res_4x32b_rpart1_0 = _mm_add_epi32(i4_res_4x32b_rpart1_0, i4_res_4x32b_rpart1_2);
            // a1c0+a2c1 + a3c2+a4c3  a2c0+a3c1 + a5c2+a5c3 a3c0+a4c1 + a5c2+a6c3
            // a4c0+a5c1 + a6c2+a7c3
            i4_res_4x32b_rpart1_1 = _mm_add_epi32(i4_res_4x32b_rpart1_1, i4_res_4x32b_rpart1_3);
            // a4c3+a5c2 + a6c1+a7c0  a5a3+a6c2 + a7c1+a8c0 a6c3+a7c2 + a8c1+a9c0
            // a7c3+a8c2+ a9c1+a10c0
            i4_res_4x32b_rpart2_0 = _mm_add_epi32(i4_res_4x32b_rpart2_0, i4_res_4x32b_rpart2_2);
            // a5c0+a6c1 + a7c2+a8c3  a6c0+a7c1 + a8c2+a9c3 a7c0+a8c1 + a9c2+a10c3
            // a8c0+a9c1 + a10c2+a11c3
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
} /* isvcd_interpolate_base_luma_dyadic */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_intra_base_sse42                        */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                    interpolation of a component to find the intra         */
/*                     resampled value                                       */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  pu1_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore              creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_intra_base_sse42(void *pv_intra_samp_ctxt, UWORD8 *pu1_out,
                                        WORD32 i4_out_stride, WORD32 i4_refarray_wd, WORD32 i4_mb_x,
                                        WORD32 i4_mb_y, WORD32 i4_chroma_flag,
                                        WORD32 i4_refarray_flag)
{
    /* --------------------------------------------------------------------- */
    /* Index Parameters                                                         */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_map_ctxt_t *ps_map_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    UWORD8 *pu1_refarray = NULL;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;
    WORD32 i4_temp_array_ht;
    WORD32 *pi4_interp_buff;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;

    WORD32 i4_x_min;
    ref_min_max_map_t *ps_x_min_max;
    WORD8 arr_y_ref_pos_luma[16] = {0};
    WORD8 arr_x_ref_pos_luma[16] = {0};
    WORD8 arr_x_ref_pos_luma_low[16] = {0};
    WORD8 arr_x_ref_pos_luma_high[16] = {0};
    WORD8 arr_phase_luma[32] = {0};
    WORD8 *pi4_y_ref_pos_luma;
    WORD8 *pi4_x_ref_pos_luma_low;
    WORD8 *pi4_x_ref_pos_luma_high;
    WORD8 *pi4_phase_luma;
    UWORD8 *pu1_refarray_temp;

    /* --------------------------------------------------------------------- */
    /* Extracting pointers from the  context                                  */
    /* --------------------------------------------------------------------- */
    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    if(0 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_buffer;
    }
    else if(1 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_cb;
    }

    /* --------------------------------------------------------------------- */
    /* LUMA    or CHROMA */
    /* --------------------------------------------------------------------- */

    if(1 == i4_chroma_flag)
        ps_map_ctxt = &(ps_lyr_ctxt->s_chroma_map_ctxt);
    else
        ps_map_ctxt = &(ps_lyr_ctxt->s_luma_map_ctxt);

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    ps_x_min_max = ps_map_ctxt->ps_x_min_max;

    i4_frm_mb_y = i4_mb_y * i4_mb_ht;
    i4_frm_mb_x = i4_mb_x * i4_mb_wd;
    /* get the min position */
    i4_x_min = ps_x_min_max[i4_mb_x].i2_min_pos;

    /* --------------------------------------------------------------------- */
    /* Projected frame level pointers                                        */
    /* --------------------------------------------------------------------- */
    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

    /* --------------------------------------------------------------------- */
    /* Pointers and Dimenstion of the temporary buffer                         */
    /* --------------------------------------------------------------------- */
    i4_temp_array_ht = i4_mb_ht;
    pi4_interp_buff = ps_ctxt->pi4_temp_interpolation_buffer;

    if(i4_chroma_flag == 0)
    {
        /* --------------------------------------------------------------------- */
        /* Loop for interpolation in vertical direction */
        /* --------------------------------------------------------------------- */
        WORD16 *pi2_interp_buff_temp;
        pi2_interp_buff_temp = (WORD16 *) pi4_interp_buff;
        {
            __m128i out_res_8x16b_0, out_res_8x16b_1;

            __m128i inp_8x16b_r0, inp_8x16b_r01_0, phs_mask_16x8b_r0, phs_mask_16x8b_r01_0,
                out_res_8x16b_r01_0;
            __m128i inp_8x16b_r1, inp_8x16b_r23_0, phs_mask_16x8b_r1, phs_mask_16x8b_r23_0,
                out_res_8x16b_r01_1;
            __m128i inp_8x16b_r2, inp_8x16b_r01_1, phs_mask_16x8b_r2, phs_mask_16x8b_r01_1,
                out_res_8x16b_r23_0;
            __m128i inp_8x16b_r3, inp_8x16b_r23_1, phs_mask_16x8b_r3, phs_mask_16x8b_r23_1,
                out_res_8x16b_r23_1;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                arr_phase_luma[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
                arr_y_ref_pos_luma[i4_y] = (WORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
            }
            pi4_y_ref_pos_luma = arr_y_ref_pos_luma;
            pi4_phase_luma = arr_phase_luma;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                pu1_refarray_temp =
                    pu1_refarray + (pi4_y_ref_pos_luma[i4_y] * i4_refarray_wd) + (i4_x_min - 1);
                inp_8x16b_r0 = _mm_loadu_si128((__m128i *) (pu1_refarray_temp - i4_refarray_wd));
                inp_8x16b_r1 = _mm_loadu_si128((__m128i *) (pu1_refarray_temp));
                inp_8x16b_r2 = _mm_loadu_si128((__m128i *) (pu1_refarray_temp + i4_refarray_wd));
                inp_8x16b_r3 =
                    _mm_loadu_si128((__m128i *) (pu1_refarray_temp + 2 * i4_refarray_wd));

                inp_8x16b_r01_0 = _mm_unpacklo_epi8(inp_8x16b_r0, inp_8x16b_r1);
                inp_8x16b_r23_0 = _mm_unpacklo_epi8(inp_8x16b_r2, inp_8x16b_r3);
                inp_8x16b_r01_1 = _mm_unpackhi_epi8(inp_8x16b_r0, inp_8x16b_r1);
                inp_8x16b_r23_1 = _mm_unpackhi_epi8(inp_8x16b_r2, inp_8x16b_r3);

                phs_mask_16x8b_r0 = _mm_set1_epi8(g_ai1_interp_filter_luma[pi4_phase_luma[i4_y]]);
                phs_mask_16x8b_r1 =
                    _mm_set1_epi8(g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 16]);
                phs_mask_16x8b_r2 =
                    _mm_set1_epi8(g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 32]);
                phs_mask_16x8b_r3 =
                    _mm_set1_epi8(g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 48]);

                phs_mask_16x8b_r01_0 = _mm_unpacklo_epi8(phs_mask_16x8b_r0, phs_mask_16x8b_r1);
                phs_mask_16x8b_r23_0 = _mm_unpacklo_epi8(phs_mask_16x8b_r2, phs_mask_16x8b_r3);
                phs_mask_16x8b_r01_1 = _mm_unpackhi_epi8(phs_mask_16x8b_r0, phs_mask_16x8b_r1);
                phs_mask_16x8b_r23_1 = _mm_unpackhi_epi8(phs_mask_16x8b_r2, phs_mask_16x8b_r3);

                out_res_8x16b_r01_0 = _mm_maddubs_epi16(inp_8x16b_r01_0, phs_mask_16x8b_r01_0);
                out_res_8x16b_r01_1 = _mm_maddubs_epi16(inp_8x16b_r01_1, phs_mask_16x8b_r01_1);
                out_res_8x16b_r23_0 = _mm_maddubs_epi16(inp_8x16b_r23_0, phs_mask_16x8b_r23_0);
                out_res_8x16b_r23_1 = _mm_maddubs_epi16(inp_8x16b_r23_1, phs_mask_16x8b_r23_1);

                out_res_8x16b_0 = _mm_add_epi16(out_res_8x16b_r01_0, out_res_8x16b_r23_0);
                out_res_8x16b_1 = _mm_add_epi16(out_res_8x16b_r01_1, out_res_8x16b_r23_1);

                _mm_storeu_si128(
                    (__m128i *) (pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1)),
                    out_res_8x16b_0);
                _mm_storeu_si128((__m128i *) (pi2_interp_buff_temp + (i4_y * i4_refarray_wd) +
                                              (i4_x_min - 1) + 8),
                                 out_res_8x16b_1);
            }
        }
        /* --------------------------------------------------------------------- */
        /* Loop for interpolation in horizontal direction                         */
        /* --------------------------------------------------------------------- */
        {
            WORD32 strt_indx = 10, strt_indx_h = 0;

            __m128i inp_8x16b_0;
            __m128i inp_8x16b_1;

            __m128i phs_mask_16x8b_0;
            __m128i phs_mask_16x8b_1;
            __m128i x_ref_pos_luma_mask_r0_0, x_ref_pos_luma_mask_r0_1, x_ref_pos_luma_mask_r1_0,
                x_ref_pos_luma_mask_r1_1, x_ref_pos_luma_mask_r2_0, x_ref_pos_luma_mask_r2_1,
                x_ref_pos_luma_mask_r3_0, x_ref_pos_luma_mask_r3_1;

            __m128i inp_8x16b_2, inp_8x16b_3;

            WORD32 i4_x2 = 0;
            WORD32 i4_mb_wd_hlf = (i4_mb_wd >> 1);
            __m128i twos = _mm_set1_epi8(2);

            strt_indx = ps_x_pos_phase[0 + i4_frm_mb_x].i2_ref_pos - 1;
            strt_indx_h = (ps_x_pos_phase[8 + i4_frm_mb_x].i2_ref_pos - strt_indx - 1);
            for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
            {
                arr_x_ref_pos_luma[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
                arr_phase_luma[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
                arr_x_ref_pos_luma[i4_x] = arr_x_ref_pos_luma[i4_x] - strt_indx - 1;
            }

            for(i4_x = 0; i4_x < i4_mb_wd_hlf; i4_x++)
            {
                i4_x2 = i4_x << 1;
                arr_x_ref_pos_luma_low[i4_x2] = (arr_x_ref_pos_luma[i4_x]) << 1;
                arr_x_ref_pos_luma_low[i4_x2 + 1] = arr_x_ref_pos_luma_low[i4_x2] + 1;
            }
            for(i4_x = i4_mb_wd_hlf; i4_x < i4_mb_wd; i4_x++)
            {
                i4_x2 = (i4_x - i4_mb_wd_hlf) << 1;
                arr_x_ref_pos_luma_high[i4_x2] = ((arr_x_ref_pos_luma[i4_x] - strt_indx_h) << 1);
                arr_x_ref_pos_luma_high[i4_x2 + 1] = arr_x_ref_pos_luma_high[i4_x2] + 1;
            }
            pi4_x_ref_pos_luma_low = arr_x_ref_pos_luma_low;
            pi4_x_ref_pos_luma_high = arr_x_ref_pos_luma_high;
            pi4_phase_luma = arr_phase_luma;

            phs_mask_16x8b_0 = _mm_loadu_si128((__m128i *) (pi4_phase_luma));
            phs_mask_16x8b_1 = _mm_loadu_si128((__m128i *) (pi4_phase_luma + 8));

            x_ref_pos_luma_mask_r0_0 = _mm_loadu_si128((__m128i *) (pi4_x_ref_pos_luma_low));
            x_ref_pos_luma_mask_r0_1 = _mm_loadu_si128((__m128i *) (pi4_x_ref_pos_luma_high));
            x_ref_pos_luma_mask_r1_0 = _mm_add_epi8(x_ref_pos_luma_mask_r0_0, twos);
            x_ref_pos_luma_mask_r1_1 = _mm_add_epi8(x_ref_pos_luma_mask_r0_1, twos);
            x_ref_pos_luma_mask_r2_0 = _mm_add_epi8(x_ref_pos_luma_mask_r1_0, twos);
            x_ref_pos_luma_mask_r2_1 = _mm_add_epi8(x_ref_pos_luma_mask_r1_1, twos);
            x_ref_pos_luma_mask_r3_0 = x_ref_pos_luma_mask_r0_0;
            x_ref_pos_luma_mask_r3_1 = x_ref_pos_luma_mask_r0_1;

            {
                __m128i ip_filt_16x8b_r0, ip_filt_8x16b_r0_0, ip_filt_8x16b_r0_1,
                    ip_filt_8x16b_r01_l_0, ip_filt_8x16b_r01_h_0;
                __m128i ip_filt_16x8b_r1, ip_filt_8x16b_r1_0, ip_filt_8x16b_r1_1,
                    ip_filt_8x16b_r23_l_0, ip_filt_8x16b_r23_h_0;
                __m128i ip_filt_16x8b_r2, ip_filt_8x16b_r2_0, ip_filt_8x16b_r2_1,
                    ip_filt_8x16b_r01_l_1, ip_filt_8x16b_r01_h_1;
                __m128i ip_filt_16x8b_r3, ip_filt_8x16b_r3_0, ip_filt_8x16b_r3_1,
                    ip_filt_8x16b_r23_l_1, ip_filt_8x16b_r23_h_1;

                __m128i inp_8x16b_r0_0, inp_8x16b_r2_0, inp_8x16b_r01_l_0, inp_8x16b_r01_h_0,
                    out_res_4x32b_r01_l_0, out_res_4x32b_r01_h_0;
                __m128i inp_8x16b_r0_1, inp_8x16b_r2_1, inp_8x16b_r23_l_0, inp_8x16b_r23_h_0,
                    out_res_4x32b_r01_l_1, out_res_4x32b_r01_h_1;
                __m128i inp_8x16b_r1_0, inp_8x16b_r3_0, inp_8x16b_r01_l_1, inp_8x16b_r01_h_1,
                    out_res_4x32b_r23_l_0, out_res_4x32b_r23_h_0;
                __m128i inp_8x16b_r1_1, inp_8x16b_r3_1, inp_8x16b_r23_l_1, inp_8x16b_r23_h_1,
                    out_res_4x32b_r23_l_1, out_res_4x32b_r23_h_1;

                __m128i out_res_4x32b_l_0;
                __m128i out_res_4x32b_l_1;
                __m128i out_res_4x32b_h_0;
                __m128i out_res_4x32b_h_1;

                __m128i out_res_8x16b_l;
                __m128i out_res_8x16b_h;

                __m128i out_res_16x8b;
                __m128i const_512 = _mm_set1_epi32(512);

                ip_filt_16x8b_r0 = _mm_loadu_si128((__m128i *) (g_ai1_interp_filter_luma));
                ip_filt_16x8b_r1 = _mm_loadu_si128((__m128i *) (g_ai1_interp_filter_luma + 16));
                ip_filt_16x8b_r2 = _mm_loadu_si128((__m128i *) (g_ai1_interp_filter_luma + 32));
                ip_filt_16x8b_r3 = _mm_loadu_si128((__m128i *) (g_ai1_interp_filter_luma + 48));

                ip_filt_8x16b_r0_0 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r0, phs_mask_16x8b_0));
                ip_filt_8x16b_r1_0 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r1, phs_mask_16x8b_0));
                ip_filt_8x16b_r2_0 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r2, phs_mask_16x8b_0));
                ip_filt_8x16b_r3_0 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r3, phs_mask_16x8b_0));

                ip_filt_8x16b_r0_1 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r0, phs_mask_16x8b_1));
                ip_filt_8x16b_r1_1 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r1, phs_mask_16x8b_1));
                ip_filt_8x16b_r2_1 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r2, phs_mask_16x8b_1));
                ip_filt_8x16b_r3_1 =
                    _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r3, phs_mask_16x8b_1));

                ip_filt_8x16b_r01_l_0 = _mm_unpacklo_epi16(ip_filt_8x16b_r0_0, ip_filt_8x16b_r1_0);
                ip_filt_8x16b_r23_l_0 = _mm_unpacklo_epi16(ip_filt_8x16b_r2_0, ip_filt_8x16b_r3_0);
                ip_filt_8x16b_r01_l_1 = _mm_unpackhi_epi16(ip_filt_8x16b_r0_0, ip_filt_8x16b_r1_0);
                ip_filt_8x16b_r23_l_1 = _mm_unpackhi_epi16(ip_filt_8x16b_r2_0, ip_filt_8x16b_r3_0);

                ip_filt_8x16b_r01_h_0 = _mm_unpacklo_epi16(ip_filt_8x16b_r0_1, ip_filt_8x16b_r1_1);
                ip_filt_8x16b_r23_h_0 = _mm_unpacklo_epi16(ip_filt_8x16b_r2_1, ip_filt_8x16b_r3_1);
                ip_filt_8x16b_r01_h_1 = _mm_unpackhi_epi16(ip_filt_8x16b_r0_1, ip_filt_8x16b_r1_1);
                ip_filt_8x16b_r23_h_1 = _mm_unpackhi_epi16(ip_filt_8x16b_r2_1, ip_filt_8x16b_r3_1);

                for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
                {
                    inp_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_interp_buff_temp + strt_indx));
                    inp_8x16b_1 = _mm_loadu_si128(
                        (__m128i *) (pi2_interp_buff_temp + strt_indx + strt_indx_h));
                    inp_8x16b_2 =
                        _mm_loadu_si128((__m128i *) (pi2_interp_buff_temp + strt_indx + 3));
                    inp_8x16b_3 = _mm_loadu_si128(
                        (__m128i *) (pi2_interp_buff_temp + strt_indx + strt_indx_h + 3));
                    pi2_interp_buff_temp += i4_refarray_wd;

                    inp_8x16b_r0_0 = _mm_shuffle_epi8(inp_8x16b_0, x_ref_pos_luma_mask_r0_0);
                    inp_8x16b_r0_1 = _mm_shuffle_epi8(inp_8x16b_1, x_ref_pos_luma_mask_r0_1);
                    inp_8x16b_r1_0 = _mm_shuffle_epi8(inp_8x16b_0, x_ref_pos_luma_mask_r1_0);
                    inp_8x16b_r1_1 = _mm_shuffle_epi8(inp_8x16b_1, x_ref_pos_luma_mask_r1_1);

                    inp_8x16b_r2_0 = _mm_shuffle_epi8(inp_8x16b_0, x_ref_pos_luma_mask_r2_0);
                    inp_8x16b_r2_1 = _mm_shuffle_epi8(inp_8x16b_1, x_ref_pos_luma_mask_r2_1);
                    inp_8x16b_r3_0 = _mm_shuffle_epi8(inp_8x16b_2, x_ref_pos_luma_mask_r3_0);
                    inp_8x16b_r3_1 = _mm_shuffle_epi8(inp_8x16b_3, x_ref_pos_luma_mask_r3_1);

                    inp_8x16b_r01_l_0 = _mm_unpacklo_epi16(inp_8x16b_r0_0, inp_8x16b_r1_0);
                    inp_8x16b_r23_l_0 = _mm_unpacklo_epi16(inp_8x16b_r2_0, inp_8x16b_r3_0);
                    inp_8x16b_r01_l_1 = _mm_unpackhi_epi16(inp_8x16b_r0_0, inp_8x16b_r1_0);
                    inp_8x16b_r23_l_1 = _mm_unpackhi_epi16(inp_8x16b_r2_0, inp_8x16b_r3_0);

                    inp_8x16b_r01_h_0 = _mm_unpacklo_epi16(inp_8x16b_r0_1, inp_8x16b_r1_1);
                    inp_8x16b_r23_h_0 = _mm_unpacklo_epi16(inp_8x16b_r2_1, inp_8x16b_r3_1);
                    inp_8x16b_r01_h_1 = _mm_unpackhi_epi16(inp_8x16b_r0_1, inp_8x16b_r1_1);
                    inp_8x16b_r23_h_1 = _mm_unpackhi_epi16(inp_8x16b_r2_1, inp_8x16b_r3_1);

                    out_res_4x32b_r01_l_0 =
                        _mm_madd_epi16(inp_8x16b_r01_l_0, ip_filt_8x16b_r01_l_0);
                    out_res_4x32b_r01_l_1 =
                        _mm_madd_epi16(inp_8x16b_r01_l_1, ip_filt_8x16b_r01_l_1);
                    out_res_4x32b_r23_l_0 =
                        _mm_madd_epi16(inp_8x16b_r23_l_0, ip_filt_8x16b_r23_l_0);
                    out_res_4x32b_r23_l_1 =
                        _mm_madd_epi16(inp_8x16b_r23_l_1, ip_filt_8x16b_r23_l_1);

                    out_res_4x32b_r01_h_0 =
                        _mm_madd_epi16(inp_8x16b_r01_h_0, ip_filt_8x16b_r01_h_0);
                    out_res_4x32b_r01_h_1 =
                        _mm_madd_epi16(inp_8x16b_r01_h_1, ip_filt_8x16b_r01_h_1);
                    out_res_4x32b_r23_h_0 =
                        _mm_madd_epi16(inp_8x16b_r23_h_0, ip_filt_8x16b_r23_h_0);
                    out_res_4x32b_r23_h_1 =
                        _mm_madd_epi16(inp_8x16b_r23_h_1, ip_filt_8x16b_r23_h_1);

                    out_res_4x32b_l_0 = _mm_add_epi32(out_res_4x32b_r01_l_0, out_res_4x32b_r23_l_0);
                    out_res_4x32b_l_1 = _mm_add_epi32(out_res_4x32b_r01_l_1, out_res_4x32b_r23_l_1);
                    out_res_4x32b_h_0 = _mm_add_epi32(out_res_4x32b_r01_h_0, out_res_4x32b_r23_h_0);
                    out_res_4x32b_h_1 = _mm_add_epi32(out_res_4x32b_r01_h_1, out_res_4x32b_r23_h_1);

                    out_res_4x32b_l_0 =
                        _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_l_0, const_512), 10);
                    out_res_4x32b_l_1 =
                        _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_l_1, const_512), 10);
                    out_res_4x32b_h_0 =
                        _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_h_0, const_512), 10);
                    out_res_4x32b_h_1 =
                        _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_h_1, const_512), 10);

                    out_res_8x16b_l = _mm_packs_epi32(out_res_4x32b_l_0, out_res_4x32b_l_1);
                    out_res_8x16b_h = _mm_packs_epi32(out_res_4x32b_h_0, out_res_4x32b_h_1);

                    out_res_16x8b = _mm_packus_epi16(out_res_8x16b_l, out_res_8x16b_h);
                    _mm_storeu_si128((__m128i *) (pu1_out + (i4_y * i4_out_stride)), out_res_16x8b);
                }
            }
        }
    }
    else
    {
        WORD16 *pi2_interp_buff_temp;
        pi2_interp_buff_temp = (WORD16 *) pi4_interp_buff;

        {
            __m128i inp_8x16b_r0, inp_8x16b_r01_0, phs_mask_16x8b_r0, phs_mask_16x8b_r01_0,
                out_res_8x16b_r01_0;
            __m128i inp_8x16b_r1, phs_mask_16x8b_r1, out_res_8x16b_r01_1;
            __m128i inp_8x16b_r01_1, phs_mask_16x8b_r01_1;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                arr_y_ref_pos_luma[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos;
                arr_phase_luma[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            }
            pi4_y_ref_pos_luma = arr_y_ref_pos_luma;
            pi4_phase_luma = arr_phase_luma;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                pu1_refarray_temp =
                    pu1_refarray + (pi4_y_ref_pos_luma[i4_y] * i4_refarray_wd) + (i4_x_min - 1);
                inp_8x16b_r0 = _mm_loadu_si128((__m128i *) (pu1_refarray_temp));
                inp_8x16b_r1 = _mm_loadu_si128((__m128i *) (pu1_refarray_temp + i4_refarray_wd));

                inp_8x16b_r01_0 = _mm_unpacklo_epi8(inp_8x16b_r0, inp_8x16b_r1);
                inp_8x16b_r01_1 = _mm_unpackhi_epi8(inp_8x16b_r0, inp_8x16b_r1);

                phs_mask_16x8b_r0 = _mm_set1_epi8(g_au1_interp_filter_chroma[pi4_phase_luma[i4_y]]);
                phs_mask_16x8b_r1 =
                    _mm_set1_epi8(g_au1_interp_filter_chroma[pi4_phase_luma[i4_y] + 16]);

                phs_mask_16x8b_r01_0 = _mm_unpacklo_epi8(phs_mask_16x8b_r0, phs_mask_16x8b_r1);
                phs_mask_16x8b_r01_1 = _mm_unpackhi_epi8(phs_mask_16x8b_r0, phs_mask_16x8b_r1);

                out_res_8x16b_r01_0 = _mm_maddubs_epi16(inp_8x16b_r01_0, phs_mask_16x8b_r01_0);
                out_res_8x16b_r01_1 = _mm_maddubs_epi16(inp_8x16b_r01_1, phs_mask_16x8b_r01_1);

                _mm_storeu_si128(
                    (__m128i *) (pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1)),
                    out_res_8x16b_r01_0);
                _mm_storeu_si128((__m128i *) (pi2_interp_buff_temp + (i4_y * i4_refarray_wd) +
                                              (i4_x_min - 1) + 8),
                                 out_res_8x16b_r01_1);
            }
        }

        {
            WORD32 strt_indx = 10;
            __m128i inp_8x16b_0, inp_8x16b_r0_0;
            __m128i phs_mask_16x8b_0;
            __m128i x_ref_pos_luma_mask_r0_0, x_ref_pos_luma_mask_r1_0;
            __m128i ip_filt_16x8b_r0, ip_filt_8x16b_r0_0, ip_filt_8x16b_r01_l_0;
            __m128i ip_filt_16x8b_r1, ip_filt_8x16b_r1_0, ip_filt_8x16b_r01_l_1;
            __m128i inp_8x16b_r1_0, inp_8x16b_r01_l_0, out_res_4x32b_r01_l_0;
            __m128i inp_8x16b_r01_l_1, out_res_4x32b_r01_l_1;

            __m128i out_res_4x32b_l_0;
            __m128i out_res_4x32b_l_1;
            __m128i out_res_8x16b_l;
            __m128i out_16x8b_r1;
            __m128i chroma_mask;
            __m128i const_512 = _mm_set1_epi32(512);

            WORD32 i4_x2 = 0;
            __m128i twos = _mm_set1_epi8(2);
            strt_indx = ps_x_pos_phase[0 + i4_frm_mb_x].i2_ref_pos;
            for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
            {
                arr_x_ref_pos_luma[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
                arr_phase_luma[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
                arr_x_ref_pos_luma[i4_x] = arr_x_ref_pos_luma[i4_x] - strt_indx;
                i4_x2 = i4_x << 1;
                arr_x_ref_pos_luma_low[i4_x2] = (arr_x_ref_pos_luma[i4_x]) << 1;
                arr_x_ref_pos_luma_low[i4_x2 + 1] = arr_x_ref_pos_luma_low[i4_x2] + 1;
            }

            pi4_x_ref_pos_luma_low = arr_x_ref_pos_luma_low;
            pi4_phase_luma = arr_phase_luma;
            phs_mask_16x8b_0 = _mm_loadu_si128((__m128i *) (pi4_phase_luma));
            x_ref_pos_luma_mask_r0_0 = _mm_loadu_si128((__m128i *) (pi4_x_ref_pos_luma_low));
            x_ref_pos_luma_mask_r1_0 = _mm_add_epi8(x_ref_pos_luma_mask_r0_0, twos);

            ip_filt_16x8b_r0 = _mm_loadu_si128((__m128i *) (g_au1_interp_filter_chroma));
            ip_filt_16x8b_r1 = _mm_loadu_si128((__m128i *) (g_au1_interp_filter_chroma + 16));

            ip_filt_8x16b_r0_0 =
                _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r0, phs_mask_16x8b_0));
            ip_filt_8x16b_r1_0 =
                _mm_cvtepi8_epi16(_mm_shuffle_epi8(ip_filt_16x8b_r1, phs_mask_16x8b_0));

            ip_filt_8x16b_r01_l_0 = _mm_unpacklo_epi16(ip_filt_8x16b_r0_0, ip_filt_8x16b_r1_0);
            ip_filt_8x16b_r01_l_1 = _mm_unpackhi_epi16(ip_filt_8x16b_r0_0, ip_filt_8x16b_r1_0);

            for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
            {
                inp_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_interp_buff_temp + strt_indx));
                pi2_interp_buff_temp += i4_refarray_wd;

                inp_8x16b_r0_0 = _mm_shuffle_epi8(inp_8x16b_0, x_ref_pos_luma_mask_r0_0);
                inp_8x16b_r1_0 = _mm_shuffle_epi8(inp_8x16b_0, x_ref_pos_luma_mask_r1_0);

                inp_8x16b_r01_l_0 = _mm_unpacklo_epi16(inp_8x16b_r0_0, inp_8x16b_r1_0);
                inp_8x16b_r01_l_1 = _mm_unpackhi_epi16(inp_8x16b_r0_0, inp_8x16b_r1_0);

                out_res_4x32b_r01_l_0 = _mm_madd_epi16(inp_8x16b_r01_l_0, ip_filt_8x16b_r01_l_0);
                out_res_4x32b_r01_l_1 = _mm_madd_epi16(inp_8x16b_r01_l_1, ip_filt_8x16b_r01_l_1);

                out_res_4x32b_l_0 =
                    _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_r01_l_0, const_512), 10);
                out_res_4x32b_l_1 =
                    _mm_srai_epi32(_mm_add_epi32(out_res_4x32b_r01_l_1, const_512), 10);

                out_res_8x16b_l = _mm_packs_epi32(out_res_4x32b_l_0, out_res_4x32b_l_1);

                chroma_mask = _mm_set1_epi16(0xFF00);
                out_16x8b_r1 = _mm_loadu_si128((__m128i *) (pu1_out + (i4_y * i4_out_stride)));
                out_16x8b_r1 = _mm_and_si128(out_16x8b_r1, chroma_mask);
                out_16x8b_r1 = _mm_add_epi8(out_res_8x16b_l, out_16x8b_r1);
                _mm_storeu_si128((__m128i *) (pu1_out + (i4_y * i4_out_stride)), out_16x8b_r1);
            }
        }
    }
    return;
} /* End of Interpolation Function */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_1_sse42                 */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                    interpolation of a component to find the intra         */
/*                     resampled value                                       */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  pu1_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation on horizontal direction        */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore              creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_1_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
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

    i4_coeff_0 = (WORD8) (8 - i4_phase_0);
    i4_coeff_1 = (WORD8) (i4_phase_0);
    i4_coeff_2 = (WORD8) (8 - i4_phase_1);
    i4_coeff_3 = (WORD8) (i4_phase_1);

    i4_c0_c1_16x8b =
        _mm_set_epi8(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0);
    i4_c2_c3_16x8b =
        _mm_set_epi8(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2);

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;

    i4_samp_16x8b_0 = _mm_loadl_epi64((__m128i *) (pu1_inp));
    i4_samp_16x8b_1 = _mm_loadl_epi64((__m128i *) (pu1_inp + i4_src_stride));
    i4_samp_16x8b_2 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1)));
    i4_samp_16x8b_3 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
    i4_samp_16x8b_4 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 2)));
    i4_samp_16x8b_5 = _mm_loadl_epi64((__m128i *) (pu1_inp + (i4_src_stride << 2) + i4_src_stride));

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

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_2_sse42                 */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                    chroma_phase_y_plus1:                                  */
/*                        ref_lyr        cur_lyr                             */
/*                            0            1                                 */
/*                            0            2                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_2_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                               WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD8 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;
    __m128i i4_samp_16x8b_0, i4_samp_16x8b_1, i4_samp_16x8b_2, i4_samp_16x8b_3, i4_samp_16x8b_4;
    __m128i i4_res_8x16b_r0, i4_res_8x16b_r1, i4_res_8x16b_r2, i4_res_8x16b_r3, i4_res_8x16b_r4,
        i4_res_8x16b_r5, i4_res_8x16b_r6, i4_res_8x16b_r7;
    __m128i i4_res_8x16b_r7_temp, i4_c0_c1_16x8b, i4_c2_c3_16x8b;
    i4_coeff_0 = (WORD8) (8 - i4_phase_0);
    i4_coeff_1 = (WORD8) (i4_phase_0);
    i4_coeff_2 = (WORD8) (8 - i4_phase_1);
    i4_coeff_3 = (WORD8) (i4_phase_1);

    i4_c0_c1_16x8b =
        _mm_set_epi8(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0);
    i4_c2_c3_16x8b =
        _mm_set_epi8(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2);

    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf + i4_src_stride;

    i4_samp_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_inp));
    i4_samp_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_inp + i4_src_stride));
    i4_samp_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1)));
    i4_samp_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
    i4_samp_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 2)));

    i4_samp_16x8b_0 = _mm_unpacklo_epi8(i4_samp_16x8b_0, i4_samp_16x8b_1);
    i4_res_8x16b_r0 = _mm_maddubs_epi16(i4_samp_16x8b_0, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp), i4_res_8x16b_r0);

    i4_res_8x16b_r1 = _mm_maddubs_epi16(i4_samp_16x8b_0, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + i4_filt_stride), i4_res_8x16b_r1);

    i4_samp_16x8b_1 = _mm_unpacklo_epi8(i4_samp_16x8b_1, i4_samp_16x8b_2);
    i4_res_8x16b_r2 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1)), i4_res_8x16b_r2);

    i4_res_8x16b_r3 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1) + i4_filt_stride),
                     i4_res_8x16b_r3);

    i4_samp_16x8b_2 = _mm_unpacklo_epi8(i4_samp_16x8b_2, i4_samp_16x8b_3);
    i4_res_8x16b_r4 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2)), i4_res_8x16b_r4);

    i4_res_8x16b_r5 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + i4_filt_stride),
                     i4_res_8x16b_r5);

    i4_samp_16x8b_3 = _mm_unpacklo_epi8(i4_samp_16x8b_3, i4_samp_16x8b_4);
    i4_res_8x16b_r6 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c0_c1_16x8b);
    _mm_storel_epi64((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1)),
                     i4_res_8x16b_r6);

    i4_res_8x16b_r7 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c2_c3_16x8b);
    i4_res_8x16b_r6 = _mm_shuffle_epi32(i4_res_8x16b_r6, 78);
    i4_res_8x16b_r7 = _mm_shuffle_epi32(i4_res_8x16b_r7, 147);
    i4_res_8x16b_r7_temp = _mm_blend_epi16(i4_res_8x16b_r6, i4_res_8x16b_r7, 252);

    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1) + 4),
                     i4_res_8x16b_r7_temp);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_3_sse42                 */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_3_sse42(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                               WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD8 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;
    __m128i i4_samp_16x8b_0, i4_samp_16x8b_1, i4_samp_16x8b_2, i4_samp_16x8b_3, i4_samp_16x8b_4;
    __m128i i4_res_8x16b_r0, i4_res_8x16b_r1, i4_res_8x16b_r2, i4_res_8x16b_r3, i4_res_8x16b_r4,
        i4_res_8x16b_r5, i4_res_8x16b_r6, i4_res_8x16b_r7;
    __m128i i4_res_8x16b_r7_temp, i4_c0_c1_16x8b, i4_c2_c3_16x8b;
    i4_coeff_0 = (WORD8) (8 - i4_phase_0);
    i4_coeff_1 = (WORD8) (i4_phase_0);
    i4_coeff_2 = (WORD8) (8 - i4_phase_1);
    i4_coeff_3 = (WORD8) (i4_phase_1);

    i4_c0_c1_16x8b =
        _mm_set_epi8(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0,
                     i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0);
    i4_c2_c3_16x8b =
        _mm_set_epi8(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2,
                     i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2);

    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf;

    i4_samp_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_inp));
    i4_samp_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_inp + i4_src_stride));
    i4_samp_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1)));
    i4_samp_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 1) + i4_src_stride));
    i4_samp_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_inp + (i4_src_stride << 2)));

    i4_samp_16x8b_0 = _mm_unpacklo_epi8(i4_samp_16x8b_0, i4_samp_16x8b_1);
    i4_res_8x16b_r0 = _mm_maddubs_epi16(i4_samp_16x8b_0, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp), i4_res_8x16b_r0);

    i4_res_8x16b_r1 = _mm_maddubs_epi16(i4_samp_16x8b_0, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + i4_filt_stride), i4_res_8x16b_r1);

    i4_samp_16x8b_1 = _mm_unpacklo_epi8(i4_samp_16x8b_1, i4_samp_16x8b_2);
    i4_res_8x16b_r2 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1)), i4_res_8x16b_r2);

    i4_res_8x16b_r3 = _mm_maddubs_epi16(i4_samp_16x8b_1, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 1) + i4_filt_stride),
                     i4_res_8x16b_r3);

    i4_samp_16x8b_2 = _mm_unpacklo_epi8(i4_samp_16x8b_2, i4_samp_16x8b_3);
    i4_res_8x16b_r4 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c0_c1_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2)), i4_res_8x16b_r4);

    i4_res_8x16b_r5 = _mm_maddubs_epi16(i4_samp_16x8b_2, i4_c2_c3_16x8b);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + i4_filt_stride),
                     i4_res_8x16b_r5);

    i4_samp_16x8b_3 = _mm_unpacklo_epi8(i4_samp_16x8b_3, i4_samp_16x8b_4);
    i4_res_8x16b_r6 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c0_c1_16x8b);
    _mm_storel_epi64((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1)),
                     i4_res_8x16b_r6);

    i4_res_8x16b_r7 = _mm_maddubs_epi16(i4_samp_16x8b_3, i4_c2_c3_16x8b);
    i4_res_8x16b_r6 = _mm_shuffle_epi32(i4_res_8x16b_r6, 78);
    i4_res_8x16b_r7 = _mm_shuffle_epi32(i4_res_8x16b_r7, 147);
    i4_res_8x16b_r7_temp = _mm_blend_epi16(i4_res_8x16b_r6, i4_res_8x16b_r7, 252);
    _mm_storeu_si128((__m128i *) (pi2_tmp + (i4_filt_stride << 2) + (i4_filt_stride << 1) + 4),
                     i4_res_8x16b_r7_temp);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_1_sse42                 */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold       */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_horz_interpol_chroma_dyadic_1_sse42(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                               WORD32 i4_out_stride, WORD32 i4_phase_0,
                                               WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
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

    __m128i i4_res_final_8x16b_r1;
    __m128i i4_res_final_8x16b_r2;
    __m128i i4_res_final_8x16b_r3;
    __m128i i4_res_final_8x16b_r4;
    __m128i i4_res_final_8x16b_r5;
    __m128i i4_res_final_8x16b_r6;
    __m128i i4_res_final_8x16b_r7;
    __m128i i4_res_final_8x16b_r8;

    __m128i out_16x8b_r1;
    __m128i out_16x8b_r2;
    __m128i out_16x8b_r3;
    __m128i out_16x8b_r4;
    __m128i out_16x8b_r5;
    __m128i out_16x8b_r6;
    __m128i out_16x8b_r7;
    __m128i out_16x8b_r8;
    __m128i i4_res_final_8x16b_r12_0, i4_res_final_8x16b_r12_1;
    __m128i i4_res_final_8x16b_r34_0, i4_res_final_8x16b_r34_1;
    __m128i i4_res_final_8x16b_r56_0, i4_res_final_8x16b_r56_1;
    __m128i i4_res_final_8x16b_r67_0, i4_res_final_8x16b_r67_1;
    __m128i chroma_mask, chroma_mask2;
    __m128i coeff_c0_c1_8x16b, coeff_c2_c3_8x16b, res_32;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;
    coeff_c0_c1_8x16b = _mm_set_epi16(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1,
                                      i4_coeff_0, i4_coeff_1, i4_coeff_0);
    coeff_c2_c3_8x16b = _mm_set_epi16(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3,
                                      i4_coeff_2, i4_coeff_3, i4_coeff_2);
    res_32 = _mm_set1_epi32(32);
    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_dst_stride = i4_out_stride;

    i4_dst_stride2 = i4_dst_stride << 1;
    i4_dst_stride4 = i4_dst_stride << 2;

    /* Horizontal interpolation */
    /* x = 0, x_phase = phase_0 */
    i4_samp_8x16b_r1_0 = _mm_loadu_si128((__m128i *) pi2_tmp);         // a0 a1 a2 a3 a4 a5 a6 a7
    i4_samp_8x16b_r2_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 6));   // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r3_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 12));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r4_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 18));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r5_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 24));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r6_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 30));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r7_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 36));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r8_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 42));  // b0 b1 b2 b3 b4 b5 b6 b7

    i4_samp_8x16b_r1_1 = _mm_srli_si128(i4_samp_8x16b_r1_0, 2);        // a1 a2 a3 a4 a5 a6 a7 0
    i4_samp_8x16b_r1_2 = _mm_srli_si128(i4_samp_8x16b_r1_0, 4);        // a2 a3 a4 a5 a6 a7 0 0

    i4_samp_8x16b_r2_1 = _mm_srli_si128(i4_samp_8x16b_r2_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r2_2 = _mm_srli_si128(i4_samp_8x16b_r2_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r3_1 = _mm_srli_si128(i4_samp_8x16b_r3_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r3_2 = _mm_srli_si128(i4_samp_8x16b_r3_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r4_1 = _mm_srli_si128(i4_samp_8x16b_r4_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r4_2 = _mm_srli_si128(i4_samp_8x16b_r4_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r5_1 = _mm_srli_si128(i4_samp_8x16b_r5_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r5_2 = _mm_srli_si128(i4_samp_8x16b_r5_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r6_1 = _mm_srli_si128(i4_samp_8x16b_r6_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r6_2 = _mm_srli_si128(i4_samp_8x16b_r6_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r7_1 = _mm_srli_si128(i4_samp_8x16b_r7_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r7_2 = _mm_srli_si128(i4_samp_8x16b_r7_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r8_1 = _mm_srli_si128(i4_samp_8x16b_r8_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r8_2 = _mm_srli_si128(i4_samp_8x16b_r8_0, 4);        // b2 b3 b4 b5 b6 b7 0 0

    i4_samp_8x16b_r1_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r1_0,
                                            i4_samp_8x16b_r1_1);  // a0 a1  a1 a2  a2 a3  a3 a4
    i4_samp_8x16b_r2_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r2_0,
                                            i4_samp_8x16b_r2_1);  // b0 b1  b1 b2  b2 b3  b3 b4
    i4_samp_8x16b_r3_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r3_0, i4_samp_8x16b_r3_1);
    i4_samp_8x16b_r4_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r4_0, i4_samp_8x16b_r4_1);
    i4_samp_8x16b_r5_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r5_0, i4_samp_8x16b_r5_1);
    i4_samp_8x16b_r6_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r6_0, i4_samp_8x16b_r6_1);
    i4_samp_8x16b_r7_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r7_0, i4_samp_8x16b_r7_1);
    i4_samp_8x16b_r8_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r8_0, i4_samp_8x16b_r8_1);

    i4_samp_8x16b_r1_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r1_1,
                                            i4_samp_8x16b_r1_2);  // a1 a2  a2 a3  a3 a4  a4 a5
    i4_samp_8x16b_r2_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r2_1,
                                            i4_samp_8x16b_r2_2);  // b1 b2  b2 b3  b3 b4  b4 b5
    i4_samp_8x16b_r3_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r3_1, i4_samp_8x16b_r3_2);
    i4_samp_8x16b_r4_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r4_1, i4_samp_8x16b_r4_2);
    i4_samp_8x16b_r5_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r5_1, i4_samp_8x16b_r5_2);
    i4_samp_8x16b_r6_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r6_1, i4_samp_8x16b_r6_2);
    i4_samp_8x16b_r7_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r7_1, i4_samp_8x16b_r7_2);
    i4_samp_8x16b_r8_1 = _mm_unpacklo_epi16(i4_samp_8x16b_r8_1, i4_samp_8x16b_r8_2);

    // a0c0+a1c1  a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
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

    i4_res_4x32b_r1_0 = _mm_add_epi32(i4_res_4x32b_r1_0, res_32);
    i4_res_4x32b_r2_0 = _mm_add_epi32(i4_res_4x32b_r2_0, res_32);
    i4_res_4x32b_r3_0 = _mm_add_epi32(i4_res_4x32b_r3_0, res_32);
    i4_res_4x32b_r4_0 = _mm_add_epi32(i4_res_4x32b_r4_0, res_32);
    i4_res_4x32b_r5_0 = _mm_add_epi32(i4_res_4x32b_r5_0, res_32);
    i4_res_4x32b_r6_0 = _mm_add_epi32(i4_res_4x32b_r6_0, res_32);
    i4_res_4x32b_r7_0 = _mm_add_epi32(i4_res_4x32b_r7_0, res_32);
    i4_res_4x32b_r8_0 = _mm_add_epi32(i4_res_4x32b_r8_0, res_32);

    i4_res_4x32b_r1_1 = _mm_add_epi32(i4_res_4x32b_r1_1, res_32);
    i4_res_4x32b_r2_1 = _mm_add_epi32(i4_res_4x32b_r2_1, res_32);
    i4_res_4x32b_r3_1 = _mm_add_epi32(i4_res_4x32b_r3_1, res_32);
    i4_res_4x32b_r4_1 = _mm_add_epi32(i4_res_4x32b_r4_1, res_32);
    i4_res_4x32b_r5_1 = _mm_add_epi32(i4_res_4x32b_r5_1, res_32);
    i4_res_4x32b_r6_1 = _mm_add_epi32(i4_res_4x32b_r6_1, res_32);
    i4_res_4x32b_r7_1 = _mm_add_epi32(i4_res_4x32b_r7_1, res_32);
    i4_res_4x32b_r8_1 = _mm_add_epi32(i4_res_4x32b_r8_1, res_32);

    i4_res_4x32b_r1_0 = _mm_srai_epi32(i4_res_4x32b_r1_0, 6);
    i4_res_4x32b_r2_0 = _mm_srai_epi32(i4_res_4x32b_r2_0, 6);
    i4_res_4x32b_r3_0 = _mm_srai_epi32(i4_res_4x32b_r3_0, 6);
    i4_res_4x32b_r4_0 = _mm_srai_epi32(i4_res_4x32b_r4_0, 6);
    i4_res_4x32b_r5_0 = _mm_srai_epi32(i4_res_4x32b_r5_0, 6);
    i4_res_4x32b_r6_0 = _mm_srai_epi32(i4_res_4x32b_r6_0, 6);
    i4_res_4x32b_r7_0 = _mm_srai_epi32(i4_res_4x32b_r7_0, 6);
    i4_res_4x32b_r8_0 = _mm_srai_epi32(i4_res_4x32b_r8_0, 6);

    i4_res_4x32b_r1_1 = _mm_srai_epi32(i4_res_4x32b_r1_1, 6);
    i4_res_4x32b_r2_1 = _mm_srai_epi32(i4_res_4x32b_r2_1, 6);
    i4_res_4x32b_r3_1 = _mm_srai_epi32(i4_res_4x32b_r3_1, 6);
    i4_res_4x32b_r4_1 = _mm_srai_epi32(i4_res_4x32b_r4_1, 6);
    i4_res_4x32b_r5_1 = _mm_srai_epi32(i4_res_4x32b_r5_1, 6);
    i4_res_4x32b_r6_1 = _mm_srai_epi32(i4_res_4x32b_r6_1, 6);
    i4_res_4x32b_r7_1 = _mm_srai_epi32(i4_res_4x32b_r7_1, 6);
    i4_res_4x32b_r8_1 = _mm_srai_epi32(i4_res_4x32b_r8_1, 6);

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
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride << 1)), out_16x8b_r3);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 3)), out_16x8b_r4);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride << 2)), out_16x8b_r5);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 5)), out_16x8b_r6);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 6)), out_16x8b_r7);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 7)), out_16x8b_r8);
    /* End of loop over x */
} /* isvcd_horz_interpol_chroma_dyadic_1_sse42 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_2_sse42                 */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_horz_interpol_chroma_dyadic_2_sse42(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                               WORD32 i4_out_stride, WORD32 i4_phase_0,
                                               WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_dst_stride, i4_dst_stride2, i4_dst_stride4;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    __m128i i4_samp_8x16b_r1_0, i4_samp_8x16b_r1_1;
    __m128i i4_samp_8x16b_r2_0, i4_samp_8x16b_r2_1;
    __m128i i4_samp_8x16b_r3_0, i4_samp_8x16b_r3_1;
    __m128i i4_samp_8x16b_r4_0, i4_samp_8x16b_r4_1;
    __m128i i4_samp_8x16b_r5_0, i4_samp_8x16b_r5_1;
    __m128i i4_samp_8x16b_r6_0, i4_samp_8x16b_r6_1;
    __m128i i4_samp_8x16b_r7_0, i4_samp_8x16b_r7_1;
    __m128i i4_samp_8x16b_r8_0, i4_samp_8x16b_r8_1;

    __m128i i4_res_4x32b_r1_0, i4_res_4x32b_r1_1;
    __m128i i4_res_4x32b_r2_0, i4_res_4x32b_r2_1;
    __m128i i4_res_4x32b_r3_0, i4_res_4x32b_r3_1;
    __m128i i4_res_4x32b_r4_0, i4_res_4x32b_r4_1;
    __m128i i4_res_4x32b_r5_0, i4_res_4x32b_r5_1;
    __m128i i4_res_4x32b_r6_0, i4_res_4x32b_r6_1;
    __m128i i4_res_4x32b_r7_0, i4_res_4x32b_r7_1;
    __m128i i4_res_4x32b_r8_0, i4_res_4x32b_r8_1;

    __m128i i4_res_final_8x16b_r1;
    __m128i i4_res_final_8x16b_r2;
    __m128i i4_res_final_8x16b_r3;
    __m128i i4_res_final_8x16b_r4;
    __m128i i4_res_final_8x16b_r5;
    __m128i i4_res_final_8x16b_r6;
    __m128i i4_res_final_8x16b_r7;
    __m128i i4_res_final_8x16b_r8;

    __m128i out_16x8b_r1;
    __m128i out_16x8b_r2;
    __m128i out_16x8b_r3;
    __m128i out_16x8b_r4;
    __m128i out_16x8b_r5;
    __m128i out_16x8b_r6;
    __m128i out_16x8b_r7;
    __m128i out_16x8b_r8;
    __m128i i4_res_final_8x16b_r12_0, i4_res_final_8x16b_r12_1;
    __m128i i4_res_final_8x16b_r34_0, i4_res_final_8x16b_r34_1;
    __m128i i4_res_final_8x16b_r56_0, i4_res_final_8x16b_r56_1;
    __m128i i4_res_final_8x16b_r67_0, i4_res_final_8x16b_r67_1;
    __m128i chroma_mask, chroma_mask2;
    __m128i coeff_c0_c1_8x16b, coeff_c2_c3_8x16b, res_32;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;
    coeff_c0_c1_8x16b = _mm_set_epi16(i4_coeff_1, i4_coeff_0, i4_coeff_1, i4_coeff_0, i4_coeff_1,
                                      i4_coeff_0, i4_coeff_1, i4_coeff_0);
    coeff_c2_c3_8x16b = _mm_set_epi16(i4_coeff_3, i4_coeff_2, i4_coeff_3, i4_coeff_2, i4_coeff_3,
                                      i4_coeff_2, i4_coeff_3, i4_coeff_2);
    res_32 = _mm_set1_epi32(32);
    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf + 1;
    i4_dst_stride = i4_out_stride;

    i4_dst_stride2 = i4_dst_stride << 1;
    i4_dst_stride4 = i4_dst_stride << 2;

    /* Horizontal interpolation */
    /* x = 0, x_phase = phase_0 */
    i4_samp_8x16b_r1_0 = _mm_loadu_si128((__m128i *) pi2_tmp);         // a0 a1 a2 a3 a4 a5 a6 a7
    i4_samp_8x16b_r2_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 6));   // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r3_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 12));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r4_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 18));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r5_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 24));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r6_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 30));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r7_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 36));  // b0 b1 b2 b3 b4 b5 b6 b7
    i4_samp_8x16b_r8_0 = _mm_loadu_si128((__m128i *) (pi2_tmp + 42));  // b0 b1 b2 b3 b4 b5 b6 b7

    i4_samp_8x16b_r1_1 = _mm_srli_si128(i4_samp_8x16b_r1_0, 2);        // a1 a2 a3 a4 a5 a6 a7 0
    i4_samp_8x16b_r2_1 = _mm_srli_si128(i4_samp_8x16b_r2_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r3_1 = _mm_srli_si128(i4_samp_8x16b_r3_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r4_1 = _mm_srli_si128(i4_samp_8x16b_r4_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r5_1 = _mm_srli_si128(i4_samp_8x16b_r5_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r6_1 = _mm_srli_si128(i4_samp_8x16b_r6_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r7_1 = _mm_srli_si128(i4_samp_8x16b_r7_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0
    i4_samp_8x16b_r8_1 = _mm_srli_si128(i4_samp_8x16b_r8_0, 2);        // b1 b2 b3 b4 b5 b6 b7 0

    i4_samp_8x16b_r1_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r1_0,
                                            i4_samp_8x16b_r1_1);  // a0 a1  a1 a2  a2 a3  a3 a4
    i4_samp_8x16b_r2_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r2_0,
                                            i4_samp_8x16b_r2_1);  // b0 b1  b1 b2  b2 b3  b3 b4
    i4_samp_8x16b_r3_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r3_0, i4_samp_8x16b_r3_1);
    i4_samp_8x16b_r4_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r4_0, i4_samp_8x16b_r4_1);
    i4_samp_8x16b_r5_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r5_0, i4_samp_8x16b_r5_1);
    i4_samp_8x16b_r6_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r6_0, i4_samp_8x16b_r6_1);
    i4_samp_8x16b_r7_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r7_0, i4_samp_8x16b_r7_1);
    i4_samp_8x16b_r8_0 = _mm_unpacklo_epi16(i4_samp_8x16b_r8_0, i4_samp_8x16b_r8_1);

    // a0c0+a1c1  a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
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
    i4_res_4x32b_r1_1 = _mm_madd_epi16(i4_samp_8x16b_r1_0, coeff_c2_c3_8x16b);
    // b1c2+b2c3  b2c2+b3c3  b3c2+b4c3  b4c2+b5c3
    i4_res_4x32b_r2_1 = _mm_madd_epi16(i4_samp_8x16b_r2_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r3_1 = _mm_madd_epi16(i4_samp_8x16b_r3_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r4_1 = _mm_madd_epi16(i4_samp_8x16b_r4_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r5_1 = _mm_madd_epi16(i4_samp_8x16b_r5_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r6_1 = _mm_madd_epi16(i4_samp_8x16b_r6_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r7_1 = _mm_madd_epi16(i4_samp_8x16b_r7_0, coeff_c2_c3_8x16b);
    i4_res_4x32b_r8_1 = _mm_madd_epi16(i4_samp_8x16b_r8_0, coeff_c2_c3_8x16b);

    i4_res_4x32b_r1_0 = _mm_add_epi32(i4_res_4x32b_r1_0, res_32);
    i4_res_4x32b_r2_0 = _mm_add_epi32(i4_res_4x32b_r2_0, res_32);
    i4_res_4x32b_r3_0 = _mm_add_epi32(i4_res_4x32b_r3_0, res_32);
    i4_res_4x32b_r4_0 = _mm_add_epi32(i4_res_4x32b_r4_0, res_32);
    i4_res_4x32b_r5_0 = _mm_add_epi32(i4_res_4x32b_r5_0, res_32);
    i4_res_4x32b_r6_0 = _mm_add_epi32(i4_res_4x32b_r6_0, res_32);
    i4_res_4x32b_r7_0 = _mm_add_epi32(i4_res_4x32b_r7_0, res_32);
    i4_res_4x32b_r8_0 = _mm_add_epi32(i4_res_4x32b_r8_0, res_32);

    i4_res_4x32b_r1_1 = _mm_add_epi32(i4_res_4x32b_r1_1, res_32);
    i4_res_4x32b_r2_1 = _mm_add_epi32(i4_res_4x32b_r2_1, res_32);
    i4_res_4x32b_r3_1 = _mm_add_epi32(i4_res_4x32b_r3_1, res_32);
    i4_res_4x32b_r4_1 = _mm_add_epi32(i4_res_4x32b_r4_1, res_32);
    i4_res_4x32b_r5_1 = _mm_add_epi32(i4_res_4x32b_r5_1, res_32);
    i4_res_4x32b_r6_1 = _mm_add_epi32(i4_res_4x32b_r6_1, res_32);
    i4_res_4x32b_r7_1 = _mm_add_epi32(i4_res_4x32b_r7_1, res_32);
    i4_res_4x32b_r8_1 = _mm_add_epi32(i4_res_4x32b_r8_1, res_32);

    i4_res_4x32b_r1_0 = _mm_srai_epi32(i4_res_4x32b_r1_0, 6);
    i4_res_4x32b_r2_0 = _mm_srai_epi32(i4_res_4x32b_r2_0, 6);
    i4_res_4x32b_r3_0 = _mm_srai_epi32(i4_res_4x32b_r3_0, 6);
    i4_res_4x32b_r4_0 = _mm_srai_epi32(i4_res_4x32b_r4_0, 6);
    i4_res_4x32b_r5_0 = _mm_srai_epi32(i4_res_4x32b_r5_0, 6);
    i4_res_4x32b_r6_0 = _mm_srai_epi32(i4_res_4x32b_r6_0, 6);
    i4_res_4x32b_r7_0 = _mm_srai_epi32(i4_res_4x32b_r7_0, 6);
    i4_res_4x32b_r8_0 = _mm_srai_epi32(i4_res_4x32b_r8_0, 6);

    i4_res_4x32b_r1_1 = _mm_srai_epi32(i4_res_4x32b_r1_1, 6);
    i4_res_4x32b_r2_1 = _mm_srai_epi32(i4_res_4x32b_r2_1, 6);
    i4_res_4x32b_r3_1 = _mm_srai_epi32(i4_res_4x32b_r3_1, 6);
    i4_res_4x32b_r4_1 = _mm_srai_epi32(i4_res_4x32b_r4_1, 6);
    i4_res_4x32b_r5_1 = _mm_srai_epi32(i4_res_4x32b_r5_1, 6);
    i4_res_4x32b_r6_1 = _mm_srai_epi32(i4_res_4x32b_r6_1, 6);
    i4_res_4x32b_r7_1 = _mm_srai_epi32(i4_res_4x32b_r7_1, 6);
    i4_res_4x32b_r8_1 = _mm_srai_epi32(i4_res_4x32b_r8_1, 6);

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
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride << 1)), out_16x8b_r3);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 3)), out_16x8b_r4);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride << 2)), out_16x8b_r5);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 5)), out_16x8b_r6);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 6)), out_16x8b_r7);
    _mm_storeu_si128((__m128i *) (pu1_out + (i4_dst_stride * 7)), out_16x8b_r8);

    /* End of loop over x */
}
