/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

#ifdef __ANDROID__
#include "log/log.h"
#include <cutils/log.h>
#endif

#include <immintrin.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_inter_pred_filters.h"

/*****************************************************************************/
/* Constant Data variables                                                   */
/*****************************************************************************/

/* coefficients for 6 tap filtering*/
//const WORD32 ih264_g_six_tap[3] ={1,-5,20};
/*****************************************************************************/
/*  Function definitions .                                                   */
/*****************************************************************************/
/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_inter_pred_luma_copy_avx2                          */
/*                                                                           */
/*  Description   : This function copies the contents of ht x wd block from  */
/*                  source to destination. (ht,wd) can be (4,4), (8,4),      */
/*                  (4,8), (8,8), (16,8), (8,16) or (16,16).                 */
/*                                                                           */
/*  Inputs        : puc_src  - pointer to source                             */
/*                  puc_dst  - pointer to destination                        */
/*                  src_strd - stride for source                             */
/*                  dst_strd - stride for destination                        */
/*                  ht       - height of the block                           */
/*                  wd       - width of the block                            */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         13 02 2015   Kaushik         Initial Version                      */
/*                      Senthoor                                             */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_inter_pred_luma_copy_avx2(UWORD8 *pu1_src,
                                      UWORD8 *pu1_dst,
                                      WORD32 src_strd,
                                      WORD32 dst_strd,
                                      WORD32 ht,
                                      WORD32 wd,
                                      UWORD8* pu1_tmp,
                                      WORD32 dydx)
{
    WORD32 src_strd2, src_strd3, src_strd4, dst_strd2, dst_strd3, dst_strd4;
    UNUSED(pu1_tmp);
    UNUSED(dydx);

    src_strd2 = src_strd << 1;
    dst_strd2 = dst_strd << 1;
    src_strd4 = src_strd << 2;
    dst_strd4 = dst_strd << 2;
    src_strd3 = src_strd2 + src_strd;
    dst_strd3 = dst_strd2 + dst_strd;
    if(wd == 4)
    {
        do
        {
            *((WORD32 *)(pu1_dst)) =  *((WORD32 *)(pu1_src));
            *((WORD32 *)(pu1_dst + dst_strd)) = *((WORD32 *)(pu1_src + src_strd));
            *((WORD32 *)(pu1_dst + dst_strd2)) = *((WORD32 *)(pu1_src + src_strd2));
            *((WORD32 *)(pu1_dst + dst_strd3)) = *((WORD32 *)(pu1_src + src_strd3));

            ht -= 4;
            pu1_src += src_strd4;
            pu1_dst += dst_strd4;
        }
        while(ht > 0);
    }
    else if(wd == 8)
    {
        __m128i y_0_16x8b, y_1_16x8b, y_2_16x8b, y_3_16x8b;
        do
        {

            y_0_16x8b = _mm_loadl_epi64((__m128i *)pu1_src);
            y_1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src + src_strd));
            y_2_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src + src_strd2));
            y_3_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src + src_strd3));

            _mm_storel_epi64((__m128i *)pu1_dst, y_0_16x8b);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd), y_1_16x8b);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd2), y_2_16x8b);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd3), y_3_16x8b);

            ht -= 4;
            pu1_src += src_strd4;
            pu1_dst += dst_strd4;
        }
        while(ht > 0);
    }
    else // wd == 16
    {
        __m256i y_0_16x8b, y_1_16x8b, y_2_16x8b, y_3_16x8b;
        WORD32 src_strd5, src_strd6, src_strd7, src_strd8;
        WORD32 dst_strd5, dst_strd6, dst_strd7, dst_strd8;

        __m256i y_4_16x8b, y_5_16x8b, y_6_16x8b, y_7_16x8b,y_0_1,y_2_3,y_4_5,y_6_7;


        src_strd5 = src_strd2 + src_strd3;
        dst_strd5 = dst_strd2 + dst_strd3;
        src_strd6 = src_strd3 << 1;
        dst_strd6 = dst_strd3 << 1;
        src_strd7 = src_strd3 + src_strd4;
        dst_strd7 = dst_strd3 + dst_strd4;
        src_strd8 = src_strd << 3;
        dst_strd8 = dst_strd << 3;

        do
        {

            y_0_1 = _mm256_loadu2_m128i((__m128i *)(pu1_src + src_strd),(__m128i *)pu1_src);
            y_2_3 = _mm256_loadu2_m128i((__m128i *)(pu1_src + src_strd3),(__m128i *)(pu1_src + src_strd2));
            y_4_5 = _mm256_loadu2_m128i((__m128i *)(pu1_src + src_strd5),(__m128i *)(pu1_src + src_strd4));
            y_6_7 = _mm256_loadu2_m128i((__m128i *)(pu1_src + src_strd7),(__m128i *)(pu1_src + src_strd6));

            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd),(__m128i *)pu1_dst,y_0_1);
            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd3),(__m128i *)(pu1_dst + dst_strd2),y_2_3);
            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd5),(__m128i *)(pu1_dst + dst_strd4),y_4_5);
            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd7),(__m128i *)(pu1_dst + dst_strd6),y_6_7);

            ht -= 8;
            pu1_src += src_strd8;
            pu1_dst += dst_strd8;
        }
        while(ht > 0);
    }
}




/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_inter_pred_chroma_avx2                           */
/*                                                                           */
/*  Description   : This function implements a four-tap 2D filter as         */
/*                  mentioned in sec. 8.4.2.2.2 titled "Chroma sample        */
/*                  "interpolation process". (ht,wd) can be (2,2), (4,2),    */
/*                  (2,4), (4,4), (8,4), (4,8) or (8,8).                     */
/*                                                                           */
/*  Inputs        : puc_src  - pointer to source                             */
/*                  puc_dst  - pointer to destination                        */
/*                  src_strd - stride for source                             */
/*                  dst_strd - stride for destination                        */
/*                  dx       - x position of destination value               */
/*                  dy       - y position of destination value               */
/*                  ht       - height of the block                           */
/*                  wd       - width of the block                            */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         13 02 2015   Kaushik         Initial Version                      */
/*                      Senthoor                                             */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_inter_pred_chroma_avx2(UWORD8 *pu1_src,
                                   UWORD8 *pu1_dst,
                                   WORD32 src_strd,
                                   WORD32 dst_strd,
                                   WORD32 dx,
                                   WORD32 dy,
                                   WORD32 ht,
                                   WORD32 wd)
{

    WORD32 i, j, A, B, C, D;
    i = 8 - dx;
    j = 8 - dy;

    A = i * j;
    B = dx * j;
    C = i * dy;
    D = dx * dy;
    if(wd == 2)
    {
        WORD32 tmp1, tmp2, tmp3, tmp4;

        do
        {
            //U
            tmp1 = A * pu1_src[0] + B * pu1_src[2] + C * pu1_src[src_strd] + D * pu1_src[src_strd + 2];
            tmp2 = A * pu1_src[2] + B * pu1_src[4] + C * pu1_src[src_strd + 2] + D * pu1_src[src_strd + 4];
            //V
            tmp3 = A * pu1_src[1] + B * pu1_src[3] + C * pu1_src[src_strd + 1] + D * pu1_src[src_strd + 3];
            tmp4 = A * pu1_src[3] + B * pu1_src[5] + C * pu1_src[src_strd + 3] + D * pu1_src[src_strd + 5];

            tmp1 = (tmp1 + 32) >> 6;
            tmp2 = (tmp2 + 32) >> 6;
            tmp3 = (tmp3 + 32) >> 6;
            tmp4 = (tmp4 + 32) >> 6;

            pu1_dst[0] = CLIP_U8(tmp1);
            pu1_dst[2] = CLIP_U8(tmp2);
            pu1_dst[1] = CLIP_U8(tmp3);
            pu1_dst[3] = CLIP_U8(tmp4);

            pu1_src += src_strd;
            pu1_dst += dst_strd;

            tmp1 = A * pu1_src[0] + B * pu1_src[2] + C * pu1_src[src_strd] + D * pu1_src[src_strd + 2];
            tmp2 = A * pu1_src[2] + B * pu1_src[4] + C * pu1_src[src_strd + 2] + D * pu1_src[src_strd + 4];
            tmp3 = A * pu1_src[1] + B * pu1_src[3] + C * pu1_src[src_strd + 1] + D * pu1_src[src_strd + 3];
            tmp4 = A * pu1_src[3] + B * pu1_src[5] + C * pu1_src[src_strd + 3] + D * pu1_src[src_strd + 5];

            tmp1 = (tmp1 + 32) >> 6;
            tmp2 = (tmp2 + 32) >> 6;
            tmp3 = (tmp3 + 32) >> 6;
            tmp4 = (tmp4 + 32) >> 6;

            pu1_dst[0] = CLIP_U8(tmp1);
            pu1_dst[2] = CLIP_U8(tmp2);
            pu1_dst[1] = CLIP_U8(tmp3);
            pu1_dst[3] = CLIP_U8(tmp4);

            ht -= 2;
            pu1_src += src_strd;
            pu1_dst += dst_strd;
        }
        while(ht > 0);

    }
    else if(wd == 4)
    {
        WORD32 AB, CD;

        __m256i coeffAB_32x8b, coeffCD_32x8b, round_add32_8x32b;
        __m256i const_shuff_32x8b;

        __m256i  src_r23_32x8b,src_r12_32x8b,res12_AB_8x32b,res12_CD_8x32b,res1_8x32b;
        __m128i res1,src_r1_16x8b_128;
        __m128i const_shuff_16x8b_128;

        AB = (B << 8) + A;
        CD = (D << 8) + C;

        coeffAB_32x8b = _mm256_set1_epi16(AB);
        coeffCD_32x8b = _mm256_set1_epi16(CD);

        round_add32_8x32b = _mm256_set1_epi16(32);
        const_shuff_16x8b_128 = _mm_setr_epi32(0x03010200, 0x05030402, 0x07050604, 0x09070806);
        const_shuff_32x8b = _mm256_setr_epi32(0x03010200, 0x05030402, 0x07050604, 0x09070806,0x03010200, 0x05030402, 0x07050604, 0x09070806);


        src_r1_16x8b_128 = _mm_loadu_si128((__m128i *)pu1_src);
        src_r1_16x8b_128 = _mm_shuffle_epi8(src_r1_16x8b_128, const_shuff_16x8b_128);
        pu1_src += src_strd;
        do
        {

            src_r23_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + src_strd),(__m128i *)pu1_src);
            src_r23_32x8b = _mm256_shuffle_epi8(src_r23_32x8b, const_shuff_32x8b);
            src_r12_32x8b = _mm256_set_m128i(_mm256_castsi256_si128(src_r23_32x8b),src_r1_16x8b_128);

            res12_AB_8x32b = _mm256_maddubs_epi16(src_r12_32x8b, coeffAB_32x8b);
            res12_CD_8x32b = _mm256_maddubs_epi16(src_r23_32x8b, coeffCD_32x8b);

            res1_8x32b = _mm256_add_epi16(res12_AB_8x32b, res12_CD_8x32b);
            res1_8x32b = _mm256_add_epi16(res1_8x32b, round_add32_8x32b);

            res1_8x32b = _mm256_srai_epi16(res1_8x32b, 6);

            res1_8x32b = _mm256_packus_epi16(res1_8x32b,res1_8x32b);
            res1_8x32b = _mm256_permute4x64_epi64(res1_8x32b, 0xD8);
            res1  = _mm256_castsi256_si128(res1_8x32b);
            _mm_storel_epi64((__m128i *)pu1_dst, res1);

            res1 = _mm_srli_si128(res1, 8);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd), res1);

            src_r1_16x8b_128 = _mm256_castsi256_si128(_mm256_permute2x128_si256(src_r23_32x8b,src_r23_32x8b,0x1));;

            ht -= 2;
            pu1_src += src_strd << 1;
            pu1_dst += dst_strd << 1;

        }
        while(ht > 0);
    }
    else // wd == 8
    {

        WORD32 AB, CD;
        __m256i src_r1lh_32x8b, src_r2lh_32x8b;

        __m256i res_lh_AB_8x32b, res_lh_CD_8x32b,res_1lh_8x32b,res_2lh_8x32b;
        __m256i res_lh_8x32b, res_l_8x32b,res_h_8x32b, res_32x8b;

        __m256i coeffAB_32x8b, coeffCD_32x8b, round_add32_8x32b;
        __m256i const_shuff_32x8b;
        __m256i zero = _mm256_setzero_si256();
        __m128i res_1,res_2;

        AB = (B << 8) + A;
        CD = (D << 8) + C;

        coeffAB_32x8b = _mm256_set1_epi16(AB);
        coeffCD_32x8b = _mm256_set1_epi16(CD);

        round_add32_8x32b = _mm256_set1_epi16(32);

        const_shuff_32x8b = _mm256_setr_epi32(0x03010200, 0x05030402, 0x07050604, 0x09070806,0x03010200, 0x05030402, 0x07050604, 0x09070806);

        src_r1lh_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + 8),(__m128i *)pu1_src);
        src_r1lh_32x8b = _mm256_shuffle_epi8(src_r1lh_32x8b, const_shuff_32x8b);
        pu1_src += src_strd;

        do
        {
            //row 1

            src_r2lh_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + 8),(__m128i *)pu1_src);
            src_r2lh_32x8b = _mm256_shuffle_epi8(src_r2lh_32x8b, const_shuff_32x8b);
            res_lh_AB_8x32b = _mm256_maddubs_epi16(src_r1lh_32x8b, coeffAB_32x8b);

            res_lh_CD_8x32b = _mm256_maddubs_epi16(src_r2lh_32x8b, coeffCD_32x8b);
            res_lh_8x32b = _mm256_add_epi16(res_lh_AB_8x32b, round_add32_8x32b);
            res_lh_8x32b = _mm256_add_epi16(res_lh_8x32b, res_lh_CD_8x32b);
            res_1lh_8x32b = _mm256_srai_epi16(res_lh_8x32b, 6);

            pu1_src += src_strd;
            //row 2
            src_r1lh_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + 8),(__m128i *)pu1_src);
            src_r1lh_32x8b = _mm256_shuffle_epi8(src_r1lh_32x8b, const_shuff_32x8b);

            res_lh_AB_8x32b = _mm256_maddubs_epi16(src_r2lh_32x8b, coeffAB_32x8b);
            res_lh_CD_8x32b = _mm256_maddubs_epi16(src_r1lh_32x8b, coeffCD_32x8b);

            res_lh_8x32b = _mm256_add_epi16(res_lh_AB_8x32b, round_add32_8x32b);
            res_lh_8x32b = _mm256_add_epi16(res_lh_8x32b, res_lh_CD_8x32b);

            res_2lh_8x32b = _mm256_srai_epi16(res_lh_8x32b, 6);

            res_1lh_8x32b =  _mm256_packus_epi16(res_1lh_8x32b, res_2lh_8x32b);
            res_1lh_8x32b = _mm256_permute4x64_epi64(res_1lh_8x32b, 0xD8);
            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd),(__m128i *)(pu1_dst),res_1lh_8x32b);

            pu1_src += src_strd;
            pu1_dst += dst_strd;
            pu1_dst += dst_strd;

            //row 3
            src_r2lh_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + 8),(__m128i *)pu1_src);
            src_r2lh_32x8b = _mm256_shuffle_epi8(src_r2lh_32x8b, const_shuff_32x8b);

            res_lh_AB_8x32b = _mm256_maddubs_epi16(src_r1lh_32x8b, coeffAB_32x8b);
            res_lh_CD_8x32b = _mm256_maddubs_epi16(src_r2lh_32x8b, coeffCD_32x8b);

            res_lh_8x32b = _mm256_add_epi16(res_lh_AB_8x32b, round_add32_8x32b);
            res_lh_8x32b = _mm256_add_epi16(res_lh_8x32b, res_lh_CD_8x32b);

            res_1lh_8x32b = _mm256_srai_epi16(res_lh_8x32b, 6);
            pu1_src += src_strd;

            //row 1
            src_r1lh_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_src + 8),(__m128i *)pu1_src);
            src_r1lh_32x8b = _mm256_shuffle_epi8(src_r1lh_32x8b, const_shuff_32x8b);

            res_lh_AB_8x32b = _mm256_maddubs_epi16(src_r2lh_32x8b, coeffAB_32x8b);
            res_lh_CD_8x32b = _mm256_maddubs_epi16(src_r1lh_32x8b, coeffCD_32x8b);

            res_lh_8x32b = _mm256_add_epi16(res_lh_AB_8x32b, round_add32_8x32b);
            res_lh_8x32b = _mm256_add_epi16(res_lh_8x32b, res_lh_CD_8x32b);

            res_2lh_8x32b = _mm256_srai_epi16(res_lh_8x32b, 6);

            res_1lh_8x32b =  _mm256_packus_epi16(res_1lh_8x32b, res_2lh_8x32b);
            res_1lh_8x32b = _mm256_permute4x64_epi64(res_1lh_8x32b, 0xD8);
            _mm256_storeu2_m128i((__m128i *)(pu1_dst + dst_strd),(__m128i *)(pu1_dst),res_1lh_8x32b);


            ht -= 4;
            pu1_src += src_strd;
            pu1_dst += dst_strd;
            pu1_dst += dst_strd;
        }
        while(ht > 0);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_inter_pred_luma_horz_qpel_vert_qpel_avx2           */
/*                                                                           */
/*  Description   : This function implements a six-tap filter vertically and */
/*                  horizontally on ht x wd block separately and averages    */
/*                  the two sets of values to calculate values at (1/4,1/4), */
/*                  (1/4, 3/4), (3/4, 1/4) or (3/4, 3/4) as mentioned in     */
/*                  sec. 8.4.2.2.1 titled "Luma sample interpolation         */
/*                  process". (ht,wd) can be (4,4), (8,4), (4,8), (8,8),     */
/*                  (16,8), (8,16) or (16,16).                               */
/*                                                                           */
/*  Inputs        : puc_src  - pointer to source                             */
/*                  puc_dst  - pointer to destination                        */
/*                  src_strd - stride for source                             */
/*                  dst_strd - stride for destination                        */
/*                  ht       - height of the block                           */
/*                  wd       - width of the block                            */
/*                  pu1_tmp  - pointer to temporary buffer                   */
/*                  dydx     - x and y reference offset for q-pel            */
/*                             calculations                                  */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         13 02 2015   Kaushik         Initial Version                      */
/*                      Senthoor                                             */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_inter_pred_luma_horz_qpel_vert_qpel_avx2(UWORD8 *pu1_src,
                                                     UWORD8 *pu1_dst,
                                                     WORD32 src_strd,
                                                     WORD32 dst_strd,
                                                     WORD32 ht,
                                                     WORD32 wd,
                                                     UWORD8* pu1_tmp,
                                                     WORD32 dydx)
{
    WORD32 ht_temp;
    UWORD8 *pu1_pred_vert,*pu1_pred_horiz;
    UWORD8 *pu1_tmp1, *pu1_tmp2;
    WORD32 x_offset, y_offset;

    __m128i coeff0_1_16x8b, coeff2_3_16x8b, coeff4_5_16x8b;
    __m128i const_val16_8x16b;
    __m256i coeff0_1_32x8b, coeff2_3_32x8b, coeff4_5_32x8b,coeff_32x8b;
    __m256i const_val16_8x32b;
    __m256i zero = _mm256_setzero_si256();
    pu1_tmp1 = pu1_tmp;

    dydx &= 0xf;
    ht_temp = ht;
    x_offset = dydx & 0x3;
    y_offset = dydx >> 2;
    pu1_tmp2 = pu1_tmp1;

    pu1_pred_vert  = pu1_src + (x_offset >> 1) - 2*src_strd;
    pu1_pred_horiz = pu1_src + (y_offset >> 1) * src_strd - 2;
    //the filter input starts from x[-2] (till x[3])

    coeff0_1_16x8b = _mm_set1_epi32(0xFB01FB01);  //c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1
    coeff2_3_16x8b = _mm_set1_epi32(0x14141414);  //c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3
    coeff4_5_16x8b = _mm_set1_epi32(0x01FB01FB);  //c4 c5 c5 c5 c4 c5 c5 c5 c4 c5 c5 c5 c4 c5 c5 c5
                                                  //c0 = c5 = 1, c1 = c4 = -5, c2 = c3 = 20
    const_val16_8x16b = _mm_set1_epi16(16);
    coeff_32x8b = _mm256_set_m128i(coeff2_3_16x8b,coeff0_1_16x8b);
    coeff0_1_32x8b = _mm256_set1_epi32(0xFB01FB01);  //c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1 c0 c1
    coeff2_3_32x8b = _mm256_set1_epi32(0x14141414);  //c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3 c2 c3
    coeff4_5_32x8b = _mm256_set1_epi32(0x01FB01FB);  //c4 c5 c5 c5 c4 c5 c5 c5 c4 c5 c5 c5 c4 c5 c5 c5
                                                     //c0 = c5 = 1, c1 = c4 = -5, c2 = c3 = 20
    const_val16_8x32b = _mm256_set1_epi16(16);

    if(wd == 4)
    {
        //vertical q-pel filter
        {
            __m128i src_r0_16x8b, src_r1_16x8b, src_r2_16x8b, src_r3_16x8b, src_r4_16x8b;
            __m128i src_r5_16x8b, src_r6_16x8b;
            __m128i src_r0r1_16x8b, src_r2r3_16x8b, src_r4r5_16x8b;

            __m128i res_r0r1_16x8b, res_t1_8x16b, res_t2_8x16b, res_t3_8x16b;

            //epilogue: Load all the pred rows except sixth  and seventh row for the
            //first and second row processing.
            src_r0_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;

            src_r1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r0_16x8b = _mm_unpacklo_epi32(src_r0_16x8b, src_r1_16x8b);

            src_r2_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r1_16x8b = _mm_unpacklo_epi32(src_r1_16x8b, src_r2_16x8b);

            src_r3_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r2_16x8b = _mm_unpacklo_epi32(src_r2_16x8b, src_r3_16x8b);

            src_r4_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r3_16x8b = _mm_unpacklo_epi32(src_r3_16x8b, src_r4_16x8b);

            //Core Loop: Process all the rows.
            do
            {
                src_r5_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
                src_r4_16x8b = _mm_unpacklo_epi32(src_r4_16x8b, src_r5_16x8b);

                src_r6_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert + src_strd));
                src_r5_16x8b = _mm_unpacklo_epi32(src_r5_16x8b, src_r6_16x8b);

                src_r0r1_16x8b = _mm_unpacklo_epi8(src_r0_16x8b, src_r1_16x8b);
                src_r2r3_16x8b = _mm_unpacklo_epi8(src_r2_16x8b, src_r3_16x8b);
                src_r4r5_16x8b = _mm_unpacklo_epi8(src_r4_16x8b, src_r5_16x8b);

                res_t1_8x16b = _mm_maddubs_epi16(src_r0r1_16x8b, coeff0_1_16x8b);
                res_t2_8x16b = _mm_maddubs_epi16(src_r2r3_16x8b, coeff2_3_16x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t2_8x16b);
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t3_8x16b, res_t1_8x16b);

                res_t1_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits.
                res_r0r1_16x8b = _mm_packus_epi16(res_t1_8x16b, res_t1_8x16b);

                _mm_storel_epi64((__m128i *)pu1_tmp1, res_r0r1_16x8b);

                src_r0_16x8b = src_r2_16x8b;
                src_r1_16x8b = src_r3_16x8b;
                src_r2_16x8b = src_r4_16x8b;
                src_r3_16x8b = src_r5_16x8b;
                src_r4_16x8b = src_r6_16x8b;

                ht_temp -= 2;
                pu1_pred_vert += src_strd << 1;
                pu1_tmp1 += 8;
            }
            while(ht_temp > 0);
        }

        //horizontal q-pel filter
        {
            __m128i res_r0r1_16x8b_128, src_r0r1_vpel_16x8b,res_16x8b;
            __m256i src_r0r1_sht_32x8b, src_r0r1_32x8b,src_r0r1_t1_32x8b;

            __m256i res_r0r1_t1_8x32b, res_r0r1_t2_8x32b, res_r0r1_t3_8x32b;
            __m256i res_r0r1_32x8b;

            //Row0 : a0 a1 a2 a3 a4 a5 a6 a7 a8 a9.....
            //Row1 : b0 b1 b2 b3 b4 b5 b6 b7 b8 b9.....

            do
            {
                src_r0r1_vpel_16x8b = _mm_loadl_epi64((__m128i *)pu1_tmp2);
                src_r0r1_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_pred_horiz + src_strd), (__m128i *)pu1_pred_horiz);

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 1);
                src_r0r1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi64(src_r0r1_32x8b, zero);                //a0 a1 a1 a2 a2 a3 a3 a4 b0 b1 b1 b2 b2 b3 b3 b4
                res_r0r1_t1_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff0_1_32x8b);    //a0*c0+a1*c1 a1*c0+a2*c1 a2*c0+a3*c1 a3*c0+a4*c1
                                                                                                //b0*c0+b1*c1 b1*c0+b2*c1 b2*c0+b3*c1 b3*c0+b4*c1

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 4);                              //a2 a3 a3 a4 a4 a5 a5 a6 a6 a7 a7 a8  0  0  0  0
                                                                                                 //b2 b3 b3 b4 b4 b5 b5 b6 b6 b7 b7 b8  0  0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi64(src_r0r1_32x8b, zero);                //a2 a3 a3 a4 a4 a5 a5 a6 b2 b3 b3 b4 b4 b5 b5 b6
                res_r0r1_t2_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff2_3_32x8b);    //a2*c2+a3*c3 a3*c2+a4*c3 a4*c2+a5*c3 a5*c2+a6*c3
                                                                                                //b2*c2+b3*c3 b3*c2+b4*c3 b4*c2+b5*c3 b5*c2+b6*c3

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 4);                              //a4 a5 a5 a6 a6 a7 a7 a8  0  0  0  0  0  0  0  0
                                                                                                //b4 b5 b5 b6 b6 b7 b7 b8  0  0  0  0  0  0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi64(src_r0r1_32x8b, zero);                //a4 a5 a5 a6 a6 a7 a7 a8 b4 b5 b5 b6 b6 b7 b7 b8
                res_r0r1_t3_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff4_5_32x8b);    //a4*c4+a5*c5 a5*c4+a6*c5 a6*c4+a7*c5 a7*c4+a8*c5
                                                                                                //b4*c4+b5*c5 b5*c4+b6*c5 b4*c6+b7*c5 b7*c4+b8*c5

                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t2_8x32b);
                res_r0r1_t3_8x32b = _mm256_add_epi16(const_val16_8x32b, res_r0r1_t3_8x32b);
                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t3_8x32b);     //a0*c0+a1*c1+a2*c2+a3*c3+a4*a4+a5*c5 + 15;
                                                                                                //a1*c0+a2*c1+a2*c2+a3*c3+a5*a4+a6*c5 + 15;
                                                                                                //a2*c0+a3*c1+a4*c2+a5*c3+a6*a4+a7*c5 + 15;
                                                                                                //a3*c0+a4*c1+a5*c2+a6*c3+a6*a4+a8*c5 + 15;
                                                                                                //b0*c0+b1*c1+b2*c2+b3*c3+b4*b4+b5*c5 + 15;
                                                                                                //b1*c0+b2*c1+b2*c2+b3*c3+b5*b4+b6*c5 + 15;
                                                                                                //b2*c0+b3*c1+b4*c2+b5*c3+b6*b4+b7*c5 + 15;
                                                                                                //b3*c0+b4*c1+b5*c2+b6*c3+b6*b4+b8*c5 + 15;
                res_r0r1_t1_8x32b = _mm256_srai_epi16(res_r0r1_t1_8x32b, 5);                       //shifting right by 5 bits.

                res_r0r1_16x8b_128 = _mm256_castsi256_si128(_mm256_packus_epi16(res_r0r1_t1_8x32b,
                                                                           res_r0r1_t1_8x32b));

                res_16x8b = _mm_avg_epu8(res_r0r1_16x8b_128,src_r0r1_vpel_16x8b);

                *((WORD32 *)(pu1_dst)) = _mm_cvtsi128_si32(res_16x8b);
                res_16x8b = _mm_srli_si128(res_16x8b, 4);
                *((WORD32 *)(pu1_dst + dst_strd)) = _mm_cvtsi128_si32(res_16x8b);

                ht -= 2;
                pu1_pred_horiz += src_strd << 1;
                pu1_tmp2 += 8;
                pu1_dst += dst_strd << 1;
            }
            while(ht > 0);
        }
    }
    else if(wd == 8)
    {
        //vertical q-pel filter
        {
            __m128i src_r0_16x8b, src_r1_16x8b, src_r2_16x8b, src_r3_16x8b;
            __m128i src_r4_16x8b, src_r5_16x8b, src_r6_16x8b;
            __m128i src_r0r1_16x8b, src_r2r3_16x8b, src_r4r5_16x8b;

            __m128i res_16x8b, res_t1_8x16b, res_t2_8x16b, res_t3_8x16b;

            //epilogue: Load all the pred rows except sixth  and seventh row for the
            //first and second row processing.
            src_r0_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;

            src_r1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r0_16x8b = _mm_unpacklo_epi64(src_r0_16x8b, src_r1_16x8b);

            src_r2_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r1_16x8b = _mm_unpacklo_epi64(src_r1_16x8b, src_r2_16x8b);

            src_r3_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r2_16x8b = _mm_unpacklo_epi64(src_r2_16x8b, src_r3_16x8b);

            src_r4_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
            pu1_pred_vert = pu1_pred_vert + src_strd;
            src_r3_16x8b = _mm_unpacklo_epi64(src_r3_16x8b, src_r4_16x8b);

            //Core Loop: Process all the rows.
            do
            {
                src_r5_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert));
                src_r4_16x8b = _mm_unpacklo_epi64(src_r4_16x8b, src_r5_16x8b);

                src_r6_16x8b = _mm_loadl_epi64((__m128i *)(pu1_pred_vert + src_strd));
                src_r5_16x8b = _mm_unpacklo_epi64(src_r5_16x8b, src_r6_16x8b);

                src_r0r1_16x8b = _mm_unpacklo_epi8(src_r0_16x8b, src_r1_16x8b);
                src_r2r3_16x8b = _mm_unpacklo_epi8(src_r2_16x8b, src_r3_16x8b);
                src_r4r5_16x8b = _mm_unpacklo_epi8(src_r4_16x8b, src_r5_16x8b);

                res_t1_8x16b = _mm_maddubs_epi16(src_r0r1_16x8b, coeff0_1_16x8b);
                res_t2_8x16b = _mm_maddubs_epi16(src_r2r3_16x8b, coeff2_3_16x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t2_8x16b);
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t3_8x16b, res_t1_8x16b);

                res_t1_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits.
                res_16x8b = _mm_packus_epi16(res_t1_8x16b, res_t1_8x16b);

                _mm_storel_epi64((__m128i *)(pu1_tmp1), res_16x8b);

                src_r0r1_16x8b = _mm_unpackhi_epi8(src_r0_16x8b, src_r1_16x8b);
                src_r2r3_16x8b = _mm_unpackhi_epi8(src_r2_16x8b, src_r3_16x8b);
                src_r4r5_16x8b = _mm_unpackhi_epi8(src_r4_16x8b, src_r5_16x8b);

                res_t1_8x16b = _mm_maddubs_epi16(src_r0r1_16x8b, coeff0_1_16x8b);
                res_t2_8x16b = _mm_maddubs_epi16(src_r2r3_16x8b, coeff2_3_16x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t2_8x16b);
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t3_8x16b, res_t1_8x16b);

                res_t1_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits.
                res_16x8b = _mm_packus_epi16(res_t1_8x16b, res_t1_8x16b);

                _mm_storel_epi64((__m128i *)(pu1_tmp1 + 8), res_16x8b);

                src_r0_16x8b = src_r2_16x8b;
                src_r1_16x8b = src_r3_16x8b;
                src_r2_16x8b = src_r4_16x8b;
                src_r3_16x8b = src_r5_16x8b;
                src_r4_16x8b = src_r6_16x8b;

                ht_temp -= 2;
                pu1_pred_vert += src_strd << 1;
                pu1_tmp1 += 16;
            }
            while(ht_temp > 0);
        }

        //horizontal q-pel filter
        {

            __m256i src_r0r1_32x8b, src_r0r1_sht_32x8b,src_r0r1_t1_32x8b,res_32x8b;
            __m128i src_r0r1_vpel_128,res_16x8b_128;
            __m256i src_r0r1_vpel_32x8b,res_r0r1_t1_8x32b, res_r0r1_t2_8x32b, res_r0r1_t3_8x32b;

            //Row0 : a0 a1 a2 a3 a4 a5 a6 a7 a8 a9.....
            //Row1 : b0 b1 b2 b3 b4 b5 b6 b7 b8 b9.....

            do
            {
                src_r0r1_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_pred_horiz + src_strd),    //a2 a3 a4 a5 a6 a7 a8....a15 0 or
                                                     (__m128i *)pu1_pred_horiz);                //a3 a4 a5 a6 a7 a8 a9....a15 0
                                                                                                //b2 b3 b4 b5 b6 b7 b8....b15 0 or
                                                                                                //b3 b4 b5 b6 b7 b8 b9....b15 0

                src_r0r1_vpel_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_tmp2 + 8),
                                                         (__m128i *)(pu1_tmp2));
                src_r0r1_vpel_32x8b = _mm256_unpacklo_epi64(src_r0r1_vpel_32x8b,zero);

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 1);                        //a1 a2 a3 a4 a5 a6 a7 a8 a9....a15 0
                                                                                              //b1 b2 b3 b4 b5 b6 b7 b8 b9....b15 0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);       //a0 a1 a1 a2 a2 a3 a3 a4 a4 a5 a5 a6 a6 a7 a7 a8

                res_r0r1_t1_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff0_1_32x8b);        //a0*c0+a1*c1 a1*c0+a2*c1 a2*c0+a3*c1 a3*c0+a4*c1
                                                                                                    //a4*c0+a5*c1 a5*c0+a6*c1 a6*c0+a7*c1 a7*c0+a8*c1

                                                                                                    //b4*c0+b5*c1 b5*c0+b6*c1 b6*c0+b7*c1 b7*c0+b8*c1

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 2);                              //a2 a3 a4 a5 a6 a7 a8 a9....a15 0 0
                                                                                                    //b2 b3 b4 b5 b6 b7 b8 b9....b15 0 0

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_sht_32x8b, 2);                      //a3 a4 a5 a6 a7 a8 a9....a15 0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);       //a2 a3 a3 a4 a4 a5 a5 a6 a6 a7 a7 a8 a8 a9 a9 a10
                                                                                                    //b2 b3 b3 b4 b4 b5 b5 b6 b6 b7 b7 b8 a8 a9 a9 a10

                res_r0r1_t2_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff2_3_32x8b);      //a2*c2+a3*c3 a3*c2+a4*c3 a4*c2+a5*c3 a5*c2+a6*c3
                                                                                                  //a6*c2+a7*c3 a7*c2+a8*c3 a8*c2+a9*c3 a9*c2+a10*c3
                                                                                                  //b2*c2+b3*c3 b3*c2+b4*c3 b2*c4+b5*c3 b5*c2+b6*c3
                                                                                                  //b6*c2+b7*c3 b7*c2+b8*c3 b8*c2+b9*c3 b9*c2+b10*c3

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 2);                            //a4 a5 a6 a7 a8 a9....a15 0  0  0  0

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_sht_32x8b, 2);                    //a5 a6 a7 a8 a9....a15 0  0  0  0  0
                                                                                                  //b5 b6 b7 b8 b9....b15 0  0  0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);       //a4 a5 a5 a6 a6 a7 a7 a8 a8 a9 a9 a10 a10 a11 a11 a12
                                                                                                    //b4 b5 b5 b6 b6 b7 b7 b8 b8 b9 b9 b10 b10 b11 b11 b12

                res_r0r1_t3_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff4_5_32x8b);      //a4*c4+a5*c5 a5*c4+a6*c5 a6*c4+a7*c5 a7*c4+a8*c5
                                                                                                  //a8*c4+a9*c5 a9*c4+a10*c5 a10*c4+a11*c5 a11*c4+a12*c5
                                                                                                  //b4*c4+b5*c5 b5*c4+b6*c5 b6*c4+b7*c5 b7*c4+b8*c5

                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t2_8x32b);
                res_r0r1_t3_8x32b = _mm256_add_epi16(const_val16_8x32b, res_r0r1_t3_8x32b);
                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t3_8x32b);
                res_r0r1_t1_8x32b = _mm256_srai_epi16(res_r0r1_t1_8x32b, 5);                      //shifting right by 5 bits.

                res_32x8b = _mm256_packus_epi16(res_r0r1_t1_8x32b, res_r0r1_t1_8x32b);

                res_32x8b = _mm256_avg_epu8(res_32x8b,src_r0r1_vpel_32x8b);

                _mm_storel_epi64((__m128i *)(pu1_dst), _mm256_castsi256_si128(res_32x8b));
                _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd), _mm256_castsi256_si128(_mm256_permute2x128_si256(res_32x8b,res_32x8b,0x1)));

                ht -= 2;
                pu1_pred_horiz += src_strd << 1;
                pu1_dst += dst_strd << 1;
                pu1_tmp2 += 16;
            }
            while(ht > 0);
        }
    }
    else // wd == 16
    {
       //vertical q-pel filter
        {
            __m256i src_r0r2_32x8b, src_r1r3_32x8b, src_r2r4_32x8b,src_32x8b;
            __m128i src_r2_16x8b,src_r4_16x8b,src_r5_16x8b,src_r6_16x8b,src_r4r5_16x8b;

            __m256i res_t1_8x32b;
            __m128i res_t0_8x16b,res_t3_8x16b,res_t1_8x16b,res_16x8b;

            //epilogue: Load all the pred rows except sixth  and seventh row for the
            //first and second row processing.
            src_r0r2_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_pred_vert + 2 * src_strd ),
                                               (__m128i *)(pu1_pred_vert));

            src_r2_16x8b = _mm_loadu_si128((__m128i *)(pu1_pred_vert + 2 * src_strd));
            src_r1r3_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_pred_vert + 3 * src_strd ),
                                               (__m128i *)(pu1_pred_vert + src_strd));
            src_r4_16x8b = _mm_loadu_si128((__m128i *)(pu1_pred_vert + 4 * src_strd));
            pu1_pred_vert =  pu1_pred_vert + (5 * src_strd ) ;


            //Core Loop: Process all the rows.
            do
            {
                src_r5_16x8b =  _mm_loadu_si128((__m128i *)(pu1_pred_vert));
                src_r6_16x8b =  _mm_loadu_si128((__m128i *)(pu1_pred_vert + src_strd));
                src_r2r4_32x8b = _mm256_set_m128i(src_r4_16x8b,src_r2_16x8b);

                src_32x8b = _mm256_unpacklo_epi8(src_r0r2_32x8b, src_r1r3_32x8b);
                src_r4r5_16x8b =  _mm_unpacklo_epi8(src_r4_16x8b, src_r5_16x8b);
                
                res_t1_8x32b = _mm256_maddubs_epi16(src_32x8b, coeff_32x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(_mm256_castsi256_si128(res_t1_8x32b), _mm256_extracti128_si256                                                (res_t1_8x32b,0x1));
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t3_8x16b);
                res_t0_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits


                src_32x8b = _mm256_unpackhi_epi8(src_r0r2_32x8b, src_r1r3_32x8b);
                src_r4r5_16x8b = _mm_unpackhi_epi8(src_r4_16x8b, src_r5_16x8b);

                res_t1_8x32b = _mm256_maddubs_epi16(src_32x8b, coeff_32x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(_mm256_castsi256_si128(res_t1_8x32b),_mm256_extracti128_si256                                                (res_t1_8x32b,0x1));
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits.

                res_16x8b = _mm_packus_epi16(res_t0_8x16b, res_t1_8x16b);
                _mm_storeu_si128((__m128i *)(pu1_tmp1), res_16x8b);

                src_32x8b = _mm256_unpacklo_epi8(src_r1r3_32x8b, src_r2r4_32x8b);
                src_r4r5_16x8b =  _mm_unpacklo_epi8(src_r5_16x8b, src_r6_16x8b);

                res_t1_8x32b = _mm256_maddubs_epi16(src_32x8b, coeff_32x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(_mm256_castsi256_si128(res_t1_8x32b), _mm256_extracti128_si256                                                (res_t1_8x32b,0x1));
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                 res_t0_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits


                src_32x8b = _mm256_unpackhi_epi8(src_r1r3_32x8b, src_r2r4_32x8b);
                src_r4r5_16x8b =  _mm_unpackhi_epi8(src_r5_16x8b, src_r6_16x8b);

                res_t1_8x32b = _mm256_maddubs_epi16(src_32x8b, coeff_32x8b);
                res_t3_8x16b = _mm_maddubs_epi16(src_r4r5_16x8b, coeff4_5_16x8b);

                res_t1_8x16b = _mm_add_epi16(_mm256_castsi256_si128(res_t1_8x32b), _mm256_extracti128_si256                                                (res_t1_8x32b,0x1));
                res_t3_8x16b = _mm_add_epi16(const_val16_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_add_epi16(res_t1_8x16b, res_t3_8x16b);
                res_t1_8x16b = _mm_srai_epi16(res_t1_8x16b, 5); //shifting right by 5 bits


                res_16x8b = _mm_packus_epi16(res_t0_8x16b, res_t1_8x16b);
                _mm_storeu_si128((__m128i *)(pu1_tmp1+16), res_16x8b);
                 src_r0r2_32x8b = src_r2r4_32x8b;
                src_r1r3_32x8b = _mm256_set_m128i(src_r5_16x8b,_mm256_extracti128_si256(src_r1r3_32x8b,0x1));
                src_r2_16x8b = src_r4_16x8b;
                src_r4_16x8b = src_r6_16x8b;

                ht_temp -= 2;
                pu1_pred_vert += src_strd << 1;
                pu1_tmp1 += 32;
            }
            while(ht_temp > 0);
        }
        //horizontal q-pel filter
        {

            __m256i src_r0r1_32x8b,src_r0r1_sht_32x8b,src_r0r1_t1_32x8b,res_32x8b;
            __m128i src_vpel_16x8b,res_16x8b_128;

            __m256i res_r0r1_t1_8x32b, res_r0r1_t2_8x32b, res_r0r1_t3_8x32b;

            //Row0 : a0 a1 a2 a3 a4 a5 a6 a7 a8 a9.....
            //Row0 :                         b0 b1 b2 b3 b4 b5 b6 b7 b8 b9.....
            //b0 is same a8. Similarly other bn pixels are same as a(n+8) pixels.

            do
            {
                src_r0r1_32x8b = _mm256_loadu2_m128i((__m128i *)(pu1_pred_horiz + 8),
                                                    (__m128i *)(pu1_pred_horiz));        //a0 a1 a2 a3 a4 a5 a6 a7 a8 a9....a15
                                                                                         //b0 b1 b2 b3 b4 b5 b6 b7 b8 b9....b15
                src_vpel_16x8b = _mm_loadu_si128((__m128i *)(pu1_tmp2));

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 1);                      //a1 a2 a3 a4 a5 a6 a7 a8 a9....a15 0
                                                                                             //b1 b2 b3 b4 b5 b6 b7 b8 b9....b15 0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);     //a0 a1 a1 a2 a2 a3 a3 a4 a4 a5 a5 a6 a6 a7 a7 a8
                                                                                                  //b0 b1 b1 b2 b2 b3 b3 b4 b4 b5 b5 b6 b6 b7 b7 b8

                res_r0r1_t1_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff0_1_32x8b);    //a0*c0+a1*c1 a1*c0+a2*c1 a2*c0+a3*c1 a3*c0+a4*c1
                                                                                                //a4*c0+a5*c1 a5*c0+a6*c1 a6*c0+a7*c1 a7*c0+a8*c1
                                                                                                //b0*c0+b1*c1 b1*c0+b2*c1 b2*c0+b3*c1 b3*c0+b4*c1
                                                                                                //b4*c0+b5*c1 b5*c0+b6*c1 b6*c0+b7*c1 b7*c0+b8*c1

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 2);                          //a2 a3 a4 a5 a6 a7 a8 a9....a15 0 0
                                                                                                //b2 b3 b4 b5 b6 b7 b8 b9....b15 0 0

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_sht_32x8b, 2);                  //a3 a4 a5 a6 a7 a8 a9....a15 0  0  0
                                                                                                //b3 b4 b5 b6 b7 b8 b9....b15 0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);     //a2 a3 a3 a4 a4 a5 a5 a6 a6 a7 a7 a8 a8 a9 a9 a10
                                                                                                //b2 b3 b3 b4 b4 b5 b5 b6 b6 b7 b7 b8 a8 a9 a9 a10

                res_r0r1_t2_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff2_3_32x8b);    //a2*c2+a3*c3 a3*c2+a4*c3 a4*c2+a5*c3 a5*c2+a6*c3
                                                                                         //a6*c2+a7*c3 a7*c2+a8*c3 a8*c2+a9*c3 a9*c2+a10*c3
                                                                                         //b2*c2+b3*c3 b3*c2+b4*c3 b2*c4+b5*c3 b5*c2+b6*c3
                                                                                         //b6*c2+b7*c3 b7*c2+b8*c3 b8*c2+b9*c3 b9*c2+b10*c3

                src_r0r1_32x8b = _mm256_srli_si256(src_r0r1_32x8b, 2);                          //a4 a5 a6 a7 a8 a9....a15 0  0  0  0
                                                                                                //b4 b5 b6 b7 b8 b9....b15 0  0  0  0

                src_r0r1_sht_32x8b = _mm256_srli_si256(src_r0r1_sht_32x8b, 2);                  //a5 a6 a7 a8 a9....a15 0  0  0  0  0
                                                                                                //b5 b6 b7 b8 b9....b15 0  0  0  0  0

                src_r0r1_t1_32x8b = _mm256_unpacklo_epi8(src_r0r1_32x8b, src_r0r1_sht_32x8b);     //a4 a5 a5 a6 a6 a7 a7 a8 a8 a9 a9 a10 a10 a11 a11 a12
                                                                                             //b4 b5 b5 b6 b6 b7 b7 b8 b8 b9 b9 b10 b10 b11 b11 b12

                res_r0r1_t3_8x32b = _mm256_maddubs_epi16(src_r0r1_t1_32x8b, coeff4_5_32x8b);    //a4*c4+a5*c5 a5*c4+a6*c5 a6*c4+a7*c5 a7*c4+a8*c5
                                                                                         //a8*c4+a9*c5 a9*c4+a10*c5 a10*c4+a11*c5 a11*c4+a12*c5
                                                                                         //b4*c4+b5*c5 b5*c4+b6*c5 b6*c4+b7*c5 b7*c4+b8*c5
                                                                                         //b8*c4+b9*c5 b9*c4+b10*c5 b10*c4+b11*c5 b11*c4+b12*c5
                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t2_8x32b);
                res_r0r1_t3_8x32b = _mm256_add_epi16(const_val16_8x32b, res_r0r1_t3_8x32b);
                res_r0r1_t1_8x32b = _mm256_add_epi16(res_r0r1_t1_8x32b, res_r0r1_t3_8x32b);
                res_r0r1_t1_8x32b = _mm256_srai_epi16(res_r0r1_t1_8x32b, 5);                    //shifting right by 5 bits.


                res_32x8b = _mm256_packus_epi16(res_r0r1_t1_8x32b, res_r0r1_t1_8x32b);
                res_16x8b_128 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(res_32x8b, 0xD8));
                res_16x8b_128 = _mm_avg_epu8(res_16x8b_128, src_vpel_16x8b);
                _mm_storeu_si128((__m128i *)(pu1_dst), res_16x8b_128);

                ht --;
                pu1_pred_horiz  += src_strd;
                pu1_dst += dst_strd;
                pu1_tmp2 += 16;
            }
            while(ht > 0);
        }
    }
}
