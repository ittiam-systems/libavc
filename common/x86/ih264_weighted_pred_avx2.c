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
/* File Includes                                                             */
/*****************************************************************************/

#include <immintrin.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_weighted_pred.h"
#include <stdint.h>
#include <string.h>

#include <stdio.h>


/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_weighted_bi_pred_luma_avx2                         */
/*                                                                           */
/*  Description   : This function performs the weighted biprediction as      */
/*                  described in sec 8.4.2.3.2 titled "Weighted sample       */
/*                  prediction process" for luma. The function gets two      */
/*                  ht x wd blocks, weights them, adds them, rounds off the  */
/*                  sum, offsets it, saturates it to unsigned 8-bit and      */
/*                  stores it in the destination block. (ht,wd) can be       */
/*                  (4,4), (8,4), (4,8), (8,8), (16,8), (8,16) or (16,16).   */
/*                                                                           */
/*  Inputs        : pu1_src1  - Pointer to source 1                          */
/*                  pu1_src2  - Pointer to source 2                          */
/*                  pu1_dst   - Pointer to destination                       */
/*                  src_strd1 - stride for source 1                          */
/*                  src_strd2 - stride for source 2                          */
/*                  dst_strd2 - stride for destination                       */
/*                  log_wd    - number of bits to be rounded off             */
/*                  wt1       - weight value for source 1                    */
/*                  wt2       - weight value for source 2                    */
/*                  ofst1     - offset value for source 1                    */
/*                  ofst2     - offset value for source 2                    */
/*                  ht        - height of the block                          */
/*                  wd        - width of the block                           */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         04 02 2015   Kaushik         Initial Version                      */
/*                      Senthoor                                             */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_weighted_bi_pred_luma_avx2(UWORD8 *pu1_src1,
                                       UWORD8 *pu1_src2,
                                       UWORD8 *pu1_dst,
                                       WORD32 src_strd1,
                                       WORD32 src_strd2,
                                       WORD32 dst_strd,
                                       WORD32 log_wd,
                                       WORD32 wt1,
                                       WORD32 wt2,
                                       WORD32 ofst1,
                                       WORD32 ofst2,
                                       WORD32 ht,
                                       WORD32 wd)
{

    __m256i wt1_8x32b, wt2_8x32b;
    __m256i ofst_8x32b, round_8x32b;
    __m256i zero;
    zero = _mm256_set1_epi8(0);

    WORD32 ofst;
    WORD32 round_val, shft;

    wt1 = (WORD16)(wt1 & 0xffff);
    wt2 = (WORD16)(wt2 & 0xffff);
    round_val = 1 << log_wd;
    shft = log_wd + 1;
    ofst1 = (WORD8)(ofst1 & 0xff);
    ofst2 = (WORD8)(ofst2 & 0xff);
    ofst = (ofst1 + ofst2 + 1) >> 1;

    wt1_8x32b = _mm256_set1_epi16(wt1);
    wt2_8x32b = _mm256_set1_epi16(wt2);
    round_8x32b = _mm256_set1_epi16(round_val);
    ofst_8x32b = _mm256_set1_epi16(ofst);


    if(wd == 4)
    {
        __m128i y1_2_16x8b, y1_3_16x8b;
        __m128i y2_2_16x8b, y2_3_16x8b;

        __m256i y1_02_32x8b,y1_13_32x8b,y2_02_32x8b,y2_13_32x8b,y1_0_32x8b,y2_0_32x8b,y1_0_8x32b,y2_1_8x32b,y2_0_8x32b;
        __m128i y1_0_16x8b_128,y2_0_16x8b_128,y1_1_16x8b_128,y1_2_16x8b_128,y1_3_16x8b_128;

        do
        {
            y1_02_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src1 + (src_strd1 << 1)), (__m128i *)(pu1_src1));
            y1_13_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src1 + src_strd1 * 3), (__m128i *)(pu1_src1 + src_strd1));

            y2_02_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src2 + (src_strd2 << 1)), (__m128i *)(pu1_src2));
            y2_13_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src2 + src_strd2 * 3), (__m128i *)(pu1_src2 + src_strd2));

            y1_02_32x8b = _mm256_unpacklo_epi64(y1_02_32x8b, zero);
            y1_13_32x8b = _mm256_unpacklo_epi64(y1_13_32x8b, zero);
            y2_02_32x8b = _mm256_unpacklo_epi64(y2_02_32x8b, zero);
            y2_13_32x8b = _mm256_unpacklo_epi64(y2_13_32x8b, zero);

            y1_0_32x8b = _mm256_unpacklo_epi32(y1_02_32x8b, y1_13_32x8b);
            y2_0_32x8b = _mm256_unpacklo_epi32(y2_02_32x8b, y2_13_32x8b);
            y1_0_16x8b_128 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(y1_0_32x8b, 0xD8));
            y2_0_16x8b_128 = _mm256_castsi256_si128(_mm256_permute4x64_epi64(y2_0_32x8b, 0xD8));

            y1_0_8x32b = _mm256_cvtepu8_epi16(y1_0_16x8b_128); // 8 to 16
            y2_0_8x32b = _mm256_cvtepu8_epi16(y2_0_16x8b_128);

            y1_0_8x32b = _mm256_mullo_epi16(y1_0_8x32b, wt1_8x32b);
            y2_0_8x32b = _mm256_mullo_epi16(y2_0_8x32b, wt2_8x32b);

            y1_0_8x32b = _mm256_adds_epi16(y1_0_8x32b, y2_0_8x32b);

            y1_0_8x32b = _mm256_srai_epi16(y1_0_8x32b, shft);

            y1_0_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0_8x32b);

            y1_0_16x8b_128 = _mm256_castsi256_si128(_mm256_packus_epi16(y1_0_8x32b, y1_0_8x32b));
            y1_2_16x8b_128 = _mm_srli_si128(y1_0_16x8b_128, 4);
            y1_1_16x8b_128 = _mm_srli_si128(y1_0_16x8b_128, 8);
            y1_3_16x8b_128 = _mm_srli_si128(y1_0_16x8b_128, 12);

            *((WORD32 *)(pu1_dst)) = _mm_cvtsi128_si32(y1_0_16x8b_128);
            *((WORD32 *)(pu1_dst + dst_strd)) = _mm_cvtsi128_si32(y1_1_16x8b_128);
            *((WORD32 *)(pu1_dst + (dst_strd << 1))) = _mm_cvtsi128_si32(y1_2_16x8b_128);
            *((WORD32 *)(pu1_dst + dst_strd * 3)) = _mm_cvtsi128_si32(y1_3_16x8b_128);

            ht -= 4;
            pu1_src1 += src_strd1 << 2;
            pu1_src2 += src_strd2 << 2;
            pu1_dst += dst_strd << 2;
        }
        while(ht > 0);
    }
    else if(wd == 8)
    {
        __m128i y1_0_16x8b_128,y2_0_16x8b_128,y1_2_16x8b_128,y1_1_16x8b_128,y1_3_16x8b_128;
        __m256i y1_02_32x8b,y1_13_32x8b,y2_02_32x8b,y2_13_32x8b,y1_0_32x8b,y2_0_32x8b,y1_0_8x32b;
        __m256i y1_1_8x32b,y2_0_8x32b,y2_1_8x32b;

        do
        {

            y1_02_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src1 + (src_strd1 << 1)), (__m128i *)(pu1_src1));
            y1_13_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src1 + src_strd1 * 3), (__m128i *)(pu1_src1 + src_strd1));

            y2_02_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src2 + (src_strd2 << 1)), (__m128i *)(pu1_src2));
            y2_13_32x8b =  _mm256_loadu2_m128i((__m128i *)(pu1_src2 + src_strd2 * 3), (__m128i *)(pu1_src2 + src_strd2));

            y1_02_32x8b = _mm256_unpacklo_epi64(y1_02_32x8b, zero);
            y1_13_32x8b = _mm256_unpacklo_epi64(y1_13_32x8b, zero);
            y2_02_32x8b = _mm256_unpacklo_epi64(y2_02_32x8b, zero);
            y2_13_32x8b = _mm256_unpacklo_epi64(y2_13_32x8b, zero);

            y1_0_32x8b = _mm256_unpacklo_epi64(y1_02_32x8b, y1_13_32x8b);
            y2_0_32x8b = _mm256_unpacklo_epi64(y2_02_32x8b, y2_13_32x8b);

            y1_0_8x32b = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(y1_0_32x8b));
            y1_0_16x8b_128 = _mm256_castsi256_si128(_mm256_permute2x128_si256(y1_0_32x8b,y1_0_32x8b,0x1));
            y1_1_8x32b = _mm256_cvtepu8_epi16(y1_0_16x8b_128);

            y2_0_8x32b = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(y2_0_32x8b));
            y2_0_16x8b_128 = _mm256_castsi256_si128(_mm256_permute2x128_si256(y2_0_32x8b,y2_0_32x8b,0x1));
            y2_1_8x32b = _mm256_cvtepu8_epi16(y2_0_16x8b_128);


            y1_0_8x32b = _mm256_mullo_epi16(y1_0_8x32b, wt1_8x32b);
            y2_0_8x32b = _mm256_mullo_epi16(y2_0_8x32b, wt2_8x32b);
            y1_1_8x32b = _mm256_mullo_epi16(y1_1_8x32b, wt1_8x32b);
            y2_1_8x32b = _mm256_mullo_epi16(y2_1_8x32b, wt2_8x32b);

            y1_0_8x32b = _mm256_adds_epi16(y1_0_8x32b, y2_0_8x32b);
            y1_1_8x32b = _mm256_adds_epi16(y1_1_8x32b, y2_1_8x32b);

            y1_0_8x32b = _mm256_srai_epi16(y1_0_8x32b, shft);
            y1_1_8x32b = _mm256_srai_epi16(y1_1_8x32b, shft);

            y1_0_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0_8x32b);
            y1_1_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_1_8x32b);

            y1_0_32x8b = _mm256_packus_epi16(y1_0_8x32b, y1_1_8x32b);
            y1_0_16x8b_128 = _mm256_castsi256_si128(y1_0_32x8b);
            y1_2_16x8b_128 = _mm256_castsi256_si128(_mm256_srli_si256(y1_0_32x8b, 8));

            y1_0_32x8b = _mm256_permute2x128_si256(y1_0_32x8b,y1_0_32x8b,1);
            y1_1_16x8b_128 = _mm256_castsi256_si128(y1_0_32x8b);
            y1_3_16x8b_128 = _mm256_castsi256_si128(_mm256_srli_si256(y1_0_32x8b, 8));

            _mm_storel_epi64((__m128i *)pu1_dst, y1_0_16x8b_128);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd), y1_1_16x8b_128);
            _mm_storel_epi64((__m128i *)(pu1_dst + (dst_strd << 1)), y1_2_16x8b_128);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd * 3), y1_3_16x8b_128);

            ht -= 4;
            pu1_src1 += src_strd1 << 2;
            pu1_src2 += src_strd2 << 2;
            pu1_dst += dst_strd << 2;

        }
        while(ht > 0);
    }
    else // wd == 16
    {
        __m256i y1_0L_8x32b, y1_0H_8x32b, y1_1L_8x32b, y1_1H_8x32b;
        __m256i y2_0L_8x32b, y2_0H_8x32b, y2_1L_8x32b, y2_1H_8x32b;

        __m256i zero_32x8b,y1_0_32x8b,y2_0_32x8b;
        zero_32x8b = _mm256_set1_epi8(0);

        do
        {

            y1_0_32x8b = _mm256_loadu_si256((__m256i *)pu1_src1);
            y2_0_32x8b = _mm256_loadu_si256((__m256i *)pu1_src2);

            y1_0L_8x32b = _mm256_unpacklo_epi8(y1_0_32x8b, zero_32x8b);
            y1_0H_8x32b = _mm256_unpackhi_epi8(y1_0_32x8b, zero_32x8b);

            y2_0L_8x32b = _mm256_unpacklo_epi8(y2_0_32x8b,zero_32x8b);
            y2_0H_8x32b = _mm256_unpackhi_epi8(y2_0_32x8b, zero_32x8b);

            y1_0L_8x32b = _mm256_mullo_epi16(y1_0L_8x32b, wt1_8x32b);
            y1_0H_8x32b = _mm256_mullo_epi16(y1_0H_8x32b, wt1_8x32b);

            y2_0L_8x32b = _mm256_mullo_epi16(y2_0L_8x32b, wt2_8x32b);
            y2_0H_8x32b = _mm256_mullo_epi16(y2_0H_8x32b, wt2_8x32b);

            y1_0L_8x32b = _mm256_adds_epi16(y1_0L_8x32b, y2_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(y1_0H_8x32b, y2_0H_8x32b);

            y1_0L_8x32b = _mm256_adds_epi16(round_8x32b, y1_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(round_8x32b, y1_0H_8x32b);

            y1_0L_8x32b = _mm256_srai_epi16(y1_0L_8x32b, shft);
            y1_0H_8x32b = _mm256_srai_epi16(y1_0H_8x32b, shft);

            y1_0L_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0H_8x32b);

            y1_0_32x8b = _mm256_packus_epi16(y1_0L_8x32b, y1_0H_8x32b);

            _mm256_storeu_si256((__m256i *)pu1_dst, y1_0_32x8b);

            ht -= 2;
            pu1_src1 += src_strd1 << 1;
            pu1_src2 += src_strd2 << 1;
            pu1_dst += dst_strd << 1;
        }
        while(ht > 0);
    }
}


/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_weighted_bi_pred_chroma_avx2                       */
/*                                                                           */
/*  Description   : This function performs the weighted biprediction as      */
/*                  described in sec 8.4.2.3.2 titled "Weighted sample       */
/*                  prediction process" for chroma. The function gets two    */
/*                  ht x wd blocks, weights them, adds them, rounds off the  */
/*                  sum, offsets it, saturates it to unsigned 8-bit and      */
/*                  stores it in the destination block. (ht,wd) can be       */
/*                  (2,2), (4,2), (2,4), (4,4), (8,4), (4,8) or (8,8).       */
/*                                                                           */
/*  Inputs        : pu1_src1  - Pointer to source 1                          */
/*                  pu1_src2  - Pointer to source 2                          */
/*                  pu1_dst   - Pointer to destination                       */
/*                  src_strd1 - stride for source 1                          */
/*                  src_strd2 - stride for source 2                          */
/*                  dst_strd2 - stride for destination                       */
/*                  log_wd    - number of bits to be rounded off             */
/*                  wt1       - weight values for u and v in source 1        */
/*                  wt2       - weight values for u and v in source 2        */
/*                  ofst1     - offset value for u and v in source 1         */
/*                  ofst2     - offset value for u and v in source 2         */
/*                  ht        - height of the block                          */
/*                  wd        - width of the block                           */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         04 02 2015   Kaushik         Initial Version                      */
/*                      Senthoor                                             */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_weighted_bi_pred_chroma_avx2(UWORD8 *pu1_src1,
                                         UWORD8 *pu1_src2,
                                         UWORD8 *pu1_dst,
                                         WORD32 src_strd1,
                                         WORD32 src_strd2,
                                         WORD32 dst_strd,
                                         WORD32 log_wd,
                                         WORD32 wt1,
                                         WORD32 wt2,
                                         WORD32 ofst1,
                                         WORD32 ofst2,
                                         WORD32 ht,
                                         WORD32 wd)
{

    __m128i y1_0_16x8b, y1_1_16x8b;
    __m128i y2_0_16x8b, y2_1_16x8b;

    __m128i wt1_8x16b, wt2_8x16b;
    __m128i ofst_8x16b, round_8x16b;

    WORD32 ofst1_u, ofst2_u, ofst_u;
    WORD32 ofst1_v, ofst2_v, ofst_v;
    WORD32 round_val, shft, ofst_val,ofst_val_256;

    round_val = 1 << log_wd;
    shft = log_wd + 1;

    ofst1_u = (WORD8)(ofst1 & 0xff);
    ofst1_v = (WORD8)(ofst1 >> 8);
    ofst2_u = (WORD8)(ofst2 & 0xff);
    ofst2_v = (WORD8)(ofst2 >> 8);

    wt1_8x16b = _mm_set1_epi32(wt1);
    wt2_8x16b = _mm_set1_epi32(wt2);

    ofst_u = (ofst1_u + ofst2_u + 1) >> 1;
    ofst_v = (ofst1_v + ofst2_v + 1) >> 1;
    ofst_val = (ofst_u & 0xffff) | (ofst_v << 16);
    ofst_val_256 = (ofst_u & 0xffff) | (ofst_v << 16);

    round_8x16b = _mm_set1_epi16(round_val);
    ofst_8x16b =  _mm_set1_epi32(ofst_val);

    if(wd == 2)
    {
        __m128i y1_0_8x16b, y2_0_8x16b;

        do
        {
            y1_0_16x8b = _mm_loadl_epi64((__m128i *)pu1_src1);                 //Loading 64 bits from diff location
            y1_1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src1 + src_strd1));

            y2_0_16x8b = _mm_loadl_epi64((__m128i *)pu1_src2);
            y2_1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src2 + src_strd2));

            y1_0_16x8b = _mm_unpacklo_epi32(y1_0_16x8b, y1_1_16x8b);
            y2_0_16x8b = _mm_unpacklo_epi32(y2_0_16x8b, y2_1_16x8b);

            y1_0_8x16b = _mm_cvtepu8_epi16(y1_0_16x8b);
            y2_0_8x16b = _mm_cvtepu8_epi16(y2_0_16x8b);

            y1_0_8x16b = _mm_mullo_epi16(y1_0_8x16b, wt1_8x16b);
            y2_0_8x16b = _mm_mullo_epi16(y2_0_8x16b, wt2_8x16b);

            y1_0_8x16b = _mm_adds_epi16(y1_0_8x16b, y2_0_8x16b);
            y1_0_8x16b = _mm_adds_epi16(round_8x16b, y1_0_8x16b);

            y1_0_8x16b = _mm_srai_epi16(y1_0_8x16b, shft);
            y1_0_8x16b = _mm_adds_epi16(ofst_8x16b, y1_0_8x16b);

            y1_0_16x8b = _mm_packus_epi16(y1_0_8x16b, y1_0_8x16b);
            y1_1_16x8b = _mm_srli_si128(y1_0_16x8b, 4);

            *((WORD32 *)(pu1_dst)) = _mm_cvtsi128_si32(y1_0_16x8b);
            *((WORD32 *)(pu1_dst + dst_strd)) = _mm_cvtsi128_si32(y1_1_16x8b);

            ht -= 2;
            pu1_src1 += src_strd1 << 1;
            pu1_src2 += src_strd2 << 1;
            pu1_dst += dst_strd << 1;
        }
        while(ht > 0);
    }
    else if(wd == 4)
    {
        __m128i y1_0_8x16b, y1_1_8x16b;
        __m128i y2_0_8x16b, y2_1_8x16b;

        do
        {
            y1_0_16x8b = _mm_loadl_epi64((__m128i *)pu1_src1);                  //Loading 64 bits from diff location
            y1_1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src1 + src_strd1));

            y2_0_16x8b = _mm_loadl_epi64((__m128i *)pu1_src2);
            y2_1_16x8b = _mm_loadl_epi64((__m128i *)(pu1_src2 + src_strd2));

            y1_0_8x16b = _mm_cvtepu8_epi16(y1_0_16x8b);
            y1_1_8x16b = _mm_cvtepu8_epi16(y1_1_16x8b);

            y2_0_8x16b = _mm_cvtepu8_epi16(y2_0_16x8b);
            y2_1_8x16b = _mm_cvtepu8_epi16(y2_1_16x8b);

            y1_0_8x16b = _mm_mullo_epi16(y1_0_8x16b, wt1_8x16b);
            y2_0_8x16b = _mm_mullo_epi16(y2_0_8x16b, wt2_8x16b);
            y1_1_8x16b = _mm_mullo_epi16(y1_1_8x16b, wt1_8x16b);
            y2_1_8x16b = _mm_mullo_epi16(y2_1_8x16b, wt2_8x16b);

            y1_0_8x16b = _mm_adds_epi16(y1_0_8x16b, y2_0_8x16b);
            y1_1_8x16b = _mm_adds_epi16(y1_1_8x16b, y2_1_8x16b);

            y1_0_8x16b = _mm_adds_epi16(round_8x16b, y1_0_8x16b);
            y1_1_8x16b = _mm_adds_epi16(round_8x16b, y1_1_8x16b);

            y1_0_8x16b = _mm_srai_epi16(y1_0_8x16b, shft);
            y1_1_8x16b = _mm_srai_epi16(y1_1_8x16b, shft);

            y1_0_8x16b = _mm_adds_epi16(ofst_8x16b, y1_0_8x16b);
            y1_1_8x16b = _mm_adds_epi16(ofst_8x16b, y1_1_8x16b);

            y1_0_16x8b = _mm_packus_epi16(y1_0_8x16b, y1_1_8x16b);
            y1_1_16x8b = _mm_srli_si128(y1_0_16x8b, 8);

            _mm_storel_epi64((__m128i *)pu1_dst, y1_0_16x8b);
            _mm_storel_epi64((__m128i *)(pu1_dst + dst_strd), y1_1_16x8b);

            ht -= 2;
            pu1_src1 += src_strd1 << 1;
            pu1_src2 += src_strd2 << 1;
            pu1_dst += dst_strd << 1;
        }
        while(ht > 0);
    }
    else // wd == 8
    {
        __m256i y1_0L_8x32b, y1_0H_8x32b, y1_1L_8x32b, y1_1H_8x32b;
        __m256i y2_0L_8x32b, y2_0H_8x32b, y2_1L_8x32b, y2_1H_8x32b;
        __m256i y1_0_32x8b,y2_0_32x8b,ofst_8x32b,round_8x32b;
        __m256i wt1_8x32b, wt2_8x32b;
        __m256i zero_32x8b;

        wt1_8x32b = _mm256_set1_epi16(wt1);
        wt2_8x32b = _mm256_set1_epi16(wt2);
        round_8x32b = _mm256_set1_epi16(round_val);
        ofst_8x32b = _mm256_set1_epi32(ofst_val_256);
        zero_32x8b = _mm256_set1_epi8(0);

        do
        {
            y1_0_32x8b = _mm256_loadu_si256((__m256i *)pu1_src1);
            y2_0_32x8b = _mm256_loadu_si256((__m256i *)pu1_src2);
            y1_0L_8x32b = _mm256_unpacklo_epi8(y1_0_32x8b, zero_32x8b);
            y1_0H_8x32b = _mm256_unpackhi_epi8(y1_0_32x8b, zero_32x8b);
            y2_0L_8x32b = _mm256_unpacklo_epi8(y2_0_32x8b, zero_32x8b);
            y2_0H_8x32b = _mm256_unpackhi_epi8(y2_0_32x8b, zero_32x8b);
            y1_0L_8x32b = _mm256_mullo_epi16(y1_0L_8x32b, wt1_8x32b);
            y1_0H_8x32b = _mm256_mullo_epi16(y1_0H_8x32b, wt1_8x32b);

            y2_0L_8x32b = _mm256_mullo_epi16(y2_0L_8x32b, wt2_8x32b);
            y2_0H_8x32b = _mm256_mullo_epi16(y2_0H_8x32b, wt2_8x32b);

            y1_0L_8x32b = _mm256_adds_epi16(y1_0L_8x32b, y2_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(y1_0H_8x32b, y2_0H_8x32b);

            y1_0L_8x32b = _mm256_adds_epi16(round_8x32b, y1_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(round_8x32b, y1_0H_8x32b);

            y1_0L_8x32b = _mm256_srai_epi16(y1_0L_8x32b, shft);
            y1_0H_8x32b = _mm256_srai_epi16(y1_0H_8x32b, shft);

            y1_0L_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0L_8x32b);
            y1_0H_8x32b = _mm256_adds_epi16(ofst_8x32b, y1_0H_8x32b);


            y1_0_32x8b = _mm256_packus_epi16(y1_0L_8x32b, y1_0H_8x32b);
            _mm256_storeu_si256((__m256i *)pu1_dst, y1_0_32x8b);

            ht -= 2;
            pu1_src1 += src_strd1 << 1;
            pu1_src2 += src_strd2 << 1;
            pu1_dst += dst_strd << 1;
        }
        while(ht > 0);
    }
}
