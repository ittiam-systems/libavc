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
/*                                                                           */
/*  File Name         : ih264_deblk_luma_avx2.c                              */
/*                                                                           */
/*  Description       : Contains function definitions for deblocking         */
/*                                                                           */
/*  List of Functions : ih264_deblk_luma_horz_bslt4_avx2()                   */
/*                      ih264_deblk_luma_vert_bslt4_avx2()                   */
/*                                                                           */
/*  Issues / Problems : None                                                 */
/*                                                                           */
/*  Revision History  :                                                      */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         12 02 2015   Naveen Kumar P  Added luma deblocking ssse3          */
/*                                      intrinsics                           */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#ifdef __ANDROID__
#include "log/log.h"
#include <cutils/log.h>
#endif

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_platform_macros.h"
#include "ih264_deblk_edge_filters.h"
#include "ih264_macros.h"


/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_deblk_luma_horz_bslt4_avx2()                       */
/*                                                                           */
/*  Description   : This function performs filtering of a luma block         */
/*                  horizontal edge when boundary strength is less than 4.   */
/*                                                                           */
/*  Inputs        : pu1_src       - pointer to the src sample q0             */
/*                  src_strd      - source stride                            */
/*                  alpha         - alpha value for the boundary             */
/*                  beta          - beta value for the boundary              */
/*                  u4_bs         - packed Boundary strength array           */
/*                  pu1_cliptab   - tc0_table                                */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : This operation is described in Sec. 8.7.2.3 under the    */
/*                  title "Filtering process for edges for bS less than 4"   */
/*                  in ITU T Rec H.264.                                      */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         12 02 2015   Naveen Kumar P  Initial version                      */
/*         15 09 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_deblk_luma_horz_bslt4_avx2(UWORD8 *pu1_src,
                                       WORD32 src_strd,
                                       WORD32 alpha,
                                       WORD32 beta,
                                       UWORD32 u4_bs,
                                       const UWORD8 *pu1_cliptab)
{

    WORD16 i16_posP2, i16_posP1, i16_posP0, i16_posQ1, i16_posQ2;
    UWORD8 *pu1_HorzPixel;
    __m256i zero = _mm256_setzero_si256();
    __m128i zero_128 = _mm_setzero_si128();
    __m128i Alpha_8x16,bs_flag_16x8b, C0_16x8, C0_8x16, C0_hi_8x16;
    __m256i Beta_8x32,in_macro_32x8,in_macro_1,in_macro_2,flag1_32x8,flag2_32x8;
    __m256i C_8x32,C0_8x32_res,temp1,temp2,temp3,temp4,res1,res2,q0p1_32x8,p0q1_32x8;
    __m128i p0_16x8,q0_16x8,temp1_128,temp2_128,flag1_16x8_128;
    __m256i const_val4_8x32,p0q0_32x8,p1q1_32x8,p2q2_32x8,q0p0_32x8;
    UWORD8 u1_Bs0, u1_Bs1, u1_Bs2, u1_Bs3;
    UWORD8 clip0, clip1, clip2, clip3;

    pu1_HorzPixel = pu1_src - (src_strd << 2);

    i16_posQ1 = src_strd;
    i16_posQ2 = X2(src_strd);
    i16_posP0 = X3(src_strd);
    i16_posP1 = X2(src_strd);
    i16_posP2 = src_strd;

    p0q0_32x8 = _mm256_loadu2_m128i((__m128i *)(pu1_src), (__m128i *)(pu1_HorzPixel + i16_posP0)); //lower -p0 higher-q0
    p1q1_32x8 = _mm256_loadu2_m128i((__m128i *)(pu1_src + i16_posQ1), (__m128i *)(pu1_HorzPixel + i16_posP1)); //l= p1, h=q1
    p2q2_32x8 = _mm256_loadu2_m128i((__m128i *)(pu1_src + i16_posQ2), (__m128i *)(pu1_HorzPixel + i16_posP2));

    u1_Bs0 = (u4_bs >> 24) & 0xff;
    u1_Bs1 = (u4_bs >> 16) & 0xff;
    u1_Bs2 = (u4_bs >> 8) & 0xff;
    u1_Bs3 = (u4_bs >> 0) & 0xff;
    clip0 = pu1_cliptab[u1_Bs0];
    clip1 = pu1_cliptab[u1_Bs1];
    clip2 = pu1_cliptab[u1_Bs2];
    clip3 = pu1_cliptab[u1_Bs3];

    Alpha_8x16 = _mm_set1_epi16(alpha);
    Beta_8x32 = _mm256_set1_epi16(beta);

    bs_flag_16x8b = _mm_set_epi8(u1_Bs3, u1_Bs3, u1_Bs3, u1_Bs3,
                                   u1_Bs2, u1_Bs2, u1_Bs2, u1_Bs2,
				                   u1_Bs1, u1_Bs1, u1_Bs1, u1_Bs1,
                                   u1_Bs0, u1_Bs0, u1_Bs0, u1_Bs0);

    C0_16x8 = _mm_set_epi8(clip3, clip3, clip3, clip3, clip2, clip2, clip2,
                           clip2, clip1, clip1, clip1, clip1, clip0, clip0,
                           clip0, clip0);

    bs_flag_16x8b = _mm_cmpeq_epi8(bs_flag_16x8b, zero_128);
    bs_flag_16x8b = _mm_xor_si128(bs_flag_16x8b, _mm_set1_epi8(0xFF)); //Invert for required mask
    C0_8x16 = _mm_unpacklo_epi8(C0_16x8, zero_128);
    C0_hi_8x16 = _mm_unpackhi_epi8(C0_16x8, zero_128);
    C0_8x32_res = _mm256_set_m128i(C0_hi_8x16,C0_8x16);

    //Cond1 (ABS(p0 - q0) < alpha)
    p0_16x8 = _mm256_castsi256_si128(p0q0_32x8);
    q0p0_32x8 = _mm256_permute2x128_si256(p0q0_32x8, p0q0_32x8, 0x1);
    q0_16x8   = _mm256_castsi256_si128(p0q0_32x8);
    temp1_128 = _mm_subs_epu8(q0_16x8, p0_16x8);
    temp2_128 = _mm_subs_epu8(p0_16x8, q0_16x8);
    temp1_128 = _mm_add_epi8(temp1_128, temp2_128);

    temp2_128 = _mm_unpacklo_epi8(temp1_128, zero_128);
    temp1_128 = _mm_unpackhi_epi8(temp1_128, zero_128);

    temp2_128 = _mm_cmpgt_epi16(Alpha_8x16, temp2_128);
    temp1_128 = _mm_cmpgt_epi16(Alpha_8x16, temp1_128);
    flag1_16x8_128 = _mm_packs_epi16(temp2_128, temp1_128);
    flag1_16x8_128 = _mm_and_si128(flag1_16x8_128, bs_flag_16x8b);

    flag1_32x8  = _mm256_set_m128i(flag1_16x8_128,flag1_16x8_128);

    //Cond2 (ABS(q1 - q0) < beta) & Cond3 (ABS(p1 - p0) < beta)
    temp1 = _mm256_subs_epu8(p0q0_32x8, p1q1_32x8);
    temp2 = _mm256_subs_epu8(p1q1_32x8, p0q0_32x8);
    temp1 = _mm256_add_epi8(temp1, temp2);

    temp2 = _mm256_unpacklo_epi8(temp1, zero);
    temp1 = _mm256_unpackhi_epi8(temp1, zero);

    temp2 = _mm256_cmpgt_epi16(Beta_8x32, temp2);
    temp1 = _mm256_cmpgt_epi16(Beta_8x32, temp1);

    flag2_32x8 = _mm256_packs_epi16(temp2, temp1);

    //!((ABS(p0 - q0) < alpha) || (ABS(q1 - q0) < beta) || (ABS(p1 - p0) < beta))
    flag1_32x8 = _mm256_and_si256(flag1_32x8, flag2_32x8);

   //(ABS(p2 - p0) < beta) & (ABS(q2 - q0) < beta)
    temp1 = _mm256_subs_epu8(p0q0_32x8, p2q2_32x8);
    temp2 = _mm256_subs_epu8(p2q2_32x8, p0q0_32x8);
    temp1 = _mm256_add_epi8(temp1, temp2);

    temp2 = _mm256_unpacklo_epi8(temp1, zero);
    temp1 = _mm256_unpackhi_epi8(temp1, zero);
    temp2 = _mm256_cmpgt_epi16(Beta_8x32, temp2);
    temp1 = _mm256_cmpgt_epi16(Beta_8x32, temp1);

    flag2_32x8 = _mm256_packs_epi16(temp2, temp1);
    flag2_32x8 = _mm256_and_si256(flag1_32x8, flag2_32x8);

    temp2 = _mm256_subs_epi16(zero, temp2);
    temp1 = _mm256_subs_epi16(zero, temp1);

    temp3 = _mm256_permute2x128_si256(temp2,temp1,0x20); // low adding
    temp4 = _mm256_permute2x128_si256(temp2,temp1,0x31); //high adding
    temp2 = _mm256_add_epi16(temp3,temp4);
    C_8x32 = _mm256_add_epi16(C0_8x32_res, temp2);       //
    const_val4_8x32 = _mm256_set1_epi16(4);

    res1 = _mm256_permute4x64_epi64(q0p0_32x8, 0xD8);
    res2 = _mm256_permute4x64_epi64(p1q1_32x8, 0xD8);

    temp3 = _mm256_subs_epi16(_mm256_unpacklo_epi8(res1, zero),
                              _mm256_unpackhi_epi8(res1, zero));
    temp4 = _mm256_subs_epi16(_mm256_unpacklo_epi8(res2, zero),
                              _mm256_unpackhi_epi8(res2, zero));

    temp1 = _mm256_slli_epi16(temp3, 2);
    temp1 = _mm256_add_epi16(temp1, temp4);
    temp1 = _mm256_add_epi16(temp1, const_val4_8x32);
    in_macro_32x8 = _mm256_srai_epi16(temp1, 3);

    in_macro_32x8 = _mm256_min_epi16(C_8x32, in_macro_32x8); //CLIP3
    C_8x32 = _mm256_subs_epi16(zero, C_8x32);
    in_macro_32x8 = _mm256_max_epi16(C_8x32, in_macro_32x8); //CLIP3

    temp3 = _mm256_unpacklo_epi8(res1, zero); //q0
    temp4 = _mm256_unpackhi_epi8(res1, zero); //p0

    temp1 = _mm256_add_epi16(temp4, in_macro_32x8);
    temp2 = _mm256_sub_epi16(temp3, in_macro_32x8);

    temp1 = _mm256_packus_epi16(temp2, temp1); // Suffle needed

    temp1 = _mm256_and_si256(temp1, flag1_32x8); //q0 p0

    temp2 = _mm256_and_si256(res1,
                          _mm256_xor_si256(flag1_32x8, _mm256_set1_epi16(0xFFFF)));

    temp1 = _mm256_add_epi8(temp1, temp2);
    temp1 = _mm256_permute4x64_epi64(temp1, 0xD8);
    _mm256_storeu2_m128i((__m128i *)(pu1_HorzPixel + i16_posP0),(__m128i *)(pu1_src),temp1);

    //if(Ap < Beta) if(Aq < Beta)
     temp1 = _mm256_avg_epu16(_mm256_unpacklo_epi8(res1, zero),
                              _mm256_unpackhi_epi8(res1, zero));

     temp2 = _mm256_slli_epi16(_mm256_unpacklo_epi8(p1q1_32x8, zero), 1);
     temp3 = _mm256_subs_epi16(_mm256_unpacklo_epi8(p2q2_32x8, zero), temp2);

     temp2 = _mm256_slli_epi16(_mm256_unpackhi_epi8(p1q1_32x8, zero), 1);
     temp2 = _mm256_subs_epi16(_mm256_unpackhi_epi8(p2q2_32x8, zero), temp2);

     temp4 = _mm256_permute2x128_si256(temp3, temp2, 0x20);     //p0 q0
     temp3 = _mm256_permute2x128_si256(temp3, temp2, 0x31);
     temp4 = _mm256_add_epi16(temp1, temp4); //p
     in_macro_1  = _mm256_srai_epi16(temp4, 1);
     temp3 = _mm256_add_epi16(temp1, temp3); //q
     in_macro_2  = _mm256_srai_epi16(temp3, 1);

     in_macro_1 = _mm256_min_epi16(C0_8x32_res, in_macro_1); //CLIP3
     C0_8x32_res = _mm256_subs_epi16(zero, C0_8x32_res);
     in_macro_1 = _mm256_max_epi16(C0_8x32_res, in_macro_1); //CLIP3

     in_macro_2 = _mm256_max_epi16(C0_8x32_res, in_macro_2); //CLIP3
     C0_8x32_res = _mm256_subs_epi16(zero, C0_8x32_res);
     in_macro_2 = _mm256_min_epi16(C0_8x32_res, in_macro_2); //CLIP3

     temp1 = _mm256_unpacklo_epi8(res2, zero);
     temp2 = _mm256_unpackhi_epi8(res2, zero);

     temp1 = _mm256_add_epi16(temp1, in_macro_1);
     temp2 = _mm256_add_epi16(temp2, in_macro_2);
     temp1 = _mm256_packus_epi16(temp1, temp2); // pl ph ql qh
     temp1 = _mm256_and_si256(temp1, flag2_32x8);
     temp2 = _mm256_and_si256(res2,_mm256_xor_si256(flag2_32x8, _mm256_set1_epi16(0xFFFF)));
     temp1 = _mm256_add_epi8(temp1, temp2);
     temp1 = _mm256_permute4x64_epi64(temp1, 0xD8);
     _mm256_storeu2_m128i((__m128i *)(pu1_src + i16_posQ1),(__m128i *)(pu1_HorzPixel + i16_posP1),temp1);

}
