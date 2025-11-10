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

#include <stdint.h>
#include <string.h>
#include <immintrin.h>



/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_deblk_chroma_vert_bslt4_avx2()                    */
/*                                                                           */
/*  Description   : This function performs filtering of a chroma block       */
/*                  vertical edge when the boundary strength is less than 4  */
/*                  in high profile.                                         */
/*                                                                           */
/*  Inputs        : pu1_src          - pointer to the src sample q0 of U     */
/*                  src_strd         - source stride                         */
/*                  alpha_cb         - alpha value for the boundary in U     */
/*                  beta_cb          - beta value for the boundary in U      */
/*                  alpha_cr         - alpha value for the boundary in V     */
/*                  beta_cr          - beta value for the boundary in V      */
/*                  u4_bs            - packed Boundary strength array        */
/*                  pu1_cliptab_cb   - tc0_table for U                       */
/*                  pu1_cliptab_cr   - tc0_table for V                       */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : This operation is described in Sec. 8.7.2.3 under the    */
/*                  title "Filtering process for edges for bS less than 4"   */
/*                  in ITU T Rec H.264 with alpha and beta values different  */
/*                  in U and V.                                              */
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

void ih264_deblk_chroma_vert_bslt4_avx2(UWORD8 *pu1_src,
                                         WORD32 src_strd,
                                         WORD32 alpha_cb,
                                         WORD32 beta_cb,
                                         WORD32 alpha_cr,
                                         WORD32 beta_cr,
                                         UWORD32 u4_bs,
                                         const UWORD8 *pu1_cliptab_cb,
                                         const UWORD8 *pu1_cliptab_cr)
{
    UWORD8 *pu1_src_uv = pu1_src; /* Pointer to the src sample q0 of plane U*/
    UWORD8 u1_Bs0, u1_Bs1, u1_Bs2, u1_Bs3;
    WORD32 alpha_cbcr = (alpha_cr << 16) + alpha_cb;
    WORD32 beta_cbcr = (beta_cr << 16) + beta_cb;
    __m128i linea, lineb, linec, lined, linee, linef, lineg, lineh;
    __m256i lineab, linecd, lineef, linegh, lineae, linebf, linecg, linedh;
    __m256i temp1, temp2, temp3, temp4;
    __m256i t1,t3, t2,t4,pq0_uv_32x8,pq1_uv_32x8,tmp1,tmp2,p0_uv_8x32,q0_uv_8x32;

    __m256i pq0_uv_8x32, pq1_uv_8x32, p1_uv_8x32,pq0_uv_8x32_1,pq0_uv_8x32_2;
    __m256i flag_bs, flag1, flag2;
    __m256i diff, diff1, alpha_cbcr_32x8, beta_cbcr_32x8, in_macro;
    __m256i zero = _mm256_setzero_si256();
    __m256i C0_uv_8x32;
    __m256i p0_uv_8x32_1, p0_uv_8x32_2, q0_uv_8x32_1, q0_uv_8x32_2,p0_uv_32x8_1,q0_uv_32x8_1;

    u1_Bs0 = (u4_bs >> 24) & 0xff;
    u1_Bs1 = (u4_bs >> 16) & 0xff;
    u1_Bs2 = (u4_bs >> 8) & 0xff;
    u1_Bs3 = (u4_bs >> 0) & 0xff;

    flag_bs = _mm256_set_epi8(u1_Bs3, u1_Bs3, u1_Bs3, u1_Bs3,u1_Bs3, u1_Bs3, u1_Bs3, u1_Bs3,
                              u1_Bs2, u1_Bs2, u1_Bs2, u1_Bs2, u1_Bs2, u1_Bs2,u1_Bs2, u1_Bs2,
                              u1_Bs1, u1_Bs1, u1_Bs1, u1_Bs1,u1_Bs1, u1_Bs1, u1_Bs1, u1_Bs1,
                              u1_Bs0, u1_Bs0, u1_Bs0, u1_Bs0,u1_Bs0, u1_Bs0, u1_Bs0, u1_Bs0);
    flag_bs = _mm256_cmpeq_epi8(flag_bs, zero); //Set flag to 1s and 0s
    flag_bs = _mm256_xor_si256(flag_bs, _mm256_set1_epi8(0xFF)); //Invert for required mask

    /* Load and transpose the pixel values */
    lineab =  _mm256_loadu2_m128i((__m128i *)(pu1_src_uv - 4 + src_strd), (__m128i *)(pu1_src_uv - 4));
    linecd =  _mm256_loadu2_m128i((__m128i *)(pu1_src_uv - 4 + 3 * src_strd), (__m128i *)(pu1_src_uv - 4 + 2 * src_strd));
    lineef =  _mm256_loadu2_m128i((__m128i *)(pu1_src_uv - 4 + 5 * src_strd), (__m128i *)(pu1_src_uv - 4 + 4 * src_strd));
    linegh =  _mm256_loadu2_m128i((__m128i *)(pu1_src_uv - 4 + 7 * src_strd), (__m128i *)(pu1_src_uv - 4 + 6 * src_strd));

    temp1 = _mm256_unpacklo_epi64(lineab, zero);  //a0 -- a7  000.. b0..b7 000
    temp2 = _mm256_unpacklo_epi64(linecd, zero);
    temp3 = _mm256_unpacklo_epi64(lineef, zero);  //e0 -- e7  000.. f0..f7 000
    temp4 = _mm256_unpacklo_epi64(linegh, zero);

    temp1 = _mm256_unpacklo_epi16(temp1, temp2);  //a0 a1 c0 c1 --  a6 a7 c6 c7 b0 b1 d0 d1.. b6 b7 d6 d7
    temp2 = _mm256_unpacklo_epi16(temp3, temp4);  //e0 e1 g0 g1                f0 f1 h0 h1

    t2 = _mm256_permute2f128_si256(temp1, temp2, 0x20);
    t3 = _mm256_permute2f128_si256(temp1, temp2, 0x31);

    tmp1 = _mm256_unpacklo_epi16(t2, t3);    //a0 a1 b0 b1 c0 c1 d0 d1 -a2 a3 b2 b3 ....  e0 e1 f0 f1 g0 g1 h0 h1  -e2 e3..
    tmp2 = _mm256_unpackhi_epi16(t2, t3);    //a4 a5 b4 b5             -a6 a7 b6 b7


    temp1 = _mm256_unpacklo_epi8(tmp1,zero); // a0 0 a1 0 b0 0 b1 0 c0 0 c1 0 d0 0 d1 0 -  e0 0 e1 0 ..   => p1
    temp2 = _mm256_unpackhi_epi8(tmp1,zero); // a2 0 a3 0                                                 => p0
    temp3 = _mm256_unpacklo_epi8(tmp2,zero); //a4 0 a5 0                                                  => q0
    temp4 = _mm256_unpackhi_epi8(tmp2,zero); //a6 0 a7 0                                                 => q1

    pq1_uv_32x8 = _mm256_packus_epi16(temp1,temp4);     // 0213
    pq0_uv_32x8 = _mm256_packus_epi16(temp2,temp3);     //0213

    diff = _mm256_subs_epi16(temp2, temp3); //Condn 1    (p0 -q0) - set (3), set(3)
    diff = _mm256_abs_epi16(diff);
    alpha_cbcr_32x8 = _mm256_set1_epi32(alpha_cbcr);
    flag1 = _mm256_cmpgt_epi16(alpha_cbcr_32x8, diff);

    diff = _mm256_subs_epi16(temp4, temp3); //Condtn 2   (q1 -q0)
    diff = _mm256_abs_epi16(diff);
    beta_cbcr_32x8 = _mm256_set1_epi32(beta_cbcr);
    flag1 = _mm256_and_si256(flag1, _mm256_cmpgt_epi16(beta_cbcr_32x8, diff));


    diff = _mm256_subs_epi16(temp1, temp2); //Condtn 3  (p1 -p0)
    diff = _mm256_abs_epi16(diff);
    flag1 = _mm256_and_si256(flag1, _mm256_cmpgt_epi16(beta_cbcr_32x8, diff));

    diff = _mm256_subs_epi16(temp3, temp2);     //(q0 -p0)
    diff = _mm256_slli_epi16(diff, 2);

    diff1 = _mm256_subs_epi16(temp1, temp4);     //(p1 -q1)
    diff = _mm256_add_epi16(diff, diff1);

    diff = _mm256_add_epi16(diff, _mm256_set1_epi16(4));
    in_macro = _mm256_srai_epi16(diff, 3);


    C0_uv_8x32 = _mm256_set_epi16(pu1_cliptab_cr[u1_Bs1], pu1_cliptab_cb[u1_Bs1],
                               pu1_cliptab_cr[u1_Bs1], pu1_cliptab_cb[u1_Bs1],
                               pu1_cliptab_cr[u1_Bs0], pu1_cliptab_cb[u1_Bs0],
                               pu1_cliptab_cr[u1_Bs0], pu1_cliptab_cb[u1_Bs0],
                               pu1_cliptab_cr[u1_Bs3], pu1_cliptab_cb[u1_Bs3],
                               pu1_cliptab_cr[u1_Bs3], pu1_cliptab_cb[u1_Bs3],
                               pu1_cliptab_cr[u1_Bs2], pu1_cliptab_cb[u1_Bs2],
                               pu1_cliptab_cr[u1_Bs2], pu1_cliptab_cb[u1_Bs2]);

    C0_uv_8x32 = _mm256_add_epi16(C0_uv_8x32, _mm256_set1_epi16(1));

    in_macro = _mm256_min_epi16(C0_uv_8x32, in_macro); //CLIP3
    C0_uv_8x32 = _mm256_subs_epi16(zero, C0_uv_8x32);
    in_macro = _mm256_max_epi16(C0_uv_8x32, in_macro);

    p0_uv_8x32_1 = _mm256_add_epi16(temp2, in_macro);
    q0_uv_8x32_1 = _mm256_sub_epi16(temp3, in_macro);


    flag1 = _mm256_and_si256(flag1, flag_bs);
    flag1 = _mm256_packs_epi16(flag1, flag1);  // 0213

    pq0_uv_8x32 = _mm256_packus_epi16(p0_uv_8x32_1,q0_uv_8x32_1); //0213

    pq0_uv_8x32_1 = _mm256_and_si256(pq0_uv_32x8,
                                 _mm256_xor_si256(flag1, _mm256_set1_epi8(0xFF)));
    pq0_uv_8x32_2 = _mm256_and_si256(pq0_uv_8x32, flag1);
    pq0_uv_32x8 = _mm256_add_epi8(pq0_uv_8x32_1, pq0_uv_8x32_2);


    t1 = _mm256_unpacklo_epi16(pq1_uv_32x8, pq0_uv_32x8);   // temp1 temp3
    t2 = _mm256_unpackhi_epi16(pq1_uv_32x8, pq0_uv_32x8);   // temp2 temp4

    t4 = _mm256_shufflelo_epi16(t2, _MM_SHUFFLE(2, 3, 0, 1));  // pshuflw
    t4 = _mm256_shufflehi_epi16(t4, _MM_SHUFFLE(2, 3, 0, 1));

    lineae = _mm256_unpacklo_epi32(t1, t4);   // temp1 temp3
    linecg = _mm256_unpackhi_epi32(t1, t4);   // temp2 temp4

    linea =  _mm256_castsi256_si128(lineae);
    lineb = _mm256_castsi256_si128(_mm256_srli_si256(lineae, 8));
    lineae =  _mm256_permute2f128_si256(lineae, lineae, 0x1);
    linee =  _mm256_castsi256_si128(lineae);
    linef = _mm256_castsi256_si128(_mm256_srli_si256(lineae, 8));


    linec =  _mm256_castsi256_si128(linecg);
    lined =  _mm256_castsi256_si128(_mm256_srli_si256(linecg, 8));
    linecg =  _mm256_permute2f128_si256(linecg, linecg, 0x1);
    lineg =   _mm256_castsi256_si128(linecg);
    lineh = _mm256_castsi256_si128(_mm256_srli_si256(linecg, 8));

    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4), linea);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + src_strd), lineb);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 2 * src_strd), linec);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 3 * src_strd), lined);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 4 * src_strd), linee);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 5 * src_strd), linef);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 6 * src_strd), lineg);
    _mm_storel_epi64((__m128i *)(pu1_src_uv - 4 + 7 * src_strd), lineh);

}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264_deblk_chroma_horz_bslt4_avx2()                    */
/*                                                                           */
/*  Description   : This function performs filtering of a chroma block       */
/*                  horizontal edge when the boundary strength is less than  */
/*                  4 in high profile.                                       */
/*                                                                           */
/*  Inputs        : pu1_src          - pointer to the src sample q0 of U     */
/*                  src_strd         - source stride                         */
/*                  alpha_cb         - alpha value for the boundary in U     */
/*                  beta_cb          - beta value for the boundary in U      */
/*                  alpha_cr         - alpha value for the boundary in V     */
/*                  beta_cr          - beta value for the boundary in V      */
/*                  u4_bs            - packed Boundary strength array        */
/*                  pu1_cliptab_cb   - tc0_table for U                       */
/*                  pu1_cliptab_cr   - tc0_table for V                       */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : This operation is described in Sec. 8.7.2.3 under the    */
/*                  title "Filtering process for edges for bS less than 4"   */
/*                  in ITU T Rec H.264 with alpha and beta values different  */
/*                  in U and V.                                              */
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
/*         12 10 2020   Priyanka Bose   AVX2 Intel Intrinsics Support        */
/*****************************************************************************/
void ih264_deblk_chroma_horz_bslt4_avx2 (UWORD8 *pu1_src,
                                         WORD32 src_strd,
                                         WORD32 alpha_cb,
                                         WORD32 beta_cb,
                                         WORD32 alpha_cr,
                                         WORD32 beta_cr,
                                         UWORD32 u4_bs,
                                         const UWORD8 *pu1_cliptab_cb,
                                         const UWORD8 *pu1_cliptab_cr)
{
    UWORD8 *pu1_src_uv = pu1_src; /* Pointer to the src sample q0 of plane U*/
    WORD16 i16_posP1, i16_posP0, i16_posQ1;
    UWORD8 u1_Bs0, u1_Bs1, u1_Bs2, u1_Bs3;

    UWORD8 *pu1_HorzPixelUV; /*! < Pointer to the first pixel of the boundary */
    WORD32 alpha_cbcr = (alpha_cr << 16) + alpha_cb;
    WORD32 beta_cbcr = (beta_cr << 16) + beta_cb;
    __m256i p0q0_uv_32x8,p1q1_uv_32x8;
    __m256i temp1,temp2,temp3,temp4;
    __m256i flag_bs, flag1, flag2;
    __m256i diff, diff1, alpha_cbcr_32x8, beta_cbcr_32x8, in_macro;
    __m256i zero = _mm256_setzero_si256();
    __m256i C0_uv_8x32;
    __m256i p0q0_uv_8x32_1, p0q0_uv_8x32_2,res1,res2,p0_uv_8x32_1,q0_uv_8x32_1;

    pu1_HorzPixelUV = pu1_src_uv - (src_strd << 1);

    i16_posQ1 = src_strd;
    i16_posP0 = src_strd;
    i16_posP1 = 0;

    u1_Bs0 = (u4_bs >> 24) & 0xff;
    u1_Bs1 = (u4_bs >> 16) & 0xff;
    u1_Bs2 = (u4_bs >> 8) & 0xff;
    u1_Bs3 = (u4_bs >> 0) & 0xff;

    flag_bs = _mm256_set_epi8(u1_Bs3, u1_Bs3, u1_Bs3, u1_Bs3,
                           u1_Bs2, u1_Bs2, u1_Bs2, u1_Bs2,
                           u1_Bs3, u1_Bs3, u1_Bs3, u1_Bs3,
                           u1_Bs2, u1_Bs2,u1_Bs2, u1_Bs2,
                           u1_Bs1, u1_Bs1, u1_Bs1, u1_Bs1,
                           u1_Bs0, u1_Bs0, u1_Bs0, u1_Bs0,
                           u1_Bs1, u1_Bs1, u1_Bs1, u1_Bs1,
                           u1_Bs0, u1_Bs0, u1_Bs0, u1_Bs0);
    flag_bs = _mm256_cmpeq_epi8(flag_bs, zero); //Set flag to 1s and 0s
    flag_bs = _mm256_xor_si256(flag_bs, _mm256_set1_epi8(0xFF)); //Invert for required mask

    p0q0_uv_32x8 = _mm256_loadu2_m128i((__m128i *)(pu1_src_uv), (__m128i *)(pu1_HorzPixelUV + i16_posP0));
    p1q1_uv_32x8 = _mm256_loadu2_m128i((__m128i *)(pu1_src_uv + i16_posQ1), (__m128i *)(pu1_HorzPixelUV + i16_posP1));

    res1 = _mm256_permute4x64_epi64(p0q0_uv_32x8,0xD8);
    res2 = _mm256_permute4x64_epi64(p1q1_uv_32x8,0xD8);

    temp3 = _mm256_unpacklo_epi8(res1, zero); //p0 l 0 h 0
    temp4 = _mm256_unpackhi_epi8(res1, zero); //q0
    temp1 = _mm256_unpacklo_epi8(res2, zero); //p1
    temp2 = _mm256_unpackhi_epi8(res2, zero); //q1

    diff = _mm256_subs_epi16(temp3, temp4); //Condn 1 //p0 l h - q0 l h
    diff = _mm256_abs_epi16(diff);
    alpha_cbcr_32x8 = _mm256_set1_epi32(alpha_cbcr);
    flag1 = _mm256_cmpgt_epi16(alpha_cbcr_32x8, diff);

    diff = _mm256_subs_epi16(temp2, temp4); //Condtn 2
    diff = _mm256_abs_epi16(diff);
    beta_cbcr_32x8 = _mm256_set1_epi32(beta_cbcr);
    flag1 = _mm256_and_si256(flag1, _mm256_cmpgt_epi16(beta_cbcr_32x8, diff));

    diff = _mm256_subs_epi16(temp1, temp3); //Condtn 3
    diff = _mm256_abs_epi16(diff);
    flag1 = _mm256_and_si256(flag1, _mm256_cmpgt_epi16(beta_cbcr_32x8, diff));

    diff = _mm256_subs_epi16(temp4, temp3);
    diff = _mm256_slli_epi16(diff, 2);
    diff1 = _mm256_subs_epi16(temp1, temp2);
    diff = _mm256_add_epi16(diff, diff1);
    diff = _mm256_add_epi16(diff, _mm256_set1_epi16(4));
    in_macro = _mm256_srai_epi16(diff, 3);

    C0_uv_8x32 = _mm256_set_epi16(
                               pu1_cliptab_cr[u1_Bs3], pu1_cliptab_cb[u1_Bs3],
                               pu1_cliptab_cr[u1_Bs3], pu1_cliptab_cb[u1_Bs3],
                               pu1_cliptab_cr[u1_Bs2], pu1_cliptab_cb[u1_Bs2],
                               pu1_cliptab_cr[u1_Bs2], pu1_cliptab_cb[u1_Bs2],
                               pu1_cliptab_cr[u1_Bs1], pu1_cliptab_cb[u1_Bs1],
                               pu1_cliptab_cr[u1_Bs1], pu1_cliptab_cb[u1_Bs1],
                               pu1_cliptab_cr[u1_Bs0], pu1_cliptab_cb[u1_Bs0],
                               pu1_cliptab_cr[u1_Bs0], pu1_cliptab_cb[u1_Bs0]);

    C0_uv_8x32 = _mm256_add_epi16(C0_uv_8x32, _mm256_set1_epi16(1));

    in_macro = _mm256_min_epi16(C0_uv_8x32, in_macro); //CLIP3
    C0_uv_8x32 = _mm256_subs_epi16(zero, C0_uv_8x32);
    in_macro = _mm256_max_epi16(C0_uv_8x32, in_macro);

    p0_uv_8x32_1 = _mm256_add_epi16(temp3, in_macro);
    q0_uv_8x32_1 = _mm256_sub_epi16(temp4, in_macro);

    p0q0_uv_8x32_2 = _mm256_packus_epi16(p0_uv_8x32_1,q0_uv_8x32_1);
    flag1 = _mm256_packs_epi16(flag1, flag1);
    flag1 = _mm256_and_si256(flag1, flag_bs); //Final flag (BS condition + other 3 conditions)

    p0q0_uv_8x32_1 = _mm256_and_si256(res1,
                                 _mm256_xor_si256(flag1, _mm256_set1_epi8(0xFF)));
    p0q0_uv_8x32_2 = _mm256_and_si256(p0q0_uv_8x32_2, flag1);
    p0q0_uv_8x32_1 = _mm256_add_epi8(p0q0_uv_8x32_1, p0q0_uv_8x32_2);
    p0q0_uv_8x32_1 = _mm256_permute4x64_epi64(p0q0_uv_8x32_1,0xD8);

    _mm256_storeu2_m128i((__m128i *)(pu1_src_uv),(__m128i *)(pu1_HorzPixelUV + i16_posP0), p0q0_uv_8x32_1);

}
