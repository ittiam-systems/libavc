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
/**
 *******************************************************************************
 * @file
 *  ih264_iquant_itrans_recon_avx2.c
 *
 * @brief
 *  Contains function definitions for inverse  quantization, inverse
 * transform and reconstruction
 *
 * @author
 *  Priyanka Bose
 *
 * @par List of Functions:
 *  - ih264_iquant_itrans_recon_4x4_avx2()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */


#include <stdio.h>

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "ih264_structs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include <immintrin.h>


/*
 ********************************************************************************
 *
 * @brief This function reconstructs a 4x4 sub block from quantized resiude and
 * prediction buffer
 *
 * @par Description:
 *  The quantized residue is first inverse quantized, then inverse transformed.
 *  This inverse transformed content is added to the prediction buffer to recon-
 *  struct the end output
 *
 * @param[in] pi2_src
 *  quantized 4x4 block
 *
 * @param[in] pu1_pred
 *  prediction 4x4 block
 *
 * @param[out] pu1_out
 *  reconstructed 4x4 block
 *
 * @param[in] src_strd
 *  quantization buffer stride
 *
 * @param[in] pred_strd,
 *  Prediction buffer stride
 *
 * @param[in] out_strd
 *  recon buffer Stride
 *
 * @param[in] pu2_scaling_list
 *  pointer to scaling list
 *
 * @param[in] pu2_norm_adjust
 *  pointer to inverse scale matrix
 *
 * @param[in] u4_qp_div_6
 *  Floor (qp/6)
 *
 * @param[in] pi4_tmp
 * temporary buffer of size 1*16
 *
 * @returns none
 *
 * @remarks none
 *
 *******************************************************************************
 */
void ih264_iquant_itrans_recon_4x4_avx2(WORD16 *pi2_src,
                                   UWORD8 *pu1_pred,
                                   UWORD8 *pu1_out,
                                   WORD32 pred_strd,
                                   WORD32 out_strd,
                                   const UWORD16 *pu2_iscal_mat,
                                   const UWORD16 *pu2_weigh_mat,
                                   UWORD32 u4_qp_div_6,
                                   WORD16 *pi2_tmp,
                                   WORD32 iq_start_idx,
                                   WORD16 *pi2_dc_ld_addr)
 {

    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    __m256i src_r0_r1, src_r2_r3;
    __m256i src_r0, src_r1, src_r2, src_r3;
    __m256i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m256i sign_reg, dequant_r0_r1, dequant_r2_r3;
    __m256i zero_8x32b = _mm256_setzero_si256();          // all bits reset to zero
    __m256i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m256i resq_r0, resq_r1, resq_r2, resq_r3;
    __m256i add_rshift = _mm256_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m256i value_32 = _mm256_set1_epi32(32);
    __m128i value_32_128 = _mm_set1_epi32(32);
    __m128i r0,r1,r2,r3,t0,t1,t2,t3,t4,t5,t6,t7,sign_reg_128, de_r0_r1,de_r2_r3,r0_r1,r2_r3,scale_r0_r1,scale_r2_r3;
    __m128i zero_8x16b_128 = _mm_setzero_si128();
    UNUSED (pi2_tmp);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/
    src_r0_r1 = _mm256_loadu_si256((__m256i *) (pi2_src)); //a00 a01 a02 a03 a10 a11 a12 a13 -- the source matrix 0th,1st row
    scalemat_r0_r1 = _mm256_loadu_si256((__m256i *) (pu2_iscal_mat)); //b00 b01 b02 b03 b10 b11 b12 b13 -- the scaling matrix 0th,1st row
    dequant_r0_r1 = _mm256_loadu_si256((__m256i *) (pu2_weigh_mat)); //q00 q01 q02 q03 q10 q11 q12 q13 -- all 16 bits

    temp0 = _mm256_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);
    temp4 = _mm256_unpacklo_epi16(temp0, zero_8x32b); // b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long
    temp5 = _mm256_unpackhi_epi16(temp0, zero_8x32b); // b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long

    src_r0 = _mm256_unpacklo_epi16(src_r0_r1, zero_8x32b); // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r1 = _mm256_unpackhi_epi16(src_r0_r1, zero_8x32b); // a10 0 a11 0 a12 0 a13 0 -- 16 bit long

    temp0 = _mm256_madd_epi16(src_r0, temp4); //a00*b00*q00 a10*b10*q10 a20*b20*q20 a30*b30 q30 -- 32 bits long
    temp1 = _mm256_madd_epi16(src_r1, temp5);

    if (u4_qp_div_6 >= 4) {
        resq_r0 = _mm256_slli_epi32(temp0, u4_qp_div_6 - 4);
        resq_r1 = _mm256_slli_epi32(temp1, u4_qp_div_6 - 4);
    } else {
        temp4 = _mm256_add_epi32(temp0, add_rshift);
        temp5 = _mm256_add_epi32(temp1, add_rshift);
        resq_r0 = _mm256_srai_epi32(temp0, 4 - u4_qp_div_6);
        resq_r1 = _mm256_srai_epi32(temp1, 4 - u4_qp_div_6);
    }

    if (iq_start_idx == 1)
        resq_r0 = _mm256_insert_epi32(resq_r0,(WORD32)pi2_dc_ld_addr[0],0);
    /* Perform Inverse transform */

    /* Perform Inverse transform */
    /*-------------------------------------------------------------*/
    /* IDCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */

    temp0 = _mm256_unpacklo_epi32(resq_r0, resq_r1);                  //a0 c0 a1 c1  b0 d0 b1 d1
    temp1 = _mm256_unpackhi_epi32(resq_r0, resq_r1);                  //a2 c2

    resq_r0 = _mm256_permute2f128_si256(temp0, temp1, 0x20);
    resq_r1 = _mm256_permute2f128_si256(temp0, temp1, 0x31);

    temp0 = _mm256_unpacklo_epi64(resq_r0, resq_r1);                  // w0 w2
    temp1 = _mm256_unpackhi_epi64(resq_r0, resq_r1);                  //  w1 w3

    resq_r0 = _mm256_permute2f128_si256(temp0, temp1, 0x20);
    resq_r1 = _mm256_permute2f128_si256(temp0, temp1, 0x31);

    r0 = _mm256_extracti128_si256(resq_r0, 0x0);
    r1 = _mm256_extracti128_si256(resq_r0, 0x1);
    r2 = _mm256_extracti128_si256(resq_r1, 0x0);
    r3 = _mm256_extracti128_si256(resq_r1, 0x1);

    //Transform starts -- horizontal transform
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    t0 = _mm_add_epi32(r0, r2);
    /* z1 = w0 - w2                                             */
    t1 = _mm_sub_epi32(r0, r2);
    /* z2 = (w1 >> 1) - w3                                      */
    t2 = _mm_srai_epi32(r1, 1);                         //(w1>>1)
    t2 = _mm_sub_epi32(t2, r3);                      //(w1>>1) - w3
    /* z3 = w1 + (w3 >> 1)                                      */
    t3 = _mm_srai_epi32(r3, 1);                         //(w3>>1) + w1
    t3 = _mm_add_epi32(t3, r1);
    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    r0 = _mm_add_epi32(t0, t3);
    /* x1 = z1 + z2                                             */
    r1 = _mm_add_epi32(t1, t2);
    /* x2 = z1 - z2                                             */
    r2 = _mm_sub_epi32(t1, t2);
    /* x3 = z0 - z3                                             */
    r3 = _mm_sub_epi32(t0, t3);

    t1 = _mm_unpacklo_epi32(r0, r1);                  //a0 a1 b0 b1
    t3 = _mm_unpacklo_epi32(r2, r3);                  //a2 a3 b2 b3
    t2 = _mm_unpackhi_epi32(r0, r1);                  //c0 c1 d0 d1
    t4 = _mm_unpackhi_epi32(r2, r3);                  //c2 c3 d2 d3
    r0 = _mm_unpacklo_epi64(t1, t3);                    //a0 a1 a2 a3
    r1 = _mm_unpackhi_epi64(t1, t3);                    //b0 b1 b2 b3
    r2 = _mm_unpacklo_epi64(t2, t4);                    //c0 c1 c2 c3
    r3 = _mm_unpackhi_epi64(t2, t4);                    //d0 d1 d2 d3

    //Transform ends -- horizontal transform

    //Load pred buffer
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0])); //p00 p01 p02 p03 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[pred_strd])); //p10 p11 p12 p13 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * pred_strd])); //p20 p21 p22 p23 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * pred_strd])); //p30 p31 p32 p33 0 0 0 0 0 0 0 0 -- all 8 bits

    pred_r0 = _mm_cvtepu8_epi32(pred_r0); //p00 p01 p02 p03 -- all 32 bits
    pred_r1 = _mm_cvtepu8_epi32(pred_r1); //p10 p11 p12 p13 -- all 32 bits ///Need to look
    pred_r2 = _mm_cvtepu8_epi32(pred_r2); //p20 p21 p22 p23 -- all 32 bits
    pred_r3 = _mm_cvtepu8_epi32(pred_r3); //p30 p31 p32 p33 -- all 32 bits

    /*--------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6      */
    /*                                                              */
    /* Add the prediction and store it back to same buffer          */
    /*--------------------------------------------------------------*/


    t0 = _mm_add_epi32(r0, r2);
    /* z1j = y0j - y2j                                                        */
    t1 = _mm_sub_epi32(r0, r2);
    /* z2j = (y1j>>1) - y3j                                                        */
    t2 = _mm_srai_epi32(r1, 1);                             //(y1j>>1)
    t2 = _mm_sub_epi32(t2, r3);
    /* z3j = y1j + (y3j>>1)                                                        */
    t3 = _mm_srai_epi32(r3, 1);                             //(y3j>>1)
    t3 = _mm_add_epi32(r1, t3);


    t4 = _mm_add_epi32(t0, t3);
    t4 = _mm_add_epi32(t4, value_32_128);
    t4 = _mm_srai_epi32(t4, 6);
    t4 = _mm_add_epi32(t4, pred_r0);

    t5 = _mm_add_epi32(t1, t2);
    t5 = _mm_add_epi32(t5, value_32_128);
    t5 = _mm_srai_epi32(t5, 6);
    t5 = _mm_add_epi32(t5, pred_r1);

    t6 = _mm_sub_epi32(t1, t2);
    t6 = _mm_add_epi32(t6, value_32_128);
    t6 = _mm_srai_epi32(t6, 6);
    t6 = _mm_add_epi32(t6, pred_r2);

    t7 = _mm_sub_epi32(t0, t3);
    t7 = _mm_add_epi32(t7, value_32_128);
    t7 = _mm_srai_epi32(t7, 6);
    t7 = _mm_add_epi32(t7, pred_r3);


    // 32-bit to 16-bit conversion
    t0 = _mm_packs_epi32(t4, t5);
    t1 = _mm_packs_epi32(t6, t7);


    //Clipping the results to 8 bits
    sign_reg_128 = _mm_cmpgt_epi16(t0, zero_8x16b_128);      // sign check
    t0 = _mm_and_si128(t0, sign_reg_128);
    sign_reg_128 = _mm_cmpgt_epi16(t1, zero_8x16b_128);
    t1 = _mm_and_si128(t1, sign_reg_128);

    r0 = _mm_packus_epi16(t0, t1);
    r1 = _mm_srli_si128(r0, 4);
    r2 = _mm_srli_si128(r1, 4);
    r3 = _mm_srli_si128(r2, 4);

    *pu4_out = _mm_cvtsi128_si32(r0);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(r1);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(r2);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(r3);
}
