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
 *  isvcd_iquant_itrans_residual_recon_sse42.c
 *
 * @brief
 *  Contains function definitions for iquant_itrans_residual_recon
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_iquant_itrans_residual_recon_4x4_sse42()
 *  - isvcd_iquant_itrans_residual_recon_8x8_sse42()
 *  - isvcd_iquant_itrans_residual_recon_4x4_dc_sse42()
 *  - isvcd_iquant_itrans_residual_recon_8x8_dc_sse42()
 *  - isvcd_iquant_itrans_residual_recon_chroma_4x4_sse42()
 *  - isvcd_iquant_itrans_residual_recon_chroma_4x4_dc_sse42()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
/* User include files */
#include <immintrin.h>
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "ih264_structs.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_4x4_sse42              */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_iquant_itrans_residual_recon_4x4_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 u4_qp_div_6, WORD16 *pi2_tmp, WORD32 iq_start_idx, WORD16 *pi2_dc_ld_addr)
{
    WORD32 i4_nnz = 0;
    WORD32 row_0, row_1, row_2, row_3;
    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i rsd_r0, rsd_r1, rsd_r2, rsd_r3;
    __m128i sign_reg, dequant_r0_r1, dequant_r2_r3;
    __m128i zero_8x16b = _mm_setzero_si128();  // all bits reset to zero
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);
    __m128i dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    __m128i dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    UNUSED(pi2_tmp);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/
    // a00 a01 a02 a03 a10 a11 a12 a13 -- the source matrix 0th,1st row
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));
    // a20 a21 a22 a23 a30 a31 a32 a33 --the source matrix 2nd,3rd row
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));
    // b00 b01 b02 b03 b10 b11 b12 b13 -- the scaling matrix 0th,1st row
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));
    // b20 b21 b22 b23 b30 b31 b32 b33 -- the scaling matrix 2nd,3rd row
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));
    // q00 q01 q02 q03 q10 q11 q12 q13 -- all 16 bits
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));
    // q20 q21 q22 q23 q30 q31 q32 q33 -- all 16 bits
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));
    // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11 b12*q12 b13*q13 -- 16 bit result
    temp0 = _mm_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);
    // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11 b12*q12 b13*q13 -- 16 bit result
    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);
    // b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);
    // b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);
    // b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);
    // b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);  // a10 0 a11 0 a12 0 a13 0 -- 16 bit long
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);  // a20 0 a21 0 a22 0 a23 0 -- 16 bit long
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);  // a30 0 a31 0 a32 0 a33 0 -- 16 bit long
    // a00*b00*q00 a10*b10*q10 a20*b20*q20 a30*b30 q30 -- 32 bits long
    temp4 = _mm_madd_epi16(src_r0, temp4);
    temp5 = _mm_madd_epi16(src_r1, temp5);
    temp6 = _mm_madd_epi16(src_r2, temp6);
    temp7 = _mm_madd_epi16(src_r3, temp7);

    if(u4_qp_div_6 >= 4)
    {
        resq_r0 = _mm_slli_epi32(temp4, u4_qp_div_6 - 4);
        resq_r1 = _mm_slli_epi32(temp5, u4_qp_div_6 - 4);
        resq_r2 = _mm_slli_epi32(temp6, u4_qp_div_6 - 4);
        resq_r3 = _mm_slli_epi32(temp7, u4_qp_div_6 - 4);
    }
    else
    {
        temp4 = _mm_add_epi32(temp4, add_rshift);
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp6 = _mm_add_epi32(temp6, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r0 = _mm_srai_epi32(temp4, 4 - u4_qp_div_6);
        resq_r1 = _mm_srai_epi32(temp5, 4 - u4_qp_div_6);
        resq_r2 = _mm_srai_epi32(temp6, 4 - u4_qp_div_6);
        resq_r3 = _mm_srai_epi32(temp7, 4 - u4_qp_div_6);
    }

    if(iq_start_idx == 1) resq_r0 = _mm_insert_epi32(resq_r0, (WORD32) pi2_dc_ld_addr[0], 0);
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
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);  // a0 b0 a1 b1
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);  // c0 d0 c1 d1
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);  // a2 b2 a3 b3
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);  // c2 d2 c3 d3
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);    // a0 b0 c0 d0
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);    // a1 b1 c1 d1
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);    // a2 b2 c2 d2
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);    // a3 b3 c3 d3
    // Transform starts -- horizontal transform
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);     //(w1>>1)
    temp2 = _mm_sub_epi32(temp2, resq_r3);  //(w1>>1) - w3
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);  //(w3>>1) + w1
    temp3 = _mm_add_epi32(temp3, resq_r1);
    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    resq_r0 = _mm_add_epi32(temp0, temp3);
    /* x1 = z1 + z2                                             */
    resq_r1 = _mm_add_epi32(temp1, temp2);
    /* x2 = z1 - z2                                             */
    resq_r2 = _mm_sub_epi32(temp1, temp2);
    /* x3 = z0 - z3                                             */
    resq_r3 = _mm_sub_epi32(temp0, temp3);
    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);  // a0 a1 b0 b1
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);  // a2 a3 b2 b3
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);  // c0 c1 d0 d1
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);  // c2 c3 d2 d3
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);    // a0 a1 a2 a3
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);    // b0 b1 b2 b3
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);    // c0 c1 c2 c3
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);    // d0 d1 d2 d3
    // Transform ends -- horizontal transform

    // Load pred buffer
    // p00 p01 p02 p03 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    // p10 p11 p12 p13 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[pred_strd]));
    // p20 p21 p22 p23 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * pred_strd]));
    // p30 p31 p32 p33 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * pred_strd]));

    pred_r0 = _mm_cvtepu8_epi32(pred_r0);  // p00 p01 p02 p03 -- all 32 bits
    pred_r1 = _mm_cvtepu8_epi32(pred_r1);  // p10 p11 p12 p13 -- all 32 bits
    pred_r2 = _mm_cvtepu8_epi32(pred_r2);  // p20 p21 p22 p23 -- all 32 bits
    pred_r3 = _mm_cvtepu8_epi32(pred_r3);  // p30 p31 p32 p33 -- all 32 bits

    // Load resd buffer
    rsd_r0 = _mm_loadl_epi64((__m128i *) (&pi2_rsd[0]));
    rsd_r1 = _mm_loadl_epi64((__m128i *) (&pi2_rsd[rsd_strd]));
    rsd_r2 = _mm_loadl_epi64((__m128i *) (&pi2_rsd[2 * rsd_strd]));
    rsd_r3 = _mm_loadl_epi64((__m128i *) (&pi2_rsd[3 * rsd_strd]));

    rsd_r0 = _mm_cvtepi16_epi32(rsd_r0);
    rsd_r1 = _mm_cvtepi16_epi32(rsd_r1);
    rsd_r2 = _mm_cvtepi16_epi32(rsd_r2);
    rsd_r3 = _mm_cvtepi16_epi32(rsd_r3);

    /*--------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6      */
    /*                                                              */
    /* Add the prediction and store it back to same buffer          */
    /*--------------------------------------------------------------*/
    /* z0j = y0j + y2j                                                        */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1j = y0j - y2j                                                        */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2j = (y1j>>1) - y3j */
    temp2 = _mm_srai_epi32(resq_r1, 1);  //(y1j>>1)
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);  //(y3j>>1)
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);

    row_0 = _mm_test_all_ones(_mm_cmpeq_epi32(temp4, zero_8x16b));  // return 1 if all zeros, else 0
    temp4 = _mm_add_epi32(temp4, rsd_r0);
    temp4 = _mm_min_epi16(dupmax_8x16b, temp4);
    temp4 = _mm_max_epi16(dupmin_8x16b, temp4);
    temp4 = _mm_add_epi32(temp4, pred_r0);
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);

    row_1 = _mm_test_all_ones(_mm_cmpeq_epi32(temp5, zero_8x16b));  // return 1 if all zeros, else 0
    temp5 = _mm_add_epi32(temp5, rsd_r1);
    temp5 = _mm_min_epi16(dupmax_8x16b, temp5);
    temp5 = _mm_max_epi16(dupmin_8x16b, temp5);
    temp5 = _mm_add_epi32(temp5, pred_r1);
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);

    row_2 = _mm_test_all_ones(_mm_cmpeq_epi32(temp6, zero_8x16b));  // return 1 if all zeros, else 0
    temp6 = _mm_add_epi32(temp6, rsd_r2);
    temp6 = _mm_min_epi16(dupmax_8x16b, temp6);
    temp6 = _mm_max_epi16(dupmin_8x16b, temp6);
    temp6 = _mm_add_epi32(temp6, pred_r2);
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);

    row_3 = _mm_test_all_ones(_mm_cmpeq_epi32(temp7, zero_8x16b));  // return 1 if all zeros, else 0
    temp7 = _mm_add_epi32(temp7, rsd_r3);
    temp7 = _mm_min_epi16(dupmax_8x16b, temp7);
    temp7 = _mm_max_epi16(dupmin_8x16b, temp7);
    temp7 = _mm_add_epi32(temp7, pred_r3);

    // 32-bit to 16-bit conversion
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);
    /*------------------------------------------------------------------*/
    // Clipping the results to 8 bits
    sign_reg = _mm_cmpgt_epi16(temp0, zero_8x16b);  // sign check
    temp0 = _mm_and_si128(temp0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, sign_reg);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    *pu4_out = _mm_cvtsi128_si32(resq_r0);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r1);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r2);
    pu1_out += out_strd;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r3);

    i4_nnz = !(row_0 && row_1 && row_2 && row_3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_8x8_sse42              */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_iquant_itrans_residual_recon_8x8_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscale_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 qp_div, WORD16 *pi2_tmp, WORD32 iq_start_idx, WORD16 *pi2_dc_ld_addr)
{
    __m128i rsd_r01_b0, rsd_r23_b0, rsd_r45_b2, rsd_r67_b2;
    __m128i rsd_r01_b1, rsd_r23_b1, rsd_r45_b3, rsd_r67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;
    __m128i src_r0;
    __m128i scalemat_r0;
    __m128i zero_8x16b = _mm_setzero_si128();  // all bits reset to zero
    __m128i value_32 = _mm_set1_epi32(32);
    __m128i add_rshift = _mm_set1_epi32((qp_div < 6) ? (1 << (5 - qp_div)) : 0);
    __m128i dequant_r0;
    __m128i predload_r;
    __m128i pred_r0_1, pred_r1_1, pred_r2_1, pred_r3_1, pred_r4_1, pred_r5_1, pred_r6_1, pred_r7_1;

    __m128i rsd_r0, rsd_r1, rsd_r2, rsd_r3, rsd_r4, rsd_r5, rsd_r6, rsd_r7;
    __m128i sign_reg;
    __m128i src_r0_1, src_r0_2;
    __m128i scalemat_r0_1, scalemat_r0_2;
    __m128i temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
    __m128i temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20;
    // To store dequantization results
    __m128i resq_r0_1, resq_r0_2, resq_r1_1, resq_r1_2, resq_r2_1, resq_r2_2, resq_r3_1, resq_r3_2,
        resq_r4_1, resq_r4_2, resq_r5_1, resq_r5_2, resq_r6_1, resq_r6_2, resq_r7_1, resq_r7_2;
    __m128i dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    __m128i dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    UNUSED(pi2_tmp);
    UNUSED(iq_start_idx);
    UNUSED(pi2_dc_ld_addr);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform. Note : DC coeff is not scaled     */
    /*************************************************************/

    // Row 0 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 -- the source matrix 0th row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src));
    // b00 b01 b02 b03 b04 b05 b06 b07 -- the scaling matrix 0th row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[0]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);

    if(qp_div >= 6)
    {
        resq_r0_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r0_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r0_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r0_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 16
    // bit long
    resq_r0_1 = _mm_packs_epi32(resq_r0_1, resq_r0_2);
    // Row 1 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 1st row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 8));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 1st row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 8));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[8]));
    // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);
    // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r1_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r1_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r1_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r1_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 16
    // bit long
    resq_r1_1 = _mm_packs_epi32(resq_r1_1, resq_r1_2);
    // Row 2 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 2nd row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 16));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 2nd row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 16));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[16]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r2_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r2_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r2_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r2_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 -- 16 bit long
    resq_r2_1 = _mm_packs_epi32(resq_r2_1, resq_r2_2);
    // Row 3 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 3rd row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 24));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 3rd row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 24));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[24]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 - 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r3_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r3_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r3_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r3_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 .... -- 16 bit long
    resq_r3_1 = _mm_packs_epi32(resq_r3_1, resq_r3_2);
    // Row 4 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 4th row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 32));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 4th row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 32));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[32]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r4_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r4_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r4_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r4_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5  -- 16 bit long
    resq_r4_1 = _mm_packs_epi32(resq_r4_1, resq_r4_2);
    // Row 5 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 5th row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 40));
    //
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 5th row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 40));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[40]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r5_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r5_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r5_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r5_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 a06*b06*q6 -- 16 bit long
    resq_r5_1 = _mm_packs_epi32(resq_r5_1, resq_r5_2);
    // Row 6 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 6th row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 48));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 6th row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 48));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[48]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r6_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r6_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r6_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r6_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 -- 16 bit long
    resq_r6_1 = _mm_packs_epi32(resq_r6_1, resq_r6_2);
    // Row 7 processing
    // a00 a01 a02 a03 a04 a05 a06 a07 a08 --the source matrix 7th row
    src_r0 = _mm_loadu_si128((__m128i *) (pi2_src + 56));
    // b00 b01 b02 b03 b04 b05 b06 b07 b08 -- the scaling matrix 7th row
    scalemat_r0 = _mm_loadu_si128((__m128i *) (pu2_iscale_mat + 56));
    // q0 q1 q2 q3 q4 q5 q6 q7 -- all 16 bits
    dequant_r0 = _mm_loadu_si128((__m128i *) (&pu2_weigh_mat[56]));
    src_r0_1 = _mm_unpacklo_epi16(src_r0, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r0_2 = _mm_unpackhi_epi16(src_r0, zero_8x16b);  // a04 0 a05 0 a06 0 a07 0 -- 16 bit long
    // b00*q0 b01*q1 b02*q2 b03*q3 b04*q4 b05*q5 b06*q6 b07*q7 -- 16 bit result
    temp10 = _mm_mullo_epi16(scalemat_r0, dequant_r0);
    // b00*q0 0 b01*q1 0 b02*q2 0 b03*q3 0 -- 16 bit long
    scalemat_r0_1 = _mm_unpacklo_epi16(temp10, zero_8x16b);
    // b04*q4 0 b05*q5 0 b06*q6 0 b07*q7 0 -- 16 bit long
    scalemat_r0_2 = _mm_unpackhi_epi16(temp10, zero_8x16b);
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 -- 32 bits long
    temp5 = _mm_madd_epi16(src_r0_1, scalemat_r0_1);
    // a04*b04*q4 a05*b05*q5 a06*b06*q6 a07*b07*q7 -- 32 bits long
    temp7 = _mm_madd_epi16(src_r0_2, scalemat_r0_2);
    if(qp_div >= 6)
    {
        resq_r7_1 = _mm_slli_epi32(temp5, qp_div - 6);
        resq_r7_2 = _mm_slli_epi32(temp7, qp_div - 6);
    }
    else
    {
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r7_1 = _mm_srai_epi32(temp5, 6 - qp_div);
        resq_r7_2 = _mm_srai_epi32(temp7, 6 - qp_div);
    }
    // a00*b00*q0 a01*b01*q1 a02*b02*q2 a03*b03*q3 a04*b04*q4 a05*b05*q5 -- 16 bit long
    resq_r7_1 = _mm_packs_epi32(resq_r7_1, resq_r7_2);

    /* Perform Inverse transform */
    /*--------------------------------------------------------------------*/
    /* IDCT [ Horizontal transformation ]                                 */
    /*--------------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3 a4 a5 a6 a7
     *  b0 b1 b2 b3 b4 b5 b6 b7
     *  c0 c1 c2 c3 c4 c5 c6 c7
     *  d0 d1 d2 d3 d4 d5 d6 d7
     */
    temp1 = _mm_unpacklo_epi16(resq_r0_1, resq_r1_1);  // a0 b0 a1 b1 a2 b2 a3 b3
    temp3 = _mm_unpacklo_epi16(resq_r2_1, resq_r3_1);  // c0 d0 c1 d1 c2 d2 c3 d3
    temp2 = _mm_unpackhi_epi16(resq_r0_1, resq_r1_1);  // a4 b4 a5 b5 a6 b6 a7 b7
    temp4 = _mm_unpackhi_epi16(resq_r2_1, resq_r3_1);  // c4 d4 c5 d5 c6 d6 c7 d7
    resq_r0_1 = _mm_unpacklo_epi32(temp1, temp3);      // a0 b0 c0 d0 a1 b1 c1 d1
    resq_r1_1 = _mm_unpackhi_epi32(temp1, temp3);      // a2 b2 c2 d2 a3 b3 c3 d3
    resq_r2_1 = _mm_unpacklo_epi32(temp2, temp4);      // a4 b4 c4 d4 a5 b5 c5 d5
    resq_r3_1 = _mm_unpackhi_epi32(temp2, temp4);      // a6 b6 c6 d6 a7 b7 c7 d7
    /*
     * e0 e1 e2 e3 e4 e5 e6 e7
     * f0 f1 f2 f3 f4 f5 f6 f7
     * g0 g1 g2 g3 g4 g5 g6 g7
     * h0 h1 h2 h3 h4 h5 h6 h7
     */
    temp1 = _mm_unpacklo_epi16(resq_r4_1, resq_r5_1);  // e0 f0 e1 f1 e2 f2 e2 f3
    temp3 = _mm_unpacklo_epi16(resq_r6_1, resq_r7_1);  // g0 h0 g1 h1 g2 h2 g3 h3
    temp2 = _mm_unpackhi_epi16(resq_r4_1, resq_r5_1);  // e4 f4 e5 f5 e6 f6 e7 f7
    temp4 = _mm_unpackhi_epi16(resq_r6_1, resq_r7_1);  // g4 h4 g5 h5 g6 h6 g7 h7
    resq_r4_1 = _mm_unpacklo_epi32(temp1, temp3);      // e0 f0 g0 h0 e1 f1 g1 h1
    resq_r5_1 = _mm_unpackhi_epi32(temp1, temp3);      // e2 f2 g2 h2 e3 f3 g3 h3
    resq_r6_1 = _mm_unpacklo_epi32(temp2, temp4);      // e4 f4 g4 h4 e5 f5 g5 h5
    resq_r7_1 = _mm_unpackhi_epi32(temp2, temp4);      // e6 f6 g6 h6 e7 f7 g7 h7
    /*
     * a0 b0 c0 d0 a1 b1 c1 d1
     * a2 b2 c2 d2 a3 b3 c3 d3
     * a4 b4 c4 d4 a5 b5 c5 d5
     * a6 b6 c6 d6 a7 b7 c7 d7
     * e0 f0 g0 h0 e1 f1 g1 h1
     * e2 f2 g2 h2 e3 f3 g3 h3
     * e4 f4 g4 h4 e5 f5 g5 h5
     * e6 f6 g6 h6 e7 f7 g7 h7
     */
    resq_r0_2 = _mm_unpacklo_epi64(resq_r0_1, resq_r4_1);  // a0 b0 c0 d0 e0 f0 g0 h0
    resq_r1_2 = _mm_unpackhi_epi64(resq_r0_1, resq_r4_1);  // a1 b1 c1 d1 e1 f1 g1 h1
    resq_r2_2 = _mm_unpacklo_epi64(resq_r1_1, resq_r5_1);  // a2 b2 c2 d2 e2 f2 g2 h2
    resq_r3_2 = _mm_unpackhi_epi64(resq_r1_1, resq_r5_1);  // a3 b3 c3 d3 e3 f3 g3 h3
    resq_r4_2 = _mm_unpacklo_epi64(resq_r2_1, resq_r6_1);  // a4 b4 c4 d4 e4 f4 g4 h4
    resq_r5_2 = _mm_unpackhi_epi64(resq_r2_1, resq_r6_1);  // a5 b5 c5 d5 e5 f5 g5 h5
    resq_r6_2 = _mm_unpacklo_epi64(resq_r3_1, resq_r7_1);  // a6 b6 c6 d6 e6 f6 g6 h6
    resq_r7_2 = _mm_unpackhi_epi64(resq_r3_1, resq_r7_1);  // a7 b7 c7 d7 e7 f7 g7 h7

    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r1_2);
    resq_r1_1 = _mm_unpacklo_epi16(resq_r1_2, sign_reg);  // a1 b1 c1 d1 -- 32 bit
    resq_r1_2 = _mm_unpackhi_epi16(resq_r1_2, sign_reg);  // e1 f1 g1 h1 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r3_2);
    resq_r3_1 = _mm_unpacklo_epi16(resq_r3_2, sign_reg);  // a3 b3 c3 d3 -- 32 bit
    resq_r3_2 = _mm_unpackhi_epi16(resq_r3_2, sign_reg);  // e3 f3 g3 h3 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r5_2);
    resq_r5_1 = _mm_unpacklo_epi16(resq_r5_2, sign_reg);  // a5 b5 c5 d5 -- 32 bit
    resq_r5_2 = _mm_unpackhi_epi16(resq_r5_2, sign_reg);  // e5 f5 g5 h5 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r7_2);
    resq_r7_1 = _mm_unpacklo_epi16(resq_r7_2, sign_reg);  // a7 b7 c7 d7 -- 32 bit
    resq_r7_2 = _mm_unpackhi_epi16(resq_r7_2, sign_reg);  // e7 f7 g7 h7 -- 32 bit
    // Transform starts -- horizontal transform
    /*------------------------------------------------------------------*/
    /* y0 = w0 + w4                                                     */
    temp1 = _mm_add_epi16(resq_r0_2, resq_r4_2);
    /* y2 = w0 - w4                                                      */
    temp3 = _mm_sub_epi16(resq_r0_2, resq_r4_2);
    /* y1 = -w3 + w5 - w7 - (w7 >> 1)                                   */
    temp2 = _mm_sub_epi32(resq_r5_1, resq_r3_1);  //-w3+w5
    temp10 = _mm_sub_epi32(resq_r5_2, resq_r3_2);
    temp4 = _mm_sub_epi32(temp2, resq_r7_1);      //-w3+w5-w7
    temp12 = _mm_sub_epi32(temp10, resq_r7_2);
    temp5 = _mm_srai_epi32(resq_r7_1, 1);         // w7>>1
    temp13 = _mm_srai_epi32(resq_r7_2, 1);
    temp2 = _mm_sub_epi32(temp4, temp5);          //-w3+w5-w7 -(w7>>1)
    temp10 = _mm_sub_epi32(temp12, temp13);
    temp2 = _mm_packs_epi32(temp2, temp10);
    /* y3 = w1 + w7 - w3 - (w3 >> 1)                                    */
    temp4 = _mm_add_epi32(resq_r1_1, resq_r7_1);  // w1+w7
    temp12 = _mm_add_epi32(resq_r1_2, resq_r7_2);
    temp4 = _mm_sub_epi32(temp4, resq_r3_1);      // w1+w7-w3
    temp12 = _mm_sub_epi32(temp12, resq_r3_2);
    temp5 = _mm_srai_epi32(resq_r3_1, 1);         // w3>>1
    temp13 = _mm_srai_epi32(resq_r3_2, 1);
    temp4 = _mm_sub_epi32(temp4, temp5);          // w1+w7-w3-(w3>>1)
    temp12 = _mm_sub_epi32(temp12, temp13);
    temp4 = _mm_packs_epi32(temp4, temp12);
    /* y4 = (w2 >> 1) - w6                                              */
    temp5 = _mm_srai_epi16(resq_r2_2, 1);     // w2>>1
    temp5 = _mm_sub_epi16(temp5, resq_r6_2);  //(w2>>1)-w6
    /* y5 = -w1 + w7 + w5 + (w5 >> 1)                                   */
    temp6 = _mm_sub_epi32(resq_r7_1, resq_r1_1);  // w7-w1
    temp14 = _mm_sub_epi32(resq_r7_2, resq_r1_2);
    temp6 = _mm_add_epi32(temp6, resq_r5_1);      // w7-w1+w5
    temp14 = _mm_add_epi32(temp14, resq_r5_2);
    temp7 = _mm_srai_epi32(resq_r5_1, 1);         // w5>>1
    temp15 = _mm_srai_epi32(resq_r5_2, 1);
    temp6 = _mm_add_epi32(temp6, temp7);          // w7-w1_w5+(w5>>1)
    temp14 = _mm_add_epi32(temp14, temp15);
    temp6 = _mm_packs_epi32(temp6, temp14);
    /* y6 = w2 + (w6 >> 1)                                              */
    temp7 = _mm_srai_epi16(resq_r6_2, 1);     // w6>>1
    temp7 = _mm_add_epi16(temp7, resq_r2_2);  //(w6>>1)+w2
    /* y7 = w3 + w5 + w1 + (w1 >> 1)                                    */
    temp8 = _mm_add_epi32(resq_r3_1, resq_r5_1);  // w3+w5
    temp16 = _mm_add_epi32(resq_r3_2, resq_r5_2);
    temp8 = _mm_add_epi32(temp8, resq_r1_1);      // w3+w5+w1
    temp16 = _mm_add_epi32(temp16, resq_r1_2);
    temp17 = _mm_srai_epi32(resq_r1_1, 1);        // w1>>1
    temp18 = _mm_srai_epi32(resq_r1_2, 1);
    temp8 = _mm_add_epi32(temp8, temp17);         // w3+w5+w1+(w1>>1)
    temp16 = _mm_add_epi32(temp16, temp18);
    temp8 = _mm_packs_epi32(temp8, temp16);
    /*------------------------------------------------------------------*/
    /*------------------------------------------------------------------*/
    /* z0 = y0 + y6                                                        */
    resq_r0_1 = _mm_add_epi16(temp1, temp7);
    /* z1 = y1 + (y7 >> 2)                                                */
    resq_r1_1 = _mm_srai_epi16(temp8, 2);
    resq_r1_1 = _mm_add_epi16(resq_r1_1, temp2);
    /* z2 = y2 + y4                                                        */
    resq_r2_1 = _mm_add_epi16(temp3, temp5);
    /* z3 = y3 + (y5 >> 2)                                                */
    resq_r3_1 = _mm_srai_epi16(temp6, 2);
    resq_r3_1 = _mm_add_epi16(resq_r3_1, temp4);
    /* z4 = y2 - y4                                                        */
    resq_r4_1 = _mm_sub_epi16(temp3, temp5);
    /* z5 = (y3 >> 2) - y5                                                 */
    resq_r5_1 = _mm_srai_epi16(temp4, 2);
    resq_r5_1 = _mm_sub_epi16(resq_r5_1, temp6);
    /* z6 = y0 - y6                                                     */
    resq_r6_1 = _mm_sub_epi16(temp1, temp7);
    /* z7 = y7 - (y1 >> 2)                                                 */
    resq_r7_1 = _mm_srai_epi16(temp2, 2);
    resq_r7_1 = _mm_sub_epi16(temp8, resq_r7_1);
    /*------------------------------------------------------------------*/
    /*------------------------------------------------------------------*/
    /* x0 = z0 + z7                                                        */
    temp1 = _mm_add_epi16(resq_r0_1, resq_r7_1);
    /* x1 = z2 + z5                                                        */
    temp2 = _mm_add_epi16(resq_r2_1, resq_r5_1);
    /* x2 = z4 + z3                                                        */
    temp3 = _mm_add_epi16(resq_r4_1, resq_r3_1);
    /* x3 = z6 + z1                                                        */
    temp4 = _mm_add_epi16(resq_r6_1, resq_r1_1);
    /* x4 = z6 - z1                                                        */
    temp5 = _mm_sub_epi16(resq_r6_1, resq_r1_1);
    /* x5 = z4 - z3                                                        */
    temp6 = _mm_sub_epi16(resq_r4_1, resq_r3_1);
    /* x6 = z2 - z5                                                        */
    temp7 = _mm_sub_epi16(resq_r2_1, resq_r5_1);
    /* x7 = z0 - z7                                                        */
    temp8 = _mm_sub_epi16(resq_r0_1, resq_r7_1);
    /*------------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 b0 c0 d0 e0 f0 g0 h0
     *  a1 b1 c1 d1 e1 f1 g1 h1
     *  a2 b2 c2 d2 e2 f2 g2 h2
     *  a3 b3 c3 d3 e3 f3 g3 h3
     */
    temp17 = _mm_unpacklo_epi16(temp1, temp2);       // a0 a1 b0 b1 c0 c1 d0 d1
    temp19 = _mm_unpacklo_epi16(temp3, temp4);       // a2 a3 b2 b3 c2 c3 d2 d3
    temp18 = _mm_unpackhi_epi16(temp1, temp2);       // e0 e1 f0 f1 g0 g1 h0 h1
    temp20 = _mm_unpackhi_epi16(temp3, temp4);       // e2 e3 f2 f3 g2 g3 h2 h3

    resq_r0_1 = _mm_unpacklo_epi32(temp17, temp19);  // a0 a1 a2 a3 b0 b1 b2 b3
    resq_r1_1 = _mm_unpackhi_epi32(temp17, temp19);  // c0 c1 c2 c3 d0 d1 d2 d3
    resq_r2_1 = _mm_unpacklo_epi32(temp18, temp20);  // e0 e1 e2 e3 f0 f1 f2 f3
    resq_r3_1 = _mm_unpackhi_epi32(temp18, temp20);  // g0 g2 g2 g3 h0 h1 h2 h3
    /*
     *  a4 b4 c4 d4 e4 f4 g4 h4
     *  a5 b5 c5 d5 e5 f5 g5 h5
     *  a6 b6 c6 d6 e6 f6 g6 h6
     *  a7 b7 c7 d7 e7 f7 g7 h7
     */
    temp17 = _mm_unpacklo_epi16(temp5, temp6);       // a4 a5 b4 b5 c4 c5 d4 d5
    temp19 = _mm_unpacklo_epi16(temp7, temp8);       // a6 a7 b6 b7 c6 c7 d6 d7
    temp18 = _mm_unpackhi_epi16(temp5, temp6);       // e4 e5 f4 f5 g4 g5 h4 h5
    temp20 = _mm_unpackhi_epi16(temp7, temp8);       // e6 e7 f6 f7 g6 g7 h6 h7

    resq_r4_1 = _mm_unpacklo_epi32(temp17, temp19);  // a4 a5 a6 a7 b4 b5 b6 b7
    resq_r5_1 = _mm_unpackhi_epi32(temp17, temp19);  // c4 c5 c6 c7 d4 d5 d6 d7
    resq_r6_1 = _mm_unpacklo_epi32(temp18, temp20);  // e4 e5 e6 e7 f4 f5 f6 f7
    resq_r7_1 = _mm_unpackhi_epi32(temp18, temp20);  // g4 g5 g6 g7 h4 h5 h6 h7
    /*  a0 a1 a2 a3 b0 b1 b2 b3
     *  c0 c1 c2 c3 d0 d1 d2 d3
     *  e0 e1 e2 e3 f0 f1 f2 f3
     *  g0 g2 g2 g3 h0 h1 h2 h3
     *  a4 a5 a6 a7 b4 b5 b6 b7
     *  c4 c5 c6 c7 d4 d5 d6 d7
     *  e4 e5 e6 e7 f4 f5 f6 f7
     *  g4 g5 g6 g7 h4 h5 h6 h7
     */
    resq_r0_2 = _mm_unpacklo_epi64(resq_r0_1, resq_r4_1);  // a0 a1 a2 a3 a4 a5 a6 a7
    resq_r1_2 = _mm_unpackhi_epi64(resq_r0_1, resq_r4_1);  // b0 b1 b2 b3 b4 b5 b6 b7
    resq_r2_2 = _mm_unpacklo_epi64(resq_r1_1, resq_r5_1);  // c0 c1 c2 c3 c4 c5 c6 c7
    resq_r3_2 = _mm_unpackhi_epi64(resq_r1_1, resq_r5_1);  // d0 d1 d2 d3 d4 d5 d6 d7
    resq_r4_2 = _mm_unpacklo_epi64(resq_r2_1, resq_r6_1);  // e0 e1 e2 e3 e4 e5 e6 e7
    resq_r5_2 = _mm_unpackhi_epi64(resq_r2_1, resq_r6_1);  // f0 f1 f2 f3 f4 f5 f6 f7
    resq_r6_2 = _mm_unpacklo_epi64(resq_r3_1, resq_r7_1);  // g0 g1 g2 g3 g4 g5 g6 g7
    resq_r7_2 = _mm_unpackhi_epi64(resq_r3_1, resq_r7_1);  // h0 h1 h2 h3 h4 h5 h6 h7

    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r1_2);
    resq_r1_1 = _mm_unpacklo_epi16(resq_r1_2, sign_reg);  // a1 b1 c1 d1 -- 32 bit
    resq_r1_2 = _mm_unpackhi_epi16(resq_r1_2, sign_reg);  // e1 f1 g1 h1 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r3_2);
    resq_r3_1 = _mm_unpacklo_epi16(resq_r3_2, sign_reg);  // a3 b3 c3 d3 -- 32 bit
    resq_r3_2 = _mm_unpackhi_epi16(resq_r3_2, sign_reg);  // e3 f3 g3 h3 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r5_2);
    resq_r5_1 = _mm_unpacklo_epi16(resq_r5_2, sign_reg);  // a5 b5 c5 d5 -- 32 bit
    resq_r5_2 = _mm_unpackhi_epi16(resq_r5_2, sign_reg);  // e5 f5 g5 h5 -- 32 bit
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, resq_r7_2);
    resq_r7_1 = _mm_unpacklo_epi16(resq_r7_2, sign_reg);  // a7 b7 c7 d7 -- 32 bit
    resq_r7_2 = _mm_unpackhi_epi16(resq_r7_2, sign_reg);  // e7 f7 g7 h7 -- 32 bit

    zero_8x16b = _mm_setzero_si128();                     // all bits reset to zero
    // Load pred buffer row 0
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r0_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 1
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r1_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 2
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r2_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 3
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r3_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 4
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[4 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r4_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 5
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bit
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[5 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r5_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 6
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[6 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r6_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);
    // Load pred buffer row 7
    // p0 p1 p2 p3 p4 p5 p6 p7 0 0 0 0 0 0 0 0 -- all 8 bits
    predload_r = _mm_loadl_epi64((__m128i *) (&pu1_pred[7 * pred_strd]));
    // p0 p1 p2 p3 p4 p5 p6 p7 -- all 16 bits
    pred_r7_1 = _mm_unpacklo_epi8(predload_r, zero_8x16b);

    rsd_r0 = _mm_loadu_si128((__m128i *) (&pi2_rsd[0]));
    rsd_r1 = _mm_loadu_si128((__m128i *) (&pi2_rsd[1 * rsd_strd]));
    rsd_r2 = _mm_loadu_si128((__m128i *) (&pi2_rsd[2 * rsd_strd]));
    rsd_r3 = _mm_loadu_si128((__m128i *) (&pi2_rsd[3 * rsd_strd]));
    rsd_r4 = _mm_loadu_si128((__m128i *) (&pi2_rsd[4 * rsd_strd]));
    rsd_r5 = _mm_loadu_si128((__m128i *) (&pi2_rsd[5 * rsd_strd]));
    rsd_r6 = _mm_loadu_si128((__m128i *) (&pi2_rsd[6 * rsd_strd]));
    rsd_r7 = _mm_loadu_si128((__m128i *) (&pi2_rsd[7 * rsd_strd]));

    /*--------------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6            */
    /*                                                                    */
    /* Add the prediction and store it back to reconstructed frame buffer */
    /* [Prediction buffer itself in this case]                            */
    /*--------------------------------------------------------------------*/

    /* y0j = w0j + w4j                                                     */
    temp1 = _mm_add_epi16(resq_r0_2, resq_r4_2);
    /* y2j = w0j - w4j                                                      */
    temp3 = _mm_sub_epi16(resq_r0_2, resq_r4_2);
    /* y1j = -w3j + w5j - w7j - (w7j >> 1)                                   */
    temp2 = _mm_sub_epi32(resq_r5_1, resq_r3_1);  //-w3+w5
    temp10 = _mm_sub_epi32(resq_r5_2, resq_r3_2);
    temp4 = _mm_sub_epi32(temp2, resq_r7_1);      //-w3+w5-w7
    temp12 = _mm_sub_epi32(temp10, resq_r7_2);
    temp5 = _mm_srai_epi32(resq_r7_1, 1);         // w7>>1
    temp13 = _mm_srai_epi32(resq_r7_2, 1);
    temp2 = _mm_sub_epi32(temp4, temp5);          //-w3+w5-w7 -(w7>>1)
    temp10 = _mm_sub_epi32(temp12, temp13);
    temp2 = _mm_packs_epi32(temp2, temp10);
    /* y3j = w1j + w7j - w3j - (w3j >> 1)                                    */
    temp4 = _mm_add_epi32(resq_r1_1, resq_r7_1);  // w1+w7
    temp12 = _mm_add_epi32(resq_r1_2, resq_r7_2);
    temp4 = _mm_sub_epi32(temp4, resq_r3_1);      // w1+w7-w3
    temp12 = _mm_sub_epi32(temp12, resq_r3_2);
    temp5 = _mm_srai_epi32(resq_r3_1, 1);         // w3>>1
    temp13 = _mm_srai_epi32(resq_r3_2, 1);
    temp4 = _mm_sub_epi32(temp4, temp5);          // w1+w7-w3-(w3>>1)
    temp12 = _mm_sub_epi32(temp12, temp13);
    temp4 = _mm_packs_epi32(temp4, temp12);
    /* y4j = (w2j >> 1) - w6j                                              */
    temp5 = _mm_srai_epi16(resq_r2_2, 1);     // w2>>1
    temp5 = _mm_sub_epi16(temp5, resq_r6_2);  //(w2>>1)-w6
    /* y5j = -w1j + w7j + w5j + (w5j >> 1)                                   */
    temp6 = _mm_sub_epi32(resq_r7_1, resq_r1_1);  // w7-w1
    temp14 = _mm_sub_epi32(resq_r7_2, resq_r1_2);
    temp6 = _mm_add_epi32(temp6, resq_r5_1);      // w7-w1+w5
    temp14 = _mm_add_epi32(temp14, resq_r5_2);
    temp7 = _mm_srai_epi32(resq_r5_1, 1);         // w5>>1
    temp15 = _mm_srai_epi32(resq_r5_2, 1);
    temp6 = _mm_add_epi32(temp6, temp7);          // w7-w1_w5+(w5>>1)
    temp14 = _mm_add_epi32(temp14, temp15);
    temp6 = _mm_packs_epi32(temp6, temp14);
    /* y6j = w2j + (w6j >> 1)                                              */
    temp7 = _mm_srai_epi16(resq_r6_2, 1);     // w6>>1
    temp7 = _mm_add_epi16(temp7, resq_r2_2);  //(w6>>1)+w2
    /* y7j = w3j + w5j + w1j + (w1j >> 1)                                    */
    temp8 = _mm_add_epi32(resq_r3_1, resq_r5_1);  // w3+w5
    temp16 = _mm_add_epi32(resq_r3_2, resq_r5_2);
    temp8 = _mm_add_epi32(temp8, resq_r1_1);      // w3+w5+w1
    temp16 = _mm_add_epi32(temp16, resq_r1_2);
    temp17 = _mm_srai_epi32(resq_r1_1, 1);        // w1>>1
    temp18 = _mm_srai_epi32(resq_r1_2, 1);
    temp8 = _mm_add_epi32(temp8, temp17);         // w3+w5+w1+(w1>>1)
    temp16 = _mm_add_epi32(temp16, temp18);
    temp8 = _mm_packs_epi32(temp8, temp16);
    /*------------------------------------------------------------------*/
    /*------------------------------------------------------------------*/
    /* z0j = y0j + y6j                                                        */
    resq_r0_1 = _mm_add_epi16(temp1, temp7);
    /* z1j = y1j + (y7j >> 2)                                                */
    resq_r1_1 = _mm_srai_epi16(temp8, 2);
    resq_r1_1 = _mm_add_epi16(resq_r1_1, temp2);
    /* z2j = y2j + y4j                                                        */
    resq_r2_1 = _mm_add_epi16(temp3, temp5);
    /* z3j = y3j + (y5j >> 2)                                                */
    resq_r3_1 = _mm_srai_epi16(temp6, 2);
    resq_r3_1 = _mm_add_epi16(resq_r3_1, temp4);
    /* z4j = y2j - y4j                                                        */
    resq_r4_1 = _mm_sub_epi16(temp3, temp5);
    /* z5j = (y3j >> 2) - y5j                                                 */
    resq_r5_1 = _mm_srai_epi16(temp4, 2);
    resq_r5_1 = _mm_sub_epi16(resq_r5_1, temp6);
    /* z6j = y0j - y6j                                                     */
    resq_r6_1 = _mm_sub_epi16(temp1, temp7);
    /* z7j = y7j - (y1j >> 2)                                                 */
    resq_r7_1 = _mm_srai_epi16(temp2, 2);
    resq_r7_1 = _mm_sub_epi16(temp8, resq_r7_1);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    /* x0j = z0j + z7j                                                        */
    temp1 = _mm_add_epi16(resq_r0_1, resq_r7_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp1);
    temp10 = _mm_unpacklo_epi16(temp1, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp1, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r0 = _mm_add_epi16(temp10, rsd_r0);
    rsd_r0 = _mm_min_epi16(dupmax_8x16b, rsd_r0);
    rsd_r0 = _mm_max_epi16(dupmin_8x16b, rsd_r0);
    temp1 = _mm_add_epi16(rsd_r0, pred_r0_1);
    /* x1j = z2j + z5j                                                        */
    temp2 = _mm_add_epi16(resq_r2_1, resq_r5_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp2);
    temp10 = _mm_unpacklo_epi16(temp2, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp2, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r1 = _mm_add_epi16(temp10, rsd_r1);
    rsd_r1 = _mm_min_epi16(dupmax_8x16b, rsd_r1);
    rsd_r1 = _mm_max_epi16(dupmin_8x16b, rsd_r1);
    temp2 = _mm_add_epi16(rsd_r1, pred_r1_1);
    /* x2j = z4j + z3j                                                        */
    temp3 = _mm_add_epi16(resq_r4_1, resq_r3_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp3);
    temp10 = _mm_unpacklo_epi16(temp3, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp3, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r2 = _mm_add_epi16(temp10, rsd_r2);
    rsd_r2 = _mm_min_epi16(dupmax_8x16b, rsd_r2);
    rsd_r2 = _mm_max_epi16(dupmin_8x16b, rsd_r2);
    temp3 = _mm_add_epi16(rsd_r2, pred_r2_1);
    /* x3j = z6j + z1j                                                        */
    temp4 = _mm_add_epi16(resq_r6_1, resq_r1_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp4);
    temp10 = _mm_unpacklo_epi16(temp4, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp4, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r3 = _mm_add_epi16(temp10, rsd_r3);
    rsd_r3 = _mm_min_epi16(dupmax_8x16b, rsd_r3);
    rsd_r3 = _mm_max_epi16(dupmin_8x16b, rsd_r3);
    temp4 = _mm_add_epi16(rsd_r3, pred_r3_1);
    /* x4j = z6j - z1j                                                        */
    temp5 = _mm_sub_epi16(resq_r6_1, resq_r1_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp5);
    temp10 = _mm_unpacklo_epi16(temp5, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp5, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r4 = _mm_add_epi16(temp10, rsd_r4);
    rsd_r4 = _mm_min_epi16(dupmax_8x16b, rsd_r4);
    rsd_r4 = _mm_max_epi16(dupmin_8x16b, rsd_r4);
    temp5 = _mm_add_epi16(rsd_r4, pred_r4_1);
    /* x5j = z4j - z3j                                                        */
    temp6 = _mm_sub_epi16(resq_r4_1, resq_r3_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp6);
    temp10 = _mm_unpacklo_epi16(temp6, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp6, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r5 = _mm_add_epi16(temp10, rsd_r5);
    rsd_r5 = _mm_min_epi16(dupmax_8x16b, rsd_r5);
    rsd_r5 = _mm_max_epi16(dupmin_8x16b, rsd_r5);
    temp6 = _mm_add_epi16(rsd_r5, pred_r5_1);
    /* x6j = z2j - z5j                                                        */
    temp7 = _mm_sub_epi16(resq_r2_1, resq_r5_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp7);
    temp10 = _mm_unpacklo_epi16(temp7, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp7, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r6 = _mm_add_epi16(temp10, rsd_r6);
    rsd_r6 = _mm_min_epi16(dupmax_8x16b, rsd_r6);
    rsd_r6 = _mm_max_epi16(dupmin_8x16b, rsd_r6);
    temp7 = _mm_add_epi16(rsd_r6, pred_r6_1);
    /* x7j = z0j - z7j                                                        */
    temp8 = _mm_sub_epi16(resq_r0_1, resq_r7_1);
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, temp8);
    temp10 = _mm_unpacklo_epi16(temp8, sign_reg);
    temp11 = _mm_unpackhi_epi16(temp8, sign_reg);
    temp10 = _mm_add_epi32(temp10, value_32);
    temp11 = _mm_add_epi32(temp11, value_32);
    temp10 = _mm_srai_epi32(temp10, 6);
    temp11 = _mm_srai_epi32(temp11, 6);
    temp10 = _mm_packs_epi32(temp10, temp11);
    rsd_r7 = _mm_add_epi16(temp10, rsd_r7);
    rsd_r7 = _mm_min_epi16(dupmax_8x16b, rsd_r7);
    rsd_r7 = _mm_max_epi16(dupmin_8x16b, rsd_r7);
    temp8 = _mm_add_epi16(rsd_r7, pred_r7_1);

    rsd_r01_b0 = _mm_unpacklo_epi64(rsd_r0, rsd_r1);
    rsd_r23_b0 = _mm_unpacklo_epi64(rsd_r2, rsd_r3);
    rsd_r45_b2 = _mm_unpacklo_epi64(rsd_r4, rsd_r5);
    rsd_r67_b2 = _mm_unpacklo_epi64(rsd_r6, rsd_r7);

    rsd_r01_b1 = _mm_unpackhi_epi64(rsd_r0, rsd_r1);
    rsd_r23_b1 = _mm_unpackhi_epi64(rsd_r2, rsd_r3);
    rsd_r45_b3 = _mm_unpackhi_epi64(rsd_r4, rsd_r5);
    rsd_r67_b3 = _mm_unpackhi_epi64(rsd_r6, rsd_r7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r23_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_45_b2 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r45_b2, zero_8x16b));  // return 1 if all zeros, else 0
    row_67_b2 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r67_b2, zero_8x16b));  // return 1 if all zeros, else 0

    row_01_b1 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r01_b1, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b1 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r23_b1, zero_8x16b));  // return 1 if all zeros, else 0
    row_45_b3 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r45_b3, zero_8x16b));  // return 1 if all zeros, else 0
    row_67_b3 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_r67_b3, zero_8x16b));  // return 1 if all zeros, else 0

    /*------------------------------------------------------------------*/
    // Clipping the results to 8 bits
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);  // sign check
    temp1 = _mm_and_si128(temp1, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp2, zero_8x16b);  // sign check
    temp2 = _mm_and_si128(temp2, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp3, zero_8x16b);  // sign check
    temp3 = _mm_and_si128(temp3, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp4, zero_8x16b);  // sign check
    temp4 = _mm_and_si128(temp4, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp5, zero_8x16b);  // sign check
    temp5 = _mm_and_si128(temp5, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp6, zero_8x16b);  // sign check
    temp6 = _mm_and_si128(temp6, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp7, zero_8x16b);  // sign check
    temp7 = _mm_and_si128(temp7, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp8, zero_8x16b);  // sign check
    temp8 = _mm_and_si128(temp8, sign_reg);

    resq_r0_2 = _mm_packus_epi16(temp1, zero_8x16b);
    resq_r1_2 = _mm_packus_epi16(temp2, zero_8x16b);
    resq_r2_2 = _mm_packus_epi16(temp3, zero_8x16b);
    resq_r3_2 = _mm_packus_epi16(temp4, zero_8x16b);
    resq_r4_2 = _mm_packus_epi16(temp5, zero_8x16b);
    resq_r5_2 = _mm_packus_epi16(temp6, zero_8x16b);
    resq_r6_2 = _mm_packus_epi16(temp7, zero_8x16b);
    resq_r7_2 = _mm_packus_epi16(temp8, zero_8x16b);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), resq_r0_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[out_strd]), resq_r1_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * out_strd]), resq_r2_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * out_strd]), resq_r3_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[4 * out_strd]), resq_r4_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[5 * out_strd]), resq_r5_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[6 * out_strd]), resq_r6_2);
    _mm_storel_epi64((__m128i *) (&pu1_out[7 * out_strd]), resq_r7_2);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_4x4_dc_sse42           */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_iquant_itrans_residual_recon_4x4_dc_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 u4_qp_div_6, WORD16 *pi2_tmp, WORD32 iq_start_idx, WORD16 *pi2_dc_ld_addr)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, out_8x16b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, out_8x16b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, out_8x16b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, out_8x16b_3, out_16x8b_3;
    __m128i rsd_8x16b_01, rsd_8x16b_23;

    __m128i i_macro_8x16b, dupmax_8x16b, dupmin_8x16b;
    __m128i zero_8x16b = _mm_setzero_si128();
    WORD32 i4_nnz, row_01, row_23;
    WORD32 q0;
    WORD16 i_macro;
    WORD16 rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    UNUSED(pi2_tmp);

    if(iq_start_idx == 0)
    {
        q0 = pi2_src[0];
        INV_QUANT(q0, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);
    }
    else
    {
        q0 = pi2_dc_ld_addr[0];  // Restoring dc value for intra case3
    }
    i_macro = ((q0 + 32) >> 6);

    i_macro_8x16b = _mm_set1_epi16(i_macro);
    dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1)));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1) + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1)));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1) + rsd_strd));

    rsd_8x16b_0 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_3);

    rsd_8x16b_0 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_0);
    rsd_8x16b_0 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_1);
    rsd_8x16b_1 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_2);
    rsd_8x16b_2 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_3);
    rsd_8x16b_3 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_3);

    rsd_8x16b_01 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);

    row_01 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01, zero_8x16b));  // return 1 if all zeros, else 0
    row_23 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);

    *((WORD32 *) (pu1_out)) = _mm_cvtsi128_si32(out_16x8b_0);
    *((WORD32 *) (pu1_out + out_strd)) = _mm_cvtsi128_si32(out_16x8b_1);
    *((WORD32 *) (pu1_out + (out_strd << 1))) = _mm_cvtsi128_si32(out_16x8b_2);
    *((WORD32 *) (pu1_out + (out_strd * 3))) = _mm_cvtsi128_si32(out_16x8b_3);

    i4_nnz = !(row_01 && row_23);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_8x8_dc_sse42           */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_iquant_itrans_residual_recon_8x8_dc_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscale_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 qp_div, WORD16 *pi2_tmp, WORD32 iq_start_idx, WORD16 *pi2_dc_ld_addr)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, out_8x16b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, out_8x16b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, out_8x16b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, out_8x16b_3, out_16x8b_3;
    __m128i pred_16x8b_4, pred_8x16b_4, rsd_8x16b_4, out_8x16b_4, out_16x8b_4;
    __m128i pred_16x8b_5, pred_8x16b_5, rsd_8x16b_5, out_8x16b_5, out_16x8b_5;
    __m128i pred_16x8b_6, pred_8x16b_6, rsd_8x16b_6, out_8x16b_6, out_16x8b_6;
    __m128i pred_16x8b_7, pred_8x16b_7, rsd_8x16b_7, out_8x16b_7, out_16x8b_7;
    __m128i rsd_8x16b_01_b0, rsd_8x16b_23_b0, rsd_8x16b_45_b2, rsd_8x16b_67_b2;
    __m128i rsd_8x16b_01_b1, rsd_8x16b_23_b1, rsd_8x16b_45_b3, rsd_8x16b_67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 pred_strd2 = (pred_strd << 1);
    WORD32 pred_strd4 = (pred_strd << 2);
    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);
    WORD32 out_strd2 = (out_strd << 1);
    WORD32 out_strd4 = (out_strd << 2);

    __m128i i_macro_8x16b, dupmax_8x16b, dupmin_8x16b;
    WORD32 q;
    WORD16 i_macro;
    WORD32 rnd_fact = (qp_div < 6) ? (1 << (5 - qp_div)) : 0;
    UNUSED(pi2_tmp);
    UNUSED(iq_start_idx);
    UNUSED(pi2_dc_ld_addr);
    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform. Note : DC coeff is not scaled     */
    /*************************************************************/
    q = pi2_src[0];
    INV_QUANT(q, pu2_iscale_mat[0], pu2_weigh_mat[0], qp_div, rnd_fact, 6);
    i_macro = (q + 32) >> 6;

    i_macro_8x16b = _mm_set1_epi16(i_macro);
    dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_0 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_3);
    rsd_8x16b_4 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_4);
    rsd_8x16b_5 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_5);
    rsd_8x16b_6 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_6);
    rsd_8x16b_7 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_7);

    rsd_8x16b_0 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_0);
    rsd_8x16b_0 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_1);
    rsd_8x16b_1 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_2);
    rsd_8x16b_2 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_3);
    rsd_8x16b_3 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_3);
    rsd_8x16b_4 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_4);
    rsd_8x16b_4 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_4);
    rsd_8x16b_5 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_5);
    rsd_8x16b_5 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_5);
    rsd_8x16b_6 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_6);
    rsd_8x16b_6 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_6);
    rsd_8x16b_7 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_7);
    rsd_8x16b_7 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_7);

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(rsd_8x16b_01_b0);  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(rsd_8x16b_23_b0);
    row_01_b1 = _mm_test_all_ones(rsd_8x16b_01_b1);
    row_23_b1 = _mm_test_all_ones(rsd_8x16b_23_b1);

    row_45_b2 = _mm_test_all_ones(rsd_8x16b_45_b2);
    row_67_b2 = _mm_test_all_ones(rsd_8x16b_67_b2);
    row_45_b3 = _mm_test_all_ones(rsd_8x16b_45_b3);
    row_67_b3 = _mm_test_all_ones(rsd_8x16b_67_b3);

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_chroma_4x4_dc_sse42    */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

void isvcd_iquant_itrans_residual_recon_chroma_4x4_dc_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 u4_qp_div_6, WORD16 *pi2_tmp, WORD16 *pi2_dc_src)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, rsd_16x8b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, rsd_16x8b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, rsd_16x8b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, rsd_16x8b_3, out_16x8b_3;

    __m128i i_macro_8x16b, dupmax_8x16b, dupmin_8x16b;
    __m128i chroma_mask, chroma_mask2;
    __m128i zero_8x16b = _mm_setzero_si128();
    WORD32 q0;
    WORD16 i_macro;
    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);

    q0 = pi2_dc_src[0];  // Restoring dc value for intra case3
    i_macro = ((q0 + 32) >> 6);

    i_macro_8x16b = _mm_set1_epi16(i_macro);
    dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1)));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1) + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1)));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1) + rsd_strd));

    rsd_8x16b_0 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_add_epi16(i_macro_8x16b, rsd_8x16b_3);

    rsd_8x16b_0 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_0);
    rsd_8x16b_0 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_1);
    rsd_8x16b_1 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_2);
    rsd_8x16b_2 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_min_epi16(dupmax_8x16b, rsd_8x16b_3);
    rsd_8x16b_3 = _mm_max_epi16(dupmin_8x16b, rsd_8x16b_3);

    rsd_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    rsd_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    rsd_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    rsd_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);

    chroma_mask = _mm_set1_epi16(0xFF00);
    chroma_mask2 = _mm_set1_epi16(0x00FF);
    out_16x8b_0 = _mm_loadu_si128((__m128i *) (&pu1_out[0]));
    out_16x8b_1 = _mm_loadu_si128((__m128i *) (&pu1_out[out_strd]));
    out_16x8b_2 = _mm_loadu_si128((__m128i *) (&pu1_out[(out_strd << 1)]));
    out_16x8b_3 = _mm_loadu_si128((__m128i *) (&pu1_out[(out_strd << 1) + out_strd]));

    out_16x8b_0 = _mm_and_si128(out_16x8b_0, chroma_mask);
    out_16x8b_1 = _mm_and_si128(out_16x8b_1, chroma_mask);
    out_16x8b_2 = _mm_and_si128(out_16x8b_2, chroma_mask);
    out_16x8b_3 = _mm_and_si128(out_16x8b_3, chroma_mask);

    rsd_16x8b_0 = _mm_packus_epi16(rsd_8x16b_0, zero_8x16b);
    rsd_16x8b_1 = _mm_packus_epi16(rsd_8x16b_1, zero_8x16b);
    rsd_16x8b_2 = _mm_packus_epi16(rsd_8x16b_2, zero_8x16b);
    rsd_16x8b_3 = _mm_packus_epi16(rsd_8x16b_3, zero_8x16b);

    rsd_8x16b_0 = _mm_and_si128(rsd_16x8b_0, chroma_mask2);
    rsd_8x16b_1 = _mm_and_si128(rsd_16x8b_1, chroma_mask2);
    rsd_8x16b_2 = _mm_and_si128(rsd_16x8b_2, chroma_mask2);
    rsd_8x16b_3 = _mm_and_si128(rsd_16x8b_3, chroma_mask2);

    out_16x8b_0 = _mm_add_epi8(rsd_8x16b_0, out_16x8b_0);
    out_16x8b_1 = _mm_add_epi8(rsd_8x16b_1, out_16x8b_1);
    out_16x8b_2 = _mm_add_epi8(rsd_8x16b_2, out_16x8b_2);
    out_16x8b_3 = _mm_add_epi8(rsd_8x16b_3, out_16x8b_3);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + (out_strd << 1)), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + (out_strd * 3)), out_16x8b_3);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_recon_chroma_4x4_sse42       */
/*                                                                           */
/*  Description   : this function computes the recon output from the         */
/*                  IQ+IT+RESD                                               */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : i4_nnz                                                   */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

void isvcd_iquant_itrans_residual_recon_chroma_4x4_sse42(
    WORD16 *pi2_src, UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out, WORD32 pred_strd,
    WORD32 rsd_strd, WORD32 out_strd, const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat,
    UWORD32 u4_qp_div_6, WORD16 *pi2_tmp, WORD16 *pi2_dc_src)
{
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i rsd_r0, rsd_r1, rsd_r2, rsd_r3;
    __m128i dequant_r0_r1, dequant_r2_r3;
    __m128i zero_8x16b = _mm_setzero_si128();  // all bits reset to zero
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);
    __m128i chroma_mask = _mm_set1_epi16(0xFF00);
    __m128i chroma_mask2 = _mm_set1_epi16(0x00FF);
    __m128i dupmax_8x16b = _mm_set1_epi16(RSD_MAX);
    __m128i dupmin_8x16b = _mm_set1_epi16(RSD_MIN);

    __m128i out_16x8b_0, out_16x8b_1, out_16x8b_2, out_16x8b_3;

    UNUSED(pi2_tmp);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/
    // a00 a01 a02 a03 a10 a11 a12 a13  -- the source matrix 0th,1st row
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));
    // a20 a21 a22 a23 a30 a31 a32 a33 --  the source matrix 2nd,3rd row
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));
    // b00 b01 b02 b03 b10 b11 b12 b13 -- the  scaling matrix 0th,1st row
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));
    // b20 b21 b22 b23 b30 b31 b32 b33 -- the  scaling matrix 2nd,3rd row
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));
    // q00 q01 q02 q03 q10 q11  q12 q13 -- all 16 bits
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));
    // q20 q21 q22 q23 q30 q31 q32 q33 -- all 16 bits
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10  b11*q11 b12*q12 b13*q13 -- 16 bit result
    temp0 = _mm_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);
    // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11 b12*q12 b13*q13 -- 16 bit result
    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);
    // b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);
    // b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);
    // b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);
    // b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);  // a00 0 a01 0 a02 0 a03 0 -- 16 bit long
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);  // a10 0 a11 0 a12 0 a13 0 -- 16 bit long
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);  // a20 0 a21 0 a22 0 a23 0 -- 16 bit long
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);  // a30 0 a31 0 a32 0 a33 0 -- 16 bit long

    // a00*b00*q00 a10*b10*q10 a20*b20*q20 a30*b30 q30 -- 32 bits long
    temp4 = _mm_madd_epi16(src_r0, temp4);
    temp5 = _mm_madd_epi16(src_r1, temp5);
    temp6 = _mm_madd_epi16(src_r2, temp6);
    temp7 = _mm_madd_epi16(src_r3, temp7);

    if(u4_qp_div_6 >= 4)
    {
        resq_r0 = _mm_slli_epi32(temp4, u4_qp_div_6 - 4);
        resq_r1 = _mm_slli_epi32(temp5, u4_qp_div_6 - 4);
        resq_r2 = _mm_slli_epi32(temp6, u4_qp_div_6 - 4);
        resq_r3 = _mm_slli_epi32(temp7, u4_qp_div_6 - 4);
    }
    else
    {
        temp4 = _mm_add_epi32(temp4, add_rshift);
        temp5 = _mm_add_epi32(temp5, add_rshift);
        temp6 = _mm_add_epi32(temp6, add_rshift);
        temp7 = _mm_add_epi32(temp7, add_rshift);
        resq_r0 = _mm_srai_epi32(temp4, 4 - u4_qp_div_6);
        resq_r1 = _mm_srai_epi32(temp5, 4 - u4_qp_div_6);
        resq_r2 = _mm_srai_epi32(temp6, 4 - u4_qp_div_6);
        resq_r3 = _mm_srai_epi32(temp7, 4 - u4_qp_div_6);
    }

    resq_r0 = _mm_insert_epi32(resq_r0, (WORD32) pi2_dc_src[0], 0);
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
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);  // a0 b0 a1 b1
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);  // c0 d0 c1 d1
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);  // a2 b2 a3 b3
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);  // c2 d2 c3 d3
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);    // a0 b0 c0 d0
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);    // a1 b1 c1 d1
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);    // a2 b2 c2 d2
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);    // a3 b3 c3 d3
    // Transform starts -- horizontal transform
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);     //(w1>>1)
    temp2 = _mm_sub_epi32(temp2, resq_r3);  //(w1>>1) - w3
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);  //(w3>>1) + w1
    temp3 = _mm_add_epi32(temp3, resq_r1);
    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    resq_r0 = _mm_add_epi32(temp0, temp3);
    /* x1 = z1 + z2                                             */
    resq_r1 = _mm_add_epi32(temp1, temp2);
    /* x2 = z1 - z2                                             */
    resq_r2 = _mm_sub_epi32(temp1, temp2);
    /* x3 = z0 - z3                                             */
    resq_r3 = _mm_sub_epi32(temp0, temp3);
    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);  // a0 a1 b0 b1
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);  // a2 a3 b2 b3
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);  // c0 c1 d0 d1
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);  // c2 c3 d2 d3
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);    // a0 a1 a2 a3
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);    // b0 b1 b2 b3
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);    // c0 c1 c2 c3
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);    // d0 d1 d2 d3
    // Transform ends -- horizontal transform

    // Load pred buffer
    // p00 p01 p02 p03 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    // p10 p11 p12 p13 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[pred_strd]));
    // p20 p21 p22 p23 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * pred_strd]));
    // p30 p31 p32 p33 0 0 0 0 0 0 0 0 -- all 8 bits
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * pred_strd]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);  // p00 p01 p02 p03 -- all 16 bits
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);  // p10 p11 p12 p13 -- all 16 bits
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);  // p20 p21 p22 p23 -- all 16 bits
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);  // p30 p31 p32 p33 -- all 16 bits

    // Load resd buffer
    rsd_r0 = _mm_loadu_si128((__m128i *) (&pi2_rsd[0]));
    rsd_r1 = _mm_loadu_si128((__m128i *) (&pi2_rsd[rsd_strd]));
    rsd_r2 = _mm_loadu_si128((__m128i *) (&pi2_rsd[2 * rsd_strd]));
    rsd_r3 = _mm_loadu_si128((__m128i *) (&pi2_rsd[3 * rsd_strd]));

    /*--------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6      */
    /*                                                              */
    /* Add the prediction and store it back to same buffer          */
    /*--------------------------------------------------------------*/
    /* z0j = y0j + y2j                                                        */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1j = y0j - y2j                                                        */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2j = (y1j>>1) - y3j */
    temp2 = _mm_srai_epi32(resq_r1, 1);  //(y1j>>1)
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);  //(y3j>>1)
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    temp4 = _mm_add_epi16(temp4, rsd_r0);
    temp4 = _mm_min_epi16(dupmax_8x16b, temp4);
    temp4 = _mm_max_epi16(dupmin_8x16b, temp4);
    temp4 = _mm_add_epi16(temp4, pred_r0);
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    temp5 = _mm_add_epi16(temp5, rsd_r1);
    temp5 = _mm_min_epi16(dupmax_8x16b, temp5);
    temp5 = _mm_max_epi16(dupmin_8x16b, temp5);
    temp5 = _mm_add_epi16(temp5, pred_r1);
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    temp6 = _mm_add_epi16(temp6, rsd_r2);
    temp6 = _mm_min_epi16(dupmax_8x16b, temp6);
    temp6 = _mm_max_epi16(dupmin_8x16b, temp6);
    temp6 = _mm_add_epi16(temp6, pred_r2);
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);
    temp7 = _mm_add_epi16(temp7, rsd_r3);
    temp7 = _mm_min_epi16(dupmax_8x16b, temp7);
    temp7 = _mm_max_epi16(dupmin_8x16b, temp7);
    temp7 = _mm_add_epi16(temp7, pred_r3);

    out_16x8b_0 = _mm_loadu_si128((__m128i *) (&pu1_out[0]));
    out_16x8b_1 = _mm_loadu_si128((__m128i *) (&pu1_out[out_strd]));
    out_16x8b_2 = _mm_loadu_si128((__m128i *) (&pu1_out[(out_strd << 1)]));
    out_16x8b_3 = _mm_loadu_si128((__m128i *) (&pu1_out[(out_strd << 1) + out_strd]));

    out_16x8b_0 = _mm_and_si128(out_16x8b_0, chroma_mask);
    out_16x8b_1 = _mm_and_si128(out_16x8b_1, chroma_mask);
    out_16x8b_2 = _mm_and_si128(out_16x8b_2, chroma_mask);
    out_16x8b_3 = _mm_and_si128(out_16x8b_3, chroma_mask);

    temp4 = _mm_packus_epi16(temp4, zero_8x16b);
    temp5 = _mm_packus_epi16(temp5, zero_8x16b);
    temp6 = _mm_packus_epi16(temp6, zero_8x16b);
    temp7 = _mm_packus_epi16(temp7, zero_8x16b);

    temp4 = _mm_and_si128(temp4, chroma_mask2);
    temp5 = _mm_and_si128(temp5, chroma_mask2);
    temp6 = _mm_and_si128(temp6, chroma_mask2);
    temp7 = _mm_and_si128(temp7, chroma_mask2);

    out_16x8b_0 = _mm_add_epi8(temp4, out_16x8b_0);
    out_16x8b_1 = _mm_add_epi8(temp5, out_16x8b_1);
    out_16x8b_2 = _mm_add_epi8(temp6, out_16x8b_2);
    out_16x8b_3 = _mm_add_epi8(temp7, out_16x8b_3);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + (out_strd << 1)), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + (out_strd * 3)), out_16x8b_3);
}
