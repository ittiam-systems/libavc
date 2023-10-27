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
 *  isvc_iquant_itrans_recon_sse42.c
 *
 * @brief
 *  Contains function definitions for inverse  quantization, inverse
 * transform and reconstruction
 *
 * @author
 *  Mohit [100664]
 *
 * @par List of Functions:
 *  - isvc_iquant_itrans_recon_4x4_sse42()
 *  - isvc_iquant_itrans_recon_chroma_4x4_sse42()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
#include <immintrin.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "isvc_structs.h"
#include "isvc_trans_quant_itrans_iquant.h"

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
 * @param[in] i4_pred_stride,
 *  Prediction buffer stride
 *
 * @param[in] i4_out_stride
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

void isvc_iquant_itrans_recon_4x4_sse42(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                        buffer_container_t *ps_res_pred, buffer_container_t *ps_res,
                                        buffer_container_t *ps_rec,
                                        iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
                                        WORD16 *pi2_tmp, WORD16 *pi2_dc_src, WORD32 i4_iq_start_idx,
                                        UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_tmp_ptr = pi2_tmp;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i sign_reg, dequant_r0_r1, dequant_r2_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);

    ASSERT(4 == i4_src_stride);
    ASSERT(0 == u1_res_accumulate);

    UNUSED(i4_src_stride);
    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/

    /* a00 a01 a02 a03 a10 a11 a12 a13 -- the source
     matrix 0th,1st row */
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));

    /* a20 a21 a22 a23 a30 a31 a32 a33 -- the
      source matrix 2nd,3rd row */
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13 -- the
     scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));

    /* b20 b21 b22 b23 b30 b31 b32 b33 --b12 b13 -- the
     the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));

    /* q00 q01 q02 q03 q10 q11
     q12 q13 -- all 16 bits */
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));

    /* q20 q21 q22 q23 q30 q31
     q32 q33 -- all 16 bits */
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    /* b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11
     b12*q12 b13*q13 -- 16 bit result */
    temp0 = _mm_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);

    /* b20*q20 b21*q21 b22*q22 b23*q23 b30*q30 b31*q31
     b32*q32 b33*q33 -- 16 bit result */
    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    /* a00 0 a01 0 a02 0 a03 0 -- 16 bit long */
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);
    /* a10 0 a11 0 a12 0 a13 0 -- 16 bit long */
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);
    /* a20 0 a21 0 a22 0 a23 0 -- 16 bit long */
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);
    /* a30 0 a31 0 a32 0 a33 0 -- 16 bit long */
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);

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

    if(i4_iq_start_idx == 1) resq_r0 = _mm_insert_epi32(resq_r0, (WORD32) pi2_dc_src[0], 0);
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

    /* a0 b0 a1 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* c0 d0 c1 d1 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* a2 b2 a3 b3 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 d2 c3 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 b0 c0 d0 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* a1 b1 c1 d1 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* a2 b2 c2 d2 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* a3 b3 c3 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform starts -- horizontal transform */
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);
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

    /* a0 a1 b0 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* a2 a3 b2 b3 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* c0 c1 d0 d1 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 c3 d2 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 a1 a2 a3 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* b0 b1 b2 b3 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* c0 c1 c2 c3 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* d0 d1 d2 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform ends -- horizontal transform */

    temp0 = _mm_packs_epi32(resq_r0, resq_r1);
    temp1 = _mm_packs_epi32(resq_r2, resq_r3);

    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[2 * 4]), temp1);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

    pred_r0 = _mm_unpacklo_epi64(pred_r0, pred_r1);
    pred_r1 = _mm_unpacklo_epi64(pred_r2, pred_r3);

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
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);

    /* 32-bit to 16-bit conversion */
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp4 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp4 = _mm_min_epi16(temp4, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp5 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp5 = _mm_min_epi16(temp5, pos_255_8x16b);

    temp0 = _mm_add_epi16(temp4, pred_r0);
    temp1 = _mm_add_epi16(temp5, pred_r1);

    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(temp0, zero_8x16b);
    temp0 = _mm_and_si128(temp0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, sign_reg);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    *pu4_out = _mm_cvtsi128_si32(resq_r0);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r1);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r2);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r3);
}

void isvc_iquant_itrans_recon_res_4x4_sse42(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                            buffer_container_t *ps_res_pred,
                                            buffer_container_t *ps_res, buffer_container_t *ps_rec,
                                            iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
                                            WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
                                            WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_tmp_ptr = pi2_tmp;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i sign_reg, dequant_r0_r1, dequant_r2_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);

    ASSERT(4 == i4_src_stride);
    ASSERT(0 == u1_res_accumulate);

    UNUSED(i4_src_stride);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/

    /* a00 a01 a02 a03 a10 a11 a12 a13 -- the source
    matrix 0th,1st row */
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));

    /* a20 a21 a22 a23 a30 a31 a32 a33 -- the
    source matrix 2nd,3rd row */
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13 -- the
    scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));

    /* b20 b21 b22 b23 b30 b31 b32 b33 --b12 b13 -- the
    the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));

    /* q00 q01 q02 q03 q10 q11
    q12 q13 -- all 16 bits */
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));

    /* q20 q21 q22 q23 q30 q31
    q32 q33 -- all 16 bits */
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    /* b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11
    b12*q12 b13*q13 -- 16 bit result */
    temp0 = _mm_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);

    /* b20*q20 b21*q21 b22*q22 b23*q23 b30*q30 b31*q31
    b32*q32 b33*q33 -- 16 bit result */
    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    /* a00 0 a01 0 a02 0 a03 0 -- 16 bit long */
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);
    /* a10 0 a11 0 a12 0 a13 0 -- 16 bit long */
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);
    /* a20 0 a21 0 a22 0 a23 0 -- 16 bit long */
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);
    /* a30 0 a31 0 a32 0 a33 0 -- 16 bit long */
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);

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

    if(i4_iq_start_idx == 1) resq_r0 = _mm_insert_epi32(resq_r0, (WORD32) pi2_dc_src[0], 0);
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

    /* a0 b0 a1 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* c0 d0 c1 d1 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* a2 b2 a3 b3 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 d2 c3 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 b0 c0 d0 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* a1 b1 c1 d1 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* a2 b2 c2 d2 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* a3 b3 c3 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform starts -- horizontal transform */
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);
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

    /* a0 a1 b0 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* a2 a3 b2 b3 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* c0 c1 d0 d1 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 c3 d2 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 a1 a2 a3 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* b0 b1 b2 b3 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* c0 c1 c2 c3 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* d0 d1 d2 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform ends -- horizontal transform */

    temp0 = _mm_packs_epi32(resq_r0, resq_r1);
    temp1 = _mm_packs_epi32(resq_r2, resq_r3);

    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[2 * 4]), temp1);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

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
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);

    /* 32-bit to 16-bit conversion */
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp0 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(temp0, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp1 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp1 = _mm_min_epi16(temp1, pos_255_8x16b);

    _mm_storel_epi64((__m128i *) (&pi2_res[0]), temp0);
    _mm_storel_epi64((__m128i *) (&pi2_res[2 * i4_res_stride]), temp1);

    temp4 = _mm_add_epi16(temp0, pred_r0);
    temp0 = _mm_srli_si128(temp0, 8);
    _mm_storel_epi64((__m128i *) (&pi2_res[i4_res_stride]), temp0);

    temp6 = _mm_add_epi16(temp1, pred_r2);
    temp1 = _mm_srli_si128(temp1, 8);
    _mm_storel_epi64((__m128i *) (&pi2_res[3 * i4_res_stride]), temp1);

    temp5 = _mm_add_epi16(temp0, pred_r1);
    temp7 = _mm_add_epi16(temp1, pred_r3);

    temp4 = _mm_cvtepi16_epi32(temp4);
    temp5 = _mm_cvtepi16_epi32(temp5);
    temp6 = _mm_cvtepi16_epi32(temp6);
    temp7 = _mm_cvtepi16_epi32(temp7);

    /* 32-bit to 16-bit conversion */
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);
    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(temp0, zero_8x16b);
    temp0 = _mm_and_si128(temp0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, sign_reg);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    *pu4_out = _mm_cvtsi128_si32(resq_r0);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r1);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r2);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r3);
}

void isvc_iquant_itrans_recon_res_4x4_with_res_acc_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_tmp_ptr = pi2_tmp;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD16 *pi2_res_pred = (WORD16 *) ps_res_pred->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_res_pred_stride = ps_res_pred->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i res_pred_r0, res_pred_r1, res_pred_r2, res_pred_r3;
    __m128i res_r0, res_r1, res_r2, res_r3;
    __m128i sign_reg, dequant_r0_r1, dequant_r2_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);

    ASSERT(4 == i4_src_stride);
    ASSERT(1 == u1_res_accumulate);

    UNUSED(i4_src_stride);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/

    /* a00 a01 a02 a03 a10 a11 a12 a13 -- the source
     matrix 0th,1st row */
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));

    /* a20 a21 a22 a23 a30 a31 a32 a33 -- the
      source matrix 2nd,3rd row */
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13 -- the
     scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));

    /* b20 b21 b22 b23 b30 b31 b32 b33 --b12 b13 -- the
     the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));

    /* q00 q01 q02 q03 q10 q11
     q12 q13 -- all 16 bits */
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));

    /* q20 q21 q22 q23 q30 q31
     q32 q33 -- all 16 bits */
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    /* b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11
     b12*q12 b13*q13 -- 16 bit result */
    temp0 = _mm_mullo_epi16(scalemat_r0_r1, dequant_r0_r1);

    /* b20*q20 b21*q21 b22*q22 b23*q23 b30*q30 b31*q31
     b32*q32 b33*q33 -- 16 bit result */
    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    /* a00 0 a01 0 a02 0 a03 0 -- 16 bit long */
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);
    /* a10 0 a11 0 a12 0 a13 0 -- 16 bit long */
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);
    /* a20 0 a21 0 a22 0 a23 0 -- 16 bit long */
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);
    /* a30 0 a31 0 a32 0 a33 0 -- 16 bit long */
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);

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

    if(i4_iq_start_idx == 1) resq_r0 = _mm_insert_epi32(resq_r0, (WORD32) pi2_dc_src[0], 0);
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

    /* a0 b0 a1 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* c0 d0 c1 d1 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* a2 b2 a3 b3 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 d2 c3 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 b0 c0 d0 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* a1 b1 c1 d1 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* a2 b2 c2 d2 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* a3 b3 c3 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform starts -- horizontal transform */
    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);
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

    /* a0 a1 b0 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* a2 a3 b2 b3 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* c0 c1 d0 d1 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 c3 d2 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 a1 a2 a3 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* b0 b1 b2 b3 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* c0 c1 c2 c3 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* d0 d1 d2 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform ends -- horizontal transform */

    temp0 = _mm_packs_epi32(resq_r0, resq_r1);
    temp1 = _mm_packs_epi32(resq_r2, resq_r3);

    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_tmp_ptr[2 * 4]), temp1);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

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
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    res_r0 = temp4;
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    res_r1 = temp5;
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    res_r2 = temp6;
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);
    res_r3 = temp7;

    /* Accumulating res */
    res_pred_r0 = _mm_loadl_epi64((__m128i *) &pi2_res_pred[0]);
    res_pred_r1 = _mm_loadl_epi64((__m128i *) &pi2_res_pred[i4_res_pred_stride]);
    res_pred_r2 = _mm_loadl_epi64((__m128i *) &pi2_res_pred[2 * i4_res_pred_stride]);
    res_pred_r3 = _mm_loadl_epi64((__m128i *) &pi2_res_pred[3 * i4_res_pred_stride]);

    res_pred_r0 = _mm_cvtepi16_epi32(res_pred_r0);
    res_pred_r1 = _mm_cvtepi16_epi32(res_pred_r1);
    res_pred_r2 = _mm_cvtepi16_epi32(res_pred_r2);
    res_pred_r3 = _mm_cvtepi16_epi32(res_pred_r3);

    temp0 = _mm_add_epi32(res_r0, res_pred_r0);
    temp1 = _mm_add_epi32(res_r1, res_pred_r1);
    temp2 = _mm_add_epi32(res_r2, res_pred_r2);
    temp3 = _mm_add_epi32(res_r3, res_pred_r3);

    temp0 = _mm_packs_epi32(temp0, temp1);
    temp1 = _mm_packs_epi32(temp2, temp3);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp0 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(temp0, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp1 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp1 = _mm_min_epi16(temp1, pos_255_8x16b);

    _mm_storel_epi64((__m128i *) (&pi2_res[0]), temp0);
    _mm_storel_epi64((__m128i *) (&pi2_res[2 * i4_res_stride]), temp1);

    temp4 = _mm_add_epi16(temp0, pred_r0);
    temp0 = _mm_srli_si128(temp0, 8);
    _mm_storel_epi64((__m128i *) (&pi2_res[i4_res_stride]), temp0);

    temp6 = _mm_add_epi16(temp1, pred_r2);
    temp1 = _mm_srli_si128(temp1, 8);
    _mm_storel_epi64((__m128i *) (&pi2_res[3 * i4_res_stride]), temp1);

    temp5 = _mm_add_epi16(temp0, pred_r1);
    temp7 = _mm_add_epi16(temp1, pred_r3);

    temp4 = _mm_cvtepi16_epi32(temp4);
    temp5 = _mm_cvtepi16_epi32(temp5);
    temp6 = _mm_cvtepi16_epi32(temp6);
    temp7 = _mm_cvtepi16_epi32(temp7);

    /* 32-bit to 16-bit conversion */
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);
    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(temp0, zero_8x16b);
    temp0 = _mm_and_si128(temp0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, sign_reg);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    *pu4_out = _mm_cvtsi128_si32(resq_r0);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r1);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r2);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(resq_r3);
}

void isvc_iquant_itrans_recon_res_chroma_4x4_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD16 *pi2_res_ptr = pi2_res;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i sign_reg, dequant_r0_r1, dequant_r2_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);
    __m128i chroma_mask = _mm_set1_epi16(0xFF);
    __m128i out_r0, out_r1, out_r2, out_r3;
    __m128i res_r0, res_r1, res_r2, res_r3;

    ASSERT(4 == i4_src_stride);
    ASSERT(0 == u1_res_accumulate);

    UNUSED(i4_src_stride);
    UNUSED(u1_res_accumulate);
    UNUSED(ps_res_pred);
    UNUSED(i4_iq_start_idx);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/
    /* a00 a01 a02 a03 a10 a11 a12 a13 -- the source
    matrix 0th,1st row */
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));

    /* a20 a21 a22 a23 a30 a31 a32 a33 -- the
    source matrix 2nd,3rd row */
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13 -- the
    scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));

    /* b20 b21 b22 b23 b30 b31 b32 b33 --b12 b13 -- the
    the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));

    /* q00 q01 q02 q03 q10 q11
    q12 q13 -- all 16 bits */
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));

    /* q20 q21 q22 q23 q30 q31
    q32 q33 -- all 16 bits */
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    temp0 = _mm_mullo_epi16(scalemat_r0_r1,
                            dequant_r0_r1);  // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11
                                             // b12*q12 b13*q13 -- 16 bit result

    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    /* a00 0 a01 0 a02 0 a03 0 -- 16 bit long */
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);
    /* a10 0 a11 0 a12 0 a13 0 -- 16 bit long */
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);
    /* a20 0 a21 0 a22 0 a23 0 -- 16 bit long */
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);
    /* a30 0 a31 0 a32 0 a33 0 -- 16 bit long */
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);

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
    /* a0 b0 a1 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* c0 d0 c1 d1 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* a2 b2 a3 b3 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 d2 c3 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 b0 c0 d0 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* a1 b1 c1 d1 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* a2 b2 c2 d2 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* a3 b3 c3 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform starts -- horizontal transform */

    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3 = w1 + (w3 >> 1)                                      */
    temp3 = _mm_srai_epi32(resq_r3, 1);
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
    /* a0 a1 b0 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* a2 a3 b2 b3 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* c0 c1 d0 d1 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 c3 d2 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 a1 a2 a3 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* b0 b1 b2 b3 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* c0 c1 c2 c3 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* d0 d1 d2 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform ends -- horizontal transform */

    temp0 = _mm_packs_epi32(resq_r0, resq_r1);
    temp1 = _mm_packs_epi32(resq_r2, resq_r3);

    _mm_storeu_si128((__m128i *) (&pi2_tmp[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_tmp[2 * 4]), temp1);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    pred_r0 = _mm_cvtepu16_epi32(pred_r0);
    pred_r1 = _mm_cvtepu16_epi32(pred_r1);
    pred_r2 = _mm_cvtepu16_epi32(pred_r2);
    pred_r3 = _mm_cvtepu16_epi32(pred_r3);

    /*--------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6      */
    /*                                                              */
    /* Add the prediction and store it back to same buffer          */
    /*--------------------------------------------------------------*/
    /* z0j = y0j + y2j                                         */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1j = y0j - y2j                                                        */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2j = (y1j>>1) - y3j */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);

    /* 32-bit to 16-bit conversion */
    temp0 = _mm_packs_epi32(temp4, temp5);
    temp1 = _mm_packs_epi32(temp6, temp7);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp0 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(temp0, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp1 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp1 = _mm_min_epi16(temp1, pos_255_8x16b);

    chroma_mask = _mm_set1_epi32(0xffff0000);
    out_r0 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[0 * i4_res_stride]));
    out_r1 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[1 * i4_res_stride]));
    out_r2 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[2 * i4_res_stride]));
    out_r3 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[3 * i4_res_stride]));

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    res_r0 = _mm_cvtepu16_epi32(temp0);
    res_r2 = _mm_cvtepu16_epi32(temp1);
    res_r1 = _mm_srli_si128(temp0, 8);
    res_r3 = _mm_srli_si128(temp1, 8);
    res_r1 = _mm_cvtepu16_epi32(res_r1);
    res_r3 = _mm_cvtepu16_epi32(res_r3);

    out_r0 = _mm_add_epi16(out_r0, res_r0);
    out_r1 = _mm_add_epi16(out_r1, res_r1);
    out_r2 = _mm_add_epi16(out_r2, res_r2);
    out_r3 = _mm_add_epi16(out_r3, res_r3);

    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[0 * i4_res_stride]), out_r0);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[1 * i4_res_stride]), out_r1);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[2 * i4_res_stride]), out_r2);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[3 * i4_res_stride]), out_r3);

    resq_r0 = _mm_add_epi16(pred_r0, res_r0);
    resq_r1 = _mm_add_epi16(pred_r1, res_r1);
    resq_r2 = _mm_add_epi16(pred_r2, res_r2);
    resq_r3 = _mm_add_epi16(pred_r3, res_r3);

    temp0 = _mm_packus_epi32(resq_r0, resq_r1);
    temp1 = _mm_packus_epi32(resq_r2, resq_r3);

    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(temp0, zero_8x16b);
    temp0 = _mm_and_si128(temp0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, sign_reg);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    resq_r0 = _mm_cvtepu8_epi16(resq_r0);
    resq_r1 = _mm_cvtepu8_epi16(resq_r1);
    resq_r2 = _mm_cvtepu8_epi16(resq_r2);
    resq_r3 = _mm_cvtepu8_epi16(resq_r3);

    chroma_mask = _mm_set1_epi16(0xff00);
    out_r0 = _mm_loadl_epi64((__m128i *) (&pu1_out[0]));
    out_r1 = _mm_loadl_epi64((__m128i *) (&pu1_out[i4_out_stride]));
    out_r2 = _mm_loadl_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]));
    out_r3 = _mm_loadl_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]));

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    out_r0 = _mm_add_epi8(out_r0, resq_r0);
    out_r1 = _mm_add_epi8(out_r1, resq_r1);
    out_r2 = _mm_add_epi8(out_r2, resq_r2);
    out_r3 = _mm_add_epi8(out_r3, resq_r3);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), out_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[i4_out_stride]), out_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]), out_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]), out_r3);
}

void isvc_iquant_itrans_recon_res_chroma_4x4_with_res_acc_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD16 *pi2_res_pred = (WORD16 *) ps_res_pred->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_res_pred_stride = ps_res_pred->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    __m128i src_r0_r1, src_r2_r3;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i res_pred_r0, res_pred_r1, res_pred_r2, res_pred_r3;
    __m128i res_r0, res_r1, res_r2, res_r3;
    __m128i dequant_r0_r1, dequant_r2_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i reg_chroma = _mm_set1_epi32(0xFFFF);
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    __m128i resq_r0, resq_r1, resq_r2, resq_r3;
    __m128i add_rshift = _mm_set1_epi32((u4_qp_div_6 < 4) ? (1 << (3 - u4_qp_div_6)) : 0);
    __m128i value_32 = _mm_set1_epi32(32);
    __m128i chroma_mask = _mm_set1_epi16(0xFF);
    __m128i out_r0, out_r1, out_r2, out_r3;
    __m128i mask_r0;

    ASSERT(4 == i4_src_stride);
    ASSERT(1 == u1_res_accumulate);

    UNUSED(i4_src_stride);
    UNUSED(u1_res_accumulate);
    UNUSED(i4_iq_start_idx);

    /*************************************************************/
    /* Dequantization of coefficients. Will be replaced by SIMD  */
    /* operations on platform                                    */
    /*************************************************************/
    /* a00 a01 a02 a03 a10 a11 a12 a13 -- the source
    matrix 0th,1st row */
    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));

    /* a20 a21 a22 a23 a30 a31 a32 a33 -- the
    source matrix 2nd,3rd row */
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13 -- the
    scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat));

    /* b20 b21 b22 b23 b30 b31 b32 b33 --b12 b13 -- the
    the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_iscal_mat + 8));

    /* q00 q01 q02 q03 q10 q11
    q12 q13 -- all 16 bits */
    dequant_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat));

    /* q20 q21 q22 q23 q30 q31
    q32 q33 -- all 16 bits */
    dequant_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_weigh_mat + 8));

    temp0 = _mm_mullo_epi16(scalemat_r0_r1,
                            dequant_r0_r1);  // b00*q00 b01*q01 b02*q02 b03*q03 b10*q10 b11*q11
                                             // b12*q12 b13*q13 -- 16 bit result

    temp1 = _mm_mullo_epi16(scalemat_r2_r3, dequant_r2_r3);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp4 = _mm_unpacklo_epi16(temp0, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp5 = _mm_unpackhi_epi16(temp0, zero_8x16b);

    /* b00*q00 0 b01*q01 0 b02*q02 0 b03*q03 0 -- 16 bit long */
    temp6 = _mm_unpacklo_epi16(temp1, zero_8x16b);

    /* b10*q10 0 b11*q11 0 b12*q12 0 b13*q13 0 -- 16 bit long */
    temp7 = _mm_unpackhi_epi16(temp1, zero_8x16b);

    /* a00 0 a01 0 a02 0 a03 0 -- 16 bit long */
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, zero_8x16b);
    /* a10 0 a11 0 a12 0 a13 0 -- 16 bit long */
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, zero_8x16b);
    /* a20 0 a21 0 a22 0 a23 0 -- 16 bit long */
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, zero_8x16b);
    /* a30 0 a31 0 a32 0 a33 0 -- 16 bit long */
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, zero_8x16b);

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
    /* a0 b0 a1 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* c0 d0 c1 d1 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* a2 b2 a3 b3 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 d2 c3 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 b0 c0 d0 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* a1 b1 c1 d1 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* a2 b2 c2 d2 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* a3 b3 c3 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform starts -- horizontal transform */

    /*------------------------------------------------------------------*/
    /* z0 = w0 + w2                                             */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1 = w0 - w2                                             */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2 = (w1 >> 1) - w3                                      */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
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
    /* a0 a1 b0 b1 */
    temp1 = _mm_unpacklo_epi32(resq_r0, resq_r1);
    /* a2 a3 b2 b3 */
    temp3 = _mm_unpacklo_epi32(resq_r2, resq_r3);
    /* c0 c1 d0 d1 */
    temp2 = _mm_unpackhi_epi32(resq_r0, resq_r1);
    /* c2 c3 d2 d3 */
    temp4 = _mm_unpackhi_epi32(resq_r2, resq_r3);
    /* a0 a1 a2 a3 */
    resq_r0 = _mm_unpacklo_epi64(temp1, temp3);
    /* b0 b1 b2 b3 */
    resq_r1 = _mm_unpackhi_epi64(temp1, temp3);
    /* c0 c1 c2 c3 */
    resq_r2 = _mm_unpacklo_epi64(temp2, temp4);
    /* d0 d1 d2 d3 */
    resq_r3 = _mm_unpackhi_epi64(temp2, temp4);
    /* Transform ends -- horizontal transform */

    temp0 = _mm_packs_epi32(resq_r0, resq_r1);
    temp1 = _mm_packs_epi32(resq_r2, resq_r3);

    _mm_storeu_si128((__m128i *) (&pi2_tmp[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_tmp[2 * 4]), temp1);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    /*--------------------------------------------------------------*/
    /* IDCT [ Vertical transformation] and Xij = (xij + 32)>>6      */
    /*                                                              */
    /* Add the prediction and store it back to same buffer          */
    /*--------------------------------------------------------------*/
    /* z0j = y0j + y2j                                         */
    temp0 = _mm_add_epi32(resq_r0, resq_r2);
    /* z1j = y0j - y2j                                                        */
    temp1 = _mm_sub_epi32(resq_r0, resq_r2);
    /* z2j = (y1j>>1) - y3j */
    temp2 = _mm_srai_epi32(resq_r1, 1);
    temp2 = _mm_sub_epi32(temp2, resq_r3);
    /* z3j = y1j + (y3j>>1) */
    temp3 = _mm_srai_epi32(resq_r3, 1);
    temp3 = _mm_add_epi32(temp3, resq_r1);

    /* x0j = z0j + z3j                                                        */
    temp4 = _mm_add_epi32(temp0, temp3);
    temp4 = _mm_add_epi32(temp4, value_32);
    temp4 = _mm_srai_epi32(temp4, 6);
    res_r0 = temp4;
    /* x1j = z1j + z2j                                                        */
    temp5 = _mm_add_epi32(temp1, temp2);
    temp5 = _mm_add_epi32(temp5, value_32);
    temp5 = _mm_srai_epi32(temp5, 6);
    res_r1 = temp5;
    /* x2j = z1j - z2j                                                        */
    temp6 = _mm_sub_epi32(temp1, temp2);
    temp6 = _mm_add_epi32(temp6, value_32);
    temp6 = _mm_srai_epi32(temp6, 6);
    res_r2 = temp6;
    /* x3j = z0j - z3j                                                        */
    temp7 = _mm_sub_epi32(temp0, temp3);
    temp7 = _mm_add_epi32(temp7, value_32);
    temp7 = _mm_srai_epi32(temp7, 6);
    res_r3 = temp7;

    res_pred_r0 = _mm_loadu_si128((__m128i *) &pi2_res_pred[0 * i4_res_pred_stride]);
    res_pred_r1 = _mm_loadu_si128((__m128i *) &pi2_res_pred[1 * i4_res_pred_stride]);
    res_pred_r2 = _mm_loadu_si128((__m128i *) &pi2_res_pred[2 * i4_res_pred_stride]);
    res_pred_r3 = _mm_loadu_si128((__m128i *) &pi2_res_pred[3 * i4_res_pred_stride]);

    res_pred_r0 = _mm_and_si128(res_pred_r0, reg_chroma);
    res_pred_r1 = _mm_and_si128(res_pred_r1, reg_chroma);
    res_pred_r2 = _mm_and_si128(res_pred_r2, reg_chroma);
    res_pred_r3 = _mm_and_si128(res_pred_r3, reg_chroma);

    temp0 = _mm_packs_epi32(res_r0, res_r1);
    temp1 = _mm_packs_epi32(res_r2, res_r3);

    res_r0 = _mm_cvtepu16_epi32(temp0);
    res_r2 = _mm_cvtepu16_epi32(temp1);
    res_r1 = _mm_srli_si128(temp0, 8);
    res_r3 = _mm_srli_si128(temp1, 8);
    res_r1 = _mm_cvtepu16_epi32(res_r1);
    res_r3 = _mm_cvtepu16_epi32(res_r3);

    res_r0 = _mm_add_epi16(res_pred_r0, res_r0);
    res_r1 = _mm_add_epi16(res_pred_r1, res_r1);
    res_r2 = _mm_add_epi16(res_pred_r2, res_r2);
    res_r3 = _mm_add_epi16(res_pred_r3, res_r3);

    temp0 = _mm_packus_epi32(res_r0, res_r1);
    temp1 = _mm_packus_epi32(res_r2, res_r3);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp0 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(temp0, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp1 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp1 = _mm_min_epi16(temp1, pos_255_8x16b);

    res_r0 = _mm_cvtepu16_epi32(temp0);
    res_r1 = _mm_srli_si128(temp0, 8);
    res_r1 = _mm_cvtepu16_epi32(res_r1);

    res_r2 = _mm_cvtepu16_epi32(temp1);
    res_r3 = _mm_srli_si128(temp1, 8);
    res_r3 = _mm_cvtepu16_epi32(res_r3);

    chroma_mask = _mm_set1_epi32(0xffff0000);
    out_r0 = _mm_loadu_si128((__m128i *) (&pi2_res[0 * i4_res_stride]));
    out_r1 = _mm_loadu_si128((__m128i *) (&pi2_res[1 * i4_res_stride]));
    out_r2 = _mm_loadu_si128((__m128i *) (&pi2_res[2 * i4_res_stride]));
    out_r3 = _mm_loadu_si128((__m128i *) (&pi2_res[3 * i4_res_stride]));

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    out_r0 = _mm_add_epi16(out_r0, res_r0);
    out_r1 = _mm_add_epi16(out_r1, res_r1);
    out_r2 = _mm_add_epi16(out_r2, res_r2);
    out_r3 = _mm_add_epi16(out_r3, res_r3);

    _mm_storeu_si128((__m128i *) (&pi2_res[0 * i4_res_stride]), out_r0);
    _mm_storeu_si128((__m128i *) (&pi2_res[1 * i4_res_stride]), out_r1);
    _mm_storeu_si128((__m128i *) (&pi2_res[2 * i4_res_stride]), out_r2);
    _mm_storeu_si128((__m128i *) (&pi2_res[3 * i4_res_stride]), out_r3);

    pred_r0 = _mm_cvtepu16_epi32(pred_r0);
    pred_r1 = _mm_cvtepu16_epi32(pred_r1);
    pred_r2 = _mm_cvtepu16_epi32(pred_r2);
    pred_r3 = _mm_cvtepu16_epi32(pred_r3);

    resq_r0 = _mm_add_epi16(pred_r0, res_r0);
    resq_r1 = _mm_add_epi16(pred_r1, res_r1);
    resq_r2 = _mm_add_epi16(pred_r2, res_r2);
    resq_r3 = _mm_add_epi16(pred_r3, res_r3);

    temp0 = _mm_packus_epi32(resq_r0, resq_r1);
    temp1 = _mm_packus_epi32(resq_r2, resq_r3);

    /* Clipping the results to 8 bits */
    mask_r0 = _mm_cmpgt_epi16(temp0, zero_8x16b);
    temp0 = _mm_and_si128(temp0, mask_r0);
    mask_r0 = _mm_cmpgt_epi16(temp1, zero_8x16b);
    temp1 = _mm_and_si128(temp1, mask_r0);

    resq_r0 = _mm_packus_epi16(temp0, temp1);
    resq_r1 = _mm_srli_si128(resq_r0, 4);
    resq_r2 = _mm_srli_si128(resq_r1, 4);
    resq_r3 = _mm_srli_si128(resq_r2, 4);

    resq_r0 = _mm_cvtepu8_epi16(resq_r0);
    resq_r1 = _mm_cvtepu8_epi16(resq_r1);
    resq_r2 = _mm_cvtepu8_epi16(resq_r2);
    resq_r3 = _mm_cvtepu8_epi16(resq_r3);

    chroma_mask = _mm_set1_epi16(0xFF00);
    out_r0 = _mm_loadl_epi64((__m128i *) (&pu1_out[0 * i4_out_stride]));
    out_r1 = _mm_loadl_epi64((__m128i *) (&pu1_out[1 * i4_out_stride]));
    out_r2 = _mm_loadl_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]));
    out_r3 = _mm_loadl_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]));

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    out_r0 = _mm_add_epi8(out_r0, resq_r0);
    out_r1 = _mm_add_epi8(out_r1, resq_r1);
    out_r2 = _mm_add_epi8(out_r2, resq_r2);
    out_r3 = _mm_add_epi8(out_r3, resq_r3);

    _mm_storel_epi64((__m128i *) (&pu1_out[0 * i4_out_stride]), out_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[1 * i4_out_stride]), out_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]), out_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]), out_r3);
}

void isvc_iquant_itrans_recon_dc_4x4_sse42(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                           buffer_container_t *ps_res_pred,
                                           buffer_container_t *ps_res, buffer_container_t *ps_rec,
                                           iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
                                           WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
                                           WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD32 *pu4_out = (UWORD32 *) pu1_out;
    WORD32 q0 = ((WORD16 *) (ps_src->pv_data))[0];
    WORD16 i_macro, rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i sign_reg;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i temp4, temp5, temp6, temp7;
    __m128i value_add;

    ASSERT(0 == u1_res_accumulate);

    UNUSED(pi2_tmp);
    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    INV_QUANT(q0, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);

    /* Restoring dc value for intra case */
    if(i4_iq_start_idx != 0)
    {
        q0 = pi2_dc_src[0];
    }

    i_macro = ((q0 + 32) >> 6);

    value_add = _mm_set1_epi16(i_macro);

    zero_8x16b = _mm_setzero_si128();

    /* Load pred buffer */

    /* p00 p01 p02 p03 0 0 0 0 -- all 8 bits */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));

    /* p10 p11 p12 p13 0 0 0 0 -- all 8 bits */
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));

    /* p20 p21 p22 p23 0 0 0 0 -- all 8 bits */
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));

    /* p30 p31 p32 p33 0 0 0 0 -- all 8 bits */
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

    pred_r0 = _mm_unpacklo_epi64(pred_r0, pred_r1);
    pred_r2 = _mm_unpacklo_epi64(pred_r2, pred_r3);

    temp4 = _mm_add_epi16(value_add, pred_r0);
    temp5 = _mm_add_epi16(value_add, pred_r2);
    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(temp4, zero_8x16b);
    temp4 = _mm_and_si128(temp4, sign_reg);
    sign_reg = _mm_cmpgt_epi16(temp5, zero_8x16b);
    temp5 = _mm_and_si128(temp5, sign_reg);

    temp4 = _mm_packus_epi16(temp4, temp5);
    temp5 = _mm_srli_si128(temp4, 4);
    temp6 = _mm_srli_si128(temp5, 4);
    temp7 = _mm_srli_si128(temp6, 4);

    *pu4_out = _mm_cvtsi128_si32(temp4);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(temp5);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(temp6);
    pu1_out += i4_out_stride;
    pu4_out = (UWORD32 *) (pu1_out);
    *(pu4_out) = _mm_cvtsi128_si32(temp7);
}

void isvc_iquant_itrans_recon_res_chroma_4x4_dc_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD16 *pi2_res_ptr = pi2_res;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    /* DC value won't be dequantized for chroma
    inverse transform */
    WORD16 q0 = pi2_dc_src[0];
    WORD16 i_macro = ((q0 + 32) >> 6);

    __m128i pred_r0, pred_r1, pred_r2, pred_r3, sign_reg;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i chroma_mask = _mm_set1_epi16(0xFF);
    __m128i value_add = _mm_set1_epi16(isvc_get_residue(i_macro, 0, 0));
    __m128i out_r0, out_r1, out_r2, out_r3;

    ASSERT(0 == u1_res_accumulate);

    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(ps_res_pred);
    UNUSED(i4_iq_start_idx);
    UNUSED(u1_res_accumulate);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));

    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));

    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));

    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    /* Mask alternate pred values from the interleaved pred buf */
    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    /* Pack the first four 16 bit values of 2 regs into a single reg*/
    pred_r0 = _mm_unpacklo_epi64(pred_r0, pred_r1);
    pred_r2 = _mm_unpacklo_epi64(pred_r2, pred_r3);

    /* Compute out pixel by adding res to pred */
    pred_r0 = _mm_add_epi16(value_add, pred_r0);
    pred_r2 = _mm_add_epi16(value_add, pred_r2);

    /* Convert res from 16 bits to 32 bits  */
    value_add = _mm_cvtepu16_epi32(value_add);

    out_r0 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[0 * i4_res_stride]));
    out_r1 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[1 * i4_res_stride]));
    out_r2 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[2 * i4_res_stride]));
    out_r3 = _mm_loadu_si128((__m128i *) (&pi2_res_ptr[3 * i4_res_stride]));

    /* Mask the loaded res in order to save the U/V res data computed in
    this function call without thrashing the U/V res data that was saved
    during an earlier function call */
    chroma_mask = _mm_set1_epi32(0xffff0000);
    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    /* Save the res in alternate locations */
    out_r0 = _mm_add_epi16(out_r0, value_add);
    out_r1 = _mm_add_epi16(out_r1, value_add);
    out_r2 = _mm_add_epi16(out_r2, value_add);
    out_r3 = _mm_add_epi16(out_r3, value_add);

    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[0 * i4_res_stride]), out_r0);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[1 * i4_res_stride]), out_r1);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[2 * i4_res_stride]), out_r2);
    _mm_storeu_si128((__m128i *) (&pi2_res_ptr[3 * i4_res_stride]), out_r3);
    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    sign_reg = _mm_cmpgt_epi16(pred_r0, zero_8x16b);
    pred_r0 = _mm_and_si128(pred_r0, sign_reg);
    sign_reg = _mm_cmpgt_epi16(pred_r2, zero_8x16b);
    pred_r2 = _mm_and_si128(pred_r2, sign_reg);

    pred_r0 = _mm_packus_epi16(pred_r0, pred_r2);
    pred_r1 = _mm_srli_si128(pred_r0, 4);
    pred_r2 = _mm_srli_si128(pred_r1, 4);
    pred_r3 = _mm_srli_si128(pred_r2, 4);

    /* p00 p01 p02 p03 -- all 16 bits */
    pred_r0 = _mm_unpacklo_epi8(pred_r0, zero_8x16b);
    /* p10 p11 p12 p13 -- all 16 bits */
    pred_r1 = _mm_unpacklo_epi8(pred_r1, zero_8x16b);
    /* p20 p21 p22 p23 -- all 16 bits */
    pred_r2 = _mm_unpacklo_epi8(pred_r2, zero_8x16b);
    /* p30 p31 p32 p33 -- all 16 bits */
    pred_r3 = _mm_unpacklo_epi8(pred_r3, zero_8x16b);

    /* Load interleaved out buffer */
    out_r0 = _mm_loadl_epi64((__m128i *) (&pu1_out[0]));
    out_r1 = _mm_loadl_epi64((__m128i *) (&pu1_out[i4_out_stride]));
    out_r2 = _mm_loadl_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]));
    out_r3 = _mm_loadl_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]));

    /* Mask the interleaved out buf in order to save the U/V out pixel computed in
    this function call without thrashing the U/V out pixel that was saved
    during an earlier function call */
    chroma_mask = _mm_set1_epi16(0xFF00);

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    /* Save the out pixels in alternate locations */
    out_r0 = _mm_add_epi8(out_r0, pred_r0);
    out_r1 = _mm_add_epi8(out_r1, pred_r1);
    out_r2 = _mm_add_epi8(out_r2, pred_r2);
    out_r3 = _mm_add_epi8(out_r3, pred_r3);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), out_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[i4_out_stride]), out_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]), out_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]), out_r3);
}

void isvc_iquant_itrans_recon_res_chroma_4x4_dc_with_res_acc_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD16 *pi2_res_pred = (WORD16 *) ps_res_pred->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_res_pred_stride = ps_res_pred->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    /* DC value won't be dequantized for chroma
    inverse transform */
    WORD16 q0 = pi2_dc_src[0];
    WORD16 i_macro = ((q0 + 32) >> 6);

    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i chroma_mask = _mm_set1_epi16(0xFF);
    __m128i reg_chroma = _mm_set_epi16(0, 0xFFFF, 0, 0xFFFF, 0, 0xFFFF, 0, 0xFFFF);
    __m128i value_add = _mm_set1_epi16(i_macro);
    __m128i out_r0, out_r1, out_r2, out_r3;
    __m128i res_r0, res_r1, res_r2, res_r3;
    __m128i res_pred_r0, res_pred_r1, res_pred_r2, res_pred_r3;
    __m128i temp0, temp1;
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));

    ASSERT(1 == u1_res_accumulate);

    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(i4_iq_start_idx);
    UNUSED(u1_res_accumulate);

    /* Load pred buffer */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));

    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));

    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));

    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));
    /* Mask alternate pred values from the interleaved pred buf */
    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    /* Pack the first four 16 bit values of 2 regs into a single reg*/
    pred_r0 = _mm_unpacklo_epi64(pred_r0, pred_r1);
    pred_r2 = _mm_unpacklo_epi64(pred_r2, pred_r3);

    /* Accumulating res */

    /* load res pred buffer */
    res_pred_r0 = _mm_loadu_si128((__m128i *) &pi2_res_pred[0 * i4_res_pred_stride]);
    res_pred_r1 = _mm_loadu_si128((__m128i *) &pi2_res_pred[1 * i4_res_pred_stride]);
    res_pred_r2 = _mm_loadu_si128((__m128i *) &pi2_res_pred[2 * i4_res_pred_stride]);
    res_pred_r3 = _mm_loadu_si128((__m128i *) &pi2_res_pred[3 * i4_res_pred_stride]);

    /* Mask res pred and retain alternate values */
    res_pred_r0 = _mm_and_si128(res_pred_r0, reg_chroma);
    res_pred_r1 = _mm_and_si128(res_pred_r1, reg_chroma);
    res_pred_r2 = _mm_and_si128(res_pred_r2, reg_chroma);
    res_pred_r3 = _mm_and_si128(res_pred_r3, reg_chroma);

    /* Convert to 32 bits */
    res_r0 = _mm_cvtepu16_epi32(value_add);
    res_r2 = _mm_cvtepu16_epi32(value_add);
    res_r1 = _mm_cvtepu16_epi32(value_add);
    res_r3 = _mm_cvtepu16_epi32(value_add);

    /* Add res pred to the res obtained from inv transform */
    res_r0 = _mm_add_epi16(res_pred_r0, res_r0);
    res_r1 = _mm_add_epi16(res_pred_r1, res_r1);
    res_r2 = _mm_add_epi16(res_pred_r2, res_r2);
    res_r3 = _mm_add_epi16(res_pred_r3, res_r3);

    /* Convert 32 bit res of the format [a0 0 a1 0 a2 0 a3 0] to
    16 bits of the format [a0 a1 a2 a3] using hadd [ao + 0,
    a1 + 0, a2 + 0, a3 + 0] To be optimized */
    temp0 = _mm_hadd_epi16(res_r0, res_r1);
    temp1 = _mm_hadd_epi16(res_r2, res_r3);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp0 = _mm_max_epi16(temp0, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(temp0, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    temp1 = _mm_max_epi16(temp1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp1 = _mm_min_epi16(temp1, pos_255_8x16b);

    /* Compute out pixel by adding res to pred */
    pred_r0 = _mm_add_epi16(temp0, pred_r0);
    pred_r2 = _mm_add_epi16(temp1, pred_r2);

    res_r0 = _mm_cvtepu16_epi32(temp0);
    res_r2 = _mm_cvtepu16_epi32(temp1);
    res_r1 = _mm_srli_si128(temp0, 8);
    res_r3 = _mm_srli_si128(temp1, 8);
    res_r1 = _mm_cvtepu16_epi32(res_r1);
    res_r3 = _mm_cvtepu16_epi32(res_r3);

    /* Load res buffer */
    out_r0 = _mm_loadu_si128((__m128i *) (&pi2_res[0 * i4_res_stride]));
    out_r1 = _mm_loadu_si128((__m128i *) (&pi2_res[1 * i4_res_stride]));
    out_r2 = _mm_loadu_si128((__m128i *) (&pi2_res[2 * i4_res_stride]));
    out_r3 = _mm_loadu_si128((__m128i *) (&pi2_res[3 * i4_res_stride]));

    /* Mask the loaded res in order to save the U/V res data computed in
    this function call without thrashing the U/V res data that was saved
    during an earlier function call */
    chroma_mask = _mm_set1_epi32(0xffff0000);

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    /* Save the res in alternate locations */
    out_r0 = _mm_add_epi16(out_r0, res_r0);
    out_r1 = _mm_add_epi16(out_r1, res_r1);
    out_r2 = _mm_add_epi16(out_r2, res_r2);
    out_r3 = _mm_add_epi16(out_r3, res_r3);

    _mm_storeu_si128((__m128i *) (&pi2_res[0 * i4_res_stride]), out_r0);
    _mm_storeu_si128((__m128i *) (&pi2_res[1 * i4_res_stride]), out_r1);
    _mm_storeu_si128((__m128i *) (&pi2_res[2 * i4_res_stride]), out_r2);
    _mm_storeu_si128((__m128i *) (&pi2_res[3 * i4_res_stride]), out_r3);
    /*------------------------------------------------------------------*/
    /* Clipping the results to 8 bits */
    pred_r0 = _mm_packus_epi16(pred_r0, pred_r2);
    pred_r1 = _mm_srli_si128(pred_r0, 4);
    pred_r2 = _mm_srli_si128(pred_r1, 4);
    pred_r3 = _mm_srli_si128(pred_r2, 4);

    /* p00 p01 p02 p03 -- all 16 bits */
    pred_r0 = _mm_unpacklo_epi8(pred_r0, zero_8x16b);
    /* p10 p11 p12 p13 -- all 16 bits */
    pred_r1 = _mm_unpacklo_epi8(pred_r1, zero_8x16b);
    /* p20 p21 p22 p23 -- all 16 bits */
    pred_r2 = _mm_unpacklo_epi8(pred_r2, zero_8x16b);
    /* p30 p31 p32 p33 -- all 16 bits */
    pred_r3 = _mm_unpacklo_epi8(pred_r3, zero_8x16b);

    /* Load interleaved out buffer */
    out_r0 = _mm_loadl_epi64((__m128i *) (&pu1_out[0]));
    out_r1 = _mm_loadl_epi64((__m128i *) (&pu1_out[i4_out_stride]));
    out_r2 = _mm_loadl_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]));
    out_r3 = _mm_loadl_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]));

    /* Mask the interleaved out buf in order to save the U/V out pixel computed in
    this function call without thrashing the U/V out pixel that was saved
    during an earlier function call */
    chroma_mask = _mm_set1_epi16(0xFF00);

    out_r0 = _mm_and_si128(out_r0, chroma_mask);
    out_r1 = _mm_and_si128(out_r1, chroma_mask);
    out_r2 = _mm_and_si128(out_r2, chroma_mask);
    out_r3 = _mm_and_si128(out_r3, chroma_mask);

    /* Save the out pixels in alternate locations */
    out_r0 = _mm_add_epi8(out_r0, pred_r0);
    out_r1 = _mm_add_epi8(out_r1, pred_r1);
    out_r2 = _mm_add_epi8(out_r2, pred_r2);
    out_r3 = _mm_add_epi8(out_r3, pred_r3);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), out_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[i4_out_stride]), out_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * i4_out_stride]), out_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * i4_out_stride]), out_r3);
}
