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

 * *******************************************************************************

 * * @file
 *  isvc_resi_trans_quant_sse42.c
 *
 * @brief
 *  Contains function
 * definitions single stage  forward transform for H.264
 *  It will calculate
 * the residue, do the cf and then do quantization
 *
 * @author
 *  Mohit
 * [100664]
 *
 * @par List of Functions:
 *  -
 * isvc_resi_trans_quant_4x4_sse42()
 *  -
 * isvc_resi_trans_quant_chroma_4x4_sse42()
 *
 * @remarks
 *  None
 *

 * *******************************************************************************

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

/**|
*******************************************************************************
*
*
* @brief
*   This function performs forward transform and quantization on a 4*4
* block
*
* @par Description:
*   The function accepts source buffer and
* estimation buffer. From these, it
*   computes the residue. This is residue
* is then transformed and quantized.
*   The transform and quantization are in
* placed computed. They use the residue
*   buffer for this.
*
* @param[in]
* pu1_src
*   Pointer to source sub-block
*
* @param[in] pu1_pred
*   Pointer
* to prediction sub-block
*
* @param[in] pi2_out
*   Pointer to residual
* sub-block
*
* @param[in] i4_src_stride
*   Source stride
*
* @param[in]
* i4_pred_stride
*   Prediction stride
*
* @param[in] dst_strd
*   Destination
* stride
*
* @param[in] u4_qbits
*    QP_BITS_h264_4x4 + floor(QP/6)
*
*
* @param[in] pu2_threshold_matrix
*   Pointer to Forward Quant Threshold
* Matrix
*
* @param[in] pu2_scale_matrix
*   Pointer to Forward Quant Scale
* Matrix
*
* @param[in] u4_round_factor
*   Quantization Round factor
*
*
* @param[out] pu1_nnz
*   Total non-zero coefficients in the current
* sub-block
*
* @returns
*
* @remarks
*
* None
*
*******************************************************************************
*/
void isvc_resi_trans_quant_4x4_sse42(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                     buffer_container_t *ps_out,
                                     buffer_container_t *ps_upsampled_res,
                                     resi_trans_quant_constants_t *ps_quant_constants,
                                     UWORD8 *pu1_nnz, WORD16 *pi2_dc_out,
                                     UWORD8 u1_use_upsampled_res)
{
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 tmp_dc, u4_zero_coeff, u4_nonzero_coeff = 0;
    WORD32 mask0, mask1;
    __m128i sum0, sum1, sum2, cmp0, cmp1;
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);
    __m128i temp_2 = _mm_set1_epi16(2);
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i temp0, temp1, temp2, temp3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i sign_reg0, sign_reg2;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i threshold_r0_r1, threshold_r2_r3;
    __m128i threshold_mask_r0_r1, threshold_mask_r2_r3;

    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_out = (WORD16 *) ps_out->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_out->i4_data_stride;

    ASSERT(0 == u1_use_upsampled_res);
    ASSERT(4 == i4_out_stride);
    UNUSED(u1_use_upsampled_res);
    UNUSED(i4_out_stride);
    UNUSED(ps_upsampled_res);

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the treshold matrix 0th,1st row */
    threshold_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the threshold matrix 2nd,3rd row */
    threshold_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix + 8));

    /* a00 a01 a02 a03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    src_r0 = _mm_loadl_epi64((__m128i *) (&pu1_src[0]));

    /* a10 a11 a12 a13 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r1 = _mm_loadl_epi64((__m128i *) (&pu1_src[i4_src_stride]));

    /* a20 a21 a22 a23 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r2 = _mm_loadl_epi64((__m128i *) (&pu1_src[2 * i4_src_stride]));

    /* a30 a31 a32 a33 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r3 = _mm_loadl_epi64((__m128i *) (&pu1_src[3 * i4_src_stride]));

    src_r0 = _mm_cvtepu8_epi16(src_r0);
    src_r1 = _mm_cvtepu8_epi16(src_r1);
    src_r2 = _mm_cvtepu8_epi16(src_r2);
    src_r3 = _mm_cvtepu8_epi16(src_r3);

    /* p00 p01 p02 p03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));

    /* p10 p11 p12 p13 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));

    /* p20 p21 p22 p23 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));

    /* p30 p31 p32 p33 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

    src_r0 = _mm_sub_epi16(src_r0, pred_r0);
    src_r1 = _mm_sub_epi16(src_r1, pred_r1);
    src_r2 = _mm_sub_epi16(src_r2, pred_r2);
    src_r3 = _mm_sub_epi16(src_r3, pred_r3);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */
    /* a0 b0 a1 b1 a2 b2 a3 b3 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* c0 d0 c1 d1 c2 d2 c3 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 b0 c0 d0 a1 b1 c1 d1 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* a2 b2 c2 d2 a3 b3 c3 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 b0 c0 d0 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* a1 b1 c1 d1 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* a2 b2 c2 d2 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* a3 b3 c3 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                             */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                             */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                             */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                             */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                                      */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                             */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                                      */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    /* a0 a1 b0 b1 c0 c1 d0 d1 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* a2 a3 b2 b3 c2 c3 d2 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 a1 a2 a3 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* b0 b1 b2 b3 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* c0 c1 c2 c3 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* d0 d1 d2 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                             */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                             */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                             */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                             */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                                      */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                             */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                                      */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    /* get the first 16 bits from the register */
    tmp_dc = _mm_extract_epi16(src_r0, 0);
    *pi2_dc_out = tmp_dc;

    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    src_r0 = _mm_unpacklo_epi64(src_r0, src_r1);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    src_r2 = _mm_unpacklo_epi64(src_r2, src_r3);
    sign_reg0 = _mm_cmpgt_epi16(zero_8x16b, src_r0);
    sign_reg2 = _mm_cmpgt_epi16(zero_8x16b, src_r2);

    sign_reg0 = _mm_mullo_epi16(temp_2, sign_reg0);
    sign_reg2 = _mm_mullo_epi16(temp_2, sign_reg2);

    sign_reg0 = _mm_add_epi16(temp_1, sign_reg0);
    sign_reg2 = _mm_add_epi16(temp_1, sign_reg2);

    src_r0 = _mm_abs_epi16(src_r0);
    src_r2 = _mm_abs_epi16(src_r2);

    threshold_mask_r0_r1 = _mm_cmpgt_epi16(threshold_r0_r1, src_r0);
    threshold_mask_r2_r3 = _mm_cmpgt_epi16(threshold_r2_r3, src_r2);

    src_r1 = _mm_srli_si128(src_r0, 8);
    src_r0 = _mm_cvtepu16_epi32(src_r0);
    src_r1 = _mm_cvtepu16_epi32(src_r1);
    src_r3 = _mm_srli_si128(src_r2, 8);
    src_r2 = _mm_cvtepu16_epi32(src_r2);
    src_r3 = _mm_cvtepu16_epi32(src_r3);

    temp0 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    scalemat_r0_r1 = _mm_srli_si128(scalemat_r0_r1, 8);
    temp2 = _mm_cvtepu16_epi32(scalemat_r2_r3);
    scalemat_r2_r3 = _mm_srli_si128(scalemat_r2_r3, 8);
    temp1 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    temp3 = _mm_cvtepu16_epi32(scalemat_r2_r3);

    temp0 = _mm_mullo_epi32(temp0, src_r0);
    temp1 = _mm_mullo_epi32(temp1, src_r1);
    temp2 = _mm_mullo_epi32(temp2, src_r2);
    temp3 = _mm_mullo_epi32(temp3, src_r3);

    temp0 = _mm_add_epi32(temp0, rnd_fact);
    temp1 = _mm_add_epi32(temp1, rnd_fact);
    temp2 = _mm_add_epi32(temp2, rnd_fact);
    temp3 = _mm_add_epi32(temp3, rnd_fact);

    temp0 = _mm_srli_epi32(temp0, u4_qbits);
    temp1 = _mm_srli_epi32(temp1, u4_qbits);
    temp2 = _mm_srli_epi32(temp2, u4_qbits);
    temp3 = _mm_srli_epi32(temp3, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);
    temp2 = _mm_packs_epi32(temp2, temp3);

    temp0 = _mm_sign_epi16(temp0, sign_reg0);
    temp2 = _mm_sign_epi16(temp2, sign_reg2);

    temp0 = _mm_andnot_si128(threshold_mask_r0_r1, temp0);
    temp2 = _mm_andnot_si128(threshold_mask_r2_r3, temp2);

    _mm_storeu_si128((__m128i *) (&pi2_out[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_out[8]), temp2);

    cmp0 = _mm_cmpeq_epi16(temp0, zero_8x16b);
    cmp1 = _mm_cmpeq_epi16(temp2, zero_8x16b);

    mask0 = _mm_movemask_epi8(cmp0);
    mask1 = _mm_movemask_epi8(cmp1);
    u4_zero_coeff = 0;

    if(mask0)
    {
        if(mask0 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp0);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }
    if(mask1)
    {
        if(mask1 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }

    /* Return total nonzero coefficients in the current sub block */
    u4_nonzero_coeff = 16 - u4_zero_coeff;
    *pu1_nnz = u4_nonzero_coeff;
}

void isvc_resi_trans_quant_4x4_with_res_pred_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_out,
    buffer_container_t *ps_upsampled_res, resi_trans_quant_constants_t *ps_quant_constants,
    UWORD8 *pu1_nnz, WORD16 *pi2_dc_out, UWORD8 u1_use_upsampled_res)
{
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 tmp_dc, u4_zero_coeff, u4_nonzero_coeff = 0;
    WORD32 mask0, mask1;
    __m128i sum0, sum1, sum2, cmp0, cmp1;
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);
    __m128i temp_2 = _mm_set1_epi16(2);
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i temp0, temp1, temp2, temp3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i sign_reg0, sign_reg2;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i upsampled_res0, upsampled_res1, upsampled_res2, upsampled_res3;
    __m128i threshold_r0_r1, threshold_r2_r3;
    __m128i threshold_mask_r0_r1, threshold_mask_r2_r3;

    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_out = (WORD16 *) ps_out->pv_data;
    WORD16 *pi2_upsampled_res = ps_upsampled_res ? (WORD16 *) ps_upsampled_res->pv_data : NULL;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_out->i4_data_stride;
    WORD32 i4_upsampled_res_stride = ps_upsampled_res ? ps_upsampled_res->i4_data_stride : 0;

    ASSERT(1 == u1_use_upsampled_res);
    ASSERT(4 == i4_out_stride);
    UNUSED(u1_use_upsampled_res);
    UNUSED(i4_out_stride);
    UNUSED(ps_upsampled_res);

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the treshold matrix 0th,1st row */
    threshold_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the threshold matrix 2nd,3rd row */
    threshold_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix + 8));

    /* a00 a01 a02 a03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    src_r0 = _mm_loadl_epi64((__m128i *) (&pu1_src[0]));

    /* a10 a11 a12 a13 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r1 = _mm_loadl_epi64((__m128i *) (&pu1_src[i4_src_stride]));

    /* a20 a21 a22 a23 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r2 = _mm_loadl_epi64((__m128i *) (&pu1_src[2 * i4_src_stride]));

    /* a30 a31 a32 a33 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r3 = _mm_loadl_epi64((__m128i *) (&pu1_src[3 * i4_src_stride]));

    src_r0 = _mm_cvtepu8_epi16(src_r0);
    src_r1 = _mm_cvtepu8_epi16(src_r1);
    src_r2 = _mm_cvtepu8_epi16(src_r2);
    src_r3 = _mm_cvtepu8_epi16(src_r3);

    /* p00 p01 p02 p03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));

    /* p10 p11 p12 p13 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));

    /* p20 p21 p22 p23 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));

    /* p30 p31 p32 p33 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_cvtepu8_epi16(pred_r0);
    pred_r1 = _mm_cvtepu8_epi16(pred_r1);
    pred_r2 = _mm_cvtepu8_epi16(pred_r2);
    pred_r3 = _mm_cvtepu8_epi16(pred_r3);

    src_r0 = _mm_sub_epi16(src_r0, pred_r0);
    src_r1 = _mm_sub_epi16(src_r1, pred_r1);
    src_r2 = _mm_sub_epi16(src_r2, pred_r2);
    src_r3 = _mm_sub_epi16(src_r3, pred_r3);

    /* load upsampled residual values and subtract from
    the previous residue */
    upsampled_res0 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[0]));

    upsampled_res1 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[i4_upsampled_res_stride]));

    upsampled_res2 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[2 * i4_upsampled_res_stride]));

    upsampled_res3 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[3 * i4_upsampled_res_stride]));

    src_r0 = _mm_sub_epi16(src_r0, upsampled_res0);
    src_r1 = _mm_sub_epi16(src_r1, upsampled_res1);
    src_r2 = _mm_sub_epi16(src_r2, upsampled_res2);
    src_r3 = _mm_sub_epi16(src_r3, upsampled_res3);

    src_r1 = _mm_unpacklo_epi16(src_r0, src_r1);
    src_r3 = _mm_unpacklo_epi16(src_r2, src_r3);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    src_r1 = _mm_max_epi16(src_r1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(src_r1, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    src_r3 = _mm_max_epi16(src_r3, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp2 = _mm_min_epi16(src_r3, pos_255_8x16b);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */
    /* a0 b0 c0 d0 a1 b1 c1 d1 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* a2 b2 c2 d2 a3 b3 c3 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 b0 c0 d0 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* a1 b1 c1 d1 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* a2 b2 c2 d2 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* a3 b3 c3 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                             */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                             */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                             */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                             */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                                      */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                             */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                                      */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    /* a0 a1 b0 b1 c0 c1 d0 d1 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* a2 a3 b2 b3 c2 c3 d2 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 a1 a2 a3 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* b0 b1 b2 b3 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* c0 c1 c2 c3 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* d0 d1 d2 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                             */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                             */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                             */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                             */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                             */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                                      */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                             */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                                      */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    /* get the first 16 bits from the register */
    tmp_dc = _mm_extract_epi16(src_r0, 0);
    *pi2_dc_out = tmp_dc;

    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    src_r0 = _mm_unpacklo_epi64(src_r0, src_r1);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    src_r2 = _mm_unpacklo_epi64(src_r2, src_r3);
    sign_reg0 = _mm_cmpgt_epi16(zero_8x16b, src_r0);
    sign_reg2 = _mm_cmpgt_epi16(zero_8x16b, src_r2);

    sign_reg0 = _mm_mullo_epi16(temp_2, sign_reg0);
    sign_reg2 = _mm_mullo_epi16(temp_2, sign_reg2);

    sign_reg0 = _mm_add_epi16(temp_1, sign_reg0);
    sign_reg2 = _mm_add_epi16(temp_1, sign_reg2);

    src_r0 = _mm_abs_epi16(src_r0);
    src_r2 = _mm_abs_epi16(src_r2);

    threshold_mask_r0_r1 = _mm_cmpgt_epi16(threshold_r0_r1, src_r0);
    threshold_mask_r2_r3 = _mm_cmpgt_epi16(threshold_r2_r3, src_r2);

    src_r1 = _mm_srli_si128(src_r0, 8);
    src_r0 = _mm_cvtepu16_epi32(src_r0);
    src_r1 = _mm_cvtepu16_epi32(src_r1);
    src_r3 = _mm_srli_si128(src_r2, 8);
    src_r2 = _mm_cvtepu16_epi32(src_r2);
    src_r3 = _mm_cvtepu16_epi32(src_r3);

    temp0 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    scalemat_r0_r1 = _mm_srli_si128(scalemat_r0_r1, 8);
    temp2 = _mm_cvtepu16_epi32(scalemat_r2_r3);
    scalemat_r2_r3 = _mm_srli_si128(scalemat_r2_r3, 8);
    temp1 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    temp3 = _mm_cvtepu16_epi32(scalemat_r2_r3);

    temp0 = _mm_mullo_epi32(temp0, src_r0);
    temp1 = _mm_mullo_epi32(temp1, src_r1);
    temp2 = _mm_mullo_epi32(temp2, src_r2);
    temp3 = _mm_mullo_epi32(temp3, src_r3);

    temp0 = _mm_add_epi32(temp0, rnd_fact);
    temp1 = _mm_add_epi32(temp1, rnd_fact);
    temp2 = _mm_add_epi32(temp2, rnd_fact);
    temp3 = _mm_add_epi32(temp3, rnd_fact);

    temp0 = _mm_srli_epi32(temp0, u4_qbits);
    temp1 = _mm_srli_epi32(temp1, u4_qbits);
    temp2 = _mm_srli_epi32(temp2, u4_qbits);
    temp3 = _mm_srli_epi32(temp3, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);
    temp2 = _mm_packs_epi32(temp2, temp3);

    temp0 = _mm_sign_epi16(temp0, sign_reg0);
    temp2 = _mm_sign_epi16(temp2, sign_reg2);

    temp0 = _mm_andnot_si128(threshold_mask_r0_r1, temp0);
    temp2 = _mm_andnot_si128(threshold_mask_r2_r3, temp2);

    _mm_storeu_si128((__m128i *) (&pi2_out[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_out[8]), temp2);

    cmp0 = _mm_cmpeq_epi16(temp0, zero_8x16b);
    cmp1 = _mm_cmpeq_epi16(temp2, zero_8x16b);

    mask0 = _mm_movemask_epi8(cmp0);
    mask1 = _mm_movemask_epi8(cmp1);
    u4_zero_coeff = 0;
    if(mask0)
    {
        if(mask0 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp0);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }
    if(mask1)
    {
        if(mask1 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }

    /* Return total nonzero coefficients in the current sub block */
    u4_nonzero_coeff = 16 - u4_zero_coeff;
    *pu1_nnz = u4_nonzero_coeff;
}

/**

 * *******************************************************************************

 * *
 * @brief
 *   This function performs forward transform and quantization on
 * a 4*4 chroma
 *block
 *
 * @par Description:
 *   The function accepts source
 * buffer and estimation buffer. From these, it
 *   computes the residue. This
 * is residue is then transformed and quantized.
 *   The transform and
 * quantization are in placed computed. They use the residue
 *   buffer for
 * this.
 *
 * @param[in] pu1_src
 *   Pointer to source sub-block
 *
 *
 * @param[in] pu1_pred
 *   Pointer to prediction sub-block
 *
 * @param[in]
 * pi2_out
 *   Pointer to residual sub-block
 *
 * @param[in] i4_src_stride
 *
 * Source stride
 *
 * @param[in] i4_pred_stride
 *   Prediction stride
 *
 *
 * @param[in] dst_strd
 *   Destination stride
 *
 * @param[in] u4_qbits
 *
 * QP_BITS_h264_4x4 + floor(QP/6)
 *
 * @param[in] pu2_threshold_matrix
 *
 * Pointer to Forward Quant Threshold Matrix
 *
 * @param[in] pu2_scale_matrix

 * *   Pointer to Forward Quant Scale Matrix
 *
 * @param[in] u4_round_factor
 *
 * Quantization Round factor
 *
 * @param[out] pu1_nnz
 *   Total non-zero
 * coefficients in the current sub-block
 *
 * @returns
 *
 * @remarks
 *
 * None
 *

 * *******************************************************************************

 */
void isvc_resi_trans_quant_chroma_4x4_sse42(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                            buffer_container_t *ps_out,
                                            buffer_container_t *ps_upsampled_res,
                                            resi_trans_quant_constants_t *ps_quant_constants,
                                            UWORD8 *pu1_nnz, WORD16 *pi2_dc_out,
                                            UWORD8 u1_use_upsampled_res)
{
    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_out = (WORD16 *) ps_out->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_out->i4_data_stride;
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 tmp_dc, u4_zero_coeff, u4_nonzero_coeff = 0;
    WORD32 mask0, mask1;
    __m128i cmp0, cmp1, sum0, sum1, sum2;
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);
    __m128i temp_2 = _mm_set1_epi16(2);
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i temp0, temp1, temp2, temp3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i sign_reg0, sign_reg2;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i threshold_r0_r1, threshold_r2_r3;
    __m128i threshold_mask_r0_r1, threshold_mask_r2_r3;
    __m128i chroma_mask = _mm_set1_epi16(0xFF);

    ASSERT(0 == u1_use_upsampled_res);
    ASSERT(4 == i4_out_stride);
    UNUSED(u1_use_upsampled_res);
    UNUSED(i4_out_stride);
    UNUSED(ps_upsampled_res);

    /* b00 b01 b02 b03 b10 b11 b12 b13
   -- the scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the treshold matrix 0th,1st row */
    threshold_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the threshold matrix 2nd,3rd row */
    threshold_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix + 8));

    /* a00 a01 a02 a03 0 0 0 0 0
    0 0 0 -- all 8 bits */
    src_r0 = _mm_loadl_epi64((__m128i *) (&pu1_src[0]));
    /* a10 a11 a12 a13 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r1 = _mm_loadl_epi64((__m128i *) (&pu1_src[i4_src_stride]));
    /* a20 a21 a22 a23 0 0 0 0 0 0 0
    0 -- all 8 bits */
    src_r2 = _mm_loadl_epi64((__m128i *) (&pu1_src[2 * i4_src_stride]));
    /* a30 a31 a32 a33 0 0 0 0 0 0 0
    0 -- all 8 bits */
    src_r3 = _mm_loadl_epi64((__m128i *) (&pu1_src[3 * i4_src_stride]));

    src_r0 = _mm_and_si128(src_r0, chroma_mask);
    src_r1 = _mm_and_si128(src_r1, chroma_mask);
    src_r2 = _mm_and_si128(src_r2, chroma_mask);
    src_r3 = _mm_and_si128(src_r3, chroma_mask);

    /* p00 p01 p02 p03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    /* p10 p11 p12 p13 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    /* p20 p21 p22 p23 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    /* p30 p31 p32 p33 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    src_r0 = _mm_sub_epi16(src_r0, pred_r0);
    src_r1 = _mm_sub_epi16(src_r1, pred_r1);
    src_r2 = _mm_sub_epi16(src_r2, pred_r2);
    src_r3 = _mm_sub_epi16(src_r3, pred_r3);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */
    /* a0 b0 a1 b1 a2 b2 a3 b3 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* c0 d0 c1 d1 c2 d2 c3 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 b0 c0 d0 a1 b1 c1 d1 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* a2 b2 c2 d2 a3 b3 c3 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 b0 c0 d0 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* a1 b1 c1 d1 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* a2 b2 c2 d2 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* a3 b3 c3 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                         */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                         */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    /* a0 a1 b0 b1 c0 c1 d0 d1 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* a2 a3 b2 b3 c2 c3 d2 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 a1 a2 a3 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* b0 b1 b2 b3 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* c0 c1 c2 c3 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* d0 d1 d2 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                  */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                  */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                  */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                  */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                  */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                           */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                  */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                           */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    /* get the first 16 bits from the register */
    tmp_dc = _mm_extract_epi16(src_r0, 0);
    *pi2_dc_out = tmp_dc;

    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    src_r0 = _mm_unpacklo_epi64(src_r0, src_r1);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    src_r2 = _mm_unpacklo_epi64(src_r2, src_r3);
    sign_reg0 = _mm_cmpgt_epi16(zero_8x16b, src_r0);
    sign_reg2 = _mm_cmpgt_epi16(zero_8x16b, src_r2);

    sign_reg0 = _mm_mullo_epi16(temp_2, sign_reg0);
    sign_reg2 = _mm_mullo_epi16(temp_2, sign_reg2);

    sign_reg0 = _mm_add_epi16(temp_1, sign_reg0);
    sign_reg2 = _mm_add_epi16(temp_1, sign_reg2);

    src_r0 = _mm_abs_epi16(src_r0);
    src_r2 = _mm_abs_epi16(src_r2);

    threshold_mask_r0_r1 = _mm_cmpgt_epi16(threshold_r0_r1, src_r0);
    threshold_mask_r2_r3 = _mm_cmpgt_epi16(threshold_r2_r3, src_r2);

    src_r1 = _mm_srli_si128(src_r0, 8);
    src_r0 = _mm_cvtepu16_epi32(src_r0);
    src_r1 = _mm_cvtepu16_epi32(src_r1);
    src_r3 = _mm_srli_si128(src_r2, 8);
    src_r2 = _mm_cvtepu16_epi32(src_r2);
    src_r3 = _mm_cvtepu16_epi32(src_r3);

    temp0 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    scalemat_r0_r1 = _mm_srli_si128(scalemat_r0_r1, 8);
    temp2 = _mm_cvtepu16_epi32(scalemat_r2_r3);
    scalemat_r2_r3 = _mm_srli_si128(scalemat_r2_r3, 8);
    temp1 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    temp3 = _mm_cvtepu16_epi32(scalemat_r2_r3);

    temp0 = _mm_mullo_epi32(temp0, src_r0);
    temp1 = _mm_mullo_epi32(temp1, src_r1);
    temp2 = _mm_mullo_epi32(temp2, src_r2);
    temp3 = _mm_mullo_epi32(temp3, src_r3);

    temp0 = _mm_add_epi32(temp0, rnd_fact);
    temp1 = _mm_add_epi32(temp1, rnd_fact);
    temp2 = _mm_add_epi32(temp2, rnd_fact);
    temp3 = _mm_add_epi32(temp3, rnd_fact);

    temp0 = _mm_srli_epi32(temp0, u4_qbits);
    temp1 = _mm_srli_epi32(temp1, u4_qbits);
    temp2 = _mm_srli_epi32(temp2, u4_qbits);
    temp3 = _mm_srli_epi32(temp3, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);
    temp2 = _mm_packs_epi32(temp2, temp3);

    temp0 = _mm_sign_epi16(temp0, sign_reg0);
    temp2 = _mm_sign_epi16(temp2, sign_reg2);

    temp0 = _mm_andnot_si128(threshold_mask_r0_r1, temp0);
    temp2 = _mm_andnot_si128(threshold_mask_r2_r3, temp2);

    _mm_storeu_si128((__m128i *) (&pi2_out[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_out[8]), temp2);

    cmp0 = _mm_cmpeq_epi16(temp0, zero_8x16b);
    cmp1 = _mm_cmpeq_epi16(temp2, zero_8x16b);

    mask0 = _mm_movemask_epi8(cmp0);
    mask1 = _mm_movemask_epi8(cmp1);
    u4_zero_coeff = 0;
    if(mask0)
    {
        if(mask0 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp0);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }
    if(mask1)
    {
        if(mask1 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }

    /* Return total nonzero coefficients in the current sub block */
    u4_nonzero_coeff = 16 - u4_zero_coeff;
    *pu1_nnz = u4_nonzero_coeff;
}

void isvc_resi_trans_quant_chroma_4x4_with_res_pred_sse42(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_out,
    buffer_container_t *ps_upsampled_res, resi_trans_quant_constants_t *ps_quant_constants,
    UWORD8 *pu1_nnz, WORD16 *pi2_dc_out, UWORD8 u1_use_upsampled_res)
{
    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_out = (WORD16 *) ps_out->pv_data;
    WORD16 *pi2_upsampled_res = ps_upsampled_res ? (WORD16 *) ps_upsampled_res->pv_data : NULL;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_out->i4_data_stride;
    WORD32 i4_upsampled_res_stride = ps_upsampled_res ? ps_upsampled_res->i4_data_stride : 0;
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 tmp_dc, u4_zero_coeff, u4_nonzero_coeff = 0;
    WORD32 mask0, mask1;
    __m128i cmp0, cmp1, sum0, sum1, sum2;
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);
    __m128i temp_2 = _mm_set1_epi16(2);
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i temp0, temp1, temp2, temp3;
    /* all bits reset to zero */
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i neg_255_8x16b = _mm_set1_epi16(-((WORD16) UINT8_MAX));
    __m128i pos_255_8x16b = _mm_set1_epi16(((WORD16) UINT8_MAX));
    __m128i sign_reg0, sign_reg2;
    __m128i scalemat_r0_r1, scalemat_r2_r3;
    __m128i upsampled_res0, upsampled_res1, upsampled_res2, upsampled_res3;
    __m128i threshold_r0_r1, threshold_r2_r3;
    __m128i threshold_mask_r0_r1, threshold_mask_r2_r3;
    __m128i chroma_mask = _mm_set1_epi16(0xFF);

    ASSERT(1 == u1_use_upsampled_res);
    ASSERT(4 == i4_out_stride);
    UNUSED(u1_use_upsampled_res);
    UNUSED(i4_out_stride);
    UNUSED(ps_upsampled_res);

    /* b00 b01 b02 b03 b10 b11 b12 b13
   -- the scaling matrix 0th,1st row */
    scalemat_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the scaling matrix 2nd,3rd row */
    scalemat_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_scale_matrix + 8));

    /* b00 b01 b02 b03 b10 b11 b12 b13
     -- the treshold matrix 0th,1st row */
    threshold_r0_r1 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix));

    /* b20 b21 b22 b23 b30 b31 b32 b33
     -- the threshold matrix 2nd,3rd row */
    threshold_r2_r3 = _mm_loadu_si128((__m128i *) (pu2_threshold_matrix + 8));

    /* a00 a01 a02 a03 0 0 0 0 0
    0 0 0 -- all 8 bits */
    src_r0 = _mm_loadl_epi64((__m128i *) (&pu1_src[0]));
    /* a10 a11 a12 a13 0 0 0 0 0 0 0
     0 -- all 8 bits */
    src_r1 = _mm_loadl_epi64((__m128i *) (&pu1_src[i4_src_stride]));
    /* a20 a21 a22 a23 0 0 0 0 0 0 0
    0 -- all 8 bits */
    src_r2 = _mm_loadl_epi64((__m128i *) (&pu1_src[2 * i4_src_stride]));
    /* a30 a31 a32 a33 0 0 0 0 0 0 0
    0 -- all 8 bits */
    src_r3 = _mm_loadl_epi64((__m128i *) (&pu1_src[3 * i4_src_stride]));

    src_r0 = _mm_and_si128(src_r0, chroma_mask);
    src_r1 = _mm_and_si128(src_r1, chroma_mask);
    src_r2 = _mm_and_si128(src_r2, chroma_mask);
    src_r3 = _mm_and_si128(src_r3, chroma_mask);

    /* p00 p01 p02 p03 0 0 0 0 0
     0 0 0 -- all 8 bits */
    pred_r0 = _mm_loadl_epi64((__m128i *) (&pu1_pred[0]));
    /* p10 p11 p12 p13 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r1 = _mm_loadl_epi64((__m128i *) (&pu1_pred[i4_pred_stride]));
    /* p20 p21 p22 p23 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r2 = _mm_loadl_epi64((__m128i *) (&pu1_pred[2 * i4_pred_stride]));
    /* p30 p31 p32 p33 0 0 0 0 0
    0 0 0 -- all 8 bits */
    pred_r3 = _mm_loadl_epi64((__m128i *) (&pu1_pred[3 * i4_pred_stride]));

    pred_r0 = _mm_and_si128(pred_r0, chroma_mask);
    pred_r1 = _mm_and_si128(pred_r1, chroma_mask);
    pred_r2 = _mm_and_si128(pred_r2, chroma_mask);
    pred_r3 = _mm_and_si128(pred_r3, chroma_mask);

    src_r0 = _mm_sub_epi16(src_r0, pred_r0);
    src_r1 = _mm_sub_epi16(src_r1, pred_r1);
    src_r2 = _mm_sub_epi16(src_r2, pred_r2);
    src_r3 = _mm_sub_epi16(src_r3, pred_r3);

    /* load upsampled residual values and subtract from
    the previous residue */
    upsampled_res0 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[0]));

    upsampled_res1 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[i4_upsampled_res_stride]));

    upsampled_res2 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[2 * i4_upsampled_res_stride]));

    upsampled_res3 = _mm_loadu_si128((__m128i *) (&pi2_upsampled_res[3 * i4_upsampled_res_stride]));

    src_r0 = _mm_sub_epi16(src_r0, upsampled_res0);
    src_r1 = _mm_sub_epi16(src_r1, upsampled_res1);
    src_r2 = _mm_sub_epi16(src_r2, upsampled_res2);
    src_r3 = _mm_sub_epi16(src_r3, upsampled_res3);

    src_r1 = _mm_unpacklo_epi16(src_r0, src_r1);
    src_r3 = _mm_unpacklo_epi16(src_r2, src_r3);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    src_r1 = _mm_max_epi16(src_r1, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp0 = _mm_min_epi16(src_r1, pos_255_8x16b);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    src_r3 = _mm_max_epi16(src_r3, neg_255_8x16b);
    /* Saturate all values > 255 to 255 and retain the rest as it is */
    temp2 = _mm_min_epi16(src_r3, pos_255_8x16b);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* a2 b2 c2 d2 a3 b3 c3 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 b0 c0 d0 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* a1 b1 c1 d1 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* a2 b2 c2 d2 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* a3 b3 c3 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                         */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                         */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    /* a0 a1 b0 b1 c0 c1 d0 d1 */
    temp0 = _mm_unpacklo_epi16(src_r0, src_r1);
    /* a2 a3 b2 b3 c2 c3 d2 d3 */
    temp2 = _mm_unpacklo_epi16(src_r2, src_r3);
    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    temp1 = _mm_unpacklo_epi32(temp0, temp2);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    temp3 = _mm_unpackhi_epi32(temp0, temp2);

    /* a0 a1 a2 a3 */
    src_r0 = _mm_unpacklo_epi64(temp1, zero_8x16b);
    /* b0 b1 b2 b3 */
    src_r1 = _mm_unpackhi_epi64(temp1, zero_8x16b);
    /* c0 c1 c2 c3 */
    src_r2 = _mm_unpacklo_epi64(temp3, zero_8x16b);
    /* d0 d1 d2 d3 */
    src_r3 = _mm_unpackhi_epi64(temp3, zero_8x16b);

    /*----------------------------------------------------------*/
    /* x0 = z0 + z3                                  */
    temp0 = _mm_add_epi16(src_r0, src_r3);
    /* x1 = z1 + z2                                  */
    temp1 = _mm_add_epi16(src_r1, src_r2);
    /* x2 = z1 - z2                                  */
    temp2 = _mm_sub_epi16(src_r1, src_r2);
    /* x3 = z0 - z3                                  */
    temp3 = _mm_sub_epi16(src_r0, src_r3);

    /* z0 = x0 + x1                                  */
    src_r0 = _mm_add_epi16(temp0, temp1);
    /* z1 = (x3 << 1) + x2                           */
    src_r1 = _mm_slli_epi16(temp3, 1);
    src_r1 = _mm_add_epi16(src_r1, temp2);
    /* z2 = x0 - x1                                  */
    src_r2 = _mm_sub_epi16(temp0, temp1);
    /* z3 = x3 - (x2 << 1)                           */
    src_r3 = _mm_slli_epi16(temp2, 1);
    src_r3 = _mm_sub_epi16(temp3, src_r3);

    /* get the first 16 bits from the register */
    tmp_dc = _mm_extract_epi16(src_r0, 0);
    *pi2_dc_out = tmp_dc;

    /* a0 a1 a2 a3 b0 b1 b2 b3 */
    src_r0 = _mm_unpacklo_epi64(src_r0, src_r1);
    /* c0 c1 c2 c3 d0 d1 d2 d3 */
    src_r2 = _mm_unpacklo_epi64(src_r2, src_r3);
    sign_reg0 = _mm_cmpgt_epi16(zero_8x16b, src_r0);
    sign_reg2 = _mm_cmpgt_epi16(zero_8x16b, src_r2);

    sign_reg0 = _mm_mullo_epi16(temp_2, sign_reg0);
    sign_reg2 = _mm_mullo_epi16(temp_2, sign_reg2);

    sign_reg0 = _mm_add_epi16(temp_1, sign_reg0);
    sign_reg2 = _mm_add_epi16(temp_1, sign_reg2);

    src_r0 = _mm_abs_epi16(src_r0);
    src_r2 = _mm_abs_epi16(src_r2);

    threshold_mask_r0_r1 = _mm_cmpgt_epi16(threshold_r0_r1, src_r0);
    threshold_mask_r2_r3 = _mm_cmpgt_epi16(threshold_r2_r3, src_r2);

    src_r1 = _mm_srli_si128(src_r0, 8);
    src_r0 = _mm_cvtepu16_epi32(src_r0);
    src_r1 = _mm_cvtepu16_epi32(src_r1);
    src_r3 = _mm_srli_si128(src_r2, 8);
    src_r2 = _mm_cvtepu16_epi32(src_r2);
    src_r3 = _mm_cvtepu16_epi32(src_r3);

    temp0 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    scalemat_r0_r1 = _mm_srli_si128(scalemat_r0_r1, 8);
    temp2 = _mm_cvtepu16_epi32(scalemat_r2_r3);
    scalemat_r2_r3 = _mm_srli_si128(scalemat_r2_r3, 8);
    temp1 = _mm_cvtepu16_epi32(scalemat_r0_r1);
    temp3 = _mm_cvtepu16_epi32(scalemat_r2_r3);

    temp0 = _mm_mullo_epi32(temp0, src_r0);
    temp1 = _mm_mullo_epi32(temp1, src_r1);
    temp2 = _mm_mullo_epi32(temp2, src_r2);
    temp3 = _mm_mullo_epi32(temp3, src_r3);

    temp0 = _mm_add_epi32(temp0, rnd_fact);
    temp1 = _mm_add_epi32(temp1, rnd_fact);
    temp2 = _mm_add_epi32(temp2, rnd_fact);
    temp3 = _mm_add_epi32(temp3, rnd_fact);

    temp0 = _mm_srli_epi32(temp0, u4_qbits);
    temp1 = _mm_srli_epi32(temp1, u4_qbits);
    temp2 = _mm_srli_epi32(temp2, u4_qbits);
    temp3 = _mm_srli_epi32(temp3, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);
    temp2 = _mm_packs_epi32(temp2, temp3);

    temp0 = _mm_sign_epi16(temp0, sign_reg0);
    temp2 = _mm_sign_epi16(temp2, sign_reg2);

    temp0 = _mm_andnot_si128(threshold_mask_r0_r1, temp0);
    temp2 = _mm_andnot_si128(threshold_mask_r2_r3, temp2);

    _mm_storeu_si128((__m128i *) (&pi2_out[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_out[8]), temp2);

    cmp0 = _mm_cmpeq_epi16(temp0, zero_8x16b);
    cmp1 = _mm_cmpeq_epi16(temp2, zero_8x16b);

    mask0 = _mm_movemask_epi8(cmp0);
    mask1 = _mm_movemask_epi8(cmp1);
    u4_zero_coeff = 0;
    if(mask0)
    {
        if(mask0 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp0);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }
    if(mask1)
    {
        if(mask1 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }

    /* Return total nonzero coefficients in the current sub block */
    u4_nonzero_coeff = 16 - u4_zero_coeff;
    *pu1_nnz = u4_nonzero_coeff;
}

/**

 * *******************************************************************************

 * *
 * @brief
 *   This function performs forward hadamard transform and
 * quantization on a 4*4
 *block
 *
 * @par Description:
 *   The function
 * accepts source buffer and estimation buffer. From these, it
 *   computes the
 * residue. This is residue is then transformed and quantized.
 *   The
 * transform and quantization are in placed computed. They use the residue
 *
 * buffer for this.
 *
 * @param[in] pu1_src
 *   Pointer to source sub-block

 * *
 * @param[in] pu1_pred
 *   Pointer to prediction sub-block
 *
 *
 * @param[in] pi2_out
 *   Pointer to residual sub-block
 *
 * @param[in]
 * i4_src_stride
 *   Source stride
 *
 * @param[in] i4_pred_stride
 *
 * Prediction stride
 *
 * @param[in] dst_strd
 *   Destination stride
 *
 *
 * @param[in] u4_qbits
 *    QP_BITS_h264_4x4 + floor(QP/6)
 *
 * @param[in]
 * pu2_threshold_matrix
 *   Pointer to Forward Quant Threshold Matrix
 *
 *
 * @param[in] pu2_scale_matrix
 *   Pointer to Forward Quant Scale Matrix
 *
 *
 * @param[in] u4_round_factor
 *   Quantization Round factor
 *
 * @param[out]
 * pu1_nnz
 *   Total non-zero coefficients in the current sub-block
 *
 *
 * @returns
 *
 * @remarks
 *   None
 *
 */

void isvc_hadamard_quant_4x4_sse42(WORD16 *pi2_src, WORD16 *pi2_dst,
                                   resi_trans_quant_constants_t *ps_quant_constants,
                                   UWORD8 *pu1_nnz)
{
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 u4_zero_coeff, u4_nonzero_coeff = 0;
    __m128i cmp0, cmp1, sum0, sum1, sum2;
    WORD32 mask0, mask1;
    __m128i src_r0_r1, src_r2_r3, sign_reg;
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i temp0, temp1, temp2, temp3;
    __m128i sign_reg0, sign_reg1, sign_reg2, sign_reg3;
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);
    __m128i scale_val = _mm_set1_epi32(pu2_scale_matrix[0]);

    UNUSED(pu2_threshold_matrix);

    src_r0_r1 = _mm_loadu_si128((__m128i *) (pi2_src));  // a00 a01 a02 a03 a10 a11 a12 a13 -- the
                                                         // source matrix 0th,1st row
    src_r2_r3 = _mm_loadu_si128((__m128i *) (pi2_src + 8));  // a20 a21 a22 a23 a30 a31 a32 a33 --
                                                             // the source matrix 2nd,3rd row
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, src_r0_r1);
    src_r0 = _mm_unpacklo_epi16(src_r0_r1, sign_reg);  // a0 a1 a2 a3
    src_r1 = _mm_unpackhi_epi16(src_r0_r1, sign_reg);  // b0 b1 b2 b3
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, src_r2_r3);
    src_r2 = _mm_unpacklo_epi16(src_r2_r3, sign_reg);  // c0 c1 c2 c3
    src_r3 = _mm_unpackhi_epi16(src_r2_r3, sign_reg);  // d0 d1 d2 d3

    /* Perform Inverse transform */
    /*-------------------------------------------------------------*/
    /* Forward DC transform [ Horizontal transformation ] */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */
    temp0 = _mm_unpacklo_epi32(src_r0, src_r1);  // a0 b0 a1 b1
    temp2 = _mm_unpacklo_epi32(src_r2, src_r3);  // c0 d0 c1 d1
    temp1 = _mm_unpackhi_epi32(src_r0, src_r1);  // a2 b2 a3 b3
    temp3 = _mm_unpackhi_epi32(src_r2, src_r3);  // c2 d2 c3 d3
    src_r0 = _mm_unpacklo_epi64(temp0, temp2);   // a0 b0 c0 d0
    src_r1 = _mm_unpackhi_epi64(temp0, temp2);   // a1 b1 c1 d1
    src_r2 = _mm_unpacklo_epi64(temp1, temp3);   // a2 b2 c2 d2
    src_r3 = _mm_unpackhi_epi64(temp1, temp3);   // a3 b3 c3 d3

    temp0 = _mm_add_epi32(src_r0, src_r3);
    temp1 = _mm_add_epi32(src_r1, src_r2);
    temp2 = _mm_sub_epi32(src_r1, src_r2);
    temp3 = _mm_sub_epi32(src_r0, src_r3);

    src_r0 = _mm_add_epi32(temp0, temp1);
    src_r1 = _mm_add_epi32(temp2, temp3);
    src_r2 = _mm_sub_epi32(temp0, temp1);
    src_r3 = _mm_sub_epi32(temp3, temp2);

    /*-------------------------------------------------------------*/
    /* Forward DC transform [ Vertical transformation ] */
    /*-------------------------------------------------------------*/
    // Matrix transpose
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */
    temp0 = _mm_unpacklo_epi32(src_r0, src_r1);  // a0 a1 b0 b1
    temp2 = _mm_unpacklo_epi32(src_r2, src_r3);  // a2 a3 b2 b3
    temp1 = _mm_unpackhi_epi32(src_r0, src_r1);  // c0 c1 d0 d1
    temp3 = _mm_unpackhi_epi32(src_r2, src_r3);  // c2 c3 d2 d3
    src_r0 = _mm_unpacklo_epi64(temp0, temp2);   // a0 a1 a2 a3
    src_r1 = _mm_unpackhi_epi64(temp0, temp2);   // b0 b1 b2 b3
    src_r2 = _mm_unpacklo_epi64(temp1, temp3);   // c0 c1 c2 c3
    src_r3 = _mm_unpackhi_epi64(temp1, temp3);   // d0 d1 d2 d3

    temp0 = _mm_add_epi32(src_r0, src_r3);
    temp1 = _mm_add_epi32(src_r1, src_r2);
    temp2 = _mm_sub_epi32(src_r1, src_r2);
    temp3 = _mm_sub_epi32(src_r0, src_r3);

    src_r0 = _mm_add_epi32(temp0, temp1);
    src_r1 = _mm_add_epi32(temp2, temp3);
    src_r2 = _mm_sub_epi32(temp0, temp1);
    src_r3 = _mm_sub_epi32(temp3, temp2);

    src_r0 = _mm_srai_epi32(src_r0, 1);
    src_r1 = _mm_srai_epi32(src_r1, 1);
    src_r2 = _mm_srai_epi32(src_r2, 1);
    src_r3 = _mm_srai_epi32(src_r3, 1);

    // Quantization
    sign_reg0 =
        _mm_cmpgt_epi32(zero_8x16b, src_r0);  // Find sign of each value for later restoration
    sign_reg1 = _mm_cmpgt_epi32(zero_8x16b, src_r1);
    sign_reg2 = _mm_cmpgt_epi32(zero_8x16b, src_r2);
    sign_reg3 = _mm_cmpgt_epi32(zero_8x16b, src_r3);

    sign_reg0 = _mm_packs_epi32(sign_reg0,
                                sign_reg1);  // Sign = -1 or 0 depending on <0 or >0 respectively
    sign_reg2 = _mm_packs_epi32(sign_reg2, sign_reg3);

    sign_reg0 = _mm_slli_epi16(sign_reg0, 1);  // Sign = -2 or 0 depending on <0 or >0 respectively
    sign_reg2 = _mm_slli_epi16(sign_reg2, 1);

    sign_reg0 =
        _mm_add_epi16(temp_1, sign_reg0);  // Sign = -1 or 1 depending on <0 or >0 respectively
    sign_reg2 = _mm_add_epi16(temp_1, sign_reg2);

    src_r0 = _mm_abs_epi32(src_r0);  // Absolute values
    src_r1 = _mm_abs_epi32(src_r1);
    src_r2 = _mm_abs_epi32(src_r2);
    src_r3 = _mm_abs_epi32(src_r3);

    temp0 = _mm_mullo_epi32(scale_val, src_r0);  // multiply by
                                                 // pu2_scale_matrix[0]
    temp1 = _mm_mullo_epi32(scale_val, src_r1);
    temp2 = _mm_mullo_epi32(scale_val, src_r2);
    temp3 = _mm_mullo_epi32(scale_val, src_r3);

    temp0 = _mm_add_epi32(temp0, rnd_fact);  // Add round factor
    temp1 = _mm_add_epi32(temp1, rnd_fact);
    temp2 = _mm_add_epi32(temp2, rnd_fact);
    temp3 = _mm_add_epi32(temp3, rnd_fact);

    temp0 = _mm_srli_epi32(temp0,
                           u4_qbits);  // RIght shift by qbits, unsigned variable,
                                       // so shift right immediate works
    temp1 = _mm_srli_epi32(temp1, u4_qbits);
    temp2 = _mm_srli_epi32(temp2, u4_qbits);
    temp3 = _mm_srli_epi32(temp3, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);  // Final values are 16-bits only.
    temp2 = _mm_packs_epi32(temp2, temp3);

    temp0 = _mm_sign_epi16(temp0, sign_reg0);  // Sign restoration
    temp2 = _mm_sign_epi16(temp2, sign_reg2);

    _mm_storeu_si128((__m128i *) (&pi2_dst[0]), temp0);
    _mm_storeu_si128((__m128i *) (&pi2_dst[8]), temp2);

    cmp0 = _mm_cmpeq_epi16(temp0, zero_8x16b);
    cmp1 = _mm_cmpeq_epi16(temp2, zero_8x16b);

    mask0 = _mm_movemask_epi8(cmp0);
    mask1 = _mm_movemask_epi8(cmp1);
    u4_zero_coeff = 0;
    if(mask0)
    {
        if(mask0 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp0);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }
    if(mask1)
    {
        if(mask1 == 0xffff)
            u4_zero_coeff += 8;
        else
        {
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            sum2 = _mm_hadd_epi16(sum1, zero_8x16b);
            u4_zero_coeff += _mm_cvtsi128_si32(sum2);
        }
    }

    /* Return total nonzero coefficients in the current sub block */
    u4_nonzero_coeff = 16 - u4_zero_coeff;
    pu1_nnz[0] = u4_nonzero_coeff;
}

/**

 * *******************************************************************************

 * *
 * @brief
 *   This function performs forward hadamard transform and
 * quantization on a 2*2
 *block for both U and V planes
 *
 * @par
 * Description:
 *   The function accepts source buffer and estimation buffer.
 * From these, it
 *   computes the residue. This is residue is then transformed
 * and quantized.
 *   The transform and quantization are in placed computed.
 * They use the residue
 *   buffer for this.
 *
 * @param[in] pu1_src
 *
 * Pointer to source sub-block
 *
 * @param[in] pu1_pred
 *   Pointer to
 * prediction sub-block
 *
 * @param[in] pi2_out
 *   Pointer to residual
 * sub-block
 *
 * @param[in] i4_src_stride
 *   Source stride
 *
 * @param[in]
 * i4_pred_stride
 *   Prediction stride
 *
 * @param[in] dst_strd
 *
 * Destination stride
 *
 * @param[in] u4_qbits
 *    QP_BITS_h264_4x4 +
 * floor(QP/6)
 *
 * @param[in] pu2_threshold_matrix
 *   Pointer to Forward
 * Quant Threshold Matrix
 *
 * @param[in] pu2_scale_matrix
 *   Pointer to
 * Forward Quant Scale Matrix
 *
 * @param[in] u4_round_factor
 *   Quantization
 * Round factor
 *
 * @param[out] pu1_nnz
 *   Total non-zero coefficients in
 * the current sub-block
 *
 * @returns
 *
 * @remarks
 *   NNZ for dc is
 * populated at 0 and 5th position of pu1_nnz
 *
 */

void isvc_hadamard_quant_2x2_uv_sse42(WORD16 *pi2_src, WORD16 *pi2_dst,
                                      resi_trans_quant_constants_t *ps_quant_constants,
                                      UWORD8 *pu1_nnz)
{
    const UWORD16 *pu2_scale_matrix = ps_quant_constants->pu2_scale_matrix;
    const UWORD16 *pu2_threshold_matrix = ps_quant_constants->pu2_threshold_matrix;
    UWORD32 u4_qbits = ps_quant_constants->u4_qbits;
    UWORD32 u4_round_factor = ps_quant_constants->u4_round_factor;
    WORD32 val, nonzero_coeff_0 = 0, nonzero_coeff_1 = 0;
    __m128i cmp, cmp0, cmp1;
    __m128i sum0, sum1;
    WORD32 mask, mask0, mask1;
    __m128i src, plane_0, plane_1, temp0, temp1, sign_reg;
    __m128i zero_8x16b = _mm_setzero_si128();
    __m128i scale_val = _mm_set1_epi32(pu2_scale_matrix[0]);
    __m128i sign_reg0, sign_reg1;
    __m128i temp_1 = _mm_set1_epi16(1);
    __m128i rnd_fact = _mm_set1_epi32(u4_round_factor);

    UNUSED(pu2_threshold_matrix);

    src = _mm_loadu_si128((__m128i *) pi2_src);  // a0 a1 a2 a3 b0 b1 b2 b3
    sign_reg = _mm_cmpgt_epi16(zero_8x16b, src);
    plane_0 = _mm_unpacklo_epi16(src, sign_reg);  // a0 a1 a2 a3 -- 32 bits
    plane_1 = _mm_unpackhi_epi16(src, sign_reg);  // b0 b1 b2 b3 -- 32 bits

    temp0 = _mm_hadd_epi32(plane_0, plane_1);  // a0+a1 a2+a3 b0+b1 b2+b3
    temp1 = _mm_hsub_epi32(plane_0, plane_1);  // a0-a1 a2-a3 b0-b1 b2-b3

    plane_0 = _mm_hadd_epi32(temp0, temp1);  // a0+a1+a2+a3 b0+b1+b2+b3 a0-a1+a2-a3 b0-b1+b2-b3
    plane_1 = _mm_hsub_epi32(temp0, temp1);  // a0+a1-a2-a3 b0+b1-b2-b3 a0-a1-a2+a3 b0-b1-b2+b3

    temp0 =
        _mm_unpacklo_epi32(plane_0, plane_1);  // a0+a1+a2+a3 a0+a1-a2-a3 b0+b1+b2+b3 b0+b1-b2-b3
    temp1 =
        _mm_unpackhi_epi32(plane_0, plane_1);  // a0-a1+a2-a3 a0-a1-a2+a3 b0-b1+b2-b3 b0-b1-b2+b3

    plane_0 = _mm_unpacklo_epi64(temp0, temp1);  // a0+a1+a2+a3 a0+a1-a2-a3 a0-a1+a2-a3 a0-a1-a2+a3
    plane_1 = _mm_unpackhi_epi64(temp0, temp1);  // b0+b1+b2+b3 b0+b1-b2-b3 b0-b1+b2-b3 b0-b1-b2+b3

    plane_0 = _mm_shuffle_epi32(plane_0, 0xd8);  // a0+a1+a2+a3 a0-a1+a2-a3 a0+a1-a2-a3 a0-a1-a2+a3
    plane_1 = _mm_shuffle_epi32(plane_1, 0xd8);  // b0+b1+b2+b3 b0-b1+b2-b3 b0+b1-b2-b3 b0-b1-b2+b3
    // Quantization
    sign_reg0 =
        _mm_cmpgt_epi32(zero_8x16b, plane_0);  // Find sign of each value for later restoration
    sign_reg1 = _mm_cmpgt_epi32(zero_8x16b, plane_1);

    sign_reg0 = _mm_packs_epi32(sign_reg0,
                                sign_reg1);    // Sign = -1 or 0 depending on <0 or >0 respectively
    sign_reg0 = _mm_slli_epi16(sign_reg0, 1);  // Sign = -2 or 0 depending on <0 or >0 respectively
    sign_reg0 =
        _mm_add_epi16(temp_1, sign_reg0);  // Sign = -1 or 1 depending on <0 or >0 respectively

    plane_0 = _mm_abs_epi32(plane_0);  // Absolute values
    plane_1 = _mm_abs_epi32(plane_1);

    temp0 = _mm_mullo_epi32(scale_val, plane_0);  // multiply by pu2_scale_matrix[0]
    temp1 = _mm_mullo_epi32(scale_val, plane_1);  // multiply by pu2_scale_matrix[0]

    temp0 = _mm_add_epi32(temp0, rnd_fact);  // Add round factor
    temp1 = _mm_add_epi32(temp1, rnd_fact);

    temp0 = _mm_srli_epi32(temp0,
                           u4_qbits);  // RIght shift by qbits, unsigned variable,
                                       // so shift right immediate works
    temp1 = _mm_srli_epi32(temp1, u4_qbits);

    temp0 = _mm_packs_epi32(temp0, temp1);     // Final values are 16-bits only.
    temp0 = _mm_sign_epi16(temp0, sign_reg0);  // Sign restoration

    _mm_storeu_si128((__m128i *) (&pi2_dst[0]), temp0);

    cmp = _mm_cmpeq_epi16(temp0, zero_8x16b);
    mask = _mm_movemask_epi8(cmp);
    mask0 = mask & 0xff;
    mask1 = mask >> 8;
    if(mask0)
    {
        if(mask0 == 0xff)
            nonzero_coeff_0 += 4;
        else
        {
            cmp0 = _mm_and_si128(temp_1, cmp);
            sum0 = _mm_hadd_epi16(cmp0, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            val = _mm_cvtsi128_si32(sum1);
            val = val & 0xffff;
            nonzero_coeff_0 += val;
        }
    }
    if(mask1)
    {
        if(mask1 == 0xff)
            nonzero_coeff_1 += 4;
        else
        {
            cmp1 = _mm_srli_si128(cmp, 8);
            cmp1 = _mm_and_si128(temp_1, cmp1);
            sum0 = _mm_hadd_epi16(cmp1, zero_8x16b);
            sum1 = _mm_hadd_epi16(sum0, zero_8x16b);
            nonzero_coeff_1 += _mm_cvtsi128_si32(sum1);
        }
    }

    pu1_nnz[0] = 4 - nonzero_coeff_0;
    pu1_nnz[1] = 4 - nonzero_coeff_1;
}
