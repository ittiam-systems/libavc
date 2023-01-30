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
 *  isvc_resi_trans_quant_neon.c
 *
 * @brief
 *  neon variants of forward transform and quantization functions
 *
 * *******************************************************************************
 */

#include <arm_neon.h>
#include <string.h>

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

void isvc_resi_trans_quant_4x4_neon(buffer_container_t *ps_src, buffer_container_t *ps_pred,
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

    uint8x8_t src0, src1, src2, src3;
    uint8x8_t pred0, pred1, pred2, pred3;
    uint8x8_t temp0_u8x8, temp1_u8x8;
    uint16x4_t temp0_u16x4, temp1_u16x4, temp2_u16x4, temp3_u16x4;
    uint16x4_t scale_mat0_16x4, scale_mat1_16x4, scale_mat2_16x4, scale_mat3_16x4;
    uint16x4_t threshold0_16x4, threshold1_16x4, threshold2_16x4, threshold3_16x4;
    uint16x4_t thresholdmask0_16x4, thresholdmask1_16x4, thresholdmask2_16x4, thresholdmask3_16x4;
    int16x4_t res0_16x4, res1_16x4, res2_16x4, res3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int16x4_t temp0_16x4, temp1_16x4, temp2_16x4, temp3_16x4;
    uint16x8_t res0_16x8, res1_16x8, res2_16x8, res3_16x8;
    uint16x8_t temp0_u16x8, temp1_u16x8;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int32x4_t tx0_32x4, tx1_32x4, tx2_32x4, tx3_32x4;

    int32x4_t rnd_factor_32x4 = vdupq_n_s32(u4_round_factor);
    int32x4_t qbits_32x4 = vdupq_n_s32(u4_qbits);
    int16x4_t zeros_16x4 = vdup_n_s16(0);

    UNUSED(ps_upsampled_res);
    UNUSED(u1_use_upsampled_res);

    threshold0_16x4 = vld1_u16(pu2_threshold_matrix);
    threshold1_16x4 = vld1_u16(pu2_threshold_matrix + 4);
    threshold2_16x4 = vld1_u16(pu2_threshold_matrix + 8);
    threshold3_16x4 = vld1_u16(pu2_threshold_matrix + 12);

    scale_mat0_16x4 = vld1_u16(pu2_scale_matrix);
    scale_mat1_16x4 = vld1_u16(pu2_scale_matrix + 4);
    scale_mat2_16x4 = vld1_u16(pu2_scale_matrix + 8);
    scale_mat3_16x4 = vld1_u16(pu2_scale_matrix + 12);

    src0 = vld1_u8(&pu1_src[0 * i4_src_stride]);
    src1 = vld1_u8(&pu1_src[1 * i4_src_stride]);
    src2 = vld1_u8(&pu1_src[2 * i4_src_stride]);
    src3 = vld1_u8(&pu1_src[3 * i4_src_stride]);

    pred0 = vld1_u8(&pu1_pred[0 * i4_pred_stride]);
    pred1 = vld1_u8(&pu1_pred[1 * i4_pred_stride]);
    pred2 = vld1_u8(&pu1_pred[2 * i4_pred_stride]);
    pred3 = vld1_u8(&pu1_pred[3 * i4_pred_stride]);

    /* calculate res = src - pred */
    res0_16x8 = vsubl_u8(src0, pred0);
    res1_16x8 = vsubl_u8(src1, pred1);
    res2_16x8 = vsubl_u8(src2, pred2);
    res3_16x8 = vsubl_u8(src3, pred3);

    res0_16x4 = vreinterpret_s16_u16(vget_low_u16(res0_16x8));
    res1_16x4 = vreinterpret_s16_u16(vget_low_u16(res1_16x8));
    res2_16x4 = vreinterpret_s16_u16(vget_low_u16(res2_16x8));
    res3_16x4 = vreinterpret_s16_u16(vget_low_u16(res3_16x8));

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    /* Matrix transpose */
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */

    xx0_16x4x2 = vtrn_s16(res0_16x4, res1_16x4);
    xx1_16x4x2 = vtrn_s16(res2_16x4, res3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(xx2_16x4, temp0_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* Matrix transpose */
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */

    xx0_16x4x2 = vtrn_s16(x0_16x4, x1_16x4);
    xx1_16x4x2 = vtrn_s16(x2_16x4, x3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Vertical Transformation */

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(temp0_16x4, xx2_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* get the first 16 bits from the register */
    *pi2_dc_out = vget_lane_s16(x0_16x4, 0);

    xx0_16x4 = vabs_s16(x0_16x4);
    xx1_16x4 = vabs_s16(x1_16x4);
    xx2_16x4 = vabs_s16(x2_16x4);
    xx3_16x4 = vabs_s16(x3_16x4);

    /* compare with zero for getting sign */
    temp0_u16x4 = vcgt_s16(x0_16x4, zeros_16x4);
    temp1_u16x4 = vcgt_s16(x1_16x4, zeros_16x4);
    temp2_u16x4 = vcgt_s16(x2_16x4, zeros_16x4);
    temp3_u16x4 = vcgt_s16(x3_16x4, zeros_16x4);

    /* compare with zero for thresholding */
    thresholdmask0_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold0_16x4), xx0_16x4);
    thresholdmask1_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold1_16x4), xx1_16x4);
    thresholdmask2_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold2_16x4), xx2_16x4);
    thresholdmask3_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold3_16x4), xx3_16x4);

    /* Multiply abs values obtained with scaling matrix */
    tx0_32x4 = vmull_s16(xx0_16x4, vreinterpret_s16_u16(scale_mat0_16x4));
    tx1_32x4 = vmull_s16(xx1_16x4, vreinterpret_s16_u16(scale_mat1_16x4));
    tx2_32x4 = vmull_s16(xx2_16x4, vreinterpret_s16_u16(scale_mat2_16x4));
    tx3_32x4 = vmull_s16(xx3_16x4, vreinterpret_s16_u16(scale_mat3_16x4));

    tx0_32x4 = vaddq_s32(tx0_32x4, rnd_factor_32x4);
    tx1_32x4 = vaddq_s32(tx1_32x4, rnd_factor_32x4);
    tx2_32x4 = vaddq_s32(tx2_32x4, rnd_factor_32x4);
    tx3_32x4 = vaddq_s32(tx3_32x4, rnd_factor_32x4);

    qbits_32x4 = vnegq_s32(qbits_32x4);

    tx0_32x4 = vshlq_s32(tx0_32x4, qbits_32x4);
    tx1_32x4 = vshlq_s32(tx1_32x4, qbits_32x4);
    tx2_32x4 = vshlq_s32(tx2_32x4, qbits_32x4);
    tx3_32x4 = vshlq_s32(tx3_32x4, qbits_32x4);

    /* Convertion to 16 bits signed */
    temp0_16x4 = vmovn_s32(tx0_32x4);
    temp1_16x4 = vmovn_s32(tx1_32x4);
    temp2_16x4 = vmovn_s32(tx2_32x4);
    temp3_16x4 = vmovn_s32(tx3_32x4);

    x0_16x4 = vneg_s16(temp0_16x4);
    x1_16x4 = vneg_s16(temp1_16x4);
    x2_16x4 = vneg_s16(temp2_16x4);
    x3_16x4 = vneg_s16(temp3_16x4);

    /* Restore sign */
    x0_16x4 = vbsl_s16(temp0_u16x4, temp0_16x4, x0_16x4);
    x1_16x4 = vbsl_s16(temp1_u16x4, temp1_16x4, x1_16x4);
    x2_16x4 = vbsl_s16(temp2_u16x4, temp2_16x4, x2_16x4);
    x3_16x4 = vbsl_s16(temp3_u16x4, temp3_16x4, x3_16x4);

    xx0_16x4 = vbsl_s16(thresholdmask0_16x4, zeros_16x4, x0_16x4);
    xx1_16x4 = vbsl_s16(thresholdmask1_16x4, zeros_16x4, x1_16x4);
    xx2_16x4 = vbsl_s16(thresholdmask2_16x4, zeros_16x4, x2_16x4);
    xx3_16x4 = vbsl_s16(thresholdmask3_16x4, zeros_16x4, x3_16x4);

    /* Store Quantized outputs */
    vst1_s16(&pi2_out[0 * i4_out_stride], xx0_16x4);
    vst1_s16(&pi2_out[1 * i4_out_stride], xx1_16x4);
    vst1_s16(&pi2_out[2 * i4_out_stride], xx2_16x4);
    vst1_s16(&pi2_out[3 * i4_out_stride], xx3_16x4);

    /* NNZ calculation */

    temp0_u16x4 = vceq_s16(xx0_16x4, zeros_16x4);
    temp1_u16x4 = vceq_s16(xx1_16x4, zeros_16x4);
    temp2_u16x4 = vceq_s16(xx2_16x4, zeros_16x4);
    temp3_u16x4 = vceq_s16(xx3_16x4, zeros_16x4);

    temp0_u16x8 = vcombine_u16(temp0_u16x4, temp2_u16x4);
    temp1_u16x8 = vcombine_u16(temp1_u16x4, temp3_u16x4);

    /* Convertion to 8 bit unsigned */
    temp0_u8x8 = vmovn_u16(temp0_u16x8);
    temp1_u8x8 = vmovn_u16(temp1_u16x8);

    temp0_u8x8 = vshr_n_u8(temp0_u8x8, 7);
    temp1_u8x8 = vshr_n_u8(temp1_u8x8, 7);

    temp0_u8x8 = vadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);

    *pu1_nnz = 16 - vget_lane_u8(temp0_u8x8, 0);
}

void isvc_resi_trans_quant_4x4_with_residual_sub_neon(
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

    uint8x8_t src0, src1, src2, src3;
    uint8x8_t pred0, pred1, pred2, pred3;
    uint8x8_t temp0_u8x8, temp1_u8x8;
    uint16x4_t temp0_u16x4, temp1_u16x4, temp2_u16x4, temp3_u16x4;
    uint16x4_t scale_mat0_16x4, scale_mat1_16x4, scale_mat2_16x4, scale_mat3_16x4;
    uint16x4_t threshold0_16x4, threshold1_16x4, threshold2_16x4, threshold3_16x4;
    uint16x4_t thresholdmask0_16x4, thresholdmask1_16x4, thresholdmask2_16x4, thresholdmask3_16x4;
    int16x4_t upres0_16x4, upres1_16x4, upres2_16x4, upres3_16x4;
    int16x4_t res0_16x4, res1_16x4, res2_16x4, res3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int16x4_t temp0_16x4, temp1_16x4, temp2_16x4, temp3_16x4;
    uint16x8_t res0_16x8, res1_16x8, res2_16x8, res3_16x8;
    uint16x8_t temp0_u16x8, temp1_u16x8;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int32x4_t tx0_32x4, tx1_32x4, tx2_32x4, tx3_32x4;

    int32x4_t rnd_factor_32x4 = vdupq_n_s32(u4_round_factor);
    int32x4_t qbits_32x4 = vdupq_n_s32(u4_qbits);
    int16x4_t zeros_16x4 = vdup_n_s16(0);
    int16x4_t pos_255_16x4 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255_16x4 = vdup_n_s16(-((WORD16) UINT8_MAX));

    UNUSED(u1_use_upsampled_res);

    threshold0_16x4 = vld1_u16(pu2_threshold_matrix);
    threshold1_16x4 = vld1_u16(pu2_threshold_matrix + 4);
    threshold2_16x4 = vld1_u16(pu2_threshold_matrix + 8);
    threshold3_16x4 = vld1_u16(pu2_threshold_matrix + 12);

    scale_mat0_16x4 = vld1_u16(pu2_scale_matrix);
    scale_mat1_16x4 = vld1_u16(pu2_scale_matrix + 4);
    scale_mat2_16x4 = vld1_u16(pu2_scale_matrix + 8);
    scale_mat3_16x4 = vld1_u16(pu2_scale_matrix + 12);

    src0 = vld1_u8(&pu1_src[0 * i4_src_stride]);
    src1 = vld1_u8(&pu1_src[1 * i4_src_stride]);
    src2 = vld1_u8(&pu1_src[2 * i4_src_stride]);
    src3 = vld1_u8(&pu1_src[3 * i4_src_stride]);

    pred0 = vld1_u8(&pu1_pred[0 * i4_pred_stride]);
    pred1 = vld1_u8(&pu1_pred[1 * i4_pred_stride]);
    pred2 = vld1_u8(&pu1_pred[2 * i4_pred_stride]);
    pred3 = vld1_u8(&pu1_pred[3 * i4_pred_stride]);

    /* calculate res = src - pred */
    res0_16x8 = vsubl_u8(src0, pred0);
    res1_16x8 = vsubl_u8(src1, pred1);
    res2_16x8 = vsubl_u8(src2, pred2);
    res3_16x8 = vsubl_u8(src3, pred3);

    res0_16x4 = vreinterpret_s16_u16(vget_low_u16(res0_16x8));
    res1_16x4 = vreinterpret_s16_u16(vget_low_u16(res1_16x8));
    res2_16x4 = vreinterpret_s16_u16(vget_low_u16(res2_16x8));
    res3_16x4 = vreinterpret_s16_u16(vget_low_u16(res3_16x8));

    /* Load upsampled res */
    upres0_16x4 = vld1_s16(&pi2_upsampled_res[0 * i4_upsampled_res_stride]);
    upres1_16x4 = vld1_s16(&pi2_upsampled_res[1 * i4_upsampled_res_stride]);
    upres2_16x4 = vld1_s16(&pi2_upsampled_res[2 * i4_upsampled_res_stride]);
    upres3_16x4 = vld1_s16(&pi2_upsampled_res[3 * i4_upsampled_res_stride]);

    /* subtract upsampled res from (src - pred) to obtain final res */
    res0_16x4 = vsub_s16(res0_16x4, upres0_16x4);
    res1_16x4 = vsub_s16(res1_16x4, upres1_16x4);
    res2_16x4 = vsub_s16(res2_16x4, upres2_16x4);
    res3_16x4 = vsub_s16(res3_16x4, upres3_16x4);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    res0_16x4 = vmax_s16(res0_16x4, neg_255_16x4);
    res1_16x4 = vmax_s16(res1_16x4, neg_255_16x4);
    res2_16x4 = vmax_s16(res2_16x4, neg_255_16x4);
    res3_16x4 = vmax_s16(res3_16x4, neg_255_16x4);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    res0_16x4 = vmin_s16(res0_16x4, pos_255_16x4);
    res1_16x4 = vmin_s16(res1_16x4, pos_255_16x4);
    res2_16x4 = vmin_s16(res2_16x4, pos_255_16x4);
    res3_16x4 = vmin_s16(res3_16x4, pos_255_16x4);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    /* Matrix transpose */
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */

    xx0_16x4x2 = vtrn_s16(res0_16x4, res1_16x4);
    xx1_16x4x2 = vtrn_s16(res2_16x4, res3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(xx2_16x4, temp0_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* Matrix transpose */
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */

    xx0_16x4x2 = vtrn_s16(x0_16x4, x1_16x4);
    xx1_16x4x2 = vtrn_s16(x2_16x4, x3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Vertical Transformation */

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(temp0_16x4, xx2_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* get the first 16 bits from the register */
    *pi2_dc_out = vget_lane_s16(x0_16x4, 0);

    xx0_16x4 = vabs_s16(x0_16x4);
    xx1_16x4 = vabs_s16(x1_16x4);
    xx2_16x4 = vabs_s16(x2_16x4);
    xx3_16x4 = vabs_s16(x3_16x4);

    /* compare with zero for getting sign */
    temp0_u16x4 = vcgt_s16(x0_16x4, zeros_16x4);
    temp1_u16x4 = vcgt_s16(x1_16x4, zeros_16x4);
    temp2_u16x4 = vcgt_s16(x2_16x4, zeros_16x4);
    temp3_u16x4 = vcgt_s16(x3_16x4, zeros_16x4);

    /* compare with zero for thresholding */
    thresholdmask0_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold0_16x4), xx0_16x4);
    thresholdmask1_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold1_16x4), xx1_16x4);
    thresholdmask2_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold2_16x4), xx2_16x4);
    thresholdmask3_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold3_16x4), xx3_16x4);

    /* Multiply abs values obtained with scaling matrix */
    tx0_32x4 = vmull_s16(xx0_16x4, vreinterpret_s16_u16(scale_mat0_16x4));
    tx1_32x4 = vmull_s16(xx1_16x4, vreinterpret_s16_u16(scale_mat1_16x4));
    tx2_32x4 = vmull_s16(xx2_16x4, vreinterpret_s16_u16(scale_mat2_16x4));
    tx3_32x4 = vmull_s16(xx3_16x4, vreinterpret_s16_u16(scale_mat3_16x4));

    tx0_32x4 = vaddq_s32(tx0_32x4, rnd_factor_32x4);
    tx1_32x4 = vaddq_s32(tx1_32x4, rnd_factor_32x4);
    tx2_32x4 = vaddq_s32(tx2_32x4, rnd_factor_32x4);
    tx3_32x4 = vaddq_s32(tx3_32x4, rnd_factor_32x4);

    qbits_32x4 = vnegq_s32(qbits_32x4);

    tx0_32x4 = vshlq_s32(tx0_32x4, qbits_32x4);
    tx1_32x4 = vshlq_s32(tx1_32x4, qbits_32x4);
    tx2_32x4 = vshlq_s32(tx2_32x4, qbits_32x4);
    tx3_32x4 = vshlq_s32(tx3_32x4, qbits_32x4);

    /* Convertion to 16 bits signed */
    temp0_16x4 = vmovn_s32(tx0_32x4);
    temp1_16x4 = vmovn_s32(tx1_32x4);
    temp2_16x4 = vmovn_s32(tx2_32x4);
    temp3_16x4 = vmovn_s32(tx3_32x4);

    x0_16x4 = vneg_s16(temp0_16x4);
    x1_16x4 = vneg_s16(temp1_16x4);
    x2_16x4 = vneg_s16(temp2_16x4);
    x3_16x4 = vneg_s16(temp3_16x4);

    /* Restore sign */
    x0_16x4 = vbsl_s16(temp0_u16x4, temp0_16x4, x0_16x4);
    x1_16x4 = vbsl_s16(temp1_u16x4, temp1_16x4, x1_16x4);
    x2_16x4 = vbsl_s16(temp2_u16x4, temp2_16x4, x2_16x4);
    x3_16x4 = vbsl_s16(temp3_u16x4, temp3_16x4, x3_16x4);

    xx0_16x4 = vbsl_s16(thresholdmask0_16x4, zeros_16x4, x0_16x4);
    xx1_16x4 = vbsl_s16(thresholdmask1_16x4, zeros_16x4, x1_16x4);
    xx2_16x4 = vbsl_s16(thresholdmask2_16x4, zeros_16x4, x2_16x4);
    xx3_16x4 = vbsl_s16(thresholdmask3_16x4, zeros_16x4, x3_16x4);

    /* Store Quantized outputs */
    vst1_s16(&pi2_out[0 * i4_out_stride], xx0_16x4);
    vst1_s16(&pi2_out[1 * i4_out_stride], xx1_16x4);
    vst1_s16(&pi2_out[2 * i4_out_stride], xx2_16x4);
    vst1_s16(&pi2_out[3 * i4_out_stride], xx3_16x4);

    /* NNZ calculation */

    temp0_u16x4 = vceq_s16(xx0_16x4, zeros_16x4);
    temp1_u16x4 = vceq_s16(xx1_16x4, zeros_16x4);
    temp2_u16x4 = vceq_s16(xx2_16x4, zeros_16x4);
    temp3_u16x4 = vceq_s16(xx3_16x4, zeros_16x4);

    temp0_u16x8 = vcombine_u16(temp0_u16x4, temp2_u16x4);
    temp1_u16x8 = vcombine_u16(temp1_u16x4, temp3_u16x4);

    /* Convertion to 8 bit unsigned */
    temp0_u8x8 = vmovn_u16(temp0_u16x8);
    temp1_u8x8 = vmovn_u16(temp1_u16x8);

    temp0_u8x8 = vshr_n_u8(temp0_u8x8, 7);
    temp1_u8x8 = vshr_n_u8(temp1_u8x8, 7);

    temp0_u8x8 = vadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);

    *pu1_nnz = 16 - vget_lane_u8(temp0_u8x8, 0);
}

void isvc_resi_trans_quant_chroma_4x4_neon(buffer_container_t *ps_src, buffer_container_t *ps_pred,
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

    uint8x8_t src0, src1, src2, src3;
    uint8x8_t pred0, pred1, pred2, pred3;
    uint8x8x2_t tmp0, tmp1, tmp2, tmp3;
    uint8x8_t temp0_u8x8, temp1_u8x8;
    uint16x4_t temp0_u16x4, temp1_u16x4, temp2_u16x4, temp3_u16x4;
    uint16x4_t scale_mat0_16x4, scale_mat1_16x4, scale_mat2_16x4, scale_mat3_16x4;
    uint16x4_t threshold0_16x4, threshold1_16x4, threshold2_16x4, threshold3_16x4;
    uint16x4_t thresholdmask0_16x4, thresholdmask1_16x4, thresholdmask2_16x4, thresholdmask3_16x4;
    int16x4_t res0_16x4, res1_16x4, res2_16x4, res3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int16x4_t temp0_16x4, temp1_16x4, temp2_16x4, temp3_16x4;
    uint16x8_t res0_16x8, res1_16x8, res2_16x8, res3_16x8;
    uint16x8_t temp0_u16x8, temp1_u16x8;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int32x4_t tx0_32x4, tx1_32x4, tx2_32x4, tx3_32x4;

    int32x4_t rnd_factor_32x4 = vdupq_n_s32(u4_round_factor);
    int32x4_t qbits_32x4 = vdupq_n_s32(u4_qbits);
    int16x4_t zeros_16x4 = vdup_n_s16(0);

    UNUSED(ps_upsampled_res);
    UNUSED(u1_use_upsampled_res);

    threshold0_16x4 = vld1_u16(pu2_threshold_matrix);
    threshold1_16x4 = vld1_u16(pu2_threshold_matrix + 4);
    threshold2_16x4 = vld1_u16(pu2_threshold_matrix + 8);
    threshold3_16x4 = vld1_u16(pu2_threshold_matrix + 12);

    scale_mat0_16x4 = vld1_u16(pu2_scale_matrix);
    scale_mat1_16x4 = vld1_u16(pu2_scale_matrix + 4);
    scale_mat2_16x4 = vld1_u16(pu2_scale_matrix + 8);
    scale_mat3_16x4 = vld1_u16(pu2_scale_matrix + 12);

    src0 = vld1_u8(&pu1_src[0 * i4_src_stride]);
    src1 = vld1_u8(&pu1_src[1 * i4_src_stride]);
    src2 = vld1_u8(&pu1_src[2 * i4_src_stride]);
    src3 = vld1_u8(&pu1_src[3 * i4_src_stride]);

    /* deinterleaving source buffer */
    tmp0 = vuzp_u8(src0, src0);
    tmp1 = vuzp_u8(src1, src1);
    tmp2 = vuzp_u8(src2, src2);
    tmp3 = vuzp_u8(src3, src3);

    src0 = tmp0.val[0];
    src1 = tmp1.val[0];
    src2 = tmp2.val[0];
    src3 = tmp3.val[0];

    pred0 = vld1_u8(&pu1_pred[0 * i4_pred_stride]);
    pred1 = vld1_u8(&pu1_pred[1 * i4_pred_stride]);
    pred2 = vld1_u8(&pu1_pred[2 * i4_pred_stride]);
    pred3 = vld1_u8(&pu1_pred[3 * i4_pred_stride]);

    /* deinterleaving pred buffer */
    tmp0 = vuzp_u8(pred0, pred0);
    tmp1 = vuzp_u8(pred1, pred1);
    tmp2 = vuzp_u8(pred2, pred2);
    tmp3 = vuzp_u8(pred3, pred3);

    pred0 = tmp0.val[0];
    pred1 = tmp1.val[0];
    pred2 = tmp2.val[0];
    pred3 = tmp3.val[0];

    /* calculate res = src - pred */
    res0_16x8 = vsubl_u8(src0, pred0);
    res1_16x8 = vsubl_u8(src1, pred1);
    res2_16x8 = vsubl_u8(src2, pred2);
    res3_16x8 = vsubl_u8(src3, pred3);

    res0_16x4 = vreinterpret_s16_u16(vget_low_u16(res0_16x8));
    res1_16x4 = vreinterpret_s16_u16(vget_low_u16(res1_16x8));
    res2_16x4 = vreinterpret_s16_u16(vget_low_u16(res2_16x8));
    res3_16x4 = vreinterpret_s16_u16(vget_low_u16(res3_16x8));

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    /* Matrix transpose */
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */

    xx0_16x4x2 = vtrn_s16(res0_16x4, res1_16x4);
    xx1_16x4x2 = vtrn_s16(res2_16x4, res3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(xx2_16x4, temp0_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* Matrix transpose */
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */

    xx0_16x4x2 = vtrn_s16(x0_16x4, x1_16x4);
    xx1_16x4x2 = vtrn_s16(x2_16x4, x3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Vertical Transformation */

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(temp0_16x4, xx2_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* get the first 16 bits from the register */
    *pi2_dc_out = vget_lane_s16(x0_16x4, 0);

    xx0_16x4 = vabs_s16(x0_16x4);
    xx1_16x4 = vabs_s16(x1_16x4);
    xx2_16x4 = vabs_s16(x2_16x4);
    xx3_16x4 = vabs_s16(x3_16x4);

    /* compare with zero for getting sign */
    temp0_u16x4 = vcgt_s16(x0_16x4, zeros_16x4);
    temp1_u16x4 = vcgt_s16(x1_16x4, zeros_16x4);
    temp2_u16x4 = vcgt_s16(x2_16x4, zeros_16x4);
    temp3_u16x4 = vcgt_s16(x3_16x4, zeros_16x4);

    /* compare with zero for thresholding */
    thresholdmask0_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold0_16x4), xx0_16x4);
    thresholdmask1_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold1_16x4), xx1_16x4);
    thresholdmask2_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold2_16x4), xx2_16x4);
    thresholdmask3_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold3_16x4), xx3_16x4);

    /* Multiply abs values obtained with scaling matrix */
    tx0_32x4 = vmull_s16(xx0_16x4, vreinterpret_s16_u16(scale_mat0_16x4));
    tx1_32x4 = vmull_s16(xx1_16x4, vreinterpret_s16_u16(scale_mat1_16x4));
    tx2_32x4 = vmull_s16(xx2_16x4, vreinterpret_s16_u16(scale_mat2_16x4));
    tx3_32x4 = vmull_s16(xx3_16x4, vreinterpret_s16_u16(scale_mat3_16x4));

    tx0_32x4 = vaddq_s32(tx0_32x4, rnd_factor_32x4);
    tx1_32x4 = vaddq_s32(tx1_32x4, rnd_factor_32x4);
    tx2_32x4 = vaddq_s32(tx2_32x4, rnd_factor_32x4);
    tx3_32x4 = vaddq_s32(tx3_32x4, rnd_factor_32x4);

    qbits_32x4 = vnegq_s32(qbits_32x4);

    tx0_32x4 = vshlq_s32(tx0_32x4, qbits_32x4);
    tx1_32x4 = vshlq_s32(tx1_32x4, qbits_32x4);
    tx2_32x4 = vshlq_s32(tx2_32x4, qbits_32x4);
    tx3_32x4 = vshlq_s32(tx3_32x4, qbits_32x4);

    /* Convertion to 16 bits signed */
    temp0_16x4 = vmovn_s32(tx0_32x4);
    temp1_16x4 = vmovn_s32(tx1_32x4);
    temp2_16x4 = vmovn_s32(tx2_32x4);
    temp3_16x4 = vmovn_s32(tx3_32x4);

    x0_16x4 = vneg_s16(temp0_16x4);
    x1_16x4 = vneg_s16(temp1_16x4);
    x2_16x4 = vneg_s16(temp2_16x4);
    x3_16x4 = vneg_s16(temp3_16x4);

    /* Restore sign */
    x0_16x4 = vbsl_s16(temp0_u16x4, temp0_16x4, x0_16x4);
    x1_16x4 = vbsl_s16(temp1_u16x4, temp1_16x4, x1_16x4);
    x2_16x4 = vbsl_s16(temp2_u16x4, temp2_16x4, x2_16x4);
    x3_16x4 = vbsl_s16(temp3_u16x4, temp3_16x4, x3_16x4);

    /* Thresholding */
    xx0_16x4 = vbsl_s16(thresholdmask0_16x4, zeros_16x4, x0_16x4);
    xx1_16x4 = vbsl_s16(thresholdmask1_16x4, zeros_16x4, x1_16x4);
    xx2_16x4 = vbsl_s16(thresholdmask2_16x4, zeros_16x4, x2_16x4);
    xx3_16x4 = vbsl_s16(thresholdmask3_16x4, zeros_16x4, x3_16x4);

    /* Store Quantized outputs */
    vst1_s16(&pi2_out[0 * i4_out_stride], xx0_16x4);
    vst1_s16(&pi2_out[1 * i4_out_stride], xx1_16x4);
    vst1_s16(&pi2_out[2 * i4_out_stride], xx2_16x4);
    vst1_s16(&pi2_out[3 * i4_out_stride], xx3_16x4);

    /* NNZ calculation */

    temp0_u16x4 = vceq_s16(xx0_16x4, zeros_16x4);
    temp1_u16x4 = vceq_s16(xx1_16x4, zeros_16x4);
    temp2_u16x4 = vceq_s16(xx2_16x4, zeros_16x4);
    temp3_u16x4 = vceq_s16(xx3_16x4, zeros_16x4);

    temp0_u16x8 = vcombine_u16(temp0_u16x4, temp2_u16x4);
    temp1_u16x8 = vcombine_u16(temp1_u16x4, temp3_u16x4);

    /* Convertion to 8 bit unsigned */
    temp0_u8x8 = vmovn_u16(temp0_u16x8);
    temp1_u8x8 = vmovn_u16(temp1_u16x8);

    temp0_u8x8 = vshr_n_u8(temp0_u8x8, 7);
    temp1_u8x8 = vshr_n_u8(temp1_u8x8, 7);

    temp0_u8x8 = vadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);

    *pu1_nnz = 16 - vget_lane_u8(temp0_u8x8, 0);
}

void isvc_resi_trans_quant_chroma_4x4_with_residual_sub_neon(
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

    uint8x8_t src0, src1, src2, src3;
    uint8x8_t pred0, pred1, pred2, pred3;
    uint8x8x2_t tmp0, tmp1, tmp2, tmp3;
    uint8x8_t temp0_u8x8, temp1_u8x8;
    uint16x4_t temp0_u16x4, temp1_u16x4, temp2_u16x4, temp3_u16x4;
    uint16x4_t scale_mat0_16x4, scale_mat1_16x4, scale_mat2_16x4, scale_mat3_16x4;
    uint16x4_t threshold0_16x4, threshold1_16x4, threshold2_16x4, threshold3_16x4;
    uint16x4_t thresholdmask0_16x4, thresholdmask1_16x4, thresholdmask2_16x4, thresholdmask3_16x4;
    int16x4_t upres0_16x4, upres1_16x4, upres2_16x4, upres3_16x4;
    int16x4_t res0_16x4, res1_16x4, res2_16x4, res3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int16x4_t temp0_16x4, temp1_16x4, temp2_16x4, temp3_16x4;
    uint16x8_t res0_16x8, res1_16x8, res2_16x8, res3_16x8;
    uint16x8_t temp0_u16x8, temp1_u16x8;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int32x4_t tx0_32x4, tx1_32x4, tx2_32x4, tx3_32x4;

    int32x4_t rnd_factor_32x4 = vdupq_n_s32(u4_round_factor);
    int32x4_t qbits_32x4 = vdupq_n_s32(u4_qbits);
    int16x4_t zeros_16x4 = vdup_n_s16(0);
    int16x4_t pos_255_16x4 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255_16x4 = vdup_n_s16(-((WORD16) UINT8_MAX));

    UNUSED(u1_use_upsampled_res);

    threshold0_16x4 = vld1_u16(pu2_threshold_matrix);
    threshold1_16x4 = vld1_u16(pu2_threshold_matrix + 4);
    threshold2_16x4 = vld1_u16(pu2_threshold_matrix + 8);
    threshold3_16x4 = vld1_u16(pu2_threshold_matrix + 12);

    scale_mat0_16x4 = vld1_u16(pu2_scale_matrix);
    scale_mat1_16x4 = vld1_u16(pu2_scale_matrix + 4);
    scale_mat2_16x4 = vld1_u16(pu2_scale_matrix + 8);
    scale_mat3_16x4 = vld1_u16(pu2_scale_matrix + 12);

    src0 = vld1_u8(&pu1_src[0 * i4_src_stride]);
    src1 = vld1_u8(&pu1_src[1 * i4_src_stride]);
    src2 = vld1_u8(&pu1_src[2 * i4_src_stride]);
    src3 = vld1_u8(&pu1_src[3 * i4_src_stride]);

    /* deinterleaving source buffer */
    tmp0 = vuzp_u8(src0, src0);
    tmp1 = vuzp_u8(src1, src1);
    tmp2 = vuzp_u8(src2, src2);
    tmp3 = vuzp_u8(src3, src3);

    src0 = tmp0.val[0];
    src1 = tmp1.val[0];
    src2 = tmp2.val[0];
    src3 = tmp3.val[0];

    pred0 = vld1_u8(&pu1_pred[0 * i4_pred_stride]);
    pred1 = vld1_u8(&pu1_pred[1 * i4_pred_stride]);
    pred2 = vld1_u8(&pu1_pred[2 * i4_pred_stride]);
    pred3 = vld1_u8(&pu1_pred[3 * i4_pred_stride]);

    /* deinterleaving pred buffer */
    tmp0 = vuzp_u8(pred0, pred0);
    tmp1 = vuzp_u8(pred1, pred1);
    tmp2 = vuzp_u8(pred2, pred2);
    tmp3 = vuzp_u8(pred3, pred3);

    pred0 = tmp0.val[0];
    pred1 = tmp1.val[0];
    pred2 = tmp2.val[0];
    pred3 = tmp3.val[0];

    /* calculate res = src - pred */
    res0_16x8 = vsubl_u8(src0, pred0);
    res1_16x8 = vsubl_u8(src1, pred1);
    res2_16x8 = vsubl_u8(src2, pred2);
    res3_16x8 = vsubl_u8(src3, pred3);

    res0_16x4 = vreinterpret_s16_u16(vget_low_u16(res0_16x8));
    res1_16x4 = vreinterpret_s16_u16(vget_low_u16(res1_16x8));
    res2_16x4 = vreinterpret_s16_u16(vget_low_u16(res2_16x8));
    res3_16x4 = vreinterpret_s16_u16(vget_low_u16(res3_16x8));

    /* Load upsampled res */
    upres0_16x4 = vld1_s16(&pi2_upsampled_res[0 * i4_upsampled_res_stride]);
    upres1_16x4 = vld1_s16(&pi2_upsampled_res[1 * i4_upsampled_res_stride]);
    upres2_16x4 = vld1_s16(&pi2_upsampled_res[2 * i4_upsampled_res_stride]);
    upres3_16x4 = vld1_s16(&pi2_upsampled_res[3 * i4_upsampled_res_stride]);

    /* subtract upsampled res from (src - pred) to obtain final res */
    res0_16x4 = vsub_s16(res0_16x4, upres0_16x4);
    res1_16x4 = vsub_s16(res1_16x4, upres1_16x4);
    res2_16x4 = vsub_s16(res2_16x4, upres2_16x4);
    res3_16x4 = vsub_s16(res3_16x4, upres3_16x4);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    res0_16x4 = vmax_s16(res0_16x4, neg_255_16x4);
    res1_16x4 = vmax_s16(res1_16x4, neg_255_16x4);
    res2_16x4 = vmax_s16(res2_16x4, neg_255_16x4);
    res3_16x4 = vmax_s16(res3_16x4, neg_255_16x4);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    res0_16x4 = vmin_s16(res0_16x4, pos_255_16x4);
    res1_16x4 = vmin_s16(res1_16x4, pos_255_16x4);
    res2_16x4 = vmin_s16(res2_16x4, pos_255_16x4);
    res3_16x4 = vmin_s16(res3_16x4, pos_255_16x4);

    /* Perform Forward transform */
    /*-------------------------------------------------------------*/
    /* DCT [ Horizontal transformation ]                          */
    /*-------------------------------------------------------------*/
    /* Matrix transpose */
    /*
     *  a0 a1 a2 a3
     *  b0 b1 b2 b3
     *  c0 c1 c2 c3
     *  d0 d1 d2 d3
     */

    xx0_16x4x2 = vtrn_s16(res0_16x4, res1_16x4);
    xx1_16x4x2 = vtrn_s16(res2_16x4, res3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(xx2_16x4, temp0_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* Matrix transpose */
    /*
     *  a0 b0 c0 d0
     *  a1 b1 c1 d1
     *  a2 b2 c2 d2
     *  a3 b3 c3 d3
     */

    xx0_16x4x2 = vtrn_s16(x0_16x4, x1_16x4);
    xx1_16x4x2 = vtrn_s16(x2_16x4, x3_16x4);
    x0_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vtrn_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Vertical Transformation */

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx3_16x4, 1);
    x1_16x4 = vadd_s16(temp0_16x4, xx2_16x4);

    x2_16x4 = vsub_s16(xx0_16x4, xx1_16x4);
    temp0_16x4 = vshl_n_s16(xx2_16x4, 1);
    x3_16x4 = vsub_s16(xx3_16x4, temp0_16x4);

    /* get the first 16 bits from the register */
    *pi2_dc_out = vget_lane_s16(x0_16x4, 0);

    xx0_16x4 = vabs_s16(x0_16x4);
    xx1_16x4 = vabs_s16(x1_16x4);
    xx2_16x4 = vabs_s16(x2_16x4);
    xx3_16x4 = vabs_s16(x3_16x4);

    /* compare with zero for getting sign */
    temp0_u16x4 = vcgt_s16(x0_16x4, zeros_16x4);
    temp1_u16x4 = vcgt_s16(x1_16x4, zeros_16x4);
    temp2_u16x4 = vcgt_s16(x2_16x4, zeros_16x4);
    temp3_u16x4 = vcgt_s16(x3_16x4, zeros_16x4);

    thresholdmask0_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold0_16x4), xx0_16x4);
    thresholdmask1_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold1_16x4), xx1_16x4);
    thresholdmask2_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold2_16x4), xx2_16x4);
    thresholdmask3_16x4 = vcgt_s16(vreinterpret_s16_u16(threshold3_16x4), xx3_16x4);

    /* Multiply abs values obtained with scaling matrix */
    tx0_32x4 = vmull_s16(xx0_16x4, vreinterpret_s16_u16(scale_mat0_16x4));
    tx1_32x4 = vmull_s16(xx1_16x4, vreinterpret_s16_u16(scale_mat1_16x4));
    tx2_32x4 = vmull_s16(xx2_16x4, vreinterpret_s16_u16(scale_mat2_16x4));
    tx3_32x4 = vmull_s16(xx3_16x4, vreinterpret_s16_u16(scale_mat3_16x4));

    tx0_32x4 = vaddq_s32(tx0_32x4, rnd_factor_32x4);
    tx1_32x4 = vaddq_s32(tx1_32x4, rnd_factor_32x4);
    tx2_32x4 = vaddq_s32(tx2_32x4, rnd_factor_32x4);
    tx3_32x4 = vaddq_s32(tx3_32x4, rnd_factor_32x4);

    qbits_32x4 = vnegq_s32(qbits_32x4);

    tx0_32x4 = vshlq_s32(tx0_32x4, qbits_32x4);
    tx1_32x4 = vshlq_s32(tx1_32x4, qbits_32x4);
    tx2_32x4 = vshlq_s32(tx2_32x4, qbits_32x4);
    tx3_32x4 = vshlq_s32(tx3_32x4, qbits_32x4);

    /* Convertion to 16 bits signed */
    temp0_16x4 = vmovn_s32(tx0_32x4);
    temp1_16x4 = vmovn_s32(tx1_32x4);
    temp2_16x4 = vmovn_s32(tx2_32x4);
    temp3_16x4 = vmovn_s32(tx3_32x4);

    x0_16x4 = vneg_s16(temp0_16x4);
    x1_16x4 = vneg_s16(temp1_16x4);
    x2_16x4 = vneg_s16(temp2_16x4);
    x3_16x4 = vneg_s16(temp3_16x4);

    /* Restore sign */
    x0_16x4 = vbsl_s16(temp0_u16x4, temp0_16x4, x0_16x4);
    x1_16x4 = vbsl_s16(temp1_u16x4, temp1_16x4, x1_16x4);
    x2_16x4 = vbsl_s16(temp2_u16x4, temp2_16x4, x2_16x4);
    x3_16x4 = vbsl_s16(temp3_u16x4, temp3_16x4, x3_16x4);

    xx0_16x4 = vbsl_s16(thresholdmask0_16x4, zeros_16x4, x0_16x4);
    xx1_16x4 = vbsl_s16(thresholdmask1_16x4, zeros_16x4, x1_16x4);
    xx2_16x4 = vbsl_s16(thresholdmask2_16x4, zeros_16x4, x2_16x4);
    xx3_16x4 = vbsl_s16(thresholdmask3_16x4, zeros_16x4, x3_16x4);

    /* Store Quantized outputs */
    vst1_s16(&pi2_out[0 * i4_out_stride], xx0_16x4);
    vst1_s16(&pi2_out[1 * i4_out_stride], xx1_16x4);
    vst1_s16(&pi2_out[2 * i4_out_stride], xx2_16x4);
    vst1_s16(&pi2_out[3 * i4_out_stride], xx3_16x4);

    /* NNZ calculation */

    temp0_u16x4 = vceq_s16(xx0_16x4, zeros_16x4);
    temp1_u16x4 = vceq_s16(xx1_16x4, zeros_16x4);
    temp2_u16x4 = vceq_s16(xx2_16x4, zeros_16x4);
    temp3_u16x4 = vceq_s16(xx3_16x4, zeros_16x4);

    temp0_u16x8 = vcombine_u16(temp0_u16x4, temp2_u16x4);
    temp1_u16x8 = vcombine_u16(temp1_u16x4, temp3_u16x4);

    /* Convertion to 8 bit unsigned */
    temp0_u8x8 = vmovn_u16(temp0_u16x8);
    temp1_u8x8 = vmovn_u16(temp1_u16x8);

    temp0_u8x8 = vshr_n_u8(temp0_u8x8, 7);
    temp1_u8x8 = vshr_n_u8(temp1_u8x8, 7);

    temp0_u8x8 = vadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);
    temp0_u8x8 = vpadd_u8(temp0_u8x8, temp1_u8x8);

    *pu1_nnz = 16 - vget_lane_u8(temp0_u8x8, 0);
}
