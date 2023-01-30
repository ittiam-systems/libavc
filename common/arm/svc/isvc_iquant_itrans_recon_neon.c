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
 *  isvc_iquant_itrans_recon_neon.c
 *
 * @brief
 *  neon variants of inverse transform and quantization functions
 *
 * *******************************************************************************
 */
#include <arm_neon.h>

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

void isvc_iquant_itrans_recon_4x4_neon(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                       buffer_container_t *ps_res_pred, buffer_container_t *ps_res,
                                       buffer_container_t *ps_rec,
                                       iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
                                       WORD16 *pi2_tmp, WORD16 *pi2_dc_src, WORD32 i4_iq_start_idx,
                                       UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd01_in, resd23_in;
    int16x8_t pred01_in, pred23_in;
    uint8x8_t pred01_un, pred23_un;

    int16x8_t pos_255_16x8 = vdupq_n_s16(((WORD16) UINT8_MAX));
    int16x8_t neg_255_16x8 = vdupq_n_s16(-((WORD16) UINT8_MAX));
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);

    WORD16 rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    int32x4_t rnd_fact = vdupq_n_s32(rnd_factor);

    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    if(i4_iq_start_idx == 1)
    {
        q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);
    }

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    resd01_in = vcombine_s16(x0_16x4, x1_16x4);
    resd23_in = vcombine_s16(x2_16x4, x3_16x4);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    resd01_in = vmaxq_s16(resd01_in, neg_255_16x8);
    resd23_in = vmaxq_s16(resd23_in, neg_255_16x8);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    resd01_in = vminq_s16(resd01_in, pos_255_16x8);
    resd23_in = vminq_s16(resd23_in, pos_255_16x8);

    /* Load pred */
    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    pred01_in = vcombine_s16(vget_low_s16(pred0), vget_low_s16(pred1));
    pred23_in = vcombine_s16(vget_low_s16(pred2), vget_low_s16(pred3));

    /* Out pixel = pred + res */
    pred01_in = vaddq_s16(pred01_in, resd01_in);
    pred23_in = vaddq_s16(pred23_in, resd23_in);

    /* Convert to 8 bit unsigned with saturation */
    pred01_un = vqmovun_s16(pred01_in);
    pred23_un = vqmovun_s16(pred23_in);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred01_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred01_un), 1);
    vst1_lane_u32((uint32_t *) (pu1_out + (i4_out_stride << 1)), vreinterpret_u32_u8(pred23_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + ((i4_out_stride << 1) + i4_out_stride)),
                  vreinterpret_u32_u8(pred23_un), 1);
}

void isvc_iquant_itrans_recon_4x4_with_res_output_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd01_in, resd23_in;
    int16x8_t pred01_in, pred23_in;
    uint8x8_t pred01_un, pred23_un;

    int16x4_t pos_255_16x4 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255_16x4 = vdup_n_s16(-((WORD16) UINT8_MAX));
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);

    WORD16 rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    int32x4_t rnd_fact = vdupq_n_s32(rnd_factor);

    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    if(i4_iq_start_idx == 1)
    {
        q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);
    }

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    x0_16x4 = vmax_s16(x0_16x4, neg_255_16x4);
    x1_16x4 = vmax_s16(x1_16x4, neg_255_16x4);
    x2_16x4 = vmax_s16(x2_16x4, neg_255_16x4);
    x3_16x4 = vmax_s16(x3_16x4, neg_255_16x4);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    x0_16x4 = vmin_s16(x0_16x4, pos_255_16x4);
    x1_16x4 = vmin_s16(x1_16x4, pos_255_16x4);
    x2_16x4 = vmin_s16(x2_16x4, pos_255_16x4);
    x3_16x4 = vmin_s16(x3_16x4, pos_255_16x4);

    vst1_s16(pi2_res, x0_16x4);
    vst1_s16(pi2_res + i4_res_stride, x1_16x4);
    vst1_s16(pi2_res + (i4_res_stride << 1), x2_16x4);
    vst1_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, x3_16x4);

    resd01_in = vcombine_s16(x0_16x4, x1_16x4);
    resd23_in = vcombine_s16(x2_16x4, x3_16x4);

    /* Load pred */
    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    pred01_in = vcombine_s16(vget_low_s16(pred0), vget_low_s16(pred1));
    pred23_in = vcombine_s16(vget_low_s16(pred2), vget_low_s16(pred3));

    /* Out pixel = pred + res */
    pred01_in = vaddq_s16(pred01_in, resd01_in);
    pred23_in = vaddq_s16(pred23_in, resd23_in);

    /* Convert to 8 bit unsigned with saturation */
    pred01_un = vqmovun_s16(pred01_in);
    pred23_un = vqmovun_s16(pred23_in);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred01_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred01_un), 1);
    vst1_lane_u32((uint32_t *) (pu1_out + (i4_out_stride << 1)), vreinterpret_u32_u8(pred23_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + ((i4_out_stride << 1) + i4_out_stride)),
                  vreinterpret_u32_u8(pred23_un), 1);
}

void isvc_iquant_itrans_recon_4x4_with_res_accumulate_neon(
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

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x4_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd01_in, resd23_in;
    int16x8_t pred01_in, pred23_in;
    uint8x8_t pred01_un, pred23_un;

    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);

    WORD16 rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    int32x4_t rnd_fact = vdupq_n_s32(rnd_factor);
    int16x4_t pos_255 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255 = vdup_n_s16(-((WORD16) UINT8_MAX));

    UNUSED(u1_res_accumulate);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    if(i4_iq_start_idx == 1)
    {
        q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);
    }

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    /* Accumulating Res */

    /* Load Res pred */
    resd0_in = vld1_s16((int16_t *) pi2_res_pred);
    resd1_in = vld1_s16((int16_t *) pi2_res_pred + i4_res_pred_stride);
    resd2_in = vld1_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 2));
    resd3_in = vld1_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 3));

    /* Add res pred with res obtained */
    resd0_in = vadd_s16(resd0_in, x0_16x4);
    resd1_in = vadd_s16(resd1_in, x1_16x4);
    resd2_in = vadd_s16(resd2_in, x2_16x4);
    resd3_in = vadd_s16(resd3_in, x3_16x4);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    resd0_in = vmax_s16(resd0_in, neg_255);
    resd1_in = vmax_s16(resd1_in, neg_255);
    resd2_in = vmax_s16(resd2_in, neg_255);
    resd3_in = vmax_s16(resd3_in, neg_255);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    resd0_in = vmin_s16(resd0_in, pos_255);
    resd1_in = vmin_s16(resd1_in, pos_255);
    resd2_in = vmin_s16(resd2_in, pos_255);
    resd3_in = vmin_s16(resd3_in, pos_255);

    vst1_s16(pi2_res, resd0_in);
    vst1_s16(pi2_res + i4_res_stride, resd1_in);
    vst1_s16(pi2_res + (i4_res_stride << 1), resd2_in);
    vst1_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resd3_in);

    resd01_in = vcombine_s16(resd0_in, resd1_in);
    resd23_in = vcombine_s16(resd2_in, resd3_in);

    /* Load pred */
    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    pred01_in = vcombine_s16(vget_low_s16(pred0), vget_low_s16(pred1));
    pred23_in = vcombine_s16(vget_low_s16(pred2), vget_low_s16(pred3));

    /* Out pixel = pred + res */
    pred01_in = vaddq_s16(pred01_in, resd01_in);
    pred23_in = vaddq_s16(pred23_in, resd23_in);

    /* Convert to 8 bit unsigned with saturation */
    pred01_un = vqmovun_s16(pred01_in);
    pred23_un = vqmovun_s16(pred23_in);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred01_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred01_un), 1);
    vst1_lane_u32((uint32_t *) (pu1_out + (i4_out_stride << 1)), vreinterpret_u32_u8(pred23_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + ((i4_out_stride << 1) + i4_out_stride)),
                  vreinterpret_u32_u8(pred23_un), 1);
}

void isvc_iquant_itrans_recon_chroma_4x4_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    WORD16 i2_rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x8_t x0_16x8, x1_16x8, x2_16x8, x3_16x8;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t rec0, rec1, rec2, rec3;
    uint8x8_t rec0_un, rec1_un, rec2_un, rec3_un;
    uint8x8_t out0, out1, out2, out3;

    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));

    int16x4_t pos_255_16x4 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255_16x4 = vdup_n_s16(-((WORD16) UINT8_MAX));
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);
    int32x4_t rnd_fact = vdupq_n_s32(i2_rnd_factor);

    UNUSED(i4_iq_start_idx);
    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    x0_16x4 = vmax_s16(x0_16x4, neg_255_16x4);
    x1_16x4 = vmax_s16(x1_16x4, neg_255_16x4);
    x2_16x4 = vmax_s16(x2_16x4, neg_255_16x4);
    x3_16x4 = vmax_s16(x3_16x4, neg_255_16x4);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    x0_16x4 = vmin_s16(x0_16x4, pos_255_16x4);
    x1_16x4 = vmin_s16(x1_16x4, pos_255_16x4);
    x2_16x4 = vmin_s16(x2_16x4, pos_255_16x4);
    x3_16x4 = vmin_s16(x3_16x4, pos_255_16x4);

    x0_16x8 = vreinterpretq_s16_s32(vmovl_s16(x0_16x4));
    x1_16x8 = vreinterpretq_s16_s32(vmovl_s16(x1_16x4));
    x2_16x8 = vreinterpretq_s16_s32(vmovl_s16(x2_16x4));
    x3_16x8 = vreinterpretq_s16_s32(vmovl_s16(x3_16x4));

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = pred + res */
    rec0 = vaddq_s16(pred0, x0_16x8);
    rec1 = vaddq_s16(pred1, x1_16x8);
    rec2 = vaddq_s16(pred2, x2_16x8);
    rec3 = vaddq_s16(pred3, x3_16x8);

    out0 = vld1_u8(pu1_out);
    out1 = vld1_u8(pu1_out + i4_out_stride);
    out2 = vld1_u8(pu1_out + i4_out_stride * 2);
    out3 = vld1_u8(pu1_out + i4_out_stride * 3);

    /* Convert to 8 bit unsigned with saturation */
    rec0_un = vqmovun_s16(rec0);
    rec1_un = vqmovun_s16(rec1);
    rec2_un = vqmovun_s16(rec2);
    rec3_un = vqmovun_s16(rec3);

    /* Store in alternate postions */
    out0 = vbsl_u8(chroma_mask_8x8, rec0_un, out0);
    out1 = vbsl_u8(chroma_mask_8x8, rec1_un, out1);
    out2 = vbsl_u8(chroma_mask_8x8, rec2_un, out2);
    out3 = vbsl_u8(chroma_mask_8x8, rec3_un, out3);

    vst1_u8((pu1_out), out0);
    vst1_u8((pu1_out + i4_out_stride), out1);
    vst1_u8((pu1_out + (i4_out_stride << 1)), out2);
    vst1_u8((pu1_out + ((i4_out_stride << 1) + i4_out_stride)), out3);
}

void isvc_iquant_itrans_recon_chroma_4x4_with_res_output_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    WORD16 i2_rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x8_t x0_16x8, x1_16x8, x2_16x8, x3_16x8;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t rec0, rec1, rec2, rec3;
    uint8x8_t rec0_un, rec1_un, rec2_un, rec3_un;
    uint8x8_t out0, out1, out2, out3;
    int16x8_t resout0, resout1, resout2, resout3;

    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));
    uint16x8_t chroma_mask_16x8 = {0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000};
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);
    int32x4_t rnd_fact = vdupq_n_s32(i2_rnd_factor);
    int16x4_t pos_255_16x4 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255_16x4 = vdup_n_s16(-((WORD16) UINT8_MAX));

    UNUSED(i4_iq_start_idx);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    x0_16x4 = vmax_s16(x0_16x4, neg_255_16x4);
    x1_16x4 = vmax_s16(x1_16x4, neg_255_16x4);
    x2_16x4 = vmax_s16(x2_16x4, neg_255_16x4);
    x3_16x4 = vmax_s16(x3_16x4, neg_255_16x4);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    x0_16x4 = vmin_s16(x0_16x4, pos_255_16x4);
    x1_16x4 = vmin_s16(x1_16x4, pos_255_16x4);
    x2_16x4 = vmin_s16(x2_16x4, pos_255_16x4);
    x3_16x4 = vmin_s16(x3_16x4, pos_255_16x4);

    resout0 = vld1q_s16(pi2_res);
    resout1 = vld1q_s16(pi2_res + i4_res_stride);
    resout2 = vld1q_s16(pi2_res + i4_res_stride * 2);
    resout3 = vld1q_s16(pi2_res + i4_res_stride * 3);

    x0_16x8 = vreinterpretq_s16_s32(vmovl_s16(x0_16x4));
    x1_16x8 = vreinterpretq_s16_s32(vmovl_s16(x1_16x4));
    x2_16x8 = vreinterpretq_s16_s32(vmovl_s16(x2_16x4));
    x3_16x8 = vreinterpretq_s16_s32(vmovl_s16(x3_16x4));

    /* Storing res in alternate positions */
    resout0 = vbslq_s16(chroma_mask_16x8, x0_16x8, resout0);
    resout1 = vbslq_s16(chroma_mask_16x8, x1_16x8, resout1);
    resout2 = vbslq_s16(chroma_mask_16x8, x2_16x8, resout2);
    resout3 = vbslq_s16(chroma_mask_16x8, x3_16x8, resout3);

    vst1q_s16(pi2_res, resout0);
    vst1q_s16(pi2_res + i4_res_stride, resout1);
    vst1q_s16(pi2_res + (i4_res_stride << 1), resout2);
    vst1q_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resout3);

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = pred + res */
    rec0 = vaddq_s16(pred0, x0_16x8);
    rec1 = vaddq_s16(pred1, x1_16x8);
    rec2 = vaddq_s16(pred2, x2_16x8);
    rec3 = vaddq_s16(pred3, x3_16x8);

    out0 = vld1_u8(pu1_out);
    out1 = vld1_u8(pu1_out + i4_out_stride);
    out2 = vld1_u8(pu1_out + i4_out_stride * 2);
    out3 = vld1_u8(pu1_out + i4_out_stride * 3);

    /* Convert to 8 bit unsigned with saturation */
    rec0_un = vqmovun_s16(rec0);
    rec1_un = vqmovun_s16(rec1);
    rec2_un = vqmovun_s16(rec2);
    rec3_un = vqmovun_s16(rec3);

    /* Store output pixels in alternate positions */
    out0 = vbsl_u8(chroma_mask_8x8, rec0_un, out0);
    out1 = vbsl_u8(chroma_mask_8x8, rec1_un, out1);
    out2 = vbsl_u8(chroma_mask_8x8, rec2_un, out2);
    out3 = vbsl_u8(chroma_mask_8x8, rec3_un, out3);

    vst1_u8((pu1_out), out0);
    vst1_u8((pu1_out + i4_out_stride), out1);
    vst1_u8((pu1_out + (i4_out_stride << 1)), out2);
    vst1_u8((pu1_out + ((i4_out_stride << 1) + i4_out_stride)), out3);
}

void isvc_iquant_itrans_recon_chroma_4x4_with_res_accumulate_neon(
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

    WORD16 i2_rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x8_t x0_16x8, x1_16x8, x2_16x8, x3_16x8;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4x2_t xx0_16x4x2, xx1_16x4x2;
    int32x2x2_t x0_32x2x2, x1_32x2x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;

    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t rec0, rec1, rec2, rec3;
    uint8x8_t rec0_un, rec1_un, rec2_un, rec3_un;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd1_in_mask, resd2_in_mask, resd3_in_mask;
    uint8x8_t out0, out1, out2, out3;
    int16x8_t resout0, resout1, resout2, resout3;
    int16x8_t pos_255 = vdupq_n_s16(((WORD16) UINT8_MAX));
    int16x8_t neg_255 = vdupq_n_s16(-((WORD16) UINT8_MAX));

    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));
    uint16x8_t chroma_mask_16x8 = {0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000};

    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);
    int32x4_t rnd_fact = vdupq_n_s32(i2_rnd_factor);

    int16x8_t resd0_in_mask = {0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000};

    UNUSED(i4_iq_start_idx);
    UNUSED(u1_res_accumulate);

    resd1_in_mask = resd0_in_mask;
    resd2_in_mask = resd0_in_mask;
    resd3_in_mask = resd0_in_mask;

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    q0_32x4 = vmull_s16(weigh0_16x4, src_16x4x2.val[0]);
    q1_32x4 = vmull_s16(weigh1_16x4, src_16x4x2.val[1]);
    q2_32x4 = vmull_s16(weigh2_16x4, src_16x4x2.val[2]);
    q3_32x4 = vmull_s16(weigh3_16x4, src_16x4x2.val[3]);

    q0_32x4 = vaddq_s32(q0_32x4, rnd_fact);
    q1_32x4 = vaddq_s32(q1_32x4, rnd_fact);
    q2_32x4 = vaddq_s32(q2_32x4, rnd_fact);
    q3_32x4 = vaddq_s32(q3_32x4, rnd_fact);

    q0_32x4 = vshlq_s32(q0_32x4, qp_div_6_32x4);
    q1_32x4 = vshlq_s32(q1_32x4, qp_div_6_32x4);
    q2_32x4 = vshlq_s32(q2_32x4, qp_div_6_32x4);
    q3_32x4 = vshlq_s32(q3_32x4, qp_div_6_32x4);

    q0_16x4 = vqshrn_n_s32(q0_32x4, 4);
    q1_16x4 = vqshrn_n_s32(q1_32x4, 4);
    q2_16x4 = vqshrn_n_s32(q2_32x4, 4);
    q3_16x4 = vqshrn_n_s32(q3_32x4, 4);

    q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);

    /* row 0 to row 3 */
    xx0_16x4x2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx1_16x4x2 = vtrn_s16(xx2_16x4, xx3_16x4);
    x0_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[0]), vreinterpret_s32_s16(xx1_16x4x2.val[0]));
    x1_32x2x2 =
        vzip_s32(vreinterpret_s32_s16(xx0_16x4x2.val[1]), vreinterpret_s32_s16(xx1_16x4x2.val[1]));

    x0_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[0]);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[0]);
    x2_16x4 = vreinterpret_s16_s32(x0_32x2x2.val[1]);
    x3_16x4 = vreinterpret_s16_s32(x1_32x2x2.val[1]);

    /* Store Horz transform output into temp */
    vst1_s16(pi2_tmp, x0_16x4);
    vst1_s16(pi2_tmp + 4, x1_16x4);
    vst1_s16(pi2_tmp + 8, x2_16x4);
    vst1_s16(pi2_tmp + 12, x3_16x4);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    resd0_in = vld1q_s16((int16_t *) pi2_res_pred);
    resd1_in = vld1q_s16((int16_t *) pi2_res_pred + i4_res_pred_stride);
    resd2_in = vld1q_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 3));

    /* Mask alternate values */
    resd0_in_mask = vbslq_s16(chroma_mask_16x8, resd0_in, resd0_in_mask);
    resd1_in_mask = vbslq_s16(chroma_mask_16x8, resd1_in, resd1_in_mask);
    resd2_in_mask = vbslq_s16(chroma_mask_16x8, resd2_in, resd2_in_mask);
    resd3_in_mask = vbslq_s16(chroma_mask_16x8, resd3_in, resd3_in_mask);

    x0_16x8 = vreinterpretq_s16_s32(vmovl_s16(x0_16x4));
    x1_16x8 = vreinterpretq_s16_s32(vmovl_s16(x1_16x4));
    x2_16x8 = vreinterpretq_s16_s32(vmovl_s16(x2_16x4));
    x3_16x8 = vreinterpretq_s16_s32(vmovl_s16(x3_16x4));

    resd0_in = vaddq_s16(resd0_in_mask, x0_16x8);
    resd1_in = vaddq_s16(resd1_in_mask, x1_16x8);
    resd2_in = vaddq_s16(resd2_in_mask, x2_16x8);
    resd3_in = vaddq_s16(resd3_in_mask, x3_16x8);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    resd0_in = vmaxq_s16(resd0_in, neg_255);
    resd1_in = vmaxq_s16(resd1_in, neg_255);
    resd2_in = vmaxq_s16(resd2_in, neg_255);
    resd3_in = vmaxq_s16(resd3_in, neg_255);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    resd0_in = vminq_s16(resd0_in, pos_255);
    resd1_in = vminq_s16(resd1_in, pos_255);
    resd2_in = vminq_s16(resd2_in, pos_255);
    resd3_in = vminq_s16(resd3_in, pos_255);

    resout0 = vld1q_s16(pi2_res);
    resout1 = vld1q_s16(pi2_res + i4_res_stride);
    resout2 = vld1q_s16(pi2_res + i4_res_stride * 2);
    resout3 = vld1q_s16(pi2_res + i4_res_stride * 3);

    /* Store res in aternate positions */
    resout0 = vbslq_s16(chroma_mask_16x8, resd0_in, resout0);
    resout1 = vbslq_s16(chroma_mask_16x8, resd1_in, resout1);
    resout2 = vbslq_s16(chroma_mask_16x8, resd2_in, resout2);
    resout3 = vbslq_s16(chroma_mask_16x8, resd3_in, resout3);

    vst1q_s16(pi2_res, resout0);
    vst1q_s16(pi2_res + i4_res_stride, resout1);
    vst1q_s16(pi2_res + (i4_res_stride << 1), resout2);
    vst1q_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resout3);

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride << 1));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (i4_pred_stride * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = pred + res */
    rec0 = vaddq_s16(pred0, resout0);
    rec1 = vaddq_s16(pred1, resout1);
    rec2 = vaddq_s16(pred2, resout2);
    rec3 = vaddq_s16(pred3, resout3);

    out0 = vld1_u8(pu1_out);
    out1 = vld1_u8(pu1_out + i4_out_stride);
    out2 = vld1_u8(pu1_out + i4_out_stride * 2);
    out3 = vld1_u8(pu1_out + i4_out_stride * 3);

    /* Convert to 8 bit unsigned with saturation */
    rec0_un = vqmovun_s16(rec0);
    rec1_un = vqmovun_s16(rec1);
    rec2_un = vqmovun_s16(rec2);
    rec3_un = vqmovun_s16(rec3);

    /* Store output pixels in alternate positions */
    out0 = vbsl_u8(chroma_mask_8x8, rec0_un, out0);
    out1 = vbsl_u8(chroma_mask_8x8, rec1_un, out1);
    out2 = vbsl_u8(chroma_mask_8x8, rec2_un, out2);
    out3 = vbsl_u8(chroma_mask_8x8, rec3_un, out3);

    vst1_u8((pu1_out), out0);
    vst1_u8((pu1_out + i4_out_stride), out1);
    vst1_u8((pu1_out + (i4_out_stride << 1)), out2);
    vst1_u8((pu1_out + ((i4_out_stride << 1) + i4_out_stride)), out3);
}

void isvc_iquant_itrans_recon_4x4_dc_neon(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                          buffer_container_t *ps_res_pred,
                                          buffer_container_t *ps_res, buffer_container_t *ps_rec,
                                          iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
                                          WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
                                          WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    WORD16 rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    WORD32 i4_iq_out_temp;
    int16x8_t temp_0;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;

    UNUSED(pi2_tmp);
    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    if(i4_iq_start_idx == 0)
    {
        i4_iq_out_temp = pi2_src[0];
        INV_QUANT(i4_iq_out_temp, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);
    }
    else
    {
        i4_iq_out_temp = pi2_dc_src[0];
    }

    temp_0 = vdupq_n_s16((i4_iq_out_temp + 32) >> 6);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = Res + pred */
    pred0 = vaddq_s16(pred0, temp_0);
    pred1 = vaddq_s16(pred1, temp_0);
    pred2 = vaddq_s16(pred2, temp_0);
    pred3 = vaddq_s16(pred3, temp_0);

    /* Convert to unsigned 8 bit with saturation */
    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred0_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred1_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride * 2), vreinterpret_u32_u8(pred2_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride * 3), vreinterpret_u32_u8(pred3_in), 0);
}

void isvc_iquant_itrans_recon_4x4_dc_with_res_output_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    WORD16 rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    WORD16 i2_it_out;
    WORD32 i4_iq_out_temp;
    int16x8_t temp_0;
    int16x4_t residue_res;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;

    UNUSED(pi2_tmp);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    if(i4_iq_start_idx == 0)
    {
        i4_iq_out_temp = pi2_src[0];
        INV_QUANT(i4_iq_out_temp, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);
    }
    else
    {
        i4_iq_out_temp = pi2_dc_src[0];
    }

    i2_it_out = ((i4_iq_out_temp + 32) >> 6);
    temp_0 = vdupq_n_s16(i2_it_out);
    residue_res = vdup_n_s16(isvc_get_residue(i2_it_out, 0, 0));

    vst1_s16(pi2_res, residue_res);
    vst1_s16(pi2_res + i4_res_stride, residue_res);
    vst1_s16(pi2_res + (i4_res_stride << 1), residue_res);
    vst1_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, residue_res);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = Res + pred */
    pred0 = vaddq_s16(pred0, temp_0);
    pred1 = vaddq_s16(pred1, temp_0);
    pred2 = vaddq_s16(pred2, temp_0);
    pred3 = vaddq_s16(pred3, temp_0);

    /* Convert to unsigned 8 bit with saturation */
    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred0_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred1_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride * 2), vreinterpret_u32_u8(pred2_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride * 3), vreinterpret_u32_u8(pred3_in), 0);
}

void isvc_iquant_itrans_recon_4x4_dc_with_res_accumulate_neon(
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
    WORD16 rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;

    WORD32 i4_iq_out_temp;
    int16x4_t temp_0;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t pred01_in, pred23_in;
    uint8x8_t pred01_un, pred23_un;

    int16x4_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd01_in, resd23_in;
    int16x4_t pos_255 = vdup_n_s16(((WORD16) UINT8_MAX));
    int16x4_t neg_255 = vdup_n_s16(-((WORD16) UINT8_MAX));

    UNUSED(pi2_tmp);
    UNUSED(u1_res_accumulate);

    if(i4_iq_start_idx == 0)
    {
        i4_iq_out_temp = pi2_src[0];
        INV_QUANT(i4_iq_out_temp, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);
    }
    else
    {
        i4_iq_out_temp = pi2_dc_src[0];
    }

    temp_0 = vdup_n_s16((i4_iq_out_temp + 32) >> 6);

    resd0_in = vld1_s16((int16_t *) pi2_res_pred);
    resd1_in = vld1_s16((int16_t *) pi2_res_pred + i4_res_pred_stride);
    resd2_in = vld1_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 2));
    resd3_in = vld1_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 3));

    /* Add res pred to the res obtained */
    resd0_in = vadd_s16(resd0_in, temp_0);
    resd1_in = vadd_s16(resd1_in, temp_0);
    resd2_in = vadd_s16(resd2_in, temp_0);
    resd3_in = vadd_s16(resd3_in, temp_0);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    resd0_in = vmax_s16(resd0_in, neg_255);
    resd1_in = vmax_s16(resd1_in, neg_255);
    resd2_in = vmax_s16(resd2_in, neg_255);
    resd3_in = vmax_s16(resd3_in, neg_255);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    resd0_in = vmin_s16(resd0_in, pos_255);
    resd1_in = vmin_s16(resd1_in, pos_255);
    resd2_in = vmin_s16(resd2_in, pos_255);
    resd3_in = vmin_s16(resd3_in, pos_255);

    vst1_s16(pi2_res, resd0_in);
    vst1_s16(pi2_res + i4_res_stride, resd1_in);
    vst1_s16(pi2_res + (i4_res_stride << 1), resd2_in);
    vst1_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resd3_in);

    resd01_in = vcombine_s16(resd0_in, resd1_in);
    resd23_in = vcombine_s16(resd2_in, resd3_in);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    pred01_in = vcombine_s16(vget_low_s16(pred0), vget_low_s16(pred1));
    pred23_in = vcombine_s16(vget_low_s16(pred2), vget_low_s16(pred3));

    /* Out pixel = Res + pred */
    pred01_in = vaddq_s16(pred01_in, resd01_in);
    pred23_in = vaddq_s16(pred23_in, resd23_in);

    /* Convert to unsigned 8 bit with saturation */
    pred01_un = vqmovun_s16(pred01_in);
    pred23_un = vqmovun_s16(pred23_in);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred01_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + i4_out_stride), vreinterpret_u32_u8(pred01_un), 1);
    vst1_lane_u32((uint32_t *) (pu1_out + (i4_out_stride << 1)), vreinterpret_u32_u8(pred23_un), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + ((i4_out_stride << 1) + i4_out_stride)),
                  vreinterpret_u32_u8(pred23_un), 1);
}

void isvc_iquant_itrans_recon_chroma_4x4_dc_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    WORD32 i4_iq_out_temp;
    int16x8_t temp_0;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    uint8x8_t i4_out_horz_8x8_r0, i4_out_horz_8x8_r1, i4_out_horz_8x8_r2, i4_out_horz_8x8_r3;
    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));

    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(i4_iq_start_idx);
    UNUSED(ps_res);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    i4_iq_out_temp = pi2_dc_src[0];
    temp_0 = vdupq_n_s16((i4_iq_out_temp + 32) >> 6);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = Res + pred */
    pred0 = vaddq_s16(pred0, temp_0);
    pred1 = vaddq_s16(pred1, temp_0);
    pred2 = vaddq_s16(pred2, temp_0);
    pred3 = vaddq_s16(pred3, temp_0);

    /* Convert to unsigned 8 bit with saturation */
    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    i4_out_horz_8x8_r0 = vld1_u8(pu1_out);
    i4_out_horz_8x8_r1 = vld1_u8(pu1_out + i4_out_stride);
    i4_out_horz_8x8_r2 = vld1_u8(pu1_out + i4_out_stride * 2);
    i4_out_horz_8x8_r3 = vld1_u8(pu1_out + i4_out_stride * 3);

    /* Store out pixels in alternate positions */
    i4_out_horz_8x8_r0 = vbsl_u8(chroma_mask_8x8, pred0_in, i4_out_horz_8x8_r0);
    i4_out_horz_8x8_r1 = vbsl_u8(chroma_mask_8x8, pred1_in, i4_out_horz_8x8_r1);
    i4_out_horz_8x8_r2 = vbsl_u8(chroma_mask_8x8, pred2_in, i4_out_horz_8x8_r2);
    i4_out_horz_8x8_r3 = vbsl_u8(chroma_mask_8x8, pred3_in, i4_out_horz_8x8_r3);

    vst1_u8((uint8_t *) (pu1_out), i4_out_horz_8x8_r0);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride), i4_out_horz_8x8_r1);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 2), i4_out_horz_8x8_r2);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 3), i4_out_horz_8x8_r3);
}

void isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_output_neon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_res_pred,
    buffer_container_t *ps_res, buffer_container_t *ps_rec,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp, WORD16 *pi2_dc_src,
    WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate)
{
    WORD16 *pi2_src = (WORD16 *) ps_src->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    UWORD8 *pu1_out = (UWORD8 *) ps_rec->pv_data;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_out_stride = ps_rec->i4_data_stride;
    const UWORD16 *pu2_iscal_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;

    WORD16 i2_it_out;
    WORD32 i4_iq_out_temp;
    int16x8_t temp_0, residue_res;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resout0, resout1, resout2, resout3;

    uint8x8_t i4_out_horz_8x8_r0, i4_out_horz_8x8_r1, i4_out_horz_8x8_r2, i4_out_horz_8x8_r3;
    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));
    uint16x8_t chroma_mask_16x8 = {0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000};

    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(i4_iq_start_idx);
    UNUSED(ps_res_pred);
    UNUSED(u1_res_accumulate);

    i4_iq_out_temp = pi2_dc_src[0];

    i2_it_out = ((i4_iq_out_temp + 32) >> 6);
    temp_0 = vdupq_n_s16(i2_it_out);
    residue_res = vdupq_n_s16(isvc_get_residue(i2_it_out, 0, 0));

    resout0 = vld1q_s16(pi2_res);
    resout1 = vld1q_s16(pi2_res + i4_res_stride);
    resout2 = vld1q_s16(pi2_res + i4_res_stride * 2);
    resout3 = vld1q_s16(pi2_res + i4_res_stride * 3);

    /* Store res in alternate positions */
    resout0 = vbslq_s16(chroma_mask_16x8, residue_res, resout0);
    resout1 = vbslq_s16(chroma_mask_16x8, residue_res, resout1);
    resout2 = vbslq_s16(chroma_mask_16x8, residue_res, resout2);
    resout3 = vbslq_s16(chroma_mask_16x8, residue_res, resout3);

    vst1q_s16(pi2_res, resout0);
    vst1q_s16(pi2_res + i4_res_stride, resout1);
    vst1q_s16(pi2_res + (i4_res_stride << 1), resout2);
    vst1q_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resout3);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = Res + pred */
    pred0 = vaddq_s16(pred0, temp_0);
    pred1 = vaddq_s16(pred1, temp_0);
    pred2 = vaddq_s16(pred2, temp_0);
    pred3 = vaddq_s16(pred3, temp_0);

    /* Convert to unsigned 8 bit with saturation */
    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    /* Store out pixels in alternate positions */
    i4_out_horz_8x8_r0 = vld1_u8(pu1_out);
    i4_out_horz_8x8_r1 = vld1_u8(pu1_out + i4_out_stride);
    i4_out_horz_8x8_r2 = vld1_u8(pu1_out + i4_out_stride * 2);
    i4_out_horz_8x8_r3 = vld1_u8(pu1_out + i4_out_stride * 3);

    i4_out_horz_8x8_r0 = vbsl_u8(chroma_mask_8x8, pred0_in, i4_out_horz_8x8_r0);
    i4_out_horz_8x8_r1 = vbsl_u8(chroma_mask_8x8, pred1_in, i4_out_horz_8x8_r1);
    i4_out_horz_8x8_r2 = vbsl_u8(chroma_mask_8x8, pred2_in, i4_out_horz_8x8_r2);
    i4_out_horz_8x8_r3 = vbsl_u8(chroma_mask_8x8, pred3_in, i4_out_horz_8x8_r3);

    vst1_u8((uint8_t *) (pu1_out), i4_out_horz_8x8_r0);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride), i4_out_horz_8x8_r1);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 2), i4_out_horz_8x8_r2);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 3), i4_out_horz_8x8_r3);
}

void isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_accumulate_neon(
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

    WORD32 i4_iq_out_temp;
    int16x8_t temp_0;
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resout0, resout1, resout2, resout3;
    int16x8_t resd1_in_mask, resd2_in_mask, resd3_in_mask;
    uint8x8_t out0, out1, out2, out3;
    int16x8_t pos_255 = vdupq_n_s16(((WORD16) UINT8_MAX));
    int16x8_t neg_255 = vdupq_n_s16(-((WORD16) UINT8_MAX));
    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));
    uint16x8_t chroma_mask_16x8 = {0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000, 0xffff, 0x0000};

    int16x8_t resd0_in_mask = {0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000};

    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(i4_iq_start_idx);
    UNUSED(u1_res_accumulate);

    resd1_in_mask = resd0_in_mask;
    resd2_in_mask = resd0_in_mask;
    resd3_in_mask = resd0_in_mask;

    i4_iq_out_temp = pi2_dc_src[0];
    temp_0 = vdupq_n_s16((i4_iq_out_temp + 32) >> 6);

    resd0_in = vld1q_s16((int16_t *) pi2_res_pred);
    resd1_in = vld1q_s16((int16_t *) pi2_res_pred + i4_res_pred_stride);
    resd2_in = vld1q_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_res_pred + (i4_res_pred_stride * 3));

    /* Mask alternate values of res pred */
    resd0_in_mask = vbslq_s16(chroma_mask_16x8, resd0_in, resd0_in_mask);
    resd1_in_mask = vbslq_s16(chroma_mask_16x8, resd1_in, resd1_in_mask);
    resd2_in_mask = vbslq_s16(chroma_mask_16x8, resd2_in, resd2_in_mask);
    resd3_in_mask = vbslq_s16(chroma_mask_16x8, resd3_in, resd3_in_mask);

    /* Add res pred to res obtained */
    resd0_in = vaddq_s16(resd0_in_mask, temp_0);
    resd1_in = vaddq_s16(resd1_in_mask, temp_0);
    resd2_in = vaddq_s16(resd2_in_mask, temp_0);
    resd3_in = vaddq_s16(resd3_in_mask, temp_0);

    /* Saturate all values < -255 to -255 and retain the rest as it is */
    resd0_in = vmaxq_s16(resd0_in, neg_255);
    resd1_in = vmaxq_s16(resd1_in, neg_255);
    resd2_in = vmaxq_s16(resd2_in, neg_255);
    resd3_in = vmaxq_s16(resd3_in, neg_255);

    /* Saturate all values > 255 to 255 and retain the rest as it is */
    resd0_in = vminq_s16(resd0_in, pos_255);
    resd1_in = vminq_s16(resd1_in, pos_255);
    resd2_in = vminq_s16(resd2_in, pos_255);
    resd3_in = vminq_s16(resd3_in, pos_255);

    resout0 = vld1q_s16(pi2_res);
    resout1 = vld1q_s16(pi2_res + i4_res_stride);
    resout2 = vld1q_s16(pi2_res + i4_res_stride * 2);
    resout3 = vld1q_s16(pi2_res + i4_res_stride * 3);

    /* Store res in alternate positions */
    resout0 = vbslq_s16(chroma_mask_16x8, resd0_in, resout0);
    resout1 = vbslq_s16(chroma_mask_16x8, resd1_in, resout1);
    resout2 = vbslq_s16(chroma_mask_16x8, resd2_in, resout2);
    resout3 = vbslq_s16(chroma_mask_16x8, resd3_in, resout3);

    vst1q_s16(pi2_res, resout0);
    vst1q_s16(pi2_res + i4_res_stride, resout1);
    vst1q_s16(pi2_res + (i4_res_stride << 1), resout2);
    vst1q_s16(pi2_res + (i4_res_stride << 1) + i4_res_stride, resout3);

    pred0_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred1_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred2_in = vld1_u8(pu1_pred);
    pu1_pred = pu1_pred + i4_pred_stride;
    pred3_in = vld1_u8(pu1_pred);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    /* Out pixel = Res + pred */
    pred0 = vaddq_s16(pred0, resout0);
    pred1 = vaddq_s16(pred1, resout1);
    pred2 = vaddq_s16(pred2, resout2);
    pred3 = vaddq_s16(pred3, resout3);

    /* Convert to unsigned 8 bit with saturation */
    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    out0 = vld1_u8(pu1_out);
    out1 = vld1_u8(pu1_out + i4_out_stride);
    out2 = vld1_u8(pu1_out + i4_out_stride * 2);
    out3 = vld1_u8(pu1_out + i4_out_stride * 3);

    /* Store out pixels in alternate positions */
    out0 = vbsl_u8(chroma_mask_8x8, pred0_in, out0);
    out1 = vbsl_u8(chroma_mask_8x8, pred1_in, out1);
    out2 = vbsl_u8(chroma_mask_8x8, pred2_in, out2);
    out3 = vbsl_u8(chroma_mask_8x8, pred3_in, out3);

    vst1_u8((uint8_t *) (pu1_out), out0);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride), out1);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 2), out2);
    vst1_u8((uint8_t *) (pu1_out + i4_out_stride * 3), out3);
}
