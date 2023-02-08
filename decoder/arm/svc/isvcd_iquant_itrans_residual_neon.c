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
 *  isvcd_iquant_itrans_residual_neonintr.c
 *
 * @brief
 *  Contains definition of functions for h264 inverse quantization inverse
 *    transformation and resd comp
 *
 * @author
 *  Kishore
 *
 *  @par List of Functions:
 *  - isvcd_iquant_itrans_residual_4x4_neonintr()
 *  - isvcd_iquant_itrans_residual_8x8_neonintr()
 *  - isvcd_iquant_itrans_residual_4x4_dc_neonintr()
 *  - isvcd_iquant_itrans_residual_8x8_dc_neonintr()
 *  - isvcd_iquant_itrans_residual_chroma_4x4_neonintr()
 *  - isvcd_iquant_itrans_residual_chroma_4x4_dc_neonintr()
 *
 * @remarks
 *
 *******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

#include <string.h>
#include <arm_neon.h>

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "ih264_structs.h"
#include "isvcd_iquant_itrans_residual.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_4x4_neonintr                 */
/*                                                                           */
/*  Description   : this function computes the resd output from the          */
/*                  IQ+IT                                                    */
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

WORD32 isvcd_iquant_itrans_residual_4x4_neonintr(WORD16 *pi2_src, WORD16 *pi2_pred, WORD16 *pi2_out,
                                                 WORD32 pred_strd, WORD32 out_strd,
                                                 const UWORD16 *pu2_iscal_mat,
                                                 const UWORD16 *pu2_weigh_mat, UWORD32 u4_qp_div_6,
                                                 WORD16 *pi2_tmp, WORD32 iq_start_idx,
                                                 WORD16 *pi2_dc_ld_addr)
{
    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    WORD16 rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    int32x4_t rnd_fact = vdupq_n_s32(rnd_factor);

    int16x4_t pred0, pred1, pred2, pred3;
    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4_t xx0_0_16x4, xx0_1_16x4, xx2_0_16x4, xx2_1_16x4;
    int32x2_t x0_32x2, x1_32x2, x2_32x2, x3_32x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4, dup_min, dup_max;
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);

    int16x8_t pred01_in, pred23_in;
    WORD32 i4_nnz;
    int16x4x2_t xx0_16x4_2, xx2_16x4_2;
    int32x2x2_t x0_32x2_2, x1_32x2_2;
    int16x8_t dup_val_1, dup_val_2, dup_abs;
    UNUSED(pi2_tmp);

    src_16x4x2 = vld4_s16(pi2_src);
    iscal_16x4x2 = vld4_s16((const int16_t *) pu2_iscal_mat);
    weigh_16x4x2 = vld4_s16((const int16_t *) pu2_weigh_mat);

    weigh0_16x4 = vmul_s16(weigh_16x4x2.val[0], iscal_16x4x2.val[0]);
    weigh1_16x4 = vmul_s16(weigh_16x4x2.val[1], iscal_16x4x2.val[1]);
    weigh2_16x4 = vmul_s16(weigh_16x4x2.val[2], iscal_16x4x2.val[2]);
    weigh3_16x4 = vmul_s16(weigh_16x4x2.val[3], iscal_16x4x2.val[3]);

    dup_min = vdup_n_s16(RSD_MIN);
    dup_max = vdup_n_s16(RSD_MAX);

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

    if(iq_start_idx == 1)
    {
        q0_16x4 = vset_lane_s16(pi2_dc_ld_addr[0], q0_16x4, 0);
    }

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);      // q1 >>1
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);      // q3 >>1

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);   // x0 =  q0 + q2
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);   // x1 =  q0 - q2
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);  // x2 =  q1>>1 - q3
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);  // x2 =  q1 + q3>>1

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);  // x0+x3
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);  // x1+x2
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);  // x1-x2
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);  // x0-x3

    xx0_16x4_2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx0_0_16x4 = xx0_16x4_2.val[0];
    xx0_1_16x4 = xx0_16x4_2.val[1];
    xx2_16x4_2 = vtrn_s16(xx2_16x4, xx3_16x4);
    xx2_0_16x4 = xx2_16x4_2.val[0];
    xx2_1_16x4 = xx2_16x4_2.val[1];
    x0_32x2_2 = vtrn_s32(vreinterpret_s32_s16(xx0_0_16x4), vreinterpret_s32_s16(xx2_0_16x4));
    x1_32x2_2 = vtrn_s32(vreinterpret_s32_s16(xx0_1_16x4), vreinterpret_s32_s16(xx2_1_16x4));
    x0_32x2 = x0_32x2_2.val[0];
    x1_32x2 = x1_32x2_2.val[0];
    x2_32x2 = x0_32x2_2.val[1];
    x3_32x2 = x1_32x2_2.val[1];

    x0_16x4 = vreinterpret_s16_s32(x0_32x2);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2);
    x2_16x4 = vreinterpret_s16_s32(x2_32x2);
    x3_16x4 = vreinterpret_s16_s32(x3_32x2);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);       // q1 >> 1
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);       // q3 >> 1

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);   // x0 =  q0 + q2
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);   // x1 =  q0 - q2
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);  // x2 =  q1>>1 - q3
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);  // x3 =  q1 + q3>>1

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);  // imacro = x0 + x3
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);  // imacro = x1 + x2
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);  // imacro = x1 - x2
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);  // imacro = x0 - x3

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    pred0 = vld1_s16((int16_t *) pi2_pred);
    pred1 = vld1_s16((int16_t *) pi2_pred + pred_strd);
    pred2 = vld1_s16((int16_t *) pi2_pred + (pred_strd * 2));
    pred3 = vld1_s16((int16_t *) pi2_pred + (pred_strd * 3));

    x0_16x4 = vadd_s16(pred0, x0_16x4);
    x1_16x4 = vadd_s16(pred1, x1_16x4);
    x2_16x4 = vadd_s16(pred2, x2_16x4);
    x3_16x4 = vadd_s16(pred3, x3_16x4);

    x0_16x4 = vmin_s16(x0_16x4, dup_max);
    x0_16x4 = vmax_s16(x0_16x4, dup_min);
    x1_16x4 = vmin_s16(x1_16x4, dup_max);
    x1_16x4 = vmax_s16(x1_16x4, dup_min);
    x2_16x4 = vmin_s16(x2_16x4, dup_max);
    x2_16x4 = vmax_s16(x2_16x4, dup_min);
    x3_16x4 = vmin_s16(x3_16x4, dup_max);
    x3_16x4 = vmax_s16(x3_16x4, dup_min);

    pred01_in = vcombine_s16(x0_16x4, x1_16x4);
    pred23_in = vcombine_s16(x2_16x4, x3_16x4);

    dup_val_1 = vabsq_s16(pred01_in);
    dup_val_2 = vabsq_s16(pred23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    vst1_s16(pi2_out, x0_16x4);
    vst1_s16(pi2_out + out_strd, x1_16x4);
    vst1_s16(pi2_out + (out_strd << 1), x2_16x4);
    vst1_s16(pi2_out + ((out_strd << 1) + out_strd), x3_16x4);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_4x4_dc_neonintr              */
/*                                                                           */
/*  Description   : this function computes the resd output from the          */
/*                  IQ+IT                                                    */
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
WORD32 isvcd_iquant_itrans_residual_4x4_dc_neonintr(WORD16 *pi2_src, WORD16 *pi2_pred,
                                                    WORD16 *pi2_out, WORD32 pred_strd,
                                                    WORD32 out_strd, const UWORD16 *pu2_iscal_mat,
                                                    const UWORD16 *pu2_weigh_mat,
                                                    UWORD32 u4_qp_div_6, WORD16 *pi2_tmp,
                                                    WORD32 iq_start_idx, WORD16 *pi2_dc_ld_addr)
{
    WORD32 i4_iq_out_temp;
    int16x4_t temp_0;
    int16x4_t pred0_in, pred1_in, pred2_in, pred3_in, dup_min, dup_max;
    int16x8_t pred01_in, pred23_in;
    WORD32 i4_nnz;
    int16x8_t dup_val_1, dup_val_2, dup_abs;

    WORD16 rnd_fact = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    UNUSED(pi2_tmp);

    if(iq_start_idx == 0)
    {
        i4_iq_out_temp = pi2_src[0];
        INV_QUANT(i4_iq_out_temp, pu2_iscal_mat[0], pu2_weigh_mat[0], u4_qp_div_6, rnd_fact, 4);
    }
    else
    {
        i4_iq_out_temp = pi2_dc_ld_addr[0];
    }

    /* inv transform and residual comp */
    pred0_in = vld1_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred1_in = vld1_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred2_in = vld1_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred3_in = vld1_s16((int16_t *) pi2_pred);

    temp_0 = vdup_n_s16((i4_iq_out_temp + 32) >> 6);
    dup_min = vdup_n_s16(RSD_MIN);
    dup_max = vdup_n_s16(RSD_MAX);

    pred0_in = vadd_s16(pred0_in, temp_0);
    pred1_in = vadd_s16(pred1_in, temp_0);
    pred2_in = vadd_s16(pred2_in, temp_0);
    pred3_in = vadd_s16(pred3_in, temp_0);

    pred0_in = vmin_s16(pred0_in, dup_max);
    pred0_in = vmax_s16(pred0_in, dup_min);
    pred1_in = vmin_s16(pred1_in, dup_max);
    pred1_in = vmax_s16(pred1_in, dup_min);
    pred2_in = vmin_s16(pred2_in, dup_max);
    pred2_in = vmax_s16(pred2_in, dup_min);
    pred3_in = vmin_s16(pred3_in, dup_max);
    pred3_in = vmax_s16(pred3_in, dup_min);

    pred01_in = vcombine_s16(pred0_in, pred1_in);
    pred23_in = vcombine_s16(pred2_in, pred3_in);

    dup_val_1 = vabsq_s16(pred01_in);
    dup_val_2 = vabsq_s16(pred23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    vst1_s16((int16_t *) (pi2_out), pred0_in);
    vst1_s16((int16_t *) (pi2_out + out_strd), pred1_in);
    vst1_s16((int16_t *) (pi2_out + (out_strd * 2)), pred2_in);
    vst1_s16((int16_t *) (pi2_out + (out_strd * 3)), pred3_in);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_chroma_4x4_neonintr          */
/*                                                                           */
/*  Description   : this function computes the resd output from the          */
/*                  IQ+IT                                                    */
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

WORD32 isvcd_iquant_itrans_residual_chroma_4x4_neonintr(
    WORD16 *pi2_src, WORD16 *pi2_pred, WORD16 *pi2_out, WORD32 pred_strd, WORD32 out_strd,
    const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat, UWORD32 u4_qp_div_6,
    WORD16 *pi2_tmp, WORD16 *pi2_dc_src)
{
    int16x4x4_t src_16x4x2;
    int16x4x4_t iscal_16x4x2;
    int16x4x4_t weigh_16x4x2;

    int16x8_t pred0, pred1, pred2, pred3;
    int16x4_t q0_16x4, q1_16x4, q2_16x4, q3_16x4;
    int32x4_t q0_32x4, q1_32x4, q2_32x4, q3_32x4;
    int16x4_t rq1_16x4, rq3_16x4;
    int16x4_t x0_16x4, x1_16x4, x2_16x4, x3_16x4;
    int16x4_t xx0_16x4, xx1_16x4, xx2_16x4, xx3_16x4;
    int16x4_t xx0_0_16x4, xx0_1_16x4, xx2_0_16x4, xx2_1_16x4;
    int32x2_t x0_32x2, x1_32x2, x2_32x2, x3_32x2;
    int16x4_t weigh0_16x4, weigh1_16x4, weigh2_16x4, weigh3_16x4;
    int32x4_t qp_div_6_32x4 = vdupq_n_s32(u4_qp_div_6);
    int16x4_t zero_16x4 = vdup_n_s16(0);
    int16x4_t x_16x4_low, x_16x4_high;
    int16x8_t x0_16x8, x1_16x8, x2_16x8, x3_16x8;
    int16x8_t dup_min, dup_max;
    WORD32 i4_nnz;
    int16x4x2_t xx0_16x4_2, xx2_16x4_2, x_16x4x2_t;
    int32x2x2_t x0_32x2_2, x1_32x2_2;
    int16x8_t dup_val_1, dup_val_2, dup_val_3, dup_val_4, dup_val_5, dup_abs;
    int16x8_t i4_out_horz_16x8_r0, i4_out_horz_16x8_r1, i4_out_horz_16x8_r2, i4_out_horz_16x8_r3;
    uint16x8_t chroma_mask_16x8 = vreinterpretq_u16_u32(vdupq_n_u32(0x0000ffff));
    int16x8_t chroma_mask_16x8_2 = vreinterpretq_s16_s32(vdupq_n_s32(0x0000ffff));

    WORD16 rnd_factor = (u4_qp_div_6 < 4) ? 1 << (3 - u4_qp_div_6) : 0;
    int32x4_t rnd_fact = vdupq_n_s32(rnd_factor);
    UNUSED(pi2_tmp);

    dup_min = vdupq_n_s16(RSD_MIN);
    dup_max = vdupq_n_s16(RSD_MAX);

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

    rq1_16x4 = vshr_n_s16(q1_16x4, 1);  // q1 >>1
    rq3_16x4 = vshr_n_s16(q3_16x4, 1);  // q3 >>1

    q0_16x4 = vset_lane_s16(pi2_dc_src[0], q0_16x4, 0);

    x0_16x4 = vadd_s16(q0_16x4, q2_16x4);   // x0 =  q0 + q2
    x1_16x4 = vsub_s16(q0_16x4, q2_16x4);   // x1 =  q0 - q2
    x2_16x4 = vsub_s16(rq1_16x4, q3_16x4);  // x2 =  q1>>1 - q3
    x3_16x4 = vadd_s16(q1_16x4, rq3_16x4);  // x2 =  q1 + q3>>1

    xx0_16x4 = vadd_s16(x0_16x4, x3_16x4);  // x0+x3
    xx1_16x4 = vadd_s16(x1_16x4, x2_16x4);  // x1+x2
    xx2_16x4 = vsub_s16(x1_16x4, x2_16x4);  // x1-x2
    xx3_16x4 = vsub_s16(x0_16x4, x3_16x4);  // x0-x3

    xx0_16x4_2 = vtrn_s16(xx0_16x4, xx1_16x4);
    xx0_0_16x4 = xx0_16x4_2.val[0];
    xx0_1_16x4 = xx0_16x4_2.val[1];
    xx2_16x4_2 = vtrn_s16(xx2_16x4, xx3_16x4);
    xx2_0_16x4 = xx2_16x4_2.val[0];
    xx2_1_16x4 = xx2_16x4_2.val[1];
    x0_32x2_2 = vtrn_s32(vreinterpret_s32_s16(xx0_0_16x4), vreinterpret_s32_s16(xx2_0_16x4));
    x1_32x2_2 = vtrn_s32(vreinterpret_s32_s16(xx0_1_16x4), vreinterpret_s32_s16(xx2_1_16x4));
    x0_32x2 = x0_32x2_2.val[0];
    x1_32x2 = x1_32x2_2.val[0];
    x2_32x2 = x0_32x2_2.val[1];
    x3_32x2 = x1_32x2_2.val[1];

    x0_16x4 = vreinterpret_s16_s32(x0_32x2);
    x1_16x4 = vreinterpret_s16_s32(x1_32x2);
    x2_16x4 = vreinterpret_s16_s32(x2_32x2);
    x3_16x4 = vreinterpret_s16_s32(x3_32x2);

    /* vertical inverse transform */
    rq1_16x4 = vshr_n_s16(x1_16x4, 1);       // q1 >> 1
    rq3_16x4 = vshr_n_s16(x3_16x4, 1);       // q3 >> 1

    xx0_16x4 = vadd_s16(x0_16x4, x2_16x4);   // x0 =  q0 + q2
    xx1_16x4 = vsub_s16(x0_16x4, x2_16x4);   // x1 =  q0 - q2
    xx2_16x4 = vsub_s16(rq1_16x4, x3_16x4);  // x2 =  q1>>1 - q3
    xx3_16x4 = vadd_s16(x1_16x4, rq3_16x4);  // x3 =  q1 + q3>>1

    x0_16x4 = vadd_s16(xx0_16x4, xx3_16x4);  // imacro = x0 + x3
    x1_16x4 = vadd_s16(xx1_16x4, xx2_16x4);  // imacro = x1 + x2
    x2_16x4 = vsub_s16(xx1_16x4, xx2_16x4);  // imacro = x1 - x2
    x3_16x4 = vsub_s16(xx0_16x4, xx3_16x4);  // imacro = x0 - x3

    x0_16x4 = vrshr_n_s16(x0_16x4, 6);
    x1_16x4 = vrshr_n_s16(x1_16x4, 6);
    x2_16x4 = vrshr_n_s16(x2_16x4, 6);
    x3_16x4 = vrshr_n_s16(x3_16x4, 6);

    x_16x4x2_t = vzip_s16(x0_16x4, zero_16x4);
    x_16x4_low = x_16x4x2_t.val[0];
    x_16x4_high = x_16x4x2_t.val[1];
    x0_16x8 = vcombine_s16(x_16x4_low, x_16x4_high);

    x_16x4x2_t = vzip_s16(x1_16x4, zero_16x4);
    x_16x4_low = x_16x4x2_t.val[0];
    x_16x4_high = x_16x4x2_t.val[1];
    x1_16x8 = vcombine_s16(x_16x4_low, x_16x4_high);

    x_16x4x2_t = vzip_s16(x2_16x4, zero_16x4);
    x_16x4_low = x_16x4x2_t.val[0];
    x_16x4_high = x_16x4x2_t.val[1];
    x2_16x8 = vcombine_s16(x_16x4_low, x_16x4_high);

    x_16x4x2_t = vzip_s16(x3_16x4, zero_16x4);
    x_16x4_low = x_16x4x2_t.val[0];
    x_16x4_high = x_16x4x2_t.val[1];
    x3_16x8 = vcombine_s16(x_16x4_low, x_16x4_high);

    pred0 = vld1q_s16((int16_t *) pi2_pred);
    pred1 = vld1q_s16((int16_t *) pi2_pred + pred_strd);
    pred2 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 2));
    pred3 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 3));

    x0_16x8 = vaddq_s16(pred0, x0_16x8);
    x1_16x8 = vaddq_s16(pred1, x1_16x8);
    x2_16x8 = vaddq_s16(pred2, x2_16x8);
    x3_16x8 = vaddq_s16(pred3, x3_16x8);

    x0_16x8 = vminq_s16(x0_16x8, dup_max);
    x0_16x8 = vmaxq_s16(x0_16x8, dup_min);
    x1_16x8 = vminq_s16(x1_16x8, dup_max);
    x1_16x8 = vmaxq_s16(x1_16x8, dup_min);
    x2_16x8 = vminq_s16(x2_16x8, dup_max);
    x2_16x8 = vmaxq_s16(x2_16x8, dup_min);
    x3_16x8 = vminq_s16(x3_16x8, dup_max);
    x3_16x8 = vmaxq_s16(x3_16x8, dup_min);

    x0_16x8 = vandq_s16(x0_16x8, chroma_mask_16x8_2);
    x1_16x8 = vandq_s16(x1_16x8, chroma_mask_16x8_2);
    x2_16x8 = vandq_s16(x2_16x8, chroma_mask_16x8_2);
    x3_16x8 = vandq_s16(x3_16x8, chroma_mask_16x8_2);

    dup_val_1 = vabsq_s16(x0_16x8);
    dup_val_2 = vabsq_s16(x1_16x8);
    dup_val_3 = vabsq_s16(x2_16x8);
    dup_val_4 = vabsq_s16(x3_16x8);
    dup_val_5 = vqaddq_s16(dup_val_1, dup_val_2);
    dup_val_1 = vqaddq_s16(dup_val_3, dup_val_4);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_5);

    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    i4_out_horz_16x8_r0 = vld1q_s16(pi2_out);
    i4_out_horz_16x8_r1 = vld1q_s16(pi2_out + out_strd);
    i4_out_horz_16x8_r2 = vld1q_s16(pi2_out + out_strd * 2);
    i4_out_horz_16x8_r3 = vld1q_s16(pi2_out + out_strd * 3);

    i4_out_horz_16x8_r0 = vbslq_s16(chroma_mask_16x8, x0_16x8, i4_out_horz_16x8_r0);
    i4_out_horz_16x8_r1 = vbslq_s16(chroma_mask_16x8, x1_16x8, i4_out_horz_16x8_r1);
    i4_out_horz_16x8_r2 = vbslq_s16(chroma_mask_16x8, x2_16x8, i4_out_horz_16x8_r2);
    i4_out_horz_16x8_r3 = vbslq_s16(chroma_mask_16x8, x3_16x8, i4_out_horz_16x8_r3);

    vst1q_s16((int16_t *) (pi2_out), i4_out_horz_16x8_r0);
    vst1q_s16((int16_t *) (pi2_out + out_strd), i4_out_horz_16x8_r1);
    vst1q_s16((int16_t *) (pi2_out + out_strd * 2), i4_out_horz_16x8_r2);
    vst1q_s16((int16_t *) (pi2_out + out_strd * 3), i4_out_horz_16x8_r3);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_chroma_4x4_dc_neonintr       */
/*                                                                           */
/*  Description   : this function computes the resd output from the          */
/*                  IQ+IT                                                    */
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

WORD32 isvcd_iquant_itrans_residual_chroma_4x4_dc_neonintr(
    WORD16 *pi2_src, WORD16 *pi2_pred, WORD16 *pi2_out, WORD32 pred_strd, WORD32 out_strd,
    const UWORD16 *pu2_iscal_mat, const UWORD16 *pu2_weigh_mat, UWORD32 u4_qp_div_6,
    WORD16 *pi2_tmp, WORD16 *pi2_dc_src)
{
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t temp_0, dup_min, dup_max;
    WORD32 i4_iq_out_temp;
    int16x8_t dup_val_1, dup_val_2, dup_val_3, dup_val_4, dup_val_5, dup_abs;
    int16x8_t i4_out_horz_16x8_r0, i4_out_horz_16x8_r1, i4_out_horz_16x8_r2, i4_out_horz_16x8_r3;
    uint16x8_t chroma_mask_16x8 = vreinterpretq_u16_u32(vdupq_n_u32(0x0000ffff));
    int16x8_t chroma_mask_16x8_2 = vreinterpretq_s16_s32(vdupq_n_s32(0x0000ffff));

    WORD32 i4_nnz;
    UNUSED(pi2_src);
    UNUSED(pu2_iscal_mat);
    UNUSED(pu2_weigh_mat);
    UNUSED(u4_qp_div_6);
    UNUSED(pi2_tmp);
    UNUSED(pi2_dc_src);

    i4_iq_out_temp = pi2_dc_src[0];
    temp_0 = vdupq_n_s16((i4_iq_out_temp + 32) >> 6);
    dup_min = vdupq_n_s16(RSD_MIN);
    dup_max = vdupq_n_s16(RSD_MAX);

    pred0 = vld1q_s16((int16_t *) pi2_pred);
    pred1 = vld1q_s16((int16_t *) pi2_pred + pred_strd);
    pred2 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 2));
    pred3 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 3));

    pred0 = vaddq_s16(pred0, temp_0);
    pred1 = vaddq_s16(pred1, temp_0);
    pred2 = vaddq_s16(pred2, temp_0);
    pred3 = vaddq_s16(pred3, temp_0);

    pred0 = vminq_s16(pred0, dup_max);
    pred0 = vmaxq_s16(pred0, dup_min);
    pred1 = vminq_s16(pred1, dup_max);
    pred1 = vmaxq_s16(pred1, dup_min);
    pred2 = vminq_s16(pred2, dup_max);
    pred2 = vmaxq_s16(pred2, dup_min);
    pred3 = vminq_s16(pred3, dup_max);
    pred3 = vmaxq_s16(pred3, dup_min);

    pred0 = vandq_s16(pred0, chroma_mask_16x8_2);
    pred1 = vandq_s16(pred1, chroma_mask_16x8_2);
    pred2 = vandq_s16(pred2, chroma_mask_16x8_2);
    pred3 = vandq_s16(pred3, chroma_mask_16x8_2);

    dup_val_1 = vabsq_s16(pred0);
    dup_val_2 = vabsq_s16(pred1);
    dup_val_3 = vabsq_s16(pred2);
    dup_val_4 = vabsq_s16(pred2);
    dup_val_5 = vqaddq_s16(dup_val_1, dup_val_2);
    dup_val_1 = vqaddq_s16(dup_val_3, dup_val_4);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_5);

    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    i4_out_horz_16x8_r0 = vld1q_s16(pi2_out);
    i4_out_horz_16x8_r1 = vld1q_s16(pi2_out + out_strd);
    i4_out_horz_16x8_r2 = vld1q_s16(pi2_out + (out_strd * 2));
    i4_out_horz_16x8_r3 = vld1q_s16(pi2_out + (out_strd * 3));

    i4_out_horz_16x8_r0 = vbslq_s16(chroma_mask_16x8, pred0, i4_out_horz_16x8_r0);
    i4_out_horz_16x8_r1 = vbslq_s16(chroma_mask_16x8, pred1, i4_out_horz_16x8_r1);
    i4_out_horz_16x8_r2 = vbslq_s16(chroma_mask_16x8, pred2, i4_out_horz_16x8_r2);
    i4_out_horz_16x8_r3 = vbslq_s16(chroma_mask_16x8, pred3, i4_out_horz_16x8_r3);

    vst1q_s16((int16_t *) (pi2_out), i4_out_horz_16x8_r0);
    vst1q_s16((int16_t *) (pi2_out + out_strd), i4_out_horz_16x8_r1);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 2)), i4_out_horz_16x8_r2);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 3)), i4_out_horz_16x8_r3);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_8x8_neonintr                 */
/*                                                                           */
/*  Description   : this function computes the resd output from the          */
/*                  IQ+IT                                                    */
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

WORD32 isvcd_iquant_itrans_residual_8x8_neonintr(WORD16 *pi2_src, WORD16 *pi2_pred, WORD16 *pi2_out,
                                                 WORD32 pred_strd, WORD32 out_strd,
                                                 const UWORD16 *pu2_iscale_mat,
                                                 const UWORD16 *pu2_weigh_mat, UWORD32 qp_div,
                                                 WORD16 *pi2_tmp, WORD32 iq_start_idx,
                                                 WORD16 *pi2_dc_ld_addr)
{
    int16x8_t iscal_16x8_0, iscal_16x8_1, iscal_16x8_2, iscal_16x8_3, iscal_16x8_4, iscal_16x8_5,
        iscal_16x8_6, iscal_16x8_7;

    int16x8_t weigh_16x8_0, weigh_16x8_1, weigh_16x8_2, weigh_16x8_3, weigh_16x8_4, weigh_16x8_5,
        weigh_16x8_6, weigh_16x8_7;

    int16x8_t src_16x8_0, src_16x8_1, src_16x8_2, src_16x8_3, src_16x8_4, src_16x8_5, src_16x8_6,
        src_16x8_7;
    int16x8_t coeff_mul_16x8_0, coeff_mul_16x8_1, coeff_mul_16x8_2, coeff_mul_16x8_3,
        coeff_mul_16x8_4, coeff_mul_16x8_5, coeff_mul_16x8_6, coeff_mul_16x8_7;

    int32x4_t quant_res_32x4_l_0, quant_res_32x4_l_1, quant_res_32x4_l_2, quant_res_32x4_l_3,
        quant_res_32x4_l_4, quant_res_32x4_l_5, quant_res_32x4_l_6, quant_res_32x4_l_7;
    int32x4_t quant_res_32x4_h_0, quant_res_32x4_h_1, quant_res_32x4_h_2, quant_res_32x4_h_3,
        quant_res_32x4_h_4, quant_res_32x4_h_5, quant_res_32x4_h_6, quant_res_32x4_h_7;
    int16x4_t quant_res_16x4_l_0, quant_res_16x4_l_1, quant_res_16x4_l_2, quant_res_16x4_l_3,
        quant_res_16x4_l_4, quant_res_16x4_l_5, quant_res_16x4_l_6, quant_res_16x4_l_7;
    int16x4_t quant_res_16x4_h_0, quant_res_16x4_h_1, quant_res_16x4_h_2, quant_res_16x4_h_3,
        quant_res_16x4_h_4, quant_res_16x4_h_5, quant_res_16x4_h_6, quant_res_16x4_h_7;

    int16x8_t quant_res_16x8_0, quant_res_16x8_1, quant_res_16x8_2, quant_res_16x8_3,
        quant_res_16x8_4, quant_res_16x8_5, quant_res_16x8_6, quant_res_16x8_7;

    int16x8_t trans_16x8_0, trans_16x8_1, trans_16x8_2, trans_16x8_3, trans_16x8_4, trans_16x8_5,
        trans_16x8_6, trans_16x8_7;
    int32x4_t trans_32x4_0, trans_32x4_1, trans_32x4_2, trans_32x4_3, trans_32x4_4, trans_32x4_5,
        trans_32x4_6, trans_32x4_7;
    int64x2_t trans_64x2_0, trans_64x2_1, trans_64x2_2, trans_64x2_3, trans_64x2_4, trans_64x2_5,
        trans_64x2_6, trans_64x2_7;
    int16x4_t trans_16x4_1_l, trans_16x4_3_l, trans_16x4_5_l, trans_16x4_7_l;
    int16x8_t rs_trans_16x8_1, rs_trans_16x8_2, rs_trans_16x8_3, rs_trans_16x8_5, rs_trans_16x8_6,
        rs_trans_16x8_7;
    int32x4_t sub_3_5_l, sub_3_5_h;
    int32x4_t add_3_5_l, add_3_5_h;
    int32x4_t sub_1_7_l, sub_1_7_h;
    int32x4_t add_1_7_l, add_1_7_h;
    int32x4_t sub_357_l, sub_357_h;
    int32x4_t add_351_l, add_351_h;
    int32x4_t add_175_l, add_175_h;
    int32x4_t sub_173_l, sub_173_h;
    int32x4_t y1_32x4_l, y1_32x4_h;
    int32x4_t y3_32x4_l, y3_32x4_h;
    int32x4_t y5_32x4_l, y5_32x4_h;
    int32x4_t y7_32x4_l, y7_32x4_h;
    int16x4_t y1_16x4_l, y3_16x4_l, y5_16x4_l, y7_16x4_l;

    int16x8_t y0_16x8, y1_16x8, y2_16x8, y3_16x8, y4_16x8, y5_16x8, y6_16x8, y7_16x8;
    int16x8_t rs_y1_16x8, rs_y3_16x8, rs_y5_16x8, rs_y7_16x8;
    int16x8_t z0_16x8, z1_16x8, z2_16x8, z3_16x8, z4_16x8, z5_16x8, z6_16x8, z7_16x8;
    int16x8_t pred0, pred1, pred2, pred3, pred4, pred5, pred6, pred7;

    int64x2_t pred0_in_64x2, pred1_in_64x2, pred2_in_64x2, pred3_in_64x2, pred4_in_64x2,
        pred5_in_64x2, pred6_in_64x2, pred7_in_64x2;

    int32x4_t qp_div_6_32x4 = vdupq_n_s32(qp_div);

    int16x8_t pred_b0_r01_in;
    int16x8_t pred_b0_r23_in;
    int16x8_t pred_b1_r01_in;
    int16x8_t pred_b1_r23_in;
    int16x8_t pred_b2_r45_in;
    int16x8_t pred_b2_r67_in;
    int16x8_t pred_b3_r45_in;
    int16x8_t pred_b3_r67_in;
    int16x8_t dup_min, dup_max;

    int16x8x2_t trans_16x8_0_1, trans_16x8_2_3, trans_16x8_4_5, trans_16x8_6_7;
    int32x4x2_t trans_32x4_0_2, trans_32x4_1_3, trans_32x4_4_6, trans_32x4_5_7;
    int16x4_t y1_16x4_h, y3_16x4_h, y5_16x4_h, y7_16x4_h;
    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3, i;
    int16x8_t dup_val_1, dup_val_2, dup_abs;

    UNUSED(pi2_tmp);
    UNUSED(iq_start_idx);
    UNUSED(pi2_dc_ld_addr);

    dup_min = vdupq_n_s16(RSD_MIN);
    dup_max = vdupq_n_s16(RSD_MAX);

    iscal_16x8_0 = vld1q_s16((const int16_t *) pu2_iscale_mat);
    iscal_16x8_1 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 8));
    iscal_16x8_2 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 16));
    iscal_16x8_3 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 24));
    iscal_16x8_4 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 32));
    iscal_16x8_5 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 40));
    iscal_16x8_6 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 48));
    iscal_16x8_7 = vld1q_s16((const int16_t *) (pu2_iscale_mat + 56));

    weigh_16x8_0 = vld1q_s16((const int16_t *) pu2_weigh_mat);
    weigh_16x8_1 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 8));
    weigh_16x8_2 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 16));
    weigh_16x8_3 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 24));
    weigh_16x8_4 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 32));
    weigh_16x8_5 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 40));
    weigh_16x8_6 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 48));
    weigh_16x8_7 = vld1q_s16((const int16_t *) (pu2_weigh_mat + 56));

    src_16x8_0 = vld1q_s16((const int16_t *) pi2_src);        // a0 a1 a2 a3 a4 a5 a6 a7
    src_16x8_1 = vld1q_s16((const int16_t *) (pi2_src + 8));  // b0 b1 b2 b3 b4 b5 b6 b7
    src_16x8_2 = vld1q_s16((const int16_t *) (pi2_src + 16));
    src_16x8_3 = vld1q_s16((const int16_t *) (pi2_src + 24));
    src_16x8_4 = vld1q_s16((const int16_t *) (pi2_src + 32));
    src_16x8_5 = vld1q_s16((const int16_t *) (pi2_src + 40));
    src_16x8_6 = vld1q_s16((const int16_t *) (pi2_src + 48));
    src_16x8_7 = vld1q_s16((const int16_t *) (pi2_src + 56));

    coeff_mul_16x8_0 = vmulq_s16(iscal_16x8_0, weigh_16x8_0);
    coeff_mul_16x8_1 = vmulq_s16(iscal_16x8_1, weigh_16x8_1);
    coeff_mul_16x8_2 = vmulq_s16(iscal_16x8_2, weigh_16x8_2);
    coeff_mul_16x8_3 = vmulq_s16(iscal_16x8_3, weigh_16x8_3);
    coeff_mul_16x8_4 = vmulq_s16(iscal_16x8_4, weigh_16x8_4);
    coeff_mul_16x8_5 = vmulq_s16(iscal_16x8_5, weigh_16x8_5);
    coeff_mul_16x8_6 = vmulq_s16(iscal_16x8_6, weigh_16x8_6);
    coeff_mul_16x8_7 = vmulq_s16(iscal_16x8_7, weigh_16x8_7);

    quant_res_32x4_l_0 = vmull_s16(vget_low_s16(coeff_mul_16x8_0), vget_low_s16(src_16x8_0));
    quant_res_32x4_l_1 = vmull_s16(vget_low_s16(coeff_mul_16x8_1), vget_low_s16(src_16x8_1));
    quant_res_32x4_l_2 = vmull_s16(vget_low_s16(coeff_mul_16x8_2), vget_low_s16(src_16x8_2));
    quant_res_32x4_l_3 = vmull_s16(vget_low_s16(coeff_mul_16x8_3), vget_low_s16(src_16x8_3));
    quant_res_32x4_l_4 = vmull_s16(vget_low_s16(coeff_mul_16x8_4), vget_low_s16(src_16x8_4));
    quant_res_32x4_l_5 = vmull_s16(vget_low_s16(coeff_mul_16x8_5), vget_low_s16(src_16x8_5));
    quant_res_32x4_l_6 = vmull_s16(vget_low_s16(coeff_mul_16x8_6), vget_low_s16(src_16x8_6));
    quant_res_32x4_l_7 = vmull_s16(vget_low_s16(coeff_mul_16x8_7), vget_low_s16(src_16x8_7));

    quant_res_32x4_h_0 = vmull_s16(vget_high_s16(coeff_mul_16x8_0), vget_high_s16(src_16x8_0));
    quant_res_32x4_h_1 = vmull_s16(vget_high_s16(coeff_mul_16x8_1), vget_high_s16(src_16x8_1));
    quant_res_32x4_h_2 = vmull_s16(vget_high_s16(coeff_mul_16x8_2), vget_high_s16(src_16x8_2));
    quant_res_32x4_h_3 = vmull_s16(vget_high_s16(coeff_mul_16x8_3), vget_high_s16(src_16x8_3));
    quant_res_32x4_h_4 = vmull_s16(vget_high_s16(coeff_mul_16x8_4), vget_high_s16(src_16x8_4));
    quant_res_32x4_h_5 = vmull_s16(vget_high_s16(coeff_mul_16x8_5), vget_high_s16(src_16x8_5));
    quant_res_32x4_h_6 = vmull_s16(vget_high_s16(coeff_mul_16x8_6), vget_high_s16(src_16x8_6));
    quant_res_32x4_h_7 = vmull_s16(vget_high_s16(coeff_mul_16x8_7), vget_high_s16(src_16x8_7));

    quant_res_32x4_l_0 = vshlq_s32(quant_res_32x4_l_0, qp_div_6_32x4);
    quant_res_32x4_l_1 = vshlq_s32(quant_res_32x4_l_1, qp_div_6_32x4);
    quant_res_32x4_l_2 = vshlq_s32(quant_res_32x4_l_2, qp_div_6_32x4);
    quant_res_32x4_l_3 = vshlq_s32(quant_res_32x4_l_3, qp_div_6_32x4);
    quant_res_32x4_l_4 = vshlq_s32(quant_res_32x4_l_4, qp_div_6_32x4);
    quant_res_32x4_l_5 = vshlq_s32(quant_res_32x4_l_5, qp_div_6_32x4);
    quant_res_32x4_l_6 = vshlq_s32(quant_res_32x4_l_6, qp_div_6_32x4);
    quant_res_32x4_l_7 = vshlq_s32(quant_res_32x4_l_7, qp_div_6_32x4);

    quant_res_32x4_h_0 = vshlq_s32(quant_res_32x4_h_0, qp_div_6_32x4);
    quant_res_32x4_h_1 = vshlq_s32(quant_res_32x4_h_1, qp_div_6_32x4);
    quant_res_32x4_h_2 = vshlq_s32(quant_res_32x4_h_2, qp_div_6_32x4);
    quant_res_32x4_h_3 = vshlq_s32(quant_res_32x4_h_3, qp_div_6_32x4);
    quant_res_32x4_h_4 = vshlq_s32(quant_res_32x4_h_4, qp_div_6_32x4);
    quant_res_32x4_h_5 = vshlq_s32(quant_res_32x4_h_5, qp_div_6_32x4);
    quant_res_32x4_h_6 = vshlq_s32(quant_res_32x4_h_6, qp_div_6_32x4);
    quant_res_32x4_h_7 = vshlq_s32(quant_res_32x4_h_7, qp_div_6_32x4);

    quant_res_16x4_l_0 = vqrshrn_n_s32(quant_res_32x4_l_0, 6);
    quant_res_16x4_l_1 = vqrshrn_n_s32(quant_res_32x4_l_1, 6);
    quant_res_16x4_l_2 = vqrshrn_n_s32(quant_res_32x4_l_2, 6);
    quant_res_16x4_l_3 = vqrshrn_n_s32(quant_res_32x4_l_3, 6);
    quant_res_16x4_l_4 = vqrshrn_n_s32(quant_res_32x4_l_4, 6);
    quant_res_16x4_l_5 = vqrshrn_n_s32(quant_res_32x4_l_5, 6);
    quant_res_16x4_l_6 = vqrshrn_n_s32(quant_res_32x4_l_6, 6);
    quant_res_16x4_l_7 = vqrshrn_n_s32(quant_res_32x4_l_7, 6);

    quant_res_16x4_h_0 = vqrshrn_n_s32(quant_res_32x4_h_0, 6);
    quant_res_16x4_h_1 = vqrshrn_n_s32(quant_res_32x4_h_1, 6);
    quant_res_16x4_h_2 = vqrshrn_n_s32(quant_res_32x4_h_2, 6);
    quant_res_16x4_h_3 = vqrshrn_n_s32(quant_res_32x4_h_3, 6);
    quant_res_16x4_h_4 = vqrshrn_n_s32(quant_res_32x4_h_4, 6);
    quant_res_16x4_h_5 = vqrshrn_n_s32(quant_res_32x4_h_5, 6);
    quant_res_16x4_h_6 = vqrshrn_n_s32(quant_res_32x4_h_6, 6);
    quant_res_16x4_h_7 = vqrshrn_n_s32(quant_res_32x4_h_7, 6);

    quant_res_16x8_0 = vcombine_s16(quant_res_16x4_l_0, quant_res_16x4_h_0);
    quant_res_16x8_1 = vcombine_s16(quant_res_16x4_l_1, quant_res_16x4_h_1);
    quant_res_16x8_2 = vcombine_s16(quant_res_16x4_l_2, quant_res_16x4_h_2);
    quant_res_16x8_3 = vcombine_s16(quant_res_16x4_l_3, quant_res_16x4_h_3);
    quant_res_16x8_4 = vcombine_s16(quant_res_16x4_l_4, quant_res_16x4_h_4);
    quant_res_16x8_5 = vcombine_s16(quant_res_16x4_l_5, quant_res_16x4_h_5);
    quant_res_16x8_6 = vcombine_s16(quant_res_16x4_l_6, quant_res_16x4_h_6);
    quant_res_16x8_7 = vcombine_s16(quant_res_16x4_l_7, quant_res_16x4_h_7);

    for(i = 0; i < 2; i++)
    {
        trans_16x8_0_1 = vtrnq_s16(quant_res_16x8_0, quant_res_16x8_1);
        trans_16x8_0 = trans_16x8_0_1.val[0];
        trans_16x8_1 = trans_16x8_0_1.val[1];

        trans_16x8_2_3 = vtrnq_s16(quant_res_16x8_2, quant_res_16x8_3);
        trans_16x8_2 = trans_16x8_2_3.val[0];
        trans_16x8_3 = trans_16x8_2_3.val[1];

        trans_16x8_4_5 = vtrnq_s16(quant_res_16x8_4, quant_res_16x8_5);
        trans_16x8_4 = trans_16x8_4_5.val[0];
        trans_16x8_5 = trans_16x8_4_5.val[1];

        trans_16x8_6_7 = vtrnq_s16(quant_res_16x8_6, quant_res_16x8_7);
        trans_16x8_6 = trans_16x8_6_7.val[0];
        trans_16x8_7 = trans_16x8_6_7.val[1];

        trans_32x4_0_2 =
            vtrnq_s32(vreinterpretq_s32_s16(trans_16x8_0), vreinterpretq_s32_s16(trans_16x8_2));
        trans_32x4_0 = trans_32x4_0_2.val[0];
        trans_32x4_2 = trans_32x4_0_2.val[1];

        trans_32x4_1_3 =
            vtrnq_s32(vreinterpretq_s32_s16(trans_16x8_1), vreinterpretq_s32_s16(trans_16x8_3));
        trans_32x4_1 = trans_32x4_1_3.val[0];
        trans_32x4_3 = trans_32x4_1_3.val[1];

        trans_32x4_4_6 =
            vtrnq_s32(vreinterpretq_s32_s16(trans_16x8_4), vreinterpretq_s32_s16(trans_16x8_6));
        trans_32x4_4 = trans_32x4_4_6.val[0];
        trans_32x4_6 = trans_32x4_4_6.val[1];

        trans_32x4_5_7 =
            vtrnq_s32(vreinterpretq_s32_s16(trans_16x8_5), vreinterpretq_s32_s16(trans_16x8_7));
        trans_32x4_5 = trans_32x4_5_7.val[0];
        trans_32x4_7 = trans_32x4_5_7.val[1];

        trans_64x2_0 = vcombine_s64(vreinterpret_s64_s32(vget_low_s32(trans_32x4_0)),
                                    vreinterpret_s64_s32(vget_low_s32(trans_32x4_4)));
        trans_64x2_4 = vcombine_s64(vreinterpret_s64_s32(vget_high_s32(trans_32x4_0)),
                                    vreinterpret_s64_s32(vget_high_s32(trans_32x4_4)));

        trans_64x2_1 = vcombine_s64(vreinterpret_s64_s32(vget_low_s32(trans_32x4_1)),
                                    vreinterpret_s64_s32(vget_low_s32(trans_32x4_5)));
        trans_64x2_5 = vcombine_s64(vreinterpret_s64_s32(vget_high_s32(trans_32x4_1)),
                                    vreinterpret_s64_s32(vget_high_s32(trans_32x4_5)));

        trans_64x2_2 = vcombine_s64(vreinterpret_s64_s32(vget_low_s32(trans_32x4_2)),
                                    vreinterpret_s64_s32(vget_low_s32(trans_32x4_6)));
        trans_64x2_6 = vcombine_s64(vreinterpret_s64_s32(vget_high_s32(trans_32x4_2)),
                                    vreinterpret_s64_s32(vget_high_s32(trans_32x4_6)));

        trans_64x2_3 = vcombine_s64(vreinterpret_s64_s32(vget_low_s32(trans_32x4_3)),
                                    vreinterpret_s64_s32(vget_low_s32(trans_32x4_7)));
        trans_64x2_7 = vcombine_s64(vreinterpret_s64_s32(vget_high_s32(trans_32x4_3)),
                                    vreinterpret_s64_s32(vget_high_s32(trans_32x4_7)));

        trans_16x8_0 = vreinterpretq_s16_s64(trans_64x2_0);
        trans_16x8_1 = vreinterpretq_s16_s64(trans_64x2_1);
        trans_16x8_2 = vreinterpretq_s16_s64(trans_64x2_2);
        trans_16x8_3 = vreinterpretq_s16_s64(trans_64x2_3);
        trans_16x8_4 = vreinterpretq_s16_s64(trans_64x2_4);
        trans_16x8_5 = vreinterpretq_s16_s64(trans_64x2_5);
        trans_16x8_6 = vreinterpretq_s16_s64(trans_64x2_6);
        trans_16x8_7 = vreinterpretq_s16_s64(trans_64x2_7);

        rs_trans_16x8_1 = vshrq_n_s16(trans_16x8_1, 1);  //(pi2_tmp_ptr[1] >> 1)
        rs_trans_16x8_2 = vshrq_n_s16(trans_16x8_2, 1);  //(pi2_tmp_ptr[2] >> 1)
        rs_trans_16x8_3 = vshrq_n_s16(trans_16x8_3, 1);  //(pi2_tmp_ptr[3] >> 1)
        rs_trans_16x8_5 = vshrq_n_s16(trans_16x8_5, 1);  //(pi2_tmp_ptr[5] >> 1)
        rs_trans_16x8_6 = vshrq_n_s16(trans_16x8_6, 1);  //(pi2_tmp_ptr[6] >> 1)
        rs_trans_16x8_7 = vshrq_n_s16(trans_16x8_7, 1);  //(pi2_tmp_ptr[7] >> 1)

        y0_16x8 = vaddq_s16(trans_16x8_0,
                            trans_16x8_4);     // i_y0 = (pi2_tmp_ptr[0] + pi2_tmp_ptr[4] );
        y2_16x8 = vsubq_s16(trans_16x8_0,
                            trans_16x8_4);     // i_y2 = (pi2_tmp_ptr[0] - pi2_tmp_ptr[4] );
        y4_16x8 = vsubq_s16(rs_trans_16x8_2,
                            trans_16x8_6);     // i_y4 = ((pi2_tmp_ptr[2] >> 1) - pi2_tmp_ptr[6] );
        y6_16x8 = vaddq_s16(trans_16x8_2,
                            rs_trans_16x8_6);  // i_y6 = (pi2_tmp_ptr[2] + (pi2_tmp_ptr[6] >> 1));

        trans_16x4_3_l = vget_low_s16(trans_16x8_3);
        trans_16x4_5_l = vget_low_s16(trans_16x8_5);

        //-w3 + w5
        sub_3_5_l = vsubl_s16(vget_low_s16(trans_16x8_5), vget_low_s16(trans_16x8_3));
        sub_3_5_h = vsubl_s16(vget_high_s16(trans_16x8_5), vget_high_s16(trans_16x8_3));

        // w3 + w5
        add_3_5_l = vaddl_s16(trans_16x4_3_l, trans_16x4_5_l);
        add_3_5_h = vaddl_s16(vget_high_s16(trans_16x8_3), vget_high_s16(trans_16x8_5));

        trans_16x4_1_l = vget_low_s16(trans_16x8_1);
        trans_16x4_7_l = vget_low_s16(trans_16x8_7);

        //-w1 + w7
        sub_1_7_l = vsubl_s16(trans_16x4_7_l, trans_16x4_1_l);
        sub_1_7_h = vsubl_s16(vget_high_s16(trans_16x8_7), vget_high_s16(trans_16x8_1));

        // w1 + w7
        add_1_7_l = vaddl_s16(trans_16x4_1_l, trans_16x4_7_l);
        add_1_7_h = vaddl_s16(vget_high_s16(trans_16x8_1), vget_high_s16(trans_16x8_7));

        //-w3 + w5 - w7
        sub_357_l = vsubw_s16(sub_3_5_l, trans_16x4_7_l);
        sub_357_h = vsubw_s16(sub_3_5_h, vget_high_s16(trans_16x8_7));

        // w3 + w5 + w1
        add_351_l = vaddw_s16(add_3_5_l, trans_16x4_1_l);
        add_351_h = vaddw_s16(add_3_5_h, vget_high_s16(trans_16x8_1));

        //-w1 + w7 + w5
        add_175_l = vaddw_s16(sub_1_7_l, trans_16x4_5_l);
        add_175_h = vaddw_s16(sub_1_7_h, vget_high_s16(trans_16x8_5));

        // w1 + w7 - w3
        sub_173_l = vsubw_s16(add_1_7_l, trans_16x4_3_l);
        sub_173_h = vsubw_s16(add_1_7_h, vget_high_s16(trans_16x8_3));

        //-w3 + w5 - w7 - (w7 >> 1)
        y1_32x4_l = vsubw_s16(sub_357_l, vget_low_s16(rs_trans_16x8_7));
        y1_32x4_h = vsubw_s16(sub_357_h, vget_high_s16(rs_trans_16x8_7));

        // w1 + w7 - w3 - (w3 >> 1)
        y3_32x4_l = vsubw_s16(sub_173_l, vget_low_s16(rs_trans_16x8_3));
        y3_32x4_h = vsubw_s16(sub_173_h, vget_high_s16(rs_trans_16x8_3));

        //-w1 + w7 + w5 + (w5 >> 1)
        y5_32x4_l = vaddw_s16(add_175_l, vget_low_s16(rs_trans_16x8_5));
        y5_32x4_h = vaddw_s16(add_175_h, vget_high_s16(rs_trans_16x8_5));

        // w3 + w5 + w1 + (w1 >> 1)
        y7_32x4_l = vaddw_s16(add_351_l, vget_low_s16(rs_trans_16x8_1));
        y7_32x4_h = vaddw_s16(add_351_h, vget_high_s16(rs_trans_16x8_1));

        y1_16x4_l = vmovn_s32(y1_32x4_l);
        y1_16x4_h = vmovn_s32(y1_32x4_h);
        y1_16x8 = vcombine_s16(y1_16x4_l, y1_16x4_h);
        y3_16x4_l = vmovn_s32(y3_32x4_l);
        y3_16x4_h = vmovn_s32(y3_32x4_h);
        y3_16x8 = vcombine_s16(y3_16x4_l, y3_16x4_h);
        y5_16x4_l = vmovn_s32(y5_32x4_l);
        y5_16x4_h = vmovn_s32(y5_32x4_h);
        y5_16x8 = vcombine_s16(y5_16x4_l, y5_16x4_h);
        y7_16x4_l = vmovn_s32(y7_32x4_l);
        y7_16x4_h = vmovn_s32(y7_32x4_h);
        y7_16x8 = vcombine_s16(y7_16x4_l, y7_16x4_h);

        rs_y1_16x8 = vshrq_n_s16(y1_16x8, 2);
        rs_y3_16x8 = vshrq_n_s16(y3_16x8, 2);
        rs_y5_16x8 = vshrq_n_s16(y5_16x8, 2);
        rs_y7_16x8 = vshrq_n_s16(y7_16x8, 2);

        z0_16x8 = vaddq_s16(y0_16x8, y6_16x8);           // z0 = y0 + y6
        z1_16x8 = vaddq_s16(y1_16x8, rs_y7_16x8);        // z1 = y1 + (y7 >> 2)
        z2_16x8 = vaddq_s16(y2_16x8, y4_16x8);           // z2 = y2 + y4
        z3_16x8 = vaddq_s16(y3_16x8, rs_y5_16x8);        // z3 = y3 + (y5 >> 2)
        z4_16x8 = vsubq_s16(y2_16x8, y4_16x8);           // z4 = y2 - y4
        z5_16x8 = vsubq_s16(rs_y3_16x8, y5_16x8);        // z5 = (y3 >> 2) - y5
        z6_16x8 = vsubq_s16(y0_16x8, y6_16x8);           // z6 = y0 - y6
        z7_16x8 = vsubq_s16(y7_16x8, rs_y1_16x8);        // z7 = y7 - (y1 >> 2)

        quant_res_16x8_0 = vaddq_s16(z0_16x8, z7_16x8);  // x0 = z0 + z7
        quant_res_16x8_1 = vaddq_s16(z2_16x8, z5_16x8);  // x1 = z2 + z5
        quant_res_16x8_2 = vaddq_s16(z4_16x8, z3_16x8);  // x2 = z4 + z3
        quant_res_16x8_3 = vaddq_s16(z6_16x8, z1_16x8);  // x3 = z6 + z1
        quant_res_16x8_4 = vsubq_s16(z6_16x8, z1_16x8);  // x4 = z6 - z1
        quant_res_16x8_5 = vsubq_s16(z4_16x8, z3_16x8);  // x5 = z4 - z3
        quant_res_16x8_6 = vsubq_s16(z2_16x8, z5_16x8);  // x6 = z2 - z5
        quant_res_16x8_7 = vsubq_s16(z0_16x8, z7_16x8);  // x7 = z0 - z7
    }

    quant_res_16x8_0 = vrshrq_n_s16(quant_res_16x8_0, 6);
    quant_res_16x8_1 = vrshrq_n_s16(quant_res_16x8_1, 6);
    quant_res_16x8_2 = vrshrq_n_s16(quant_res_16x8_2, 6);
    quant_res_16x8_3 = vrshrq_n_s16(quant_res_16x8_3, 6);
    quant_res_16x8_4 = vrshrq_n_s16(quant_res_16x8_4, 6);
    quant_res_16x8_5 = vrshrq_n_s16(quant_res_16x8_5, 6);
    quant_res_16x8_6 = vrshrq_n_s16(quant_res_16x8_6, 6);
    quant_res_16x8_7 = vrshrq_n_s16(quant_res_16x8_7, 6);

    pred0 = vld1q_s16((int16_t *) pi2_pred);
    pred1 = vld1q_s16((int16_t *) pi2_pred + pred_strd);
    pred2 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 2));
    pred3 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 3));
    pred4 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 4));
    pred5 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 5));
    pred6 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 6));
    pred7 = vld1q_s16((int16_t *) pi2_pred + (pred_strd * 7));

    quant_res_16x8_0 = vaddq_s16(pred0, quant_res_16x8_0);
    quant_res_16x8_1 = vaddq_s16(pred1, quant_res_16x8_1);
    quant_res_16x8_2 = vaddq_s16(pred2, quant_res_16x8_2);
    quant_res_16x8_3 = vaddq_s16(pred3, quant_res_16x8_3);
    quant_res_16x8_4 = vaddq_s16(pred4, quant_res_16x8_4);
    quant_res_16x8_5 = vaddq_s16(pred5, quant_res_16x8_5);
    quant_res_16x8_6 = vaddq_s16(pred6, quant_res_16x8_6);
    quant_res_16x8_7 = vaddq_s16(pred7, quant_res_16x8_7);

    quant_res_16x8_0 = vminq_s16(quant_res_16x8_0, dup_max);
    quant_res_16x8_0 = vmaxq_s16(quant_res_16x8_0, dup_min);
    quant_res_16x8_1 = vminq_s16(quant_res_16x8_1, dup_max);
    quant_res_16x8_1 = vmaxq_s16(quant_res_16x8_1, dup_min);
    quant_res_16x8_2 = vminq_s16(quant_res_16x8_2, dup_max);
    quant_res_16x8_2 = vmaxq_s16(quant_res_16x8_2, dup_min);
    quant_res_16x8_3 = vminq_s16(quant_res_16x8_3, dup_max);
    quant_res_16x8_3 = vmaxq_s16(quant_res_16x8_3, dup_min);
    quant_res_16x8_4 = vminq_s16(quant_res_16x8_4, dup_max);
    quant_res_16x8_4 = vmaxq_s16(quant_res_16x8_4, dup_min);
    quant_res_16x8_5 = vminq_s16(quant_res_16x8_5, dup_max);
    quant_res_16x8_5 = vmaxq_s16(quant_res_16x8_5, dup_min);
    quant_res_16x8_6 = vminq_s16(quant_res_16x8_6, dup_max);
    quant_res_16x8_6 = vmaxq_s16(quant_res_16x8_6, dup_min);
    quant_res_16x8_7 = vminq_s16(quant_res_16x8_7, dup_max);
    quant_res_16x8_7 = vmaxq_s16(quant_res_16x8_7, dup_min);

    pred0_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_0);
    pred1_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_1);
    pred2_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_2);
    pred3_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_3);
    pred4_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_4);
    pred5_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_5);
    pred6_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_6);
    pred7_in_64x2 = vreinterpretq_s64_s16(quant_res_16x8_7);

    pred_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred0_in_64x2), vget_low_s64(pred1_in_64x2)));
    pred_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred2_in_64x2), vget_low_s64(pred3_in_64x2)));
    pred_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred0_in_64x2), vget_high_s64(pred1_in_64x2)));
    pred_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred2_in_64x2), vget_high_s64(pred3_in_64x2)));
    pred_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred4_in_64x2), vget_low_s64(pred5_in_64x2)));
    pred_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred6_in_64x2), vget_low_s64(pred7_in_64x2)));
    pred_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred4_in_64x2), vget_high_s64(pred5_in_64x2)));
    pred_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred6_in_64x2), vget_high_s64(pred7_in_64x2)));

    dup_val_1 = vabsq_s16(pred_b0_r01_in);
    dup_val_2 = vabsq_s16(pred_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b1_r01_in);
    dup_val_2 = vabsq_s16(pred_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b2_r45_in);
    dup_val_2 = vabsq_s16(pred_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b3_r45_in);
    dup_val_2 = vabsq_s16(pred_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    vst1q_s16(pi2_out, quant_res_16x8_0);
    vst1q_s16(pi2_out + out_strd, quant_res_16x8_1);
    vst1q_s16(pi2_out + out_strd * 2, quant_res_16x8_2);
    vst1q_s16(pi2_out + out_strd * 3, quant_res_16x8_3);
    vst1q_s16(pi2_out + out_strd * 4, quant_res_16x8_4);
    vst1q_s16(pi2_out + out_strd * 5, quant_res_16x8_5);
    vst1q_s16(pi2_out + out_strd * 6, quant_res_16x8_6);
    vst1q_s16(pi2_out + out_strd * 7, quant_res_16x8_7);

    return nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_iquant_itrans_residual_8x8_dc_neonintr              */
/*                                                                           */
/*  Description   : this function computes the resd dc output from the       */
/*                  IQ+IT                                                    */
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

WORD32 isvcd_iquant_itrans_residual_8x8_dc_neonintr(WORD16 *pi2_src, WORD16 *pi2_pred,
                                                    WORD16 *pi2_out, WORD32 pred_strd,
                                                    WORD32 out_strd, const UWORD16 *pu2_iscale_mat,
                                                    const UWORD16 *pu2_weigh_mat, UWORD32 qp_div,
                                                    WORD16 *pi2_tmp, WORD32 iq_start_idx,
                                                    WORD16 *pi2_dc_ld_addr)
{
    WORD32 i4_iq_out_temp;
    int16x8_t temp_0, dup_min, dup_max;
    int16x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred4_in, pred5_in, pred6_in, pred7_in;
    WORD32 rnd_fact = (qp_div < 6) ? (1 << (5 - qp_div)) : 0;

    int64x2_t pred0_in_64x2, pred1_in_64x2, pred2_in_64x2, pred3_in_64x2, pred4_in_64x2,
        pred5_in_64x2, pred6_in_64x2, pred7_in_64x2;

    int16x8_t pred_b0_r01_in;
    int16x8_t pred_b0_r23_in;
    int16x8_t pred_b1_r01_in;
    int16x8_t pred_b1_r23_in;
    int16x8_t pred_b2_r45_in;
    int16x8_t pred_b2_r67_in;
    int16x8_t pred_b3_r45_in;
    int16x8_t pred_b3_r67_in;

    int16x8_t dup_val_1, dup_val_2, dup_abs;
    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;

    UNUSED(pi2_tmp);
    UNUSED(iq_start_idx);
    UNUSED(pi2_dc_ld_addr);

    i4_iq_out_temp = pi2_src[0];
    INV_QUANT(i4_iq_out_temp, pu2_iscale_mat[0], pu2_weigh_mat[0], qp_div, rnd_fact, 6);

    /* inv transform and residual comp */
    pred0_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred1_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred2_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred3_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred4_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred5_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred6_in = vld1q_s16((int16_t *) pi2_pred);
    pi2_pred = pi2_pred + pred_strd;
    pred7_in = vld1q_s16((int16_t *) pi2_pred);

    temp_0 = vdupq_n_s16((i4_iq_out_temp + 32) >> 6);
    dup_min = vdupq_n_s16(RSD_MIN);
    dup_max = vdupq_n_s16(RSD_MAX);

    pred0_in = vaddq_s16(pred0_in, temp_0);
    pred1_in = vaddq_s16(pred1_in, temp_0);
    pred2_in = vaddq_s16(pred2_in, temp_0);
    pred3_in = vaddq_s16(pred3_in, temp_0);
    pred4_in = vaddq_s16(pred4_in, temp_0);
    pred5_in = vaddq_s16(pred5_in, temp_0);
    pred6_in = vaddq_s16(pred6_in, temp_0);
    pred7_in = vaddq_s16(pred7_in, temp_0);

    pred0_in = vminq_s16(pred0_in, dup_max);
    pred0_in = vmaxq_s16(pred0_in, dup_min);
    pred1_in = vminq_s16(pred1_in, dup_max);
    pred1_in = vmaxq_s16(pred1_in, dup_min);
    pred2_in = vminq_s16(pred2_in, dup_max);
    pred2_in = vmaxq_s16(pred2_in, dup_min);
    pred3_in = vminq_s16(pred3_in, dup_max);
    pred3_in = vmaxq_s16(pred3_in, dup_min);
    pred4_in = vminq_s16(pred4_in, dup_max);
    pred4_in = vmaxq_s16(pred4_in, dup_min);
    pred5_in = vminq_s16(pred5_in, dup_max);
    pred5_in = vmaxq_s16(pred5_in, dup_min);
    pred6_in = vminq_s16(pred6_in, dup_max);
    pred6_in = vmaxq_s16(pred6_in, dup_min);
    pred7_in = vminq_s16(pred7_in, dup_max);
    pred7_in = vmaxq_s16(pred7_in, dup_min);

    pred0_in_64x2 = vreinterpretq_s64_s16(pred0_in);
    pred1_in_64x2 = vreinterpretq_s64_s16(pred1_in);
    pred2_in_64x2 = vreinterpretq_s64_s16(pred2_in);
    pred3_in_64x2 = vreinterpretq_s64_s16(pred3_in);
    pred4_in_64x2 = vreinterpretq_s64_s16(pred4_in);
    pred5_in_64x2 = vreinterpretq_s64_s16(pred5_in);
    pred6_in_64x2 = vreinterpretq_s64_s16(pred6_in);
    pred7_in_64x2 = vreinterpretq_s64_s16(pred7_in);

    pred_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred0_in_64x2), vget_low_s64(pred1_in_64x2)));
    pred_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred2_in_64x2), vget_low_s64(pred3_in_64x2)));
    pred_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred0_in_64x2), vget_high_s64(pred1_in_64x2)));
    pred_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred2_in_64x2), vget_high_s64(pred3_in_64x2)));
    pred_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred4_in_64x2), vget_low_s64(pred5_in_64x2)));
    pred_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(pred6_in_64x2), vget_low_s64(pred7_in_64x2)));
    pred_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred4_in_64x2), vget_high_s64(pred5_in_64x2)));
    pred_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(pred6_in_64x2), vget_high_s64(pred7_in_64x2)));

    dup_val_1 = vabsq_s16(pred_b0_r01_in);
    dup_val_2 = vabsq_s16(pred_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b1_r01_in);
    dup_val_2 = vabsq_s16(pred_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b2_r45_in);
    dup_val_2 = vabsq_s16(pred_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(pred_b3_r45_in);
    dup_val_2 = vabsq_s16(pred_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    vst1q_s16((int16_t *) (pi2_out), pred0_in);
    vst1q_s16((int16_t *) (pi2_out + out_strd), pred1_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 2)), pred2_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 3)), pred3_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 4)), pred4_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 5)), pred5_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 6)), pred6_in);
    vst1q_s16((int16_t *) (pi2_out + (out_strd * 7)), pred7_in);

    return nnz;
}
