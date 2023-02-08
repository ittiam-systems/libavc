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
*  isvcd_pred_residual_recon_neonintr.c
*
* @brief
*  Contains definition of functions for h264 inverse quantization inverse
*    transformation and resd comp
*
* @author
*  Kishore
*
*  @par List of Functions:
*  - isvcd_pred_residual_recon_16x16_neonintr()
*  - isvcd_pred_residual_recon_8x8_neonintr()
*  - isvcd_pred_residual_recon_4x4_neonintr()
*  - isvcd_pred_residual_recon_chroma_4x4_neonintr()
*  - isvcd_pred_residual_recon_chroma_8x8_neonintr()
*  - isvcd_residual_luma_4x4_neonintr()
*  - isvcd_residual_luma_8x8_neonintr()
*  - isvcd_residual_luma_16x16_neonintr()
*  - isvcd_residual_chroma_cb_cr_8x8_neonintr()
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
#include "isvcd_pred_residual_recon.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_4x4_neonintr                    */
/*                                                                           */
/*  Description   : this function computes the recon data from the           */
/*                   pred and residual buffer                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_4x4_neonintr(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                              WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd01_in, resd23_in;
    WORD32 i4_nnz;
    int16x8_t dup_val_1, dup_val_2, dup_abs;

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 1));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));

    resd01_in = vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(resd0_in)),
                                                   vget_low_s64(vreinterpretq_s64_s16(resd1_in))));

    resd23_in = vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(resd2_in)),
                                                   vget_low_s64(vreinterpretq_s64_s16(resd3_in))));

    dup_val_1 = vabsq_s16(resd01_in);
    dup_val_2 = vabsq_s16(resd23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    vst1_lane_u32((uint32_t *) (pu1_out), vreinterpret_u32_u8(pred0_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + out_strd), vreinterpret_u32_u8(pred1_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + out_strd * 2), vreinterpret_u32_u8(pred2_in), 0);
    vst1_lane_u32((uint32_t *) (pu1_out + out_strd * 3), vreinterpret_u32_u8(pred3_in), 0);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_8x8_neonintr                    */
/*                                                                           */
/*  Description   : this function computes the recon data from the           */
/*                   pred and residual buffer                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_8x8_neonintr(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                              WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    uint8x8_t pred4_in, pred5_in, pred6_in, pred7_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t pred4, pred5, pred6, pred7;
    int16x8_t resd4_in, resd5_in, resd6_in, resd7_in;
    int16x8_t dup_val_1, dup_val_2, dup_abs;
    int64x2_t resd0_in_64x2, resd1_in_64x2, resd2_in_64x2, resd3_in_64x2, resd4_in_64x2,
        resd5_in_64x2, resd6_in_64x2, resd7_in_64x2;

    int16x8_t resd_b0_r01_in;
    int16x8_t resd_b0_r23_in;
    int16x8_t resd_b1_r01_in;
    int16x8_t resd_b1_r23_in;
    int16x8_t resd_b2_r45_in;
    int16x8_t resd_b2_r67_in;
    int16x8_t resd_b3_r45_in;
    int16x8_t resd_b3_r67_in;

    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 3));
    pred4_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 4));
    pred5_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 5));
    pred6_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 6));
    pred7_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 7));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));

    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    vst1_u8((uint8_t *) (pu1_out), pred0_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd), pred1_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 2), pred2_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 3), pred3_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 4), pred4_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 5), pred5_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 6), pred6_in);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 7), pred7_in);

    return nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_16x16_neonintr                  */
/*                                                                           */
/*  Description   : this function computes the recon data from the           */
/*                   pred and residual buffer                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_16x16_neonintr(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                                WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    uint8x8_t pred4_in, pred5_in, pred6_in, pred7_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t pred4, pred5, pred6, pred7;
    int16x8_t resd4_in, resd5_in, resd6_in, resd7_in;
    int16x8_t dup_val_1, dup_val_2, dup_abs;
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;

    int64x2_t resd0_in_64x2, resd1_in_64x2, resd2_in_64x2, resd3_in_64x2, resd4_in_64x2,
        resd5_in_64x2, resd6_in_64x2, resd7_in_64x2;

    int16x8_t resd_b0_r01_in;
    int16x8_t resd_b0_r23_in;
    int16x8_t resd_b1_r01_in;
    int16x8_t resd_b1_r23_in;
    int16x8_t resd_b2_r45_in;
    int16x8_t resd_b2_r67_in;
    int16x8_t resd_b3_r45_in;
    int16x8_t resd_b3_r67_in;

    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;

    /* First row of 8,  first 8x8 elements */
    pred0_in = vld1_u8((uint8_t *) pu1_pred_ptr);
    pred1_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 3));
    pred4_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 4));
    pred5_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 5));
    pred6_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 6));
    pred7_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 7));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));

    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    vst1_u8((uint8_t *) (pu1_out_ptr), pred0_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd), pred1_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 2), pred2_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 3), pred3_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 4), pred4_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 5), pred5_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 6), pred6_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 7), pred7_in);

    /* first row of 8, sec 8x8 elements */
    pu1_out_ptr = pu1_out_ptr + 8;
    pi2_rsd_ptr = pi2_rsd_ptr + 8;
    pu1_pred_ptr = pu1_pred_ptr + 8;

    pred0_in = vld1_u8((uint8_t *) pu1_pred_ptr);
    pred1_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 3));
    pred4_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 4));
    pred5_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 5));
    pred6_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 6));
    pred7_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 7));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));

    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 2 | (nnz_b1 << 3) | (nnz_b2 << 6) | (nnz_b3 << 7));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    vst1_u8((uint8_t *) (pu1_out_ptr), pred0_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd), pred1_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 2), pred2_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 3), pred3_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 4), pred4_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 5), pred5_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 6), pred6_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 7), pred7_in);

    pu1_out_ptr = pu1_out + (8 * out_strd);
    pi2_rsd_ptr = pi2_rsd + (8 * rsd_strd);
    pu1_pred_ptr = pu1_pred + (8 * pred_strd);

    /*Sec row of 8, first 8x8*/
    pred0_in = vld1_u8((uint8_t *) pu1_pred_ptr);
    pred1_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 3));
    pred4_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 4));
    pred5_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 5));
    pred6_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 6));
    pred7_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 7));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 8 | (nnz_b1 << 9) | (nnz_b2 << 12) | (nnz_b3 << 13));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    vst1_u8((uint8_t *) (pu1_out_ptr), pred0_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd), pred1_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 2), pred2_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 3), pred3_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 4), pred4_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 5), pred5_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 6), pred6_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 7), pred7_in);

    /*Sec row of 8, Sec 8x8*/
    pu1_out_ptr = pu1_out_ptr + 8;
    pi2_rsd_ptr = pi2_rsd_ptr + 8;
    pu1_pred_ptr = pu1_pred_ptr + 8;

    pred0_in = vld1_u8((uint8_t *) pu1_pred_ptr);
    pred1_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 3));
    pred4_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 4));
    pred5_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 5));
    pred6_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 6));
    pred7_in = vld1_u8((uint8_t *) pu1_pred_ptr + (pred_strd * 7));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 10 | (nnz_b1 << 11) | (nnz_b2 << 14) | (nnz_b3 << 15));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    vst1_u8((uint8_t *) (pu1_out_ptr), pred0_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd), pred1_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 2), pred2_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 3), pred3_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 4), pred4_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 5), pred5_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 6), pred6_in);
    vst1_u8((uint8_t *) (pu1_out_ptr + out_strd * 7), pred7_in);

    return nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_8x8_neonintr             */
/*                                                                           */
/*  Description   : this function computes the recon data from the           */
/*                   pred and residual buffer                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
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

void isvcd_pred_residual_recon_chroma_8x8_neonintr(UWORD8 *pu1_pred, WORD16 *pi2_rsd,
                                                   UWORD8 *pu1_out, WORD32 pred_strd,
                                                   WORD32 rsd_strd, WORD32 out_strd)
{
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    uint8x8_t pred4_in, pred5_in, pred6_in, pred7_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t pred4, pred5, pred6, pred7;
    int16x8_t resd4_in, resd5_in, resd6_in, resd7_in;

    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    uint8x16_t pred0_inp_full, pred1_inp_full, pred2_inp_full, pred3_inp_full;
    uint8x16_t pred4_inp_full, pred5_inp_full, pred6_inp_full, pred7_inp_full;
    uint8x8_t i4_out_horz_8x8_r0, i4_out_horz_8x8_r1, i4_out_horz_8x8_r2, i4_out_horz_8x8_r3;
    uint8x8_t i4_out_horz_8x8_r4, i4_out_horz_8x8_r5, i4_out_horz_8x8_r6, i4_out_horz_8x8_r7;
    uint8x8_t chroma_mask_8x8;

    pred0_inp_full = vld1q_u8((uint8_t *) pu1_pred);
    pred1_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd));
    pred2_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 2));
    pred3_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 3));
    pred4_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 4));
    pred5_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 5));
    pred6_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 6));
    pred7_inp_full = vld1q_u8((uint8_t *) pu1_pred + (pred_strd * 7));

    pred0_in = vget_low_u8(pred0_inp_full);
    pred1_in = vget_low_u8(pred1_inp_full);
    pred2_in = vget_low_u8(pred2_inp_full);
    pred3_in = vget_low_u8(pred3_inp_full);
    pred4_in = vget_low_u8(pred4_inp_full);
    pred5_in = vget_low_u8(pred5_inp_full);
    pred6_in = vget_low_u8(pred6_inp_full);
    pred7_in = vget_low_u8(pred7_inp_full);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 7));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));

    i4_out_horz_8x8_r0 = vld1_u8(pu1_out);
    i4_out_horz_8x8_r1 = vld1_u8(pu1_out + out_strd);
    i4_out_horz_8x8_r2 = vld1_u8(pu1_out + out_strd * 2);
    i4_out_horz_8x8_r3 = vld1_u8(pu1_out + out_strd * 3);
    i4_out_horz_8x8_r4 = vld1_u8(pu1_out + out_strd * 4);
    i4_out_horz_8x8_r5 = vld1_u8(pu1_out + out_strd * 5);
    i4_out_horz_8x8_r6 = vld1_u8(pu1_out + out_strd * 6);
    i4_out_horz_8x8_r7 = vld1_u8(pu1_out + out_strd * 7);

    i4_out_horz_8x8_r0 = vbsl_u8(chroma_mask_8x8, pred0_in, i4_out_horz_8x8_r0);
    i4_out_horz_8x8_r1 = vbsl_u8(chroma_mask_8x8, pred1_in, i4_out_horz_8x8_r1);
    i4_out_horz_8x8_r2 = vbsl_u8(chroma_mask_8x8, pred2_in, i4_out_horz_8x8_r2);
    i4_out_horz_8x8_r3 = vbsl_u8(chroma_mask_8x8, pred3_in, i4_out_horz_8x8_r3);
    i4_out_horz_8x8_r4 = vbsl_u8(chroma_mask_8x8, pred4_in, i4_out_horz_8x8_r4);
    i4_out_horz_8x8_r5 = vbsl_u8(chroma_mask_8x8, pred5_in, i4_out_horz_8x8_r5);
    i4_out_horz_8x8_r6 = vbsl_u8(chroma_mask_8x8, pred6_in, i4_out_horz_8x8_r6);
    i4_out_horz_8x8_r7 = vbsl_u8(chroma_mask_8x8, pred7_in, i4_out_horz_8x8_r7);

    vst1_u8((uint8_t *) (pu1_out), i4_out_horz_8x8_r0);
    vst1_u8((uint8_t *) (pu1_out + out_strd), i4_out_horz_8x8_r1);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 2), i4_out_horz_8x8_r2);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 3), i4_out_horz_8x8_r3);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 4), i4_out_horz_8x8_r4);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 5), i4_out_horz_8x8_r5);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 6), i4_out_horz_8x8_r6);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 7), i4_out_horz_8x8_r7);

    /* for the next 4 elements interleaved format */
    pred0_in = vget_high_u8(pred0_inp_full);
    pred1_in = vget_high_u8(pred1_inp_full);
    pred2_in = vget_high_u8(pred2_inp_full);
    pred3_in = vget_high_u8(pred3_inp_full);
    pred4_in = vget_high_u8(pred4_inp_full);
    pred5_in = vget_high_u8(pred5_inp_full);
    pred6_in = vget_high_u8(pred6_inp_full);
    pred7_in = vget_high_u8(pred7_inp_full);

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));
    pred4 = vreinterpretq_s16_u16(vmovl_u8(pred4_in));
    pred5 = vreinterpretq_s16_u16(vmovl_u8(pred5_in));
    pred6 = vreinterpretq_s16_u16(vmovl_u8(pred6_in));
    pred7 = vreinterpretq_s16_u16(vmovl_u8(pred7_in));

    pi2_rsd = pi2_rsd + 8;
    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 7));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);
    pred4 = vaddq_s16(pred4, resd4_in);
    pred5 = vaddq_s16(pred5, resd5_in);
    pred6 = vaddq_s16(pred6, resd6_in);
    pred7 = vaddq_s16(pred7, resd7_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);
    pred4_in = vqmovun_s16(pred4);
    pred5_in = vqmovun_s16(pred5);
    pred6_in = vqmovun_s16(pred6);
    pred7_in = vqmovun_s16(pred7);

    pu1_out = pu1_out + 8;
    i4_out_horz_8x8_r0 = vld1_u8(pu1_out);
    i4_out_horz_8x8_r1 = vld1_u8(pu1_out + out_strd);
    i4_out_horz_8x8_r2 = vld1_u8(pu1_out + out_strd * 2);
    i4_out_horz_8x8_r3 = vld1_u8(pu1_out + out_strd * 3);
    i4_out_horz_8x8_r4 = vld1_u8(pu1_out + out_strd * 4);
    i4_out_horz_8x8_r5 = vld1_u8(pu1_out + out_strd * 5);
    i4_out_horz_8x8_r6 = vld1_u8(pu1_out + out_strd * 6);
    i4_out_horz_8x8_r7 = vld1_u8(pu1_out + out_strd * 7);

    i4_out_horz_8x8_r0 = vbsl_u8(chroma_mask_8x8, pred0_in, i4_out_horz_8x8_r0);
    i4_out_horz_8x8_r1 = vbsl_u8(chroma_mask_8x8, pred1_in, i4_out_horz_8x8_r1);
    i4_out_horz_8x8_r2 = vbsl_u8(chroma_mask_8x8, pred2_in, i4_out_horz_8x8_r2);
    i4_out_horz_8x8_r3 = vbsl_u8(chroma_mask_8x8, pred3_in, i4_out_horz_8x8_r3);
    i4_out_horz_8x8_r4 = vbsl_u8(chroma_mask_8x8, pred4_in, i4_out_horz_8x8_r4);
    i4_out_horz_8x8_r5 = vbsl_u8(chroma_mask_8x8, pred5_in, i4_out_horz_8x8_r5);
    i4_out_horz_8x8_r6 = vbsl_u8(chroma_mask_8x8, pred6_in, i4_out_horz_8x8_r6);
    i4_out_horz_8x8_r7 = vbsl_u8(chroma_mask_8x8, pred7_in, i4_out_horz_8x8_r7);

    vst1_u8((uint8_t *) (pu1_out), i4_out_horz_8x8_r0);
    vst1_u8((uint8_t *) (pu1_out + out_strd), i4_out_horz_8x8_r1);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 2), i4_out_horz_8x8_r2);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 3), i4_out_horz_8x8_r3);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 4), i4_out_horz_8x8_r4);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 5), i4_out_horz_8x8_r5);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 6), i4_out_horz_8x8_r6);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 7), i4_out_horz_8x8_r7);

    pu1_out = pu1_out_ptr;
    pu1_pred = pu1_pred_ptr;
    pi2_rsd = pi2_rsd_ptr;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_4x4_neonintr             */
/*                                                                           */
/*  Description   : this function computes the recon data from the           */
/*                   pred and residual buffer                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
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

void isvcd_pred_residual_recon_chroma_4x4_neonintr(UWORD8 *pu1_pred, WORD16 *pi2_rsd,
                                                   UWORD8 *pu1_out, WORD32 pred_strd,
                                                   WORD32 rsd_strd, WORD32 out_strd)
{
    uint8x8_t pred0_in, pred1_in, pred2_in, pred3_in;
    int16x8_t pred0, pred1, pred2, pred3;
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;

    uint8x8_t i4_out_horz_8x8_r0, i4_out_horz_8x8_r1, i4_out_horz_8x8_r2, i4_out_horz_8x8_r3;
    uint8x8_t chroma_mask_8x8 = vreinterpret_u8_u16(vdup_n_u16(0x00ff));

    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;

    pred0_in = vld1_u8((uint8_t *) pu1_pred);
    pred1_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd));
    pred2_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 2));
    pred3_in = vld1_u8((uint8_t *) pu1_pred + (pred_strd * 3));

    pred0 = vreinterpretq_s16_u16(vmovl_u8(pred0_in));
    pred1 = vreinterpretq_s16_u16(vmovl_u8(pred1_in));
    pred2 = vreinterpretq_s16_u16(vmovl_u8(pred2_in));
    pred3 = vreinterpretq_s16_u16(vmovl_u8(pred3_in));

    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));

    pred0 = vaddq_s16(pred0, resd0_in);
    pred1 = vaddq_s16(pred1, resd1_in);
    pred2 = vaddq_s16(pred2, resd2_in);
    pred3 = vaddq_s16(pred3, resd3_in);

    pred0_in = vqmovun_s16(pred0);
    pred1_in = vqmovun_s16(pred1);
    pred2_in = vqmovun_s16(pred2);
    pred3_in = vqmovun_s16(pred3);

    i4_out_horz_8x8_r0 = vld1_u8(pu1_out);
    i4_out_horz_8x8_r1 = vld1_u8(pu1_out + out_strd);
    i4_out_horz_8x8_r2 = vld1_u8(pu1_out + out_strd * 2);
    i4_out_horz_8x8_r3 = vld1_u8(pu1_out + out_strd * 3);

    i4_out_horz_8x8_r0 = vbsl_u8(chroma_mask_8x8, pred0_in, i4_out_horz_8x8_r0);
    i4_out_horz_8x8_r1 = vbsl_u8(chroma_mask_8x8, pred1_in, i4_out_horz_8x8_r1);
    i4_out_horz_8x8_r2 = vbsl_u8(chroma_mask_8x8, pred2_in, i4_out_horz_8x8_r2);
    i4_out_horz_8x8_r3 = vbsl_u8(chroma_mask_8x8, pred3_in, i4_out_horz_8x8_r3);

    vst1_u8((uint8_t *) (pu1_out), i4_out_horz_8x8_r0);
    vst1_u8((uint8_t *) (pu1_out + out_strd), i4_out_horz_8x8_r1);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 2), i4_out_horz_8x8_r2);
    vst1_u8((uint8_t *) (pu1_out + out_strd * 3), i4_out_horz_8x8_r3);

    pu1_out = pu1_out_ptr;
    pu1_pred = pu1_pred_ptr;
    pi2_rsd = pi2_rsd_ptr;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_4x4_neonintr                          */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_4x4_neonintr(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd01_in, resd23_in;
    WORD32 i4_nnz;
    int16x8_t dup_val_1, dup_val_2, dup_abs;
    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));

    resd01_in = vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(resd0_in)),
                                                   vget_low_s64(vreinterpretq_s64_s16(resd1_in))));

    resd23_in = vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(resd2_in)),
                                                   vget_low_s64(vreinterpretq_s64_s16(resd3_in))));

    dup_val_1 = vabsq_s16(resd01_in);
    dup_val_2 = vabsq_s16(resd23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    i4_nnz = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_8x8_neonintr                          */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_8x8_neonintr(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd4_in, resd5_in, resd6_in, resd7_in;

    int64x2_t resd0_in_64x2, resd1_in_64x2, resd2_in_64x2, resd3_in_64x2, resd4_in_64x2,
        resd5_in_64x2, resd6_in_64x2, resd7_in_64x2;

    int16x8_t resd_b0_r01_in;
    int16x8_t resd_b0_r23_in;
    int16x8_t resd_b1_r01_in;
    int16x8_t resd_b1_r23_in;
    int16x8_t resd_b2_r45_in;
    int16x8_t resd_b2_r67_in;
    int16x8_t resd_b3_r45_in;
    int16x8_t resd_b3_r67_in;

    int16x8_t dup_val_1, dup_val_2, dup_abs;
    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;

    resd0_in = vld1q_s16((int16_t *) pi2_rsd);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    return nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_16x16_neonintr                        */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_16x16_neonintr(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    int16x8_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8_t resd4_in, resd5_in, resd6_in, resd7_in;

    WORD16 *pi2_rsd_ptr = pi2_rsd;
    int64x2_t resd0_in_64x2, resd1_in_64x2, resd2_in_64x2, resd3_in_64x2, resd4_in_64x2,
        resd5_in_64x2, resd6_in_64x2, resd7_in_64x2;

    int16x8_t resd_b0_r01_in;
    int16x8_t resd_b0_r23_in;
    int16x8_t resd_b1_r01_in;
    int16x8_t resd_b1_r23_in;
    int16x8_t resd_b2_r45_in;
    int16x8_t resd_b2_r67_in;
    int16x8_t resd_b3_r45_in;
    int16x8_t resd_b3_r67_in;

    int16x8_t dup_val_1, dup_val_2, dup_abs;
    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;

    /* First row of 8,  first 8x8 elements */
    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = (nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 4) | (nnz_b3 << 5));

    /* first row of 8, sec 8x8 elements */
    pi2_rsd_ptr = pi2_rsd_ptr + 8;

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 2 | (nnz_b1 << 3) | (nnz_b2 << 6) | (nnz_b3 << 7));

    pi2_rsd_ptr = pi2_rsd + (8 * rsd_strd);
    /*Sec row of 8, first 8x8*/
    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 8 | (nnz_b1 << 9) | (nnz_b2 << 12) | (nnz_b3 << 13));

    /*Sec row of 8, Sec 8x8*/
    pi2_rsd_ptr = pi2_rsd_ptr + 8;

    resd0_in = vld1q_s16((int16_t *) pi2_rsd_ptr);
    resd1_in = vld1q_s16((int16_t *) pi2_rsd_ptr + rsd_strd);
    resd2_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 2));
    resd3_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 3));
    resd4_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 4));
    resd5_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 5));
    resd6_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 6));
    resd7_in = vld1q_s16((int16_t *) pi2_rsd_ptr + (rsd_strd * 7));

    resd0_in_64x2 = vreinterpretq_s64_s16(resd0_in);
    resd1_in_64x2 = vreinterpretq_s64_s16(resd1_in);
    resd2_in_64x2 = vreinterpretq_s64_s16(resd2_in);
    resd3_in_64x2 = vreinterpretq_s64_s16(resd3_in);
    resd4_in_64x2 = vreinterpretq_s64_s16(resd4_in);
    resd5_in_64x2 = vreinterpretq_s64_s16(resd5_in);
    resd6_in_64x2 = vreinterpretq_s64_s16(resd6_in);
    resd7_in_64x2 = vreinterpretq_s64_s16(resd7_in);

    resd_b0_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_in_64x2), vget_low_s64(resd1_in_64x2)));
    resd_b0_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_in_64x2), vget_low_s64(resd3_in_64x2)));
    resd_b1_r01_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_in_64x2), vget_high_s64(resd1_in_64x2)));
    resd_b1_r23_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_in_64x2), vget_high_s64(resd3_in_64x2)));
    resd_b2_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_in_64x2), vget_low_s64(resd5_in_64x2)));
    resd_b2_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_in_64x2), vget_low_s64(resd7_in_64x2)));
    resd_b3_r45_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_in_64x2), vget_high_s64(resd5_in_64x2)));
    resd_b3_r67_in = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_in_64x2), vget_high_s64(resd7_in_64x2)));
    dup_val_1 = vabsq_s16(resd_b0_r01_in);
    dup_val_2 = vabsq_s16(resd_b0_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_in);
    dup_val_2 = vabsq_s16(resd_b1_r23_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_in);
    dup_val_2 = vabsq_s16(resd_b2_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_in);
    dup_val_2 = vabsq_s16(resd_b3_r67_in);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= (nnz_b0 << 10 | (nnz_b1 << 11) | (nnz_b2 << 14) | (nnz_b3 << 15));
    return nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_chroma_cb_cr_8x8_neonintr                  */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_chroma_cb_cr_8x8_neonintr(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    int16x8x2_t resd0_in, resd1_in, resd2_in, resd3_in;
    int16x8x2_t resd4_in, resd5_in, resd6_in, resd7_in;

    int64x2_t resd0_cr_64x2, resd1_cr_64x2, resd2_cr_64x2, resd3_cr_64x2, resd4_cr_64x2,
        resd5_cr_64x2, resd6_cr_64x2, resd7_cr_64x2;

    int16x8_t resd_b0_r01_cr;
    int16x8_t resd_b0_r23_cr;
    int16x8_t resd_b1_r01_cr;
    int16x8_t resd_b1_r23_cr;
    int16x8_t resd_b2_r45_cr;
    int16x8_t resd_b2_r67_cr;
    int16x8_t resd_b3_r45_cr;
    int16x8_t resd_b3_r67_cr;

    int64x2_t resd0_cb_64x2, resd1_cb_64x2, resd2_cb_64x2, resd3_cb_64x2, resd4_cb_64x2,
        resd5_cb_64x2, resd6_cb_64x2, resd7_cb_64x2;

    int16x8_t resd_b0_r01_cb;
    int16x8_t resd_b0_r23_cb;
    int16x8_t resd_b1_r01_cb;
    int16x8_t resd_b1_r23_cb;
    int16x8_t resd_b2_r45_cb;
    int16x8_t resd_b2_r67_cb;
    int16x8_t resd_b3_r45_cb;
    int16x8_t resd_b3_r67_cb;

    WORD32 nnz, nnz_b0, nnz_b1, nnz_b2, nnz_b3;
    int16x8_t dup_val_1, dup_val_2, dup_abs;

    resd0_in = vld2q_s16((int16_t *) pi2_rsd);
    resd1_in = vld2q_s16((int16_t *) pi2_rsd + rsd_strd);
    resd2_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 2));
    resd3_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 3));
    resd4_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 4));
    resd5_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 5));
    resd6_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 6));
    resd7_in = vld2q_s16((int16_t *) pi2_rsd + (rsd_strd * 7));

    resd0_cb_64x2 = vreinterpretq_s64_s16(resd0_in.val[0]);
    resd1_cb_64x2 = vreinterpretq_s64_s16(resd1_in.val[0]);
    resd2_cb_64x2 = vreinterpretq_s64_s16(resd2_in.val[0]);
    resd3_cb_64x2 = vreinterpretq_s64_s16(resd3_in.val[0]);
    resd4_cb_64x2 = vreinterpretq_s64_s16(resd4_in.val[0]);
    resd5_cb_64x2 = vreinterpretq_s64_s16(resd5_in.val[0]);
    resd6_cb_64x2 = vreinterpretq_s64_s16(resd6_in.val[0]);
    resd7_cb_64x2 = vreinterpretq_s64_s16(resd7_in.val[0]);

    resd_b0_r01_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_cb_64x2), vget_low_s64(resd1_cb_64x2)));
    resd_b0_r23_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_cb_64x2), vget_low_s64(resd3_cb_64x2)));
    resd_b1_r01_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_cb_64x2), vget_high_s64(resd1_cb_64x2)));
    resd_b1_r23_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_cb_64x2), vget_high_s64(resd3_cb_64x2)));
    resd_b2_r45_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_cb_64x2), vget_low_s64(resd5_cb_64x2)));
    resd_b2_r67_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_cb_64x2), vget_low_s64(resd7_cb_64x2)));
    resd_b3_r45_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_cb_64x2), vget_high_s64(resd5_cb_64x2)));
    resd_b3_r67_cb = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_cb_64x2), vget_high_s64(resd7_cb_64x2)));

    resd0_cr_64x2 = vreinterpretq_s64_s16(resd0_in.val[1]);
    resd1_cr_64x2 = vreinterpretq_s64_s16(resd1_in.val[1]);
    resd2_cr_64x2 = vreinterpretq_s64_s16(resd2_in.val[1]);
    resd3_cr_64x2 = vreinterpretq_s64_s16(resd3_in.val[1]);
    resd4_cr_64x2 = vreinterpretq_s64_s16(resd4_in.val[1]);
    resd5_cr_64x2 = vreinterpretq_s64_s16(resd5_in.val[1]);
    resd6_cr_64x2 = vreinterpretq_s64_s16(resd6_in.val[1]);
    resd7_cr_64x2 = vreinterpretq_s64_s16(resd7_in.val[1]);

    resd_b0_r01_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd0_cr_64x2), vget_low_s64(resd1_cr_64x2)));
    resd_b0_r23_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd2_cr_64x2), vget_low_s64(resd3_cr_64x2)));
    resd_b1_r01_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd0_cr_64x2), vget_high_s64(resd1_cr_64x2)));
    resd_b1_r23_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd2_cr_64x2), vget_high_s64(resd3_cr_64x2)));
    resd_b2_r45_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd4_cr_64x2), vget_low_s64(resd5_cr_64x2)));
    resd_b2_r67_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_low_s64(resd6_cr_64x2), vget_low_s64(resd7_cr_64x2)));
    resd_b3_r45_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd4_cr_64x2), vget_high_s64(resd5_cr_64x2)));
    resd_b3_r67_cr = vreinterpretq_s16_s64(
        vcombine_s64(vget_high_s64(resd6_cr_64x2), vget_high_s64(resd7_cr_64x2)));

    dup_val_1 = vabsq_s16(resd_b0_r01_cr);
    dup_val_2 = vabsq_s16(resd_b0_r23_cr);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_cr);
    dup_val_2 = vabsq_s16(resd_b1_r23_cr);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_cr);
    dup_val_2 = vabsq_s16(resd_b2_r67_cr);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_cr);
    dup_val_2 = vabsq_s16(resd_b3_r67_cr);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz = ((nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 2) | (nnz_b3 << 3)) << 4);

    dup_val_1 = vabsq_s16(resd_b0_r01_cb);
    dup_val_2 = vabsq_s16(resd_b0_r23_cb);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b0 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b1_r01_cb);
    dup_val_2 = vabsq_s16(resd_b1_r23_cb);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b2_r45_cb);
    dup_val_2 = vabsq_s16(resd_b2_r67_cb);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b2 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    dup_val_1 = vabsq_s16(resd_b3_r45_cb);
    dup_val_2 = vabsq_s16(resd_b3_r67_cb);
    dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
    nnz_b3 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] || dup_abs[5] ||
             dup_abs[6] || dup_abs[7];

    nnz |= ((nnz_b0 | (nnz_b1 << 1) | (nnz_b2 << 2) | (nnz_b3 << 3)));
    return nnz;
}
