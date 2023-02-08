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
 *  isvcd_intra_resamp_neonintr.c
 *
 * @brief
 *  Contains routines that resample for SVC resampling
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_interpolate_base_luma_dyadic_neonintr()
 *  - isvcd_interpolate_intra_base_neonintr()
 *  - isvcd_horz_interpol_chroma_dyadic_1_neonintr()
 *  - isvcd_horz_interpol_chroma_dyadic_2_neonintr()
 *  - isvcd_vert_interpol_chroma_dyadic_1_neonintr()
 *  - isvcd_vert_interpol_chroma_dyadic_2_neonintr()
 *  - isvcd_vert_interpol_chroma_dyadic_3_neonintr()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
#include <assert.h>
#include <string.h>
#include <arm_neon.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"
#include "ih264_debug.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_base_luma_dyadic_neonintr               */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  intra resampling for dyadic scaling ratios               */
/*  Inputs        : pu1_inp_buf : ptr to the 12x12 reference sample buffer   */
/*                    pi2_tmp_filt_buf : ptr to the 12x16 buffer to hold the */
/*                        vertically interpolated data                       */
/*                  pu1_out_buf : output buffer pointer                      */
/*                  i4_out_stride : output buffer stride                     */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         05 21 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_base_luma_dyadic_neonintr(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                                 UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_y;
    WORD16 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp, *pu1_out;
    WORD16 *pi2_tmp;

    int16x4_t i4_rslt_vert_16x4_1, i4_rslt_vert_16x4_2;
    uint8x8_t i4_samp_vert_8x8_0, i4_samp_vert_8x8_1, i4_samp_vert_8x8_2, i4_samp_vert_8x8_3;
    int16x8_t i4_rslt_vert_16x8_0, i4_rslt_vert_16x8_2;
    /* Horizontal interpolation */
    int32x4_t const_512_32x4 = vdupq_n_s32(512);
    int32x4_t i4_rslt_horz_r0_1, i4_rslt_horz_r1_1, i4_rslt_horz_r0_2, i4_rslt_horz_r1_2;
    uint16x4_t i4_rslt_horz_r0_1_tmp, i4_rslt_horz_r1_1_tmp, i4_rslt_horz_r0_2_tmp,
        i4_rslt_horz_r1_2_tmp;
    uint16x8_t rslt_16x8_t_1, rslt_16x8_t_2;
    int32x4x2_t i4_rslt_horz_32x4x2_t;
    int16x4_t i4_samp_horz_16x4_0, i4_samp_horz_16x4_1, i4_samp_horz_16x4_2, i4_samp_horz_16x4_3,
        i4_samp_horz_16x4_4;
    int16x4_t i4_samp_horz_16x4_5, i4_samp_horz_16x4_6, i4_samp_horz_16x4_7, i4_samp_horz_16x4_8;
    int16_t i4_coeff_c0 = -3;
    int16_t i4_coeff_c1 = 28;
    int16_t i4_coeff_c2 = 8;
    int16_t i4_coeff_c3 = -1;
    int32x4_t i4_rslt_horz_r0_1_tmp32, i4_rslt_horz_r1_1_tmp32, i4_rslt_horz_r0_2_tmp32,
        i4_rslt_horz_r1_2_tmp32;

    /* Filter coefficient values for phase 4 */
    i4_coeff_0 = -3;
    i4_coeff_1 = 28;
    i4_coeff_2 = 8;
    i4_coeff_3 = -1;
    i4_filt_stride = 12;
    i4_src_stride = DYADIC_REF_W_Y;

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    pu1_out = pu1_out_buf;

    /* Vertical interpolation */
    // First 64 bits
    i4_samp_vert_8x8_0 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_1 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_2 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;

    i4_rslt_vert_16x8_0 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_3);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_2);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_1);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_0);

    vst1q_s16(pi2_tmp, i4_rslt_vert_16x8_0);
    pi2_tmp += i4_filt_stride;

    for(i4_y = 1; i4_y < 15; i4_y += 2)
    {
        i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
        i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
        i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
        i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);

        i4_rslt_vert_16x8_0 =
            vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_0);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_1);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_2);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_3);

        i4_rslt_vert_16x8_2 =
            vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_3);
        i4_rslt_vert_16x8_2 = vmlaq_n_s16(
            i4_rslt_vert_16x8_2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_2);
        i4_rslt_vert_16x8_2 = vmlaq_n_s16(
            i4_rslt_vert_16x8_2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_1);
        i4_rslt_vert_16x8_2 = vmlaq_n_s16(
            i4_rslt_vert_16x8_2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_0);

        /* Storing the results */
        vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_0));
        pi2_tmp += i4_filt_stride;
        vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_2));
        pi2_tmp += i4_filt_stride;
        pu1_inp += i4_src_stride;
    } /*End of Loop over y*/

    /* y = 15, y_phase = 4 */
    i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
    i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
    i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
    i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);

    i4_rslt_vert_16x8_0 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_0);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_1);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_2);
    i4_rslt_vert_16x8_0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_3);

    vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_0));
    /* End of loop over x */

    // Remaining 32 bits
    pu1_inp = pu1_inp_buf + 8;
    pi2_tmp = pi2_tmp_filt_buf + 8;

    i4_samp_vert_8x8_0 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_1 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_2 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);
    pu1_inp += i4_src_stride;

    i4_rslt_vert_16x4_1 =
        vmul_n_s16(vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_3);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))), i4_coeff_2);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))), i4_coeff_1);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))), i4_coeff_0);

    vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
    pi2_tmp += i4_filt_stride;

    for(i4_y = 1; i4_y < 15; i4_y += 2)
    {
        i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
        i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
        i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
        i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);

        i4_rslt_vert_16x4_1 = vmul_n_s16(
            vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_0);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))),
            i4_coeff_1);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))),
            i4_coeff_2);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))),
            i4_coeff_3);

        i4_rslt_vert_16x4_2 = vmul_n_s16(
            vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_3);
        i4_rslt_vert_16x4_2 = vmla_n_s16(
            i4_rslt_vert_16x4_2, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))),
            i4_coeff_2);
        i4_rslt_vert_16x4_2 = vmla_n_s16(
            i4_rslt_vert_16x4_2, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))),
            i4_coeff_1);
        i4_rslt_vert_16x4_2 = vmla_n_s16(
            i4_rslt_vert_16x4_2, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))),
            i4_coeff_0);

        vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
        pi2_tmp += i4_filt_stride;
        vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_2));
        pi2_tmp += i4_filt_stride;
        pu1_inp += i4_src_stride;
    }

    i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
    i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
    i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
    i4_samp_vert_8x8_3 = vld1_u8((const uint8_t *) pu1_inp);

    i4_rslt_vert_16x4_1 =
        vmul_n_s16(vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_0);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))), i4_coeff_1);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))), i4_coeff_2);
    i4_rslt_vert_16x4_1 =
        vmla_n_s16(i4_rslt_vert_16x4_1,
                   vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))), i4_coeff_3);

    vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
    /* Reinitializing the ptrs */
    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 16; i4_y++)
    {
        i4_samp_horz_16x4_0 = vld1_s16(pi2_tmp);
        i4_samp_horz_16x4_1 = vld1_s16(pi2_tmp + 1);
        i4_samp_horz_16x4_2 = vld1_s16(pi2_tmp + 2);
        i4_samp_horz_16x4_3 = vld1_s16(pi2_tmp + 3);
        i4_samp_horz_16x4_4 = vld1_s16(pi2_tmp + 4);
        i4_samp_horz_16x4_5 = vld1_s16(pi2_tmp + 5);
        i4_samp_horz_16x4_6 = vld1_s16(pi2_tmp + 6);
        i4_samp_horz_16x4_7 = vld1_s16(pi2_tmp + 7);
        i4_samp_horz_16x4_8 = vld1_s16(pi2_tmp + 8);

        i4_rslt_horz_r0_1 = vmull_n_s16(i4_samp_horz_16x4_0, i4_coeff_c3);  // a0c3 a1c3  a2c3  a3c3
        i4_rslt_horz_r0_1 = vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_1,
                                        i4_coeff_c2);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
        i4_rslt_horz_r0_1 = vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_2, i4_coeff_c1);
        i4_rslt_horz_r0_1 = vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_3, i4_coeff_c0);

        i4_rslt_horz_r1_1 = vmull_n_s16(i4_samp_horz_16x4_1, i4_coeff_c0);  // a0c0 a1c0  a2c0  a3c0
        i4_rslt_horz_r1_1 = vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_2,
                                        i4_coeff_c1);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
        i4_rslt_horz_r1_1 = vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_3, i4_coeff_c2);
        i4_rslt_horz_r1_1 = vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_4, i4_coeff_c3);

        i4_rslt_horz_r0_2 = vmull_n_s16(i4_samp_horz_16x4_4, i4_coeff_c3);  // a0c3 a1c3  a2c3  a3c3
        i4_rslt_horz_r0_2 = vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_5,
                                        i4_coeff_c2);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
        i4_rslt_horz_r0_2 = vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_6, i4_coeff_c1);
        i4_rslt_horz_r0_2 = vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_7, i4_coeff_c0);

        i4_rslt_horz_r1_2 = vmull_n_s16(i4_samp_horz_16x4_5, i4_coeff_c0);  // a0c0 a1c0  a2c0  a3c0
        i4_rslt_horz_r1_2 = vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_6,
                                        i4_coeff_c1);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1
        i4_rslt_horz_r1_2 = vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_7, i4_coeff_c2);
        i4_rslt_horz_r1_2 = vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_8, i4_coeff_c3);

        i4_rslt_horz_32x4x2_t = vzipq_s32(i4_rslt_horz_r0_1, i4_rslt_horz_r1_1);
        i4_rslt_horz_r0_1_tmp32 = i4_rslt_horz_32x4x2_t.val[0];  // 0 to 3
        i4_rslt_horz_r1_1_tmp32 = i4_rslt_horz_32x4x2_t.val[1];  // 4 to 7

        i4_rslt_horz_32x4x2_t = vzipq_s32(i4_rslt_horz_r0_2, i4_rslt_horz_r1_2);
        i4_rslt_horz_r0_2_tmp32 = i4_rslt_horz_32x4x2_t.val[0];  // 8 to 11
        i4_rslt_horz_r1_2_tmp32 = i4_rslt_horz_32x4x2_t.val[1];  // 12 to 15

        i4_rslt_horz_r0_1 = vaddq_s32(i4_rslt_horz_r0_1_tmp32, const_512_32x4);
        i4_rslt_horz_r1_1 = vaddq_s32(i4_rslt_horz_r1_1_tmp32, const_512_32x4);
        i4_rslt_horz_r0_2 = vaddq_s32(i4_rslt_horz_r0_2_tmp32, const_512_32x4);
        i4_rslt_horz_r1_2 = vaddq_s32(i4_rslt_horz_r1_2_tmp32, const_512_32x4);

        i4_rslt_horz_r0_1_tmp = vqshrun_n_s32(i4_rslt_horz_r0_1, 10);
        i4_rslt_horz_r1_1_tmp = vqshrun_n_s32(i4_rslt_horz_r1_1, 10);

        i4_rslt_horz_r0_2_tmp = vqshrun_n_s32(i4_rslt_horz_r0_2, 10);
        i4_rslt_horz_r1_2_tmp = vqshrun_n_s32(i4_rslt_horz_r1_2, 10);

        rslt_16x8_t_1 = vcombine_u16(i4_rslt_horz_r0_1_tmp, i4_rslt_horz_r1_1_tmp);  // 0 to 7
        rslt_16x8_t_2 = vcombine_u16(i4_rslt_horz_r0_2_tmp, i4_rslt_horz_r1_2_tmp);  // 8 to 15

        vst1_u8(pu1_out, vqmovn_u16(rslt_16x8_t_1));
        vst1_u8(pu1_out + 8, vqmovn_u16(rslt_16x8_t_2));

        pu1_out += i4_out_stride;
        pi2_tmp += i4_filt_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_intra_base_neonintr                     */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                    interpolation of a component to find the intra         */
/*                     resampled value                                       */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  pu1_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2009   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_intra_base_neonintr(void *pv_intra_samp_ctxt, UWORD8 *pu1_out,
                                           WORD32 i4_out_stride, WORD32 i4_refarray_wd,
                                           WORD32 i4_mb_x, WORD32 i4_mb_y, WORD32 i4_chroma_flag,
                                           WORD32 i4_refarray_flag)
{
    /* --------------------------------------------------------------------- */
    /* Index Parameters                                                      */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_map_ctxt_t *ps_map_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    UWORD8 *pu1_refarray = NULL;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;
    WORD32 i4_temp_array_ht;
    WORD32 *pi4_interp_buff;

    UWORD8 arr_y_ref_pos_luma[16] = {0};
    UWORD8 arr_x_ref_pos_luma[16] = {0};
    UWORD8 arr_x_ref_pos_luma_low[16] = {0};
    UWORD8 arr_x_ref_pos_luma_high[16] = {0};
    UWORD8 arr_phase_luma[16] = {0};
    UWORD8 *pi4_y_ref_pos_luma;
    UWORD8 *pi4_x_ref_pos_luma_low;
    UWORD8 *pi4_x_ref_pos_luma_high;
    UWORD8 *pi4_phase_luma;
    WORD16 *pi2_interp_buff_temp;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD32 i4_x_min;
    ref_min_max_map_t *ps_x_min_max;
    UWORD8 *pu1_refarray_temp;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    if(0 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_buffer;
    }
    else if(1 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_cb;
    }

    /* --------------------------------------------------------------------- */
    /* LUMA or CHROMA                                                        */
    /* --------------------------------------------------------------------- */

    if(1 == i4_chroma_flag)
        ps_map_ctxt = &(ps_lyr_ctxt->s_chroma_map_ctxt);
    else
        ps_map_ctxt = &(ps_lyr_ctxt->s_luma_map_ctxt);

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    ps_x_min_max = ps_map_ctxt->ps_x_min_max;

    i4_frm_mb_y = i4_mb_y * i4_mb_ht;
    i4_frm_mb_x = i4_mb_x * i4_mb_wd;
    /* get the min and max positions */
    i4_x_min = ps_x_min_max[i4_mb_x].i2_min_pos;

    /* --------------------------------------------------------------------- */
    /* Projected frame level pointers                                        */
    /* --------------------------------------------------------------------- */
    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

    /* --------------------------------------------------------------------- */
    /* Pointers and Dimenstion of the temporary buffer                       */
    /* --------------------------------------------------------------------- */
    i4_temp_array_ht = i4_mb_ht;
    pi4_interp_buff = ps_ctxt->pi4_temp_interpolation_buffer;
    pi2_interp_buff_temp = (WORD16 *) pi4_interp_buff;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation in vertical direction                          */
    /* --------------------------------------------------------------------- */
    if(i4_chroma_flag == 0)
    {
        {
            uint8x8_t inp_8x8_r0, inp_8x8_r0_1;
            uint8x8_t inp_8x8_r1, inp_8x8_r1_1;
            uint8x8_t inp_8x8_r2, inp_8x8_r2_1;
            uint8x8_t inp_8x8_r3, inp_8x8_r3_1;
            int16x8_t out_res_16x8_r0_0, out_res_16x8_r0_1;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                arr_phase_luma[i4_y] = (UWORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
                arr_y_ref_pos_luma[i4_y] = (UWORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
            }
            pi4_y_ref_pos_luma = arr_y_ref_pos_luma;
            pi4_phase_luma = arr_phase_luma;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                pu1_refarray_temp =
                    pu1_refarray + (pi4_y_ref_pos_luma[i4_y] * i4_refarray_wd) + (i4_x_min - 1);
                inp_8x8_r0 = vld1_u8((pu1_refarray_temp - i4_refarray_wd));
                inp_8x8_r1 = vld1_u8((pu1_refarray_temp));
                inp_8x8_r2 = vld1_u8((pu1_refarray_temp + i4_refarray_wd));
                inp_8x8_r3 = vld1_u8((pu1_refarray_temp + 2 * i4_refarray_wd));

                inp_8x8_r0_1 = vld1_u8((pu1_refarray_temp + 8 - i4_refarray_wd));
                inp_8x8_r1_1 = vld1_u8((pu1_refarray_temp + 8));
                inp_8x8_r2_1 = vld1_u8((pu1_refarray_temp + 8 + i4_refarray_wd));
                inp_8x8_r3_1 = vld1_u8((pu1_refarray_temp + 8 + 2 * i4_refarray_wd));

                out_res_16x8_r0_0 = vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r0)),
                                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y]]);
                out_res_16x8_r0_0 =
                    vmlaq_n_s16(out_res_16x8_r0_0, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r1)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 16]);
                out_res_16x8_r0_0 =
                    vmlaq_n_s16(out_res_16x8_r0_0, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r2)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 32]);
                out_res_16x8_r0_0 =
                    vmlaq_n_s16(out_res_16x8_r0_0, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r3)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 48]);

                out_res_16x8_r0_1 = vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r0_1)),
                                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y]]);
                out_res_16x8_r0_1 =
                    vmlaq_n_s16(out_res_16x8_r0_1, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r1_1)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 16]);
                out_res_16x8_r0_1 =
                    vmlaq_n_s16(out_res_16x8_r0_1, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r2_1)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 32]);
                out_res_16x8_r0_1 =
                    vmlaq_n_s16(out_res_16x8_r0_1, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r3_1)),
                                g_ai1_interp_filter_luma[pi4_phase_luma[i4_y] + 48]);

                vst1q_s16((pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1)),
                          out_res_16x8_r0_0);
                vst1q_s16((pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1) + 8),
                          out_res_16x8_r0_1);
            }
        }
        /*Horizontal Interpolation*/
        {
            WORD32 strt_indx = 10;

            uint8x16_t phs_mask_8x8_0;
            uint8x16_t x_ref_pos_luma_mask_r0_0;
            uint8x16_t x_ref_pos_luma_mask_r0_1;
            uint8x16_t x_ref_pos_luma_mask_r1_0;
            uint8x16_t x_ref_pos_luma_mask_r1_1;
            uint8x16_t x_ref_pos_luma_mask_r2_0;
            uint8x16_t x_ref_pos_luma_mask_r2_1;
            uint8x16_t x_ref_pos_luma_mask_r3_0;
            uint8x16_t x_ref_pos_luma_mask_r3_1;

            WORD32 strt_indx_h = 0, i4_x2 = 0;
            WORD32 i4_mb_wd_hlf = (i4_mb_wd >> 1);
            uint8x16_t twos = vdupq_n_u8(2);
            strt_indx = ps_x_pos_phase[0 + i4_frm_mb_x].i2_ref_pos - 1;
            strt_indx_h = (ps_x_pos_phase[8 + i4_frm_mb_x].i2_ref_pos - strt_indx - 1);
            for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
            {
                arr_x_ref_pos_luma[i4_x] = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
                arr_phase_luma[i4_x] = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
                arr_x_ref_pos_luma[i4_x] = arr_x_ref_pos_luma[i4_x] - strt_indx - 1;
            }

            for(i4_x = 0; i4_x < i4_mb_wd_hlf; i4_x++)
            {
                i4_x2 = i4_x << 1;
                arr_x_ref_pos_luma_low[i4_x2] = (arr_x_ref_pos_luma[i4_x]) << 1;
                arr_x_ref_pos_luma_low[i4_x2 + 1] = arr_x_ref_pos_luma_low[i4_x2] + 1;
            }
            for(i4_x = i4_mb_wd_hlf; i4_x < i4_mb_wd; i4_x++)
            {
                i4_x2 = (i4_x - i4_mb_wd_hlf) << 1;
                arr_x_ref_pos_luma_high[i4_x2] = ((arr_x_ref_pos_luma[i4_x] - strt_indx_h) << 1);
                arr_x_ref_pos_luma_high[i4_x2 + 1] = arr_x_ref_pos_luma_high[i4_x2] + 1;
            }
            pi4_x_ref_pos_luma_low = arr_x_ref_pos_luma_low;
            pi4_x_ref_pos_luma_high = arr_x_ref_pos_luma_high;
            pi4_phase_luma = arr_phase_luma;

            phs_mask_8x8_0 = vld1q_u8((const uint8_t *) pi4_phase_luma);

            x_ref_pos_luma_mask_r0_0 = vld1q_u8(pi4_x_ref_pos_luma_low);
            x_ref_pos_luma_mask_r0_1 = vld1q_u8(pi4_x_ref_pos_luma_high);
            x_ref_pos_luma_mask_r1_0 = vaddq_u8(x_ref_pos_luma_mask_r0_0, twos);
            x_ref_pos_luma_mask_r1_1 = vaddq_u8(x_ref_pos_luma_mask_r0_1, twos);
            x_ref_pos_luma_mask_r2_0 = vaddq_u8(x_ref_pos_luma_mask_r1_0, twos);
            x_ref_pos_luma_mask_r2_1 = vaddq_u8(x_ref_pos_luma_mask_r1_1, twos);
            x_ref_pos_luma_mask_r3_0 = x_ref_pos_luma_mask_r0_0;
            x_ref_pos_luma_mask_r3_1 = x_ref_pos_luma_mask_r0_1;

            {
                int8x16_t ip_filt_8x16_r0;
                int8x16_t ip_filt_8x16_r1;
                int8x16_t ip_filt_8x16_r2;
                int8x16_t ip_filt_8x16_r3;

                int16x8_t ip_filt_16x8_r0_0, ip_filt_16x8_r0_1;
                int16x8_t ip_filt_16x8_r1_0, ip_filt_16x8_r1_1;
                int16x8_t ip_filt_16x8_r2_0, ip_filt_16x8_r2_1;
                int16x8_t ip_filt_16x8_r3_0, ip_filt_16x8_r3_1;

                int16x8_t inp_16x8_0;
                int16x8_t inp_16x8_1;
                int16x8_t inp_16x8_2;
                int16x8_t inp_16x8_3;

                int16x8_t inp_16x8_r0_0, inp_16x8_r2_0;
                int16x8_t inp_16x8_r0_1, inp_16x8_r2_1;
                int16x8_t inp_16x8_r1_0, inp_16x8_r3_0;
                int16x8_t inp_16x8_r1_1, inp_16x8_r3_1;

                int16x4_t inp_16x4_r0_0, inp_16x4_r2_0;
                int16x4_t inp_16x4_r0_1, inp_16x4_r2_1;
                int16x4_t inp_16x4_r1_0, inp_16x4_r3_0;
                int16x4_t inp_16x4_r1_1, inp_16x4_r3_1;

                int32x4_t out_res_32x4_r0_l_0;
                int32x4_t out_res_32x4_r0_l_1;
                int32x4_t out_res_32x4_r0_h_0;
                int32x4_t out_res_32x4_r0_h_1;

                uint16x4_t out_res_16x4_r0_l_0;
                uint16x4_t out_res_16x4_r0_l_1;
                uint16x4_t out_res_16x4_r0_h_0;
                uint16x4_t out_res_16x4_r0_h_1;

                uint8x8_t out_res_8x8_r0_l, out_res_8x8_r0_h;
                uint8x8x2_t u1_temp_8x8x2_t;
                uint8x8_t u1_temp_8x8_t0, u1_temp_8x8_t1;

                ip_filt_8x16_r0 = vld1q_s8((g_ai1_interp_filter_luma));
                ip_filt_8x16_r1 = vld1q_s8((g_ai1_interp_filter_luma + 16));
                ip_filt_8x16_r2 = vld1q_s8((g_ai1_interp_filter_luma + 32));
                ip_filt_8x16_r3 = vld1q_s8((g_ai1_interp_filter_luma + 48));

                u1_temp_8x8x2_t.val[0] = vreinterpret_u8_s8(vget_low_s8(ip_filt_8x16_r0));
                u1_temp_8x8x2_t.val[1] = vreinterpret_u8_s8(vget_high_s8(ip_filt_8x16_r0));
                u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r0 = vcombine_s8(vreinterpret_s8_u8(u1_temp_8x8_t0),
                                              vreinterpret_s8_u8(u1_temp_8x8_t1));

                u1_temp_8x8x2_t.val[0] = vreinterpret_u8_s8(vget_low_s8(ip_filt_8x16_r1));
                u1_temp_8x8x2_t.val[1] = vreinterpret_u8_s8(vget_high_s8(ip_filt_8x16_r1));
                u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r1 = vcombine_s8(vreinterpret_s8_u8(u1_temp_8x8_t0),
                                              vreinterpret_s8_u8(u1_temp_8x8_t1));
                u1_temp_8x8x2_t.val[0] = vreinterpret_u8_s8(vget_low_s8(ip_filt_8x16_r2));
                u1_temp_8x8x2_t.val[1] = vreinterpret_u8_s8(vget_high_s8(ip_filt_8x16_r2));
                u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r2 = vcombine_s8(vreinterpret_s8_u8(u1_temp_8x8_t0),
                                              vreinterpret_s8_u8(u1_temp_8x8_t1));
                u1_temp_8x8x2_t.val[0] = vreinterpret_u8_s8(vget_low_s8(ip_filt_8x16_r3));
                u1_temp_8x8x2_t.val[1] = vreinterpret_u8_s8(vget_high_s8(ip_filt_8x16_r3));
                u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r3 = vcombine_s8(vreinterpret_s8_u8(u1_temp_8x8_t0),
                                              vreinterpret_s8_u8(u1_temp_8x8_t1));
                ip_filt_16x8_r0_0 = vmovl_s8(vget_low_s8(ip_filt_8x16_r0));
                ip_filt_16x8_r1_0 = vmovl_s8(vget_low_s8(ip_filt_8x16_r1));
                ip_filt_16x8_r2_0 = vmovl_s8(vget_low_s8(ip_filt_8x16_r2));
                ip_filt_16x8_r3_0 = vmovl_s8(vget_low_s8(ip_filt_8x16_r3));
                ip_filt_16x8_r0_1 = vmovl_s8(vget_high_s8(ip_filt_8x16_r0));
                ip_filt_16x8_r1_1 = vmovl_s8(vget_high_s8(ip_filt_8x16_r1));
                ip_filt_16x8_r2_1 = vmovl_s8(vget_high_s8(ip_filt_8x16_r2));
                ip_filt_16x8_r3_1 = vmovl_s8(vget_high_s8(ip_filt_8x16_r3));

                for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
                {
                    inp_16x8_0 = vld1q_s16((pi2_interp_buff_temp + strt_indx));
                    inp_16x8_1 = vld1q_s16((pi2_interp_buff_temp + strt_indx + strt_indx_h));
                    inp_16x8_2 = vld1q_s16((pi2_interp_buff_temp + strt_indx + 3));
                    inp_16x8_3 = vld1q_s16((pi2_interp_buff_temp + strt_indx + strt_indx_h + 3));
                    pi2_interp_buff_temp += i4_refarray_wd;
                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r0_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r0_0));
                    inp_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r0_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r0_1));
                    inp_16x8_r0_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r1_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r1_0));
                    inp_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r1_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r1_1));
                    inp_16x8_r1_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r2_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r2_0));
                    inp_16x8_r2_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_1)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r2_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r2_1));
                    inp_16x8_r2_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_2)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_2)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r3_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r3_0));
                    inp_16x8_r3_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_3)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_3)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r3_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r3_1));
                    inp_16x8_r3_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    inp_16x4_r0_0 = vget_low_s16(inp_16x8_r0_0);
                    inp_16x4_r0_1 = vget_low_s16(inp_16x8_r0_1);
                    inp_16x4_r1_0 = vget_low_s16(inp_16x8_r1_0);
                    inp_16x4_r1_1 = vget_low_s16(inp_16x8_r1_1);

                    inp_16x4_r2_0 = vget_low_s16(inp_16x8_r2_0);
                    inp_16x4_r2_1 = vget_low_s16(inp_16x8_r2_1);
                    inp_16x4_r3_0 = vget_low_s16(inp_16x8_r3_0);
                    inp_16x4_r3_1 = vget_low_s16(inp_16x8_r3_1);

                    out_res_32x4_r0_l_0 = vmull_s16(inp_16x4_r0_0, vget_low_s16(ip_filt_16x8_r0_0));
                    out_res_32x4_r0_l_0 = vmlal_s16(out_res_32x4_r0_l_0, inp_16x4_r1_0,
                                                    vget_low_s16(ip_filt_16x8_r1_0));
                    out_res_32x4_r0_l_0 = vmlal_s16(out_res_32x4_r0_l_0, inp_16x4_r2_0,
                                                    vget_low_s16(ip_filt_16x8_r2_0));
                    out_res_32x4_r0_l_0 = vmlal_s16(out_res_32x4_r0_l_0, inp_16x4_r3_0,
                                                    vget_low_s16(ip_filt_16x8_r3_0));
                    out_res_32x4_r0_l_1 =
                        vmull_s16(vget_high_s16(inp_16x8_r0_0), vget_high_s16(ip_filt_16x8_r0_0));
                    out_res_32x4_r0_l_1 =
                        vmlal_s16(out_res_32x4_r0_l_1, vget_high_s16(inp_16x8_r1_0),
                                  vget_high_s16(ip_filt_16x8_r1_0));
                    out_res_32x4_r0_l_1 =
                        vmlal_s16(out_res_32x4_r0_l_1, vget_high_s16(inp_16x8_r2_0),
                                  vget_high_s16(ip_filt_16x8_r2_0));
                    out_res_32x4_r0_l_1 =
                        vmlal_s16(out_res_32x4_r0_l_1, vget_high_s16(inp_16x8_r3_0),
                                  vget_high_s16(ip_filt_16x8_r3_0));

                    out_res_32x4_r0_h_0 = vmull_s16(inp_16x4_r0_1, vget_low_s16(ip_filt_16x8_r0_1));
                    out_res_32x4_r0_h_0 = vmlal_s16(out_res_32x4_r0_h_0, inp_16x4_r1_1,
                                                    vget_low_s16(ip_filt_16x8_r1_1));
                    out_res_32x4_r0_h_0 = vmlal_s16(out_res_32x4_r0_h_0, inp_16x4_r2_1,
                                                    vget_low_s16(ip_filt_16x8_r2_1));
                    out_res_32x4_r0_h_0 = vmlal_s16(out_res_32x4_r0_h_0, inp_16x4_r3_1,
                                                    vget_low_s16(ip_filt_16x8_r3_1));

                    out_res_32x4_r0_h_1 =
                        vmull_s16(vget_high_s16(inp_16x8_r0_1), vget_high_s16(ip_filt_16x8_r0_1));
                    out_res_32x4_r0_h_1 =
                        vmlal_s16(out_res_32x4_r0_h_1, vget_high_s16(inp_16x8_r1_1),
                                  vget_high_s16(ip_filt_16x8_r1_1));
                    out_res_32x4_r0_h_1 =
                        vmlal_s16(out_res_32x4_r0_h_1, vget_high_s16(inp_16x8_r2_1),
                                  vget_high_s16(ip_filt_16x8_r2_1));
                    out_res_32x4_r0_h_1 =
                        vmlal_s16(out_res_32x4_r0_h_1, vget_high_s16(inp_16x8_r3_1),
                                  vget_high_s16(ip_filt_16x8_r3_1));

                    out_res_16x4_r0_l_0 = vqrshrun_n_s32(out_res_32x4_r0_l_0, 10);
                    out_res_16x4_r0_l_1 = vqrshrun_n_s32(out_res_32x4_r0_l_1, 10);
                    out_res_16x4_r0_h_0 = vqrshrun_n_s32(out_res_32x4_r0_h_0, 10);
                    out_res_16x4_r0_h_1 = vqrshrun_n_s32(out_res_32x4_r0_h_1, 10);

                    out_res_8x8_r0_l =
                        vqmovn_u16(vcombine_u16(out_res_16x4_r0_l_0, out_res_16x4_r0_l_1));
                    out_res_8x8_r0_h =
                        vqmovn_u16(vcombine_u16(out_res_16x4_r0_h_0, out_res_16x4_r0_h_1));
                    vst1q_u8((pu1_out + (i4_y * i4_out_stride)),
                             vcombine_u8(out_res_8x8_r0_l, out_res_8x8_r0_h));
                }
            }
        }
    }
    else
    {
        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            arr_y_ref_pos_luma[i4_y] = (UWORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos;
            arr_phase_luma[i4_y] = (UWORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
        }
        pi4_y_ref_pos_luma = arr_y_ref_pos_luma;
        pi4_phase_luma = arr_phase_luma;

        {
            uint8x8_t inp_8x8_r0, inp_8x8_r0_1;
            uint8x8_t inp_8x8_r1, inp_8x8_r1_1;
            int16x8_t out_res_16x8_r0_0, out_res_16x8_r0_1;

            for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
            {
                pu1_refarray_temp =
                    pu1_refarray + (pi4_y_ref_pos_luma[i4_y] * i4_refarray_wd) + (i4_x_min - 1);
                inp_8x8_r0 = vld1_u8((pu1_refarray_temp));
                inp_8x8_r1 = vld1_u8((pu1_refarray_temp + i4_refarray_wd));

                inp_8x8_r0_1 = vld1_u8((pu1_refarray_temp + 8));
                inp_8x8_r1_1 = vld1_u8((pu1_refarray_temp + 8 + i4_refarray_wd));

                out_res_16x8_r0_0 = vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r0)),
                                                g_au1_interp_filter_chroma[pi4_phase_luma[i4_y]]);
                out_res_16x8_r0_0 =
                    vmlaq_n_s16(out_res_16x8_r0_0, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r1)),
                                g_au1_interp_filter_chroma[pi4_phase_luma[i4_y] + 16]);

                out_res_16x8_r0_1 = vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r0_1)),
                                                g_au1_interp_filter_chroma[pi4_phase_luma[i4_y]]);
                out_res_16x8_r0_1 =
                    vmlaq_n_s16(out_res_16x8_r0_1, vreinterpretq_s16_u16(vmovl_u8(inp_8x8_r1_1)),
                                g_au1_interp_filter_chroma[pi4_phase_luma[i4_y] + 16]);

                vst1q_s16((pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1)),
                          out_res_16x8_r0_0);
                vst1q_s16((pi2_interp_buff_temp + (i4_y * i4_refarray_wd) + (i4_x_min - 1) + 8),
                          out_res_16x8_r0_1);
            }
        }

        {
            WORD32 strt_indx = 10;

            uint8x16_t phs_mask_8x8_0;
            uint8x16_t x_ref_pos_luma_mask_r0_0;
            uint8x16_t x_ref_pos_luma_mask_r1_0;

            WORD32 i4_x2 = 0;
            uint8x16_t twos = vdupq_n_u8(2);
            strt_indx = ps_x_pos_phase[0 + i4_frm_mb_x].i2_ref_pos;

            for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
            {
                arr_x_ref_pos_luma[i4_x] = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
                arr_phase_luma[i4_x] = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
                arr_x_ref_pos_luma[i4_x] = arr_x_ref_pos_luma[i4_x] - strt_indx;
                i4_x2 = i4_x << 1;
                arr_x_ref_pos_luma_low[i4_x2] = (arr_x_ref_pos_luma[i4_x]) << 1;
                arr_x_ref_pos_luma_low[i4_x2 + 1] = arr_x_ref_pos_luma_low[i4_x2] + 1;
            }

            pi4_x_ref_pos_luma_low = arr_x_ref_pos_luma_low;
            pi4_phase_luma = arr_phase_luma;

            phs_mask_8x8_0 = vld1q_u8(pi4_phase_luma);
            x_ref_pos_luma_mask_r0_0 = vld1q_u8(pi4_x_ref_pos_luma_low);
            x_ref_pos_luma_mask_r1_0 = vaddq_u8(x_ref_pos_luma_mask_r0_0, twos);

            {
                uint8x16_t ip_filt_8x16_r0;
                uint8x16_t ip_filt_8x16_r1;
                int16x8_t ip_filt_16x8_r0_0;
                int16x8_t ip_filt_16x8_r1_0;
                int16x8_t inp_16x8_0;
                int16x8_t inp_16x8_r0_0;
                int16x8_t inp_16x8_r1_0;
                int16x4_t inp_16x4_r0_0;
                int16x4_t inp_16x4_r1_0;
                int32x4_t out_res_32x4_r0_l_0;
                int32x4_t out_res_32x4_r0_l_1;
                uint16x4_t out_res_16x4_r0_l_0;
                uint16x4_t out_res_16x4_r0_l_1;
                uint16x8_t out_res_16x8_r0_l;
                uint8x16_t out_8x16_r0;
                uint8x8x2_t u1_incr_8x8x2_t;
                uint8x8_t u1_incr_8x8_t0, u1_incr_8x8_t1;
                uint8x8x2_t u1_temp_8x8x2_t;
                uint8x8_t u1_temp_8x8_t0, u1_temp_8x8_t1;
                uint8x16_t chroma_mask_8x16 = vreinterpretq_u8_u16(vdupq_n_u16(0x00ff));

                ip_filt_8x16_r0 = vld1q_u8((g_au1_interp_filter_chroma));
                ip_filt_8x16_r1 = vld1q_u8((g_au1_interp_filter_chroma + 16));

                u1_incr_8x8x2_t.val[0] = vget_low_u8(ip_filt_8x16_r0);
                u1_incr_8x8x2_t.val[1] = vget_high_u8(ip_filt_8x16_r0);
                u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                u1_incr_8x8x2_t.val[0] = vget_low_u8(ip_filt_8x16_r1);
                u1_incr_8x8x2_t.val[1] = vget_high_u8(ip_filt_8x16_r1);
                u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(phs_mask_8x8_0));
                u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(phs_mask_8x8_0));
                ip_filt_8x16_r1 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                ip_filt_16x8_r0_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(ip_filt_8x16_r0)));
                ip_filt_16x8_r1_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(ip_filt_8x16_r1)));

                for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
                {
                    inp_16x8_0 = vld1q_s16((pi2_interp_buff_temp + strt_indx));
                    pi2_interp_buff_temp += i4_refarray_wd;
                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r0_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r0_0));
                    inp_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(inp_16x8_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_luma_mask_r1_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_luma_mask_r1_0));
                    inp_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));
                    inp_16x4_r0_0 = vget_low_s16(inp_16x8_r0_0);
                    inp_16x4_r1_0 = vget_low_s16(inp_16x8_r1_0);

                    out_res_32x4_r0_l_0 = vmull_s16(inp_16x4_r0_0, vget_low_s16(ip_filt_16x8_r0_0));
                    out_res_32x4_r0_l_0 = vmlal_s16(out_res_32x4_r0_l_0, inp_16x4_r1_0,
                                                    vget_low_s16(ip_filt_16x8_r1_0));
                    out_res_32x4_r0_l_1 =
                        vmull_s16(vget_high_s16(inp_16x8_r0_0), vget_high_s16(ip_filt_16x8_r0_0));
                    out_res_32x4_r0_l_1 =
                        vmlal_s16(out_res_32x4_r0_l_1, vget_high_s16(inp_16x8_r1_0),
                                  vget_high_s16(ip_filt_16x8_r1_0));

                    out_res_16x4_r0_l_0 = vqrshrun_n_s32(out_res_32x4_r0_l_0, 10);
                    out_res_16x4_r0_l_1 = vqrshrun_n_s32(out_res_32x4_r0_l_1, 10);
                    out_res_16x8_r0_l = vcombine_u16(out_res_16x4_r0_l_0, out_res_16x4_r0_l_1);
                    out_8x16_r0 = vld1q_u8(pu1_out + (i4_y * i4_out_stride));
                    out_8x16_r0 = vbslq_u8(chroma_mask_8x16,
                                           vreinterpretq_u8_u16(out_res_16x8_r0_l), out_8x16_r0);
                    vst1q_u8((pu1_out + (i4_y * i4_out_stride)), out_8x16_r0);
                }
            }
        }
    }
    return;
} /* End of Interpolation Function */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_1_neonintr              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                    interpolation of a component to find the intra         */
/*                     resampled value                                       */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  pu1_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation on horizontal direction        */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2009   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/

void isvcd_horz_interpol_chroma_dyadic_1_neonintr(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                                  WORD32 i4_out_stride, WORD32 i4_phase_0,
                                                  WORD32 i4_phase_1)
{
    WORD32 i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_dst_stride;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    int16x8_t i4_samp_horz_16x8_r0_0, i4_samp_horz_16x8_r0_1, i4_samp_horz_16x8_r0_2;
    int16x8_t i4_samp_horz_16x8_r1_0, i4_samp_horz_16x8_r1_1, i4_samp_horz_16x8_r1_2;
    int16x8_t i4_rslt_horz_r0_1, i4_rslt_horz_r0_2;
    int16x8_t i4_rslt_horz_r1_1, i4_rslt_horz_r1_2;

    int16x8_t final_horz_16x8_r0_1;
    int16x8_t final_horz_16x8_r1_1;

    uint8x16_t i4_out_horz_8x16_r0, i4_out_horz_8x16_r1;
    uint8x16_t chroma_mask_8x16 = vreinterpretq_u8_u16(vdupq_n_u16(0x00ff));

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_dst_stride = i4_out_stride;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y += 2)
    {
        i4_samp_horz_16x8_r0_0 = vld1q_s16(pi2_tmp);      // a0 a1 a2 a3 a4 a5 a6 a7
        i4_samp_horz_16x8_r0_1 = vld1q_s16(pi2_tmp + 1);  // a1 a2 a3 a4
        i4_samp_horz_16x8_r0_2 = vld1q_s16(pi2_tmp + 2);  // a2 a3 a4 a5

        i4_samp_horz_16x8_r1_0 = vld1q_s16(pi2_tmp + i4_filt_stride);
        i4_samp_horz_16x8_r1_1 = vld1q_s16(pi2_tmp + i4_filt_stride + 1);
        i4_samp_horz_16x8_r1_2 = vld1q_s16(pi2_tmp + (i4_filt_stride + 2));

        i4_rslt_horz_r0_1 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_0, i4_coeff_0);  // a0c0 a1c0  a2c0 a3c0

        i4_rslt_horz_r0_2 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_1, i4_coeff_2);  // a1c2 a2c2  a3c2 a4c2
        i4_rslt_horz_r0_1 = vmlaq_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x8_r0_1,
                                        i4_coeff_1);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1

        i4_rslt_horz_r0_2 = vmlaq_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x8_r0_2,
                                        i4_coeff_3);  // a1c2+a2c3  a2c2+a3c3 a3c2+a4c3 a4c2+a5c3

        i4_rslt_horz_r1_1 = vmulq_n_s16(i4_samp_horz_16x8_r1_0, i4_coeff_0);
        i4_rslt_horz_r1_2 = vmulq_n_s16(i4_samp_horz_16x8_r1_1, i4_coeff_2);

        i4_rslt_horz_r1_1 = vmlaq_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x8_r1_1, i4_coeff_1);
        i4_rslt_horz_r1_2 = vmlaq_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x8_r1_2, i4_coeff_3);

        final_horz_16x8_r0_1 = vzipq_s16(i4_rslt_horz_r0_1, i4_rslt_horz_r0_2).val[0];
        final_horz_16x8_r1_1 = vzipq_s16(i4_rslt_horz_r1_1, i4_rslt_horz_r1_2).val[0];

        final_horz_16x8_r0_1 = vrshrq_n_s16(final_horz_16x8_r0_1, 6);

        final_horz_16x8_r1_1 = vrshrq_n_s16(final_horz_16x8_r1_1, 6);

        i4_out_horz_8x16_r0 = vld1q_u8(pu1_out);
        i4_out_horz_8x16_r1 = vld1q_u8(pu1_out + i4_dst_stride);

        i4_out_horz_8x16_r0 = vbslq_u8(chroma_mask_8x16, vreinterpretq_u8_s16(final_horz_16x8_r0_1),
                                       i4_out_horz_8x16_r0);
        i4_out_horz_8x16_r1 = vbslq_u8(chroma_mask_8x16, vreinterpretq_u8_s16(final_horz_16x8_r1_1),
                                       i4_out_horz_8x16_r1);

        vst1q_u8(pu1_out, i4_out_horz_8x16_r0);
        vst1q_u8(pu1_out + i4_dst_stride, i4_out_horz_8x16_r1);

        /* Incrementing ptr */
        pi2_tmp += (i4_filt_stride << 1);
        pu1_out += (i4_dst_stride << 1);

    } /* End of loop over y */
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_2_neonintr              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                    chroma_phase_y_plus1:                                  */
/*                        ref_lyr        cur_lyr                             */
/*                            0            1                                 */
/*                            0            2                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_horz_interpol_chroma_dyadic_2_neonintr(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                                  WORD32 i4_out_stride, WORD32 i4_phase_0,
                                                  WORD32 i4_phase_1)
{
    WORD32 i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_dst_stride;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    int16x8_t i4_samp_horz_16x8_r0_0, i4_samp_horz_16x8_r0_1;
    int16x8_t i4_samp_horz_16x8_r1_0, i4_samp_horz_16x8_r1_1;
    int16x8_t i4_rslt_horz_r0_1, i4_rslt_horz_r0_2;
    int16x8_t i4_rslt_horz_r1_1, i4_rslt_horz_r1_2;

    int16x8_t final_horz_16x8_r0_1;
    int16x8_t final_horz_16x8_r1_1;

    uint8x16_t i4_out_horz_8x16_r0, i4_out_horz_8x16_r1;
    uint8x16_t chroma_mask_8x16 = vreinterpretq_u8_u16(vdupq_n_u16(0x00ff));

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf + 1;
    i4_filt_stride = 6;
    i4_dst_stride = i4_out_stride;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y += 2)
    {
        i4_samp_horz_16x8_r0_0 = vld1q_s16(pi2_tmp);      // a0 a1 a2 a3 a4 a5 a6 a7
        i4_samp_horz_16x8_r0_1 = vld1q_s16(pi2_tmp + 1);  // a1 a2 a3 a4

        i4_samp_horz_16x8_r1_0 = vld1q_s16(pi2_tmp + i4_filt_stride);
        i4_samp_horz_16x8_r1_1 = vld1q_s16(pi2_tmp + i4_filt_stride + 1);

        i4_rslt_horz_r0_1 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_0, i4_coeff_0);  // a0c0 a1c0  a2c0 a3c0

        i4_rslt_horz_r0_2 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_0, i4_coeff_2);  // a1c2 a2c2  a3c2 a4c2
        i4_rslt_horz_r0_1 = vmlaq_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x8_r0_1,
                                        i4_coeff_1);  // a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1

        i4_rslt_horz_r0_2 = vmlaq_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x8_r0_1,
                                        i4_coeff_3);  // a1c2+a2c3  a2c2+a3c3 a3c2+a4c3 a4c2+a5c3

        i4_rslt_horz_r1_1 = vmulq_n_s16(i4_samp_horz_16x8_r1_0, i4_coeff_0);
        i4_rslt_horz_r1_2 = vmulq_n_s16(i4_samp_horz_16x8_r1_0, i4_coeff_2);

        i4_rslt_horz_r1_1 = vmlaq_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x8_r1_1, i4_coeff_1);
        i4_rslt_horz_r1_2 = vmlaq_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x8_r1_1, i4_coeff_3);

        final_horz_16x8_r0_1 = vzipq_s16(i4_rslt_horz_r0_1, i4_rslt_horz_r0_2).val[0];
        final_horz_16x8_r1_1 = vzipq_s16(i4_rslt_horz_r1_1, i4_rslt_horz_r1_2).val[0];

        final_horz_16x8_r0_1 = vrshrq_n_s16(final_horz_16x8_r0_1, 6);

        final_horz_16x8_r1_1 = vrshrq_n_s16(final_horz_16x8_r1_1, 6);

        i4_out_horz_8x16_r0 = vld1q_u8(pu1_out);
        i4_out_horz_8x16_r1 = vld1q_u8(pu1_out + i4_dst_stride);

        i4_out_horz_8x16_r0 = vbslq_u8(chroma_mask_8x16, vreinterpretq_u8_s16(final_horz_16x8_r0_1),
                                       i4_out_horz_8x16_r0);
        i4_out_horz_8x16_r1 = vbslq_u8(chroma_mask_8x16, vreinterpretq_u8_s16(final_horz_16x8_r1_1),
                                       i4_out_horz_8x16_r1);

        vst1q_u8(pu1_out, i4_out_horz_8x16_r0);
        vst1q_u8(pu1_out + i4_dst_stride, i4_out_horz_8x16_r1);

        /* Incrementing ptr */
        pi2_tmp += (i4_filt_stride << 1);
        pu1_out += (i4_dst_stride << 1);

    } /* End of loop over y */
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_1_neonintr              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold       */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_1_neonintr(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                                  WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    uint8x8_t i4_samp_vert_8x8_r0, i4_samp_vert_8x8_r1, i4_samp_vert_8x8_r2;
    uint8x8_t i4_samp_vert_8x8_r3, i4_samp_vert_8x8_r4, i4_samp_vert_8x8_r5;
    int16x8_t i4_rslt_vert_16x8_r0, i4_rslt_vert_16x8_r1, i4_rslt_vert_16x8_r2,
        i4_rslt_vert_16x8_r3;
    int16x8_t i4_rslt_vert_16x8_r4, i4_rslt_vert_16x8_r5, i4_rslt_vert_16x8_r6,
        i4_rslt_vert_16x8_r7;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_src_stride = DYADIC_REF_W_C;

    /* Vertical interpolation */
    i4_samp_vert_8x8_r0 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r1 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r2 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r3 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r4 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r5 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;

    i4_rslt_vert_16x8_r0 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r0)), i4_coeff_0);
    i4_rslt_vert_16x8_r0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_1);
    vst1q_s16(pi2_tmp, i4_rslt_vert_16x8_r0);

    i4_rslt_vert_16x8_r1 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_2);
    i4_rslt_vert_16x8_r1 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r1, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 6, i4_rslt_vert_16x8_r1);

    i4_rslt_vert_16x8_r2 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_0);
    i4_rslt_vert_16x8_r2 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 12, i4_rslt_vert_16x8_r2);

    i4_rslt_vert_16x8_r3 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_2);
    i4_rslt_vert_16x8_r3 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r3, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 18, i4_rslt_vert_16x8_r3);

    i4_rslt_vert_16x8_r4 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_0);
    i4_rslt_vert_16x8_r4 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r4, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 24, i4_rslt_vert_16x8_r4);

    i4_rslt_vert_16x8_r5 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_2);
    i4_rslt_vert_16x8_r5 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r5, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 30, i4_rslt_vert_16x8_r5);

    i4_rslt_vert_16x8_r6 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_0);
    i4_rslt_vert_16x8_r6 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r6, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 36, i4_rslt_vert_16x8_r6);

    i4_rslt_vert_16x8_r7 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_2);
    i4_rslt_vert_16x8_r7 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r7, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r5)), i4_coeff_3);
    vst1_s16(pi2_tmp + 42, vget_low_s16(i4_rslt_vert_16x8_r7));
    vst1q_lane_s16(pi2_tmp + 46, i4_rslt_vert_16x8_r7, 4);
    vst1q_lane_s16(pi2_tmp + 47, i4_rslt_vert_16x8_r7, 5);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_2_neonintr              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_2_neonintr(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                                  WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    uint8x8_t i4_samp_vert_8x8_r0, i4_samp_vert_8x8_r1, i4_samp_vert_8x8_r2, i4_samp_vert_8x8_r3;
    uint8x8_t i4_samp_vert_8x8_r4;
    int16x8_t i4_rslt_vert_16x8_r0, i4_rslt_vert_16x8_r1, i4_rslt_vert_16x8_r2,
        i4_rslt_vert_16x8_r3;
    int16x8_t i4_rslt_vert_16x8_r4, i4_rslt_vert_16x8_r5, i4_rslt_vert_16x8_r6,
        i4_rslt_vert_16x8_r7;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pi2_tmp = pi2_tmp_filt_buf;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf + i4_src_stride;

    /* Vertical interpolation */
    i4_samp_vert_8x8_r0 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r1 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r2 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r3 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r4 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;

    /* since y_phase = phase_0 for y = 0 */
    i4_rslt_vert_16x8_r0 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r0)), i4_coeff_0);
    i4_rslt_vert_16x8_r0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_1);
    vst1q_s16(pi2_tmp, i4_rslt_vert_16x8_r0);

    i4_rslt_vert_16x8_r1 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r0)), i4_coeff_2);
    i4_rslt_vert_16x8_r1 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r1, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 6, i4_rslt_vert_16x8_r1);

    i4_rslt_vert_16x8_r2 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_0);
    i4_rslt_vert_16x8_r2 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 12, i4_rslt_vert_16x8_r2);

    i4_rslt_vert_16x8_r3 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_2);
    i4_rslt_vert_16x8_r3 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r3, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 18, i4_rslt_vert_16x8_r3);

    i4_rslt_vert_16x8_r4 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_0);
    i4_rslt_vert_16x8_r4 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r4, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 24, i4_rslt_vert_16x8_r4);

    i4_rslt_vert_16x8_r5 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_2);
    i4_rslt_vert_16x8_r5 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r5, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 30, i4_rslt_vert_16x8_r5);

    i4_rslt_vert_16x8_r6 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_0);
    i4_rslt_vert_16x8_r6 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r6, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 36, i4_rslt_vert_16x8_r6);

    i4_rslt_vert_16x8_r7 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_2);
    i4_rslt_vert_16x8_r7 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r7, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_3);
    vst1_s16(pi2_tmp + 42, vget_low_s16(i4_rslt_vert_16x8_r7));

    vst1q_lane_s16(pi2_tmp + 46, i4_rslt_vert_16x8_r7, 4);
    vst1q_lane_s16(pi2_tmp + 47, i4_rslt_vert_16x8_r7, 5);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_3_neonintr              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         21 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_3_neonintr(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                                  WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    uint8x8_t i4_samp_vert_8x8_r0, i4_samp_vert_8x8_r1, i4_samp_vert_8x8_r2;
    uint8x8_t i4_samp_vert_8x8_r3, i4_samp_vert_8x8_r4;
    int16x8_t i4_rslt_vert_16x8_r0, i4_rslt_vert_16x8_r1, i4_rslt_vert_16x8_r2,
        i4_rslt_vert_16x8_r3;
    int16x8_t i4_rslt_vert_16x8_r4, i4_rslt_vert_16x8_r5, i4_rslt_vert_16x8_r6,
        i4_rslt_vert_16x8_r7;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pi2_tmp = pi2_tmp_filt_buf;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf;

    /* Vertical interpolation */
    /* y = 0, y_phase = phase_0 */
    i4_samp_vert_8x8_r0 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r1 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r2 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r3 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;
    i4_samp_vert_8x8_r4 = vld1_u8(pu1_inp);
    pu1_inp += i4_src_stride;

    /* since y_phase = phase_0 for y = 0 */
    i4_rslt_vert_16x8_r0 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r0)), i4_coeff_0);
    i4_rslt_vert_16x8_r0 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_1);
    vst1q_s16(pi2_tmp, i4_rslt_vert_16x8_r0);

    i4_rslt_vert_16x8_r1 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r0)), i4_coeff_2);
    i4_rslt_vert_16x8_r1 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r1, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 6, i4_rslt_vert_16x8_r1);

    i4_rslt_vert_16x8_r2 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_0);
    i4_rslt_vert_16x8_r2 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r2, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 12, i4_rslt_vert_16x8_r2);

    i4_rslt_vert_16x8_r3 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r1)), i4_coeff_2);
    i4_rslt_vert_16x8_r3 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r3, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 18, i4_rslt_vert_16x8_r3);

    i4_rslt_vert_16x8_r4 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_0);
    i4_rslt_vert_16x8_r4 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r4, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 24, i4_rslt_vert_16x8_r4);

    i4_rslt_vert_16x8_r5 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r2)), i4_coeff_2);
    i4_rslt_vert_16x8_r5 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r5, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_3);
    vst1q_s16(pi2_tmp + 30, i4_rslt_vert_16x8_r5);

    i4_rslt_vert_16x8_r6 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_0);
    i4_rslt_vert_16x8_r6 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r6, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_1);
    vst1q_s16(pi2_tmp + 36, i4_rslt_vert_16x8_r6);

    i4_rslt_vert_16x8_r7 =
        vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r3)), i4_coeff_2);
    i4_rslt_vert_16x8_r7 = vmlaq_n_s16(
        i4_rslt_vert_16x8_r7, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_r4)), i4_coeff_3);
    vst1_s16(pi2_tmp + 42, vget_low_s16(i4_rslt_vert_16x8_r7));

    vst1q_lane_s16(pi2_tmp + 46, i4_rslt_vert_16x8_r7, 4);
    vst1q_lane_s16(pi2_tmp + 47, i4_rslt_vert_16x8_r7, 5);
}
