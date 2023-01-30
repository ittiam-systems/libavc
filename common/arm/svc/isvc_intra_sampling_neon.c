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
 *  isvc_intra_sampling_neon.c
 *
 * @brief
 *  neon variants of intra sampling functions used by IBL mode
 *
 * *******************************************************************************
 */

#include <arm_neon.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "isvc_intra_resample.h"

void isvc_interpolate_base_luma_dyadic_neon(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                            UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_y;
    WORD16 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp = pu1_inp_buf;
    UWORD8 *pu1_out = pu1_out_buf;
    WORD16 *pi2_tmp = pi2_tmp_filt_buf;

    int16x4_t i4_rslt_vert_16x4_1, i4_rslt_vert_16x4_2;
    uint8x8_t i4_samp_vert_8x8_0, i4_samp_vert_8x8_1, i4_samp_vert_8x8_2, i4_samp_vert_8x8_3;
    int16x8_t i4_rslt_vert_16x8_0, i4_rslt_vert_16x8_2;

    /* Horizontal interpolation */
    int32x4_t i4_rslt_horz_r0_1, i4_rslt_horz_r1_1, i4_rslt_horz_r0_2, i4_rslt_horz_r1_2;
    uint16x4_t i4_rslt_horz_r0_1_tmp, i4_rslt_horz_r1_1_tmp, i4_rslt_horz_r0_2_tmp,
        i4_rslt_horz_r1_2_tmp;
    uint16x8_t rslt_16x8_t_1, rslt_16x8_t_2;

    int16x4_t i4_samp_horz_16x4_0, i4_samp_horz_16x4_1, i4_samp_horz_16x4_2, i4_samp_horz_16x4_3,
        i4_samp_horz_16x4_4;
    int16x4_t i4_samp_horz_16x4_5, i4_samp_horz_16x4_6, i4_samp_horz_16x4_7, i4_samp_horz_16x4_8;
    int16_t i4_coeff_c0 = -3;
    int16_t i4_coeff_c1 = 28;
    int16_t i4_coeff_c2 = 8;
    int16_t i4_coeff_c3 = -1;
    int32x4x2_t i4_rslt_horz_r0_tmp32, i4_rslt_horz_r1_tmp32;
    int32x4_t const_512_32x4 = vdupq_n_s32(512);

    /* Filter coefficient values for phase 4 */
    i4_coeff_0 = -3;
    i4_coeff_1 = 28;
    i4_coeff_2 = 8;
    i4_coeff_3 = -1;

    i4_filt_stride = 12;
    i4_src_stride = DYADIC_REF_W_Y;

    /* Vertical interpolation */
    {
        /* First 64 bits*/
        i4_samp_vert_8x8_0 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_1 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_2 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);
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
            i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);

            i4_rslt_vert_16x8_0 =
                vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_0);
            i4_rslt_vert_16x8_0 =
                vmlaq_n_s16(i4_rslt_vert_16x8_0,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_1);
            i4_rslt_vert_16x8_0 =
                vmlaq_n_s16(i4_rslt_vert_16x8_0,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_2);
            i4_rslt_vert_16x8_0 =
                vmlaq_n_s16(i4_rslt_vert_16x8_0,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_3);

            i4_rslt_vert_16x8_2 =
                vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_3);
            i4_rslt_vert_16x8_2 =
                vmlaq_n_s16(i4_rslt_vert_16x8_2,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_2);
            i4_rslt_vert_16x8_2 =
                vmlaq_n_s16(i4_rslt_vert_16x8_2,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_1);
            i4_rslt_vert_16x8_2 =
                vmlaq_n_s16(i4_rslt_vert_16x8_2,
                            vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_0);

            vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_0));
            pi2_tmp += i4_filt_stride;
            vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_2));
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;
        }

        /* y = 15, y_phase = 4 */
        i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
        i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
        i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
        i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);

        i4_rslt_vert_16x8_0 =
            vmulq_n_s16(vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_0)), i4_coeff_0);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_1)), i4_coeff_1);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_2)), i4_coeff_2);
        i4_rslt_vert_16x8_0 = vmlaq_n_s16(
            i4_rslt_vert_16x8_0, vreinterpretq_s16_u16(vmovl_u8(i4_samp_vert_8x8_3)), i4_coeff_3);

        vst1q_s16(pi2_tmp, (i4_rslt_vert_16x8_0));
    }

    {
        /* Remaining 32 bits */
        pu1_inp = pu1_inp_buf + 8;
        pi2_tmp = pi2_tmp_filt_buf + 8;

        i4_samp_vert_8x8_0 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_1 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_2 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;
        i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);
        pu1_inp += i4_src_stride;

        i4_rslt_vert_16x4_1 = vmul_n_s16(
            vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_3);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))),
            i4_coeff_2);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))),
            i4_coeff_1);
        i4_rslt_vert_16x4_1 = vmla_n_s16(
            i4_rslt_vert_16x4_1, vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))),
            i4_coeff_0);

        vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 15; i4_y += 2)
        {
            i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
            i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
            i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
            i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);

            i4_rslt_vert_16x4_1 = vmul_n_s16(
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_0);
            i4_rslt_vert_16x4_1 = vmla_n_s16(
                i4_rslt_vert_16x4_1,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))), i4_coeff_1);
            i4_rslt_vert_16x4_1 = vmla_n_s16(
                i4_rslt_vert_16x4_1,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))), i4_coeff_2);
            i4_rslt_vert_16x4_1 = vmla_n_s16(
                i4_rslt_vert_16x4_1,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))), i4_coeff_3);

            i4_rslt_vert_16x4_2 = vmul_n_s16(
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_0))), i4_coeff_3);
            i4_rslt_vert_16x4_2 = vmla_n_s16(
                i4_rslt_vert_16x4_2,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_1))), i4_coeff_2);
            i4_rslt_vert_16x4_2 = vmla_n_s16(
                i4_rslt_vert_16x4_2,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_2))), i4_coeff_1);
            i4_rslt_vert_16x4_2 = vmla_n_s16(
                i4_rslt_vert_16x4_2,
                vreinterpret_s16_u16(vget_low_u16(vmovl_u8(i4_samp_vert_8x8_3))), i4_coeff_0);

            vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
            pi2_tmp += i4_filt_stride;
            vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_2));
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;
        }

        i4_samp_vert_8x8_0 = i4_samp_vert_8x8_1;
        i4_samp_vert_8x8_1 = i4_samp_vert_8x8_2;
        i4_samp_vert_8x8_2 = i4_samp_vert_8x8_3;
        i4_samp_vert_8x8_3 = vld1_u8((const UWORD8 *) pu1_inp);

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

        vst1_s16(pi2_tmp, (i4_rslt_vert_16x4_1));
        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;
    }

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

        i4_rslt_horz_r0_1 =
            vmull_n_s16(i4_samp_horz_16x4_0, i4_coeff_c3); /* a0c3 a1c3  a2c3  a3c3 */
        i4_rslt_horz_r0_1 =
            vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_1,
                        i4_coeff_c2); /* a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1 */
        i4_rslt_horz_r0_1 = vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_2, i4_coeff_c1);
        i4_rslt_horz_r0_1 = vmlal_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x4_3, i4_coeff_c0);
        /* i4_rslt_horz_r0_1 : contains res at even pos:0,2,4,6 */

        i4_rslt_horz_r1_1 =
            vmull_n_s16(i4_samp_horz_16x4_1, i4_coeff_c0); /* a0c0 a1c0  a2c0  a3c0 */
        i4_rslt_horz_r1_1 =
            vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_2,
                        i4_coeff_c1); /* a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1 */
        i4_rslt_horz_r1_1 = vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_3, i4_coeff_c2);
        i4_rslt_horz_r1_1 = vmlal_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x4_4, i4_coeff_c3);
        /* i4_rslt_horz_r1_1 : contains res at odd pos:1,3,5,7 */

        i4_rslt_horz_r0_2 =
            vmull_n_s16(i4_samp_horz_16x4_4, i4_coeff_c3); /* a0c3 a1c3  a2c3  a3c3 */
        i4_rslt_horz_r0_2 =
            vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_5,
                        i4_coeff_c2); /* a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1 */
        i4_rslt_horz_r0_2 = vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_6, i4_coeff_c1);
        i4_rslt_horz_r0_2 = vmlal_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x4_7, i4_coeff_c0);
        /* i4_rslt_horz_r0_1 : contains res at even pos:8,10,12,14 */

        i4_rslt_horz_r1_2 =
            vmull_n_s16(i4_samp_horz_16x4_5, i4_coeff_c0); /* a0c0 a1c0  a2c0  a3c0 */
        i4_rslt_horz_r1_2 =
            vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_6,
                        i4_coeff_c1); /* a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1 */
        i4_rslt_horz_r1_2 = vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_7, i4_coeff_c2);
        i4_rslt_horz_r1_2 = vmlal_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x4_8, i4_coeff_c3);
        /* i4_rslt_horz_r1_1 : contains res at odd pos:1,3,5,7 */

        i4_rslt_horz_r0_tmp32 = vzipq_s32(i4_rslt_horz_r0_1, i4_rslt_horz_r1_1);
        i4_rslt_horz_r1_tmp32 = vzipq_s32(i4_rslt_horz_r0_2, i4_rslt_horz_r1_2);

        i4_rslt_horz_r0_1 = vaddq_s32(i4_rslt_horz_r0_tmp32.val[0], const_512_32x4);
        i4_rslt_horz_r1_1 = vaddq_s32(i4_rslt_horz_r0_tmp32.val[1], const_512_32x4);
        i4_rslt_horz_r0_2 = vaddq_s32(i4_rslt_horz_r1_tmp32.val[0], const_512_32x4);
        i4_rslt_horz_r1_2 = vaddq_s32(i4_rslt_horz_r1_tmp32.val[1], const_512_32x4);

        i4_rslt_horz_r0_1_tmp = vqshrun_n_s32(i4_rslt_horz_r0_1, 10);
        i4_rslt_horz_r1_1_tmp = vqshrun_n_s32(i4_rslt_horz_r1_1, 10);

        i4_rslt_horz_r0_2_tmp = vqshrun_n_s32(i4_rslt_horz_r0_2, 10);
        i4_rslt_horz_r1_2_tmp = vqshrun_n_s32(i4_rslt_horz_r1_2, 10);

        rslt_16x8_t_1 = vcombine_u16(i4_rslt_horz_r0_1_tmp, i4_rslt_horz_r1_1_tmp);
        rslt_16x8_t_2 = vcombine_u16(i4_rslt_horz_r0_2_tmp, i4_rslt_horz_r1_2_tmp);

        vst1_u8(pu1_out, vqmovn_u16(rslt_16x8_t_1));
        vst1_u8(pu1_out + 8, vqmovn_u16(rslt_16x8_t_2));

        pu1_out += i4_out_stride;
        pi2_tmp += i4_filt_stride;
    }
}

void isvc_horz_interpol_chroma_dyadic_neon(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                           WORD32 i4_out_stride, WORD32 i4_phase_0,
                                           WORD32 i4_phase_1)
{
    WORD32 i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    UWORD8 *pu1_out = pu1_out_buf;
    WORD16 *pi2_tmp = pi2_tmp_filt_buf;
    WORD32 i4_filt_stride = 6;
    WORD32 i4_dst_stride = i4_out_stride;

    int16x8_t i4_samp_horz_16x8_r0_0, i4_samp_horz_16x8_r0_1, i4_samp_horz_16x8_r0_2;
    int16x8_t i4_samp_horz_16x8_r1_0, i4_samp_horz_16x8_r1_1, i4_samp_horz_16x8_r1_2;
    int16x8_t i4_rslt_horz_r0_1, i4_rslt_horz_r0_2;
    int16x8_t i4_rslt_horz_r1_1, i4_rslt_horz_r1_2;

    int16x8x2_t temp_horz_16x8_r0;
    int16x8x2_t temp_horz_16x8_r1;
    int16x8_t final_horz_16x8_r0_1;
    int16x8_t final_horz_16x8_r1_1;

    uint8x16_t i4_out_horz_8x16_r0, i4_out_horz_8x16_r1;
    uint8x16_t chroma_mask_8x16 = vreinterpretq_u8_u16(vdupq_n_u16(0x00ff));

    i4_coeff_0 = 16 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 16 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y += 2)
    {
        i4_samp_horz_16x8_r0_0 = vld1q_s16(pi2_tmp);     /* a0 a1 a2 a3 a4 a5 a6 a7 */
        i4_samp_horz_16x8_r0_1 = vld1q_s16(pi2_tmp + 1); /* a1 a2 a3 a4 */
        i4_samp_horz_16x8_r0_2 = vld1q_s16(pi2_tmp + 2); /* a2 a3 a4 a5 */

        i4_samp_horz_16x8_r1_0 = vld1q_s16(pi2_tmp + i4_filt_stride);
        i4_samp_horz_16x8_r1_1 = vld1q_s16(pi2_tmp + i4_filt_stride + 1);
        i4_samp_horz_16x8_r1_2 = vld1q_s16(pi2_tmp + (i4_filt_stride + 2));

        i4_rslt_horz_r0_1 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_0, i4_coeff_0); /* a0c0 a1c0  a2c0  a3c0 */
        i4_rslt_horz_r0_2 =
            vmulq_n_s16(i4_samp_horz_16x8_r0_1, i4_coeff_2); /* a1c2 a2c2  a3c2 a4c2 */

        i4_rslt_horz_r0_1 = vmlaq_n_s16(i4_rslt_horz_r0_1, i4_samp_horz_16x8_r0_1,
                                        i4_coeff_1); /* a0c0+a1c1 a1c0+a2c1  a2c0+a3c1  a3c0+a4c1 */
        i4_rslt_horz_r0_2 = vmlaq_n_s16(i4_rslt_horz_r0_2, i4_samp_horz_16x8_r0_2,
                                        i4_coeff_3); /* a1c2+a2c3  a2c2+a3c3 a3c2+a4c3 a4c2+a5c3 */

        i4_rslt_horz_r1_1 = vmulq_n_s16(i4_samp_horz_16x8_r1_0, i4_coeff_0);
        i4_rslt_horz_r1_2 = vmulq_n_s16(i4_samp_horz_16x8_r1_1, i4_coeff_2);

        i4_rslt_horz_r1_1 = vmlaq_n_s16(i4_rslt_horz_r1_1, i4_samp_horz_16x8_r1_1, i4_coeff_1);
        i4_rslt_horz_r1_2 = vmlaq_n_s16(i4_rslt_horz_r1_2, i4_samp_horz_16x8_r1_2, i4_coeff_3);

        temp_horz_16x8_r0 = vzipq_s16(i4_rslt_horz_r0_1, i4_rslt_horz_r0_2);
        temp_horz_16x8_r1 = vzipq_s16(i4_rslt_horz_r1_1, i4_rslt_horz_r1_2);

        final_horz_16x8_r0_1 = temp_horz_16x8_r0.val[0];
        final_horz_16x8_r1_1 = temp_horz_16x8_r1.val[0];

        final_horz_16x8_r0_1 = vrshrq_n_s16(final_horz_16x8_r0_1, 8);
        final_horz_16x8_r1_1 = vrshrq_n_s16(final_horz_16x8_r1_1, 8);

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
    }
}

void isvc_vert_interpol_chroma_dyadic_neon(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                           WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_src_stride = DYADIC_REF_W_C;
    UWORD8 *pu1_inp = pu1_inp_buf;
    WORD16 *pi2_tmp = pi2_tmp_filt_buf;

    uint8x8_t i4_samp_vert_8x8_r0, i4_samp_vert_8x8_r1, i4_samp_vert_8x8_r2, i4_samp_vert_8x8_r3,
        i4_samp_vert_8x8_r4, i4_samp_vert_8x8_r5;

    int16x8_t i4_rslt_vert_16x8_r0, i4_rslt_vert_16x8_r1, i4_rslt_vert_16x8_r2,
        i4_rslt_vert_16x8_r3, i4_rslt_vert_16x8_r4, i4_rslt_vert_16x8_r5, i4_rslt_vert_16x8_r6,
        i4_rslt_vert_16x8_r7;

    i4_coeff_0 = 16 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 16 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

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
