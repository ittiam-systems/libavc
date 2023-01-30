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
*
* @file
*  isvce_svc_residual_pred_neon.c
*
* @brief
*  Contains functions
* used for SVC residual
* prediction
*
*******************************************************************************
*/
#include <arm_neon.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_size_defs.h"
#include "isvc_macros.h"
#include "isvc_structs.h"

void isvce_luma_residual_sampler_2x_neon(coordinates_t *ps_ref_array_positions,
                                         coordinates_t *ps_ref_array_phases,
                                         buffer_container_t *ps_inp, buffer_container_t *ps_out,
                                         buffer_container_t *ps_scratch, UWORD32 u4_ref_nnz,
                                         UWORD8 u1_ref_tx_size)
{
    WORD16 *pi2_inp_data = (WORD16 *) ps_inp->pv_data;
    WORD16 *pi2_out_res = (WORD16 *) ps_out->pv_data;
    WORD32 i4_inp_data_stride = ps_inp->i4_data_stride;
    WORD32 i4_out_res_stride = ps_out->i4_data_stride;
    WORD16 *pi2_refarray_buffer = (WORD16 *) ps_scratch->pv_data;
    WORD32 i4_blk_ctr;

    UNUSED(ps_ref_array_positions);
    UNUSED(ps_ref_array_phases);

    /* For 2x scaling, offsets always point to TL pixel outside MB */
    /* Hence, refTransBlkIdc will be different and since phase */
    /* for first refArray pos for horiz filtering samples > 8, */
    /* first row and first column from the refArray is never used */
    pi2_inp_data += 1 + i4_inp_data_stride;

    if((u1_ref_tx_size) && (0 != u4_ref_nnz))
    {
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        WORD32 i4_i, i4_j;

        /* ----------- Horizontal Interpolation ---------------- */
        int16x8_t i2_coeff_add_16x8_r0;
        int16x8_t i2_coeff_16x8_r0_0, i2_coeff_16x8_r0_1;
        int16x8_t i2_coeff_16x8_sl_r0_0, i2_coeff_16x8_sl_r0_1;
        int16x8_t result_16x8_r0_0, result_16x8_r0_1;

        int16x8_t i2_coeff_add_16x8_r1;
        int16x8_t i2_coeff_16x8_r1_0, i2_coeff_16x8_r1_1;
        int16x8_t i2_coeff_16x8_sl_r1_0, i2_coeff_16x8_sl_r1_1;
        int16x8_t result_16x8_r1_0, result_16x8_r1_1;
        int16x8x2_t final_result_16x8x2_r0, final_result_16x8x2_r1;

        pi2_ref_data_byte = pi2_inp_data;

        /* ----------- Horizontal Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

        for(i4_i = 0; i4_i < BLK8x8SIZE; i4_i += 2)
        {
            i2_coeff_16x8_r0_0 = vld1q_s16(pi2_ref_data_byte);
            i2_coeff_16x8_r0_1 = vld1q_s16((pi2_ref_data_byte + 1));

            i2_coeff_16x8_r1_0 = vld1q_s16(pi2_ref_data_byte + i4_inp_data_stride);
            i2_coeff_16x8_r1_1 = vld1q_s16((pi2_ref_data_byte + i4_inp_data_stride + 1));

            i2_coeff_add_16x8_r0 = vaddq_s16(i2_coeff_16x8_r0_0, i2_coeff_16x8_r0_1);
            i2_coeff_16x8_sl_r0_0 = vshlq_n_s16(i2_coeff_16x8_r0_0, 1);
            i2_coeff_16x8_sl_r0_1 = vshlq_n_s16(i2_coeff_16x8_r0_1, 1);

            i2_coeff_add_16x8_r1 = vaddq_s16(i2_coeff_16x8_r1_0, i2_coeff_16x8_r1_1);
            i2_coeff_16x8_sl_r1_0 = vshlq_n_s16(i2_coeff_16x8_r1_0, 1);
            i2_coeff_16x8_sl_r1_1 = vshlq_n_s16(i2_coeff_16x8_r1_1, 1);

            result_16x8_r0_0 = vaddq_s16(i2_coeff_16x8_sl_r0_0, i2_coeff_add_16x8_r0);
            result_16x8_r0_1 = vaddq_s16(i2_coeff_16x8_sl_r0_1, i2_coeff_add_16x8_r0);

            result_16x8_r1_0 = vaddq_s16(i2_coeff_16x8_sl_r1_0, i2_coeff_add_16x8_r1);
            result_16x8_r1_1 = vaddq_s16(i2_coeff_16x8_sl_r1_1, i2_coeff_add_16x8_r1);

            final_result_16x8x2_r0 = vzipq_s16(result_16x8_r0_0, result_16x8_r0_1);
            final_result_16x8x2_r1 = vzipq_s16(result_16x8_r1_0, result_16x8_r1_1);

            vst1q_s32(pi4_ref_array + 1, vmovl_s16(vget_low_s16(final_result_16x8x2_r0.val[0])));
            vst1q_s32(pi4_ref_array + 5, vmovl_s16(vget_high_s16(final_result_16x8x2_r0.val[0])));
            vst1q_s32(pi4_ref_array + 9, vmovl_s16(vget_low_s16(final_result_16x8x2_r0.val[1])));
            vst1q_s32(pi4_ref_array + 13, vmovl_s16(vget_high_s16(final_result_16x8x2_r0.val[1])));

            pi4_ref_array[0] = pi2_ref_data_byte[0] << 2;
            pi4_ref_array[15] = pi2_ref_data_byte[7] << 2;
            pi4_ref_array += 16;
            pi2_ref_data_byte += i4_inp_data_stride;

            vst1q_s32(pi4_ref_array + 1, vmovl_s16(vget_low_s16(final_result_16x8x2_r1.val[0])));
            vst1q_s32(pi4_ref_array + 5, vmovl_s16(vget_high_s16(final_result_16x8x2_r1.val[0])));
            vst1q_s32(pi4_ref_array + 9, vmovl_s16(vget_low_s16(final_result_16x8x2_r1.val[1])));
            vst1q_s32(pi4_ref_array + 13, vmovl_s16(vget_high_s16(final_result_16x8x2_r1.val[1])));

            pi4_ref_array[0] = pi2_ref_data_byte[0] << 2;
            pi4_ref_array[15] = pi2_ref_data_byte[7] << 2;
            pi4_ref_array += 16;
            /* vertical loop updates */
            pi2_ref_data_byte = pi2_inp_data + ((i4_i + 2) * i4_inp_data_stride);
        }

        /* ----------- Vertical Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;
        {
            WORD32 *pi4_ref_array_temp;
            WORD16 *pi2_out;
            int32x4_t i4_horz_samp_32x4_r1_1, i4_horz_samp_32x4_r1_2, i4_horz_samp_32x4_r1_3,
                i4_horz_samp_32x4_r1_4;
            int32x4_t i4_horz_samp_32x4_r2_1, i4_horz_samp_32x4_r2_2, i4_horz_samp_32x4_r2_3,
                i4_horz_samp_32x4_r2_4;

            int32x4_t i4_horz_res_32x4_r1_1, i4_horz_res_32x4_r1_2, i4_horz_res_32x4_r1_3,
                i4_horz_res_32x4_r1_4;
            int32x4_t i4_horz_res_32x4_r2_1, i4_horz_res_32x4_r2_2, i4_horz_res_32x4_r2_3,
                i4_horz_res_32x4_r2_4;
            int32x4_t i4_horz_res_32x4_r3_1, i4_horz_res_32x4_r3_2, i4_horz_res_32x4_r3_3,
                i4_horz_res_32x4_r3_4;
            int32x4_t horz_add_32x4_r2_1, horz_add_32x4_r2_2, horz_add_32x4_r2_3,
                horz_add_32x4_r2_4;

            int16x8_t comb_horz_16x8_1, comb_horz_16x8_2, comb_horz_16x8_3, comb_horz_16x8_4;
            pi4_ref_array_temp = pi4_ref_array;
            pi2_out = pi2_out_res;

            i4_horz_samp_32x4_r1_1 = vld1q_s32(pi4_ref_array_temp);
            i4_horz_samp_32x4_r1_2 = vld1q_s32(pi4_ref_array_temp + 4);
            i4_horz_samp_32x4_r1_3 = vld1q_s32(pi4_ref_array_temp + 8);
            i4_horz_samp_32x4_r1_4 = vld1q_s32(pi4_ref_array_temp + 12);

            /* populate the first inter sample */
            i4_horz_res_32x4_r1_1 = vrshrq_n_s32(i4_horz_samp_32x4_r1_1, 2);
            i4_horz_res_32x4_r1_2 = vrshrq_n_s32(i4_horz_samp_32x4_r1_2, 2);
            i4_horz_res_32x4_r1_3 = vrshrq_n_s32(i4_horz_samp_32x4_r1_3, 2);
            i4_horz_res_32x4_r1_4 = vrshrq_n_s32(i4_horz_samp_32x4_r1_4, 2);

            comb_horz_16x8_1 =
                vcombine_s16(vmovn_s32(i4_horz_res_32x4_r1_1), vmovn_s32(i4_horz_res_32x4_r1_2));
            comb_horz_16x8_2 =
                vcombine_s16(vmovn_s32(i4_horz_res_32x4_r1_3), vmovn_s32(i4_horz_res_32x4_r1_4));
            vst1q_s16(pi2_out, comb_horz_16x8_1);
            vst1q_s16(pi2_out + 8, comb_horz_16x8_2);

            pi2_out += i4_out_res_stride;

            for(i4_j = 0; i4_j < 14; i4_j += 2)
            {
                pi4_ref_array_temp += MB_SIZE;
                i4_horz_samp_32x4_r2_1 = vld1q_s32(pi4_ref_array_temp);
                i4_horz_samp_32x4_r2_2 = vld1q_s32(pi4_ref_array_temp + 4);
                i4_horz_samp_32x4_r2_3 = vld1q_s32(pi4_ref_array_temp + 8);
                i4_horz_samp_32x4_r2_4 = vld1q_s32(pi4_ref_array_temp + 12);

                horz_add_32x4_r2_1 = vaddq_s32(i4_horz_samp_32x4_r1_1, i4_horz_samp_32x4_r2_1);
                horz_add_32x4_r2_2 = vaddq_s32(i4_horz_samp_32x4_r1_2, i4_horz_samp_32x4_r2_2);
                horz_add_32x4_r2_3 = vaddq_s32(i4_horz_samp_32x4_r1_3, i4_horz_samp_32x4_r2_3);
                horz_add_32x4_r2_4 = vaddq_s32(i4_horz_samp_32x4_r1_4, i4_horz_samp_32x4_r2_4);

                i4_horz_res_32x4_r2_1 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r1_1, 1), horz_add_32x4_r2_1);
                i4_horz_res_32x4_r2_2 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r1_2, 1), horz_add_32x4_r2_2);
                i4_horz_res_32x4_r2_3 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r1_3, 1), horz_add_32x4_r2_3);
                i4_horz_res_32x4_r2_4 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r1_4, 1), horz_add_32x4_r2_4);

                i4_horz_res_32x4_r3_1 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r2_1, 1), horz_add_32x4_r2_1);
                i4_horz_res_32x4_r3_2 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r2_2, 1), horz_add_32x4_r2_2);
                i4_horz_res_32x4_r3_3 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r2_3, 1), horz_add_32x4_r2_3);
                i4_horz_res_32x4_r3_4 =
                    vaddq_s32(vshlq_n_s32(i4_horz_samp_32x4_r2_4, 1), horz_add_32x4_r2_4);

                i4_horz_res_32x4_r2_1 = vrshrq_n_s32(i4_horz_res_32x4_r2_1, 4);
                i4_horz_res_32x4_r2_2 = vrshrq_n_s32(i4_horz_res_32x4_r2_2, 4);
                i4_horz_res_32x4_r2_3 = vrshrq_n_s32(i4_horz_res_32x4_r2_3, 4);
                i4_horz_res_32x4_r2_4 = vrshrq_n_s32(i4_horz_res_32x4_r2_4, 4);

                i4_horz_res_32x4_r3_1 = vrshrq_n_s32(i4_horz_res_32x4_r3_1, 4);
                i4_horz_res_32x4_r3_2 = vrshrq_n_s32(i4_horz_res_32x4_r3_2, 4);
                i4_horz_res_32x4_r3_3 = vrshrq_n_s32(i4_horz_res_32x4_r3_3, 4);
                i4_horz_res_32x4_r3_4 = vrshrq_n_s32(i4_horz_res_32x4_r3_4, 4);

                comb_horz_16x8_1 = vcombine_s16(vmovn_s32(i4_horz_res_32x4_r2_1),
                                                vmovn_s32(i4_horz_res_32x4_r2_2));
                comb_horz_16x8_2 = vcombine_s16(vmovn_s32(i4_horz_res_32x4_r2_3),
                                                vmovn_s32(i4_horz_res_32x4_r2_4));

                comb_horz_16x8_3 = vcombine_s16(vmovn_s32(i4_horz_res_32x4_r3_1),
                                                vmovn_s32(i4_horz_res_32x4_r3_2));
                comb_horz_16x8_4 = vcombine_s16(vmovn_s32(i4_horz_res_32x4_r3_3),
                                                vmovn_s32(i4_horz_res_32x4_r3_4));

                /* populate 2 samples based on current coeffs */
                vst1q_s16(pi2_out, comb_horz_16x8_1);
                vst1q_s16(pi2_out + 8, comb_horz_16x8_2);
                pi2_out += i4_out_res_stride;

                vst1q_s16(pi2_out, comb_horz_16x8_3);
                vst1q_s16(pi2_out + 8, comb_horz_16x8_4);
                pi2_out += i4_out_res_stride;

                /* store the coeff 2 to coeff 1 */
                /* (used in next iteration)     */
                i4_horz_samp_32x4_r1_1 = i4_horz_samp_32x4_r2_1;
                i4_horz_samp_32x4_r1_2 = i4_horz_samp_32x4_r2_2;
                i4_horz_samp_32x4_r1_3 = i4_horz_samp_32x4_r2_3;
                i4_horz_samp_32x4_r1_4 = i4_horz_samp_32x4_r2_4;
            }

            /* populate the first inter sample */
            i4_horz_res_32x4_r1_1 = vrshrq_n_s32(i4_horz_samp_32x4_r1_1, 2);
            i4_horz_res_32x4_r1_2 = vrshrq_n_s32(i4_horz_samp_32x4_r1_2, 2);
            i4_horz_res_32x4_r1_3 = vrshrq_n_s32(i4_horz_samp_32x4_r1_3, 2);
            i4_horz_res_32x4_r1_4 = vrshrq_n_s32(i4_horz_samp_32x4_r1_4, 2);

            comb_horz_16x8_1 =
                vcombine_s16(vmovn_s32(i4_horz_res_32x4_r1_1), vmovn_s32(i4_horz_res_32x4_r1_2));
            comb_horz_16x8_2 =
                vcombine_s16(vmovn_s32(i4_horz_res_32x4_r1_3), vmovn_s32(i4_horz_res_32x4_r1_4));
            vst1q_s16(pi2_out, comb_horz_16x8_1);
            vst1q_s16(pi2_out + 8, comb_horz_16x8_2);

            /* horizontal loop updates */
            pi4_ref_array++;
            pi2_out_res++;
        }
    }
    else
    {
        /* ----------------------------------------------------------------- */
        /* LOOP over number of blocks                                        */
        /* ----------------------------------------------------------------- */
        for(i4_blk_ctr = 0; i4_blk_ctr < 4; i4_blk_ctr++)
        {
            /* if reference layer is not coded then no processing */
            if(0 != (u4_ref_nnz & 0x1))
            {
                int16x8_t i2_coeff1_16x8_r0_0, i2_coeff1_16x8_r0_1;
                int16x8_t i2_coeff1_16x8_r1_0, i2_coeff1_16x8_r1_1;
                int16x8_t i2_coeff1_16x8_r2_0, i2_coeff1_16x8_r2_1;
                int16x8_t i2_coeff1_16x8_r3_0, i2_coeff1_16x8_r3_1;
                int16x8_t i2_add_16x8_r0_0;
                int16x8_t i2_add_16x8_r1_0;
                int16x8_t i2_add_16x8_r2_0;
                int16x8_t i2_add_16x8_r3_0;
                int16x8_t i2_res_16x8_r0_0, i2_res_16x8_r0_1;
                int16x8_t i2_res_16x8_r1_0, i2_res_16x8_r1_1;
                int16x8_t i2_res_16x8_r2_0, i2_res_16x8_r2_1;
                int16x8_t i2_res_16x8_r3_0, i2_res_16x8_r3_1;
                int16x4_t i4_horz_samp_16x4_r0_1, i4_horz_samp_16x4_r0_2;
                int16x4_t i4_horz_samp_16x4_r1_1, i4_horz_samp_16x4_r1_2;
                int16x4_t i4_horz_samp_16x4_r2_1, i4_horz_samp_16x4_r2_2;
                int16x4_t i4_horz_samp_16x4_r3_1, i4_horz_samp_16x4_r3_2;
                int32x4_t i4_horz_samp_32x4_r0_1, i4_horz_samp_32x4_r0_2;
                int32x4_t i4_horz_samp_32x4_r1_1, i4_horz_samp_32x4_r1_2;
                int32x4_t i4_horz_samp_32x4_r2_1, i4_horz_samp_32x4_r2_2;
                int32x4_t i4_horz_samp_32x4_r3_1, i4_horz_samp_32x4_r3_2;
                int32x4_t i4_horz_add_32x4_r1_1, i4_horz_add_32x4_r1_2;
                int32x4_t i4_horz_add_32x4_r2_1, i4_horz_add_32x4_r2_2;
                int32x4_t i4_horz_add_32x4_r3_1, i4_horz_add_32x4_r3_2;
                int16x4_t i4_horz_res_16x4_r0_1, i4_horz_res_16x4_r0_2;
                int16x4_t i4_horz_res_16x4_r1_1, i4_horz_res_16x4_r1_2;
                int16x4_t i4_horz_res_16x4_r2_1, i4_horz_res_16x4_r2_2;
                int16x4_t i4_horz_res_16x4_r3_1, i4_horz_res_16x4_r3_2;
                int16x4_t i4_horz_res_16x4_r4_1, i4_horz_res_16x4_r4_2;
                int16x4_t i4_horz_res_16x4_r5_1, i4_horz_res_16x4_r5_2;
                int16x4_t i4_horz_res_16x4_r6_1, i4_horz_res_16x4_r6_2;
                int16x4_t i4_horz_res_16x4_r7_1, i4_horz_res_16x4_r7_2;
                int32x4_t i4_horz_res_32x4_r1_1, i4_horz_res_32x4_r1_2;
                int32x4_t i4_horz_res_32x4_r2_1, i4_horz_res_32x4_r2_2;
                int32x4_t i4_horz_res_32x4_r3_1, i4_horz_res_32x4_r3_2;
                int32x4_t i4_horz_res_32x4_r4_1, i4_horz_res_32x4_r4_2;
                int32x4_t i4_horz_res_32x4_r5_1, i4_horz_res_32x4_r5_2;
                int32x4_t i4_horz_res_32x4_r6_1, i4_horz_res_32x4_r6_2;
                int16x8x2_t ti2_res_16x8x2_r0, ti2_res_16x8x2_r1;
                int16x8x2_t ti2_res_16x8x2_r2, ti2_res_16x8x2_r3;

                i2_coeff1_16x8_r0_0 = vld1q_s16(pi2_inp_data);
                i2_coeff1_16x8_r1_0 = vld1q_s16(pi2_inp_data + i4_inp_data_stride);
                i2_coeff1_16x8_r2_0 = vld1q_s16(pi2_inp_data + (i4_inp_data_stride << 1));
                i2_coeff1_16x8_r3_0 =
                    vld1q_s16(pi2_inp_data + (i4_inp_data_stride << 1) + i4_inp_data_stride);

                i2_coeff1_16x8_r0_1 = vextq_s16(i2_coeff1_16x8_r0_0, i2_coeff1_16x8_r0_0, 1);
                i2_coeff1_16x8_r1_1 = vextq_s16(i2_coeff1_16x8_r1_0, i2_coeff1_16x8_r1_0, 1);
                i2_coeff1_16x8_r2_1 = vextq_s16(i2_coeff1_16x8_r2_0, i2_coeff1_16x8_r2_0, 1);
                i2_coeff1_16x8_r3_1 = vextq_s16(i2_coeff1_16x8_r3_0, i2_coeff1_16x8_r3_0, 1);

                i2_add_16x8_r0_0 = vaddq_s16(i2_coeff1_16x8_r0_1, i2_coeff1_16x8_r0_0);
                i2_add_16x8_r1_0 = vaddq_s16(i2_coeff1_16x8_r1_1, i2_coeff1_16x8_r1_0);
                i2_add_16x8_r2_0 = vaddq_s16(i2_coeff1_16x8_r2_1, i2_coeff1_16x8_r2_0);
                i2_add_16x8_r3_0 = vaddq_s16(i2_coeff1_16x8_r3_1, i2_coeff1_16x8_r3_0);

                i2_coeff1_16x8_r0_0 = vshlq_n_s16(i2_coeff1_16x8_r0_0, 1);
                i2_coeff1_16x8_r1_0 = vshlq_n_s16(i2_coeff1_16x8_r1_0, 1);
                i2_coeff1_16x8_r2_0 = vshlq_n_s16(i2_coeff1_16x8_r2_0, 1);
                i2_coeff1_16x8_r3_0 = vshlq_n_s16(i2_coeff1_16x8_r3_0, 1);

                i2_coeff1_16x8_r0_1 = vshlq_n_s16(i2_coeff1_16x8_r0_1, 1);
                i2_coeff1_16x8_r1_1 = vshlq_n_s16(i2_coeff1_16x8_r1_1, 1);
                i2_coeff1_16x8_r2_1 = vshlq_n_s16(i2_coeff1_16x8_r2_1, 1);
                i2_coeff1_16x8_r3_1 = vshlq_n_s16(i2_coeff1_16x8_r3_1, 1);

                i2_res_16x8_r0_0 = vaddq_s16(i2_coeff1_16x8_r0_0, i2_add_16x8_r0_0);
                i2_res_16x8_r1_0 = vaddq_s16(i2_coeff1_16x8_r1_0, i2_add_16x8_r1_0);
                i2_res_16x8_r2_0 = vaddq_s16(i2_coeff1_16x8_r2_0, i2_add_16x8_r2_0);
                i2_res_16x8_r3_0 = vaddq_s16(i2_coeff1_16x8_r3_0, i2_add_16x8_r3_0);

                i2_res_16x8_r0_1 = vaddq_s16(i2_coeff1_16x8_r0_1, i2_add_16x8_r0_0);
                i2_res_16x8_r1_1 = vaddq_s16(i2_coeff1_16x8_r1_1, i2_add_16x8_r1_0);
                i2_res_16x8_r2_1 = vaddq_s16(i2_coeff1_16x8_r2_1, i2_add_16x8_r2_0);
                i2_res_16x8_r3_1 = vaddq_s16(i2_coeff1_16x8_r3_1, i2_add_16x8_r3_0);

                ti2_res_16x8x2_r0 = vzipq_s16(i2_res_16x8_r0_0, i2_res_16x8_r0_1);
                ti2_res_16x8x2_r1 = vzipq_s16(i2_res_16x8_r1_0, i2_res_16x8_r1_1);
                ti2_res_16x8x2_r2 = vzipq_s16(i2_res_16x8_r2_0, i2_res_16x8_r2_1);
                ti2_res_16x8x2_r3 = vzipq_s16(i2_res_16x8_r3_0, i2_res_16x8_r3_1);

                i2_coeff1_16x8_r0_0 = vshlq_n_s16(i2_coeff1_16x8_r0_0, 1);
                i2_coeff1_16x8_r1_0 = vshlq_n_s16(i2_coeff1_16x8_r1_0, 1);
                i2_coeff1_16x8_r2_0 = vshlq_n_s16(i2_coeff1_16x8_r2_0, 1);
                i2_coeff1_16x8_r3_0 = vshlq_n_s16(i2_coeff1_16x8_r3_0, 1);

                vst1q_s16(pi2_refarray_buffer + 1, ti2_res_16x8x2_r0.val[0]);
                vst1q_lane_s16(pi2_refarray_buffer, i2_coeff1_16x8_r0_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 7, i2_coeff1_16x8_r0_0, 3);

                vst1q_s16(pi2_refarray_buffer + 9, ti2_res_16x8x2_r1.val[0]);
                vst1q_lane_s16(pi2_refarray_buffer + 8, i2_coeff1_16x8_r1_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 15, i2_coeff1_16x8_r1_0, 3);

                vst1q_s16(pi2_refarray_buffer + 17, ti2_res_16x8x2_r2.val[0]);
                vst1q_lane_s16(pi2_refarray_buffer + 16, i2_coeff1_16x8_r2_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 23, i2_coeff1_16x8_r2_0, 3);

                vst1q_s16(pi2_refarray_buffer + 25, ti2_res_16x8x2_r3.val[0]);
                vst1q_lane_s16(pi2_refarray_buffer + 24, i2_coeff1_16x8_r3_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 31, i2_coeff1_16x8_r3_0, 3);

                i4_horz_samp_16x4_r0_1 = vld1_s16(pi2_refarray_buffer);
                i4_horz_samp_16x4_r0_2 = vld1_s16(pi2_refarray_buffer + 4);

                i4_horz_samp_16x4_r1_1 = vld1_s16(pi2_refarray_buffer + 8);
                i4_horz_samp_16x4_r1_2 = vld1_s16(pi2_refarray_buffer + 12);

                i4_horz_samp_16x4_r2_1 = vld1_s16(pi2_refarray_buffer + 16);
                i4_horz_samp_16x4_r2_2 = vld1_s16(pi2_refarray_buffer + 20);

                i4_horz_samp_16x4_r3_1 = vld1_s16(pi2_refarray_buffer + 24);
                i4_horz_samp_16x4_r3_2 = vld1_s16(pi2_refarray_buffer + 28);

                i4_horz_res_16x4_r0_1 = vrshr_n_s16(i4_horz_samp_16x4_r0_1, 2);
                i4_horz_res_16x4_r0_2 = vrshr_n_s16(i4_horz_samp_16x4_r0_2, 2);

                i4_horz_add_32x4_r1_1 = vaddl_s16(i4_horz_samp_16x4_r0_1, i4_horz_samp_16x4_r1_1);
                i4_horz_add_32x4_r1_2 = vaddl_s16(i4_horz_samp_16x4_r0_2, i4_horz_samp_16x4_r1_2);

                i4_horz_add_32x4_r2_1 = vaddl_s16(i4_horz_samp_16x4_r1_1, i4_horz_samp_16x4_r2_1);
                i4_horz_add_32x4_r2_2 = vaddl_s16(i4_horz_samp_16x4_r1_2, i4_horz_samp_16x4_r2_2);

                i4_horz_add_32x4_r3_1 = vaddl_s16(i4_horz_samp_16x4_r2_1, i4_horz_samp_16x4_r3_1);
                i4_horz_add_32x4_r3_2 = vaddl_s16(i4_horz_samp_16x4_r2_2, i4_horz_samp_16x4_r3_2);

                i4_horz_samp_32x4_r0_1 = vshll_n_s16(i4_horz_samp_16x4_r0_1, 1);
                i4_horz_samp_32x4_r0_2 = vshll_n_s16(i4_horz_samp_16x4_r0_2, 1);

                i4_horz_samp_32x4_r1_1 = vshll_n_s16(i4_horz_samp_16x4_r1_1, 1);
                i4_horz_samp_32x4_r1_2 = vshll_n_s16(i4_horz_samp_16x4_r1_2, 1);

                i4_horz_samp_32x4_r2_1 = vshll_n_s16(i4_horz_samp_16x4_r2_1, 1);
                i4_horz_samp_32x4_r2_2 = vshll_n_s16(i4_horz_samp_16x4_r2_2, 1);

                i4_horz_samp_32x4_r3_1 = vshll_n_s16(i4_horz_samp_16x4_r3_1, 1);
                i4_horz_samp_32x4_r3_2 = vshll_n_s16(i4_horz_samp_16x4_r3_2, 1);

                i4_horz_res_32x4_r1_1 = vaddq_s32(i4_horz_samp_32x4_r0_1, i4_horz_add_32x4_r1_1);
                i4_horz_res_32x4_r1_2 = vaddq_s32(i4_horz_samp_32x4_r0_2, i4_horz_add_32x4_r1_2);

                i4_horz_res_32x4_r2_1 = vaddq_s32(i4_horz_samp_32x4_r1_1, i4_horz_add_32x4_r1_1);
                i4_horz_res_32x4_r2_2 = vaddq_s32(i4_horz_samp_32x4_r1_2, i4_horz_add_32x4_r1_2);

                i4_horz_res_32x4_r3_1 = vaddq_s32(i4_horz_samp_32x4_r1_1, i4_horz_add_32x4_r2_1);
                i4_horz_res_32x4_r3_2 = vaddq_s32(i4_horz_samp_32x4_r1_2, i4_horz_add_32x4_r2_2);

                i4_horz_res_32x4_r4_1 = vaddq_s32(i4_horz_samp_32x4_r2_1, i4_horz_add_32x4_r2_1);
                i4_horz_res_32x4_r4_2 = vaddq_s32(i4_horz_samp_32x4_r2_2, i4_horz_add_32x4_r2_2);

                i4_horz_res_32x4_r5_1 = vaddq_s32(i4_horz_samp_32x4_r2_1, i4_horz_add_32x4_r3_1);
                i4_horz_res_32x4_r5_2 = vaddq_s32(i4_horz_samp_32x4_r2_2, i4_horz_add_32x4_r3_2);

                i4_horz_res_32x4_r6_1 = vaddq_s32(i4_horz_samp_32x4_r3_1, i4_horz_add_32x4_r3_1);
                i4_horz_res_32x4_r6_2 = vaddq_s32(i4_horz_samp_32x4_r3_2, i4_horz_add_32x4_r3_2);

                i4_horz_res_16x4_r1_1 = vqrshrn_n_s32(i4_horz_res_32x4_r1_1, 4);
                i4_horz_res_16x4_r1_2 = vqrshrn_n_s32(i4_horz_res_32x4_r1_2, 4);

                i4_horz_res_16x4_r2_1 = vqrshrn_n_s32(i4_horz_res_32x4_r2_1, 4);
                i4_horz_res_16x4_r2_2 = vqrshrn_n_s32(i4_horz_res_32x4_r2_2, 4);

                i4_horz_res_16x4_r3_1 = vqrshrn_n_s32(i4_horz_res_32x4_r3_1, 4);
                i4_horz_res_16x4_r3_2 = vqrshrn_n_s32(i4_horz_res_32x4_r3_2, 4);

                i4_horz_res_16x4_r4_1 = vqrshrn_n_s32(i4_horz_res_32x4_r4_1, 4);
                i4_horz_res_16x4_r4_2 = vqrshrn_n_s32(i4_horz_res_32x4_r4_2, 4);

                i4_horz_res_16x4_r5_1 = vqrshrn_n_s32(i4_horz_res_32x4_r5_1, 4);
                i4_horz_res_16x4_r5_2 = vqrshrn_n_s32(i4_horz_res_32x4_r5_2, 4);

                i4_horz_res_16x4_r6_1 = vqrshrn_n_s32(i4_horz_res_32x4_r6_1, 4);
                i4_horz_res_16x4_r6_2 = vqrshrn_n_s32(i4_horz_res_32x4_r6_2, 4);

                i4_horz_res_16x4_r7_1 = vrshr_n_s16(i4_horz_samp_16x4_r3_1, 2);
                i4_horz_res_16x4_r7_2 = vrshr_n_s16(i4_horz_samp_16x4_r3_2, 2);

                vst1_s16(pi2_out_res, i4_horz_res_16x4_r0_1);
                vst1_s16(pi2_out_res + 4, i4_horz_res_16x4_r0_2);

                vst1_s16(pi2_out_res + i4_out_res_stride, i4_horz_res_16x4_r1_1);
                vst1_s16(pi2_out_res + i4_out_res_stride + 4, i4_horz_res_16x4_r1_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride << 1), i4_horz_res_16x4_r2_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride << 1) + 4, i4_horz_res_16x4_r2_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride * 3), i4_horz_res_16x4_r3_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride * 3) + 4, i4_horz_res_16x4_r3_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride << 2), i4_horz_res_16x4_r4_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride << 2) + 4, i4_horz_res_16x4_r4_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride * 5), i4_horz_res_16x4_r5_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride * 5) + 4, i4_horz_res_16x4_r5_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride * 6), i4_horz_res_16x4_r6_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride * 6) + 4, i4_horz_res_16x4_r6_2);

                vst1_s16(pi2_out_res + (i4_out_res_stride * 7), i4_horz_res_16x4_r7_1);
                vst1_s16(pi2_out_res + (i4_out_res_stride * 7) + 4, i4_horz_res_16x4_r7_2);

                pi2_out_res += BLK8x8SIZE;
            }
            else
            {
                pi2_out_res += BLK8x8SIZE;
            }

            /* Block level loop updates */
            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= SUB_BLK_WIDTH_4x4;
                pi2_inp_data += (i4_inp_data_stride * SUB_BLK_HEIGHT_4x4);
                pi2_out_res -= MB_SIZE;
                pi2_out_res += (i4_out_res_stride * BLK8x8SIZE);
                u4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += SUB_BLK_HEIGHT_4x4;
            }
            u4_ref_nnz >>= 1;
        }
        /* The above loop iterates over all the blocks */
    }
}

UWORD32 isvce_get_sad_with_residual_pred_neon(buffer_container_t *ps_src,
                                              buffer_container_t *ps_pred,
                                              buffer_container_t *ps_res, UWORD32 u4_mb_wd,
                                              UWORD32 u4_mb_ht)
{
    UWORD32 i, j, u4_sad = 0;
    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD8 *pu1_pred = (UWORD8 *) ps_pred->pv_data;
    WORD16 *pi2_res = (WORD16 *) ps_res->pv_data;
    WORD32 i4_src_stride = ps_src->i4_data_stride;
    WORD32 i4_pred_stride = ps_pred->i4_data_stride;
    WORD32 i4_res_stride = ps_res->i4_data_stride;
    UWORD32 u4_num_rows_per_loop = 8;
    UWORD32 u4_ht_by_8 = u4_mb_ht / u4_num_rows_per_loop;
    uint8x8_t src0, src1, src2, src3;
    uint8x8_t src4, src5, src6, src7;
    uint8x8_t pred0, pred1, pred2, pred3;
    uint8x8_t pred4, pred5, pred6, pred7;
    int16x8_t res0_16x8, res1_16x8, res2_16x8, res3_16x8, res4_16x8, res5_16x8, res6_16x8,
        res7_16x8;
    uint16x8_t res0_u16x8, res1_u16x8, res2_u16x8, res3_u16x8, res4_u16x8, res5_u16x8, res6_u16x8,
        res7_u16x8;
    int16x8_t respred0_16x8, respred1_16x8, respred2_16x8, respred3_16x8, respred4_16x8,
        respred5_16x8, respred6_16x8, respred7_16x8;
    int16x8_t temp0_16x8, temp1_16x8, temp2_16x8, temp3_16x8, temp4_16x8, temp5_16x8, temp6_16x8,
        temp7_16x8;
    int32x4_t temp0_32x4;
    int32x2_t temp0_32x2;

    if((u4_mb_wd == 16) && (u4_mb_ht % 8 == 0))
    {
        for(i = 0; i < u4_ht_by_8; i++)
        {
            /* This loop processes 4 rows of 16 bytes each iteration */
            /* So, 8 rows are processed across two iterations */
            for(j = 0; j < 2; j++)
            {
                src0 = vld1_u8(pu1_src);
                src1 = vld1_u8(pu1_src + 8);

                pu1_src += i4_src_stride;

                src2 = vld1_u8(pu1_src);
                src3 = vld1_u8(pu1_src + 8);

                pu1_src += i4_src_stride;

                src4 = vld1_u8(pu1_src);
                src5 = vld1_u8(pu1_src + 8);

                pu1_src += i4_src_stride;

                src6 = vld1_u8(pu1_src);
                src7 = vld1_u8(pu1_src + 8);

                pu1_src += i4_src_stride;

                pred0 = vld1_u8(pu1_pred);
                pred1 = vld1_u8(pu1_pred + 8);

                pu1_pred += i4_pred_stride;

                pred2 = vld1_u8(pu1_pred);
                pred3 = vld1_u8(pu1_pred + 8);

                pu1_pred += i4_pred_stride;

                pred4 = vld1_u8(pu1_pred);
                pred5 = vld1_u8(pu1_pred + 8);

                pu1_pred += i4_pred_stride;

                pred6 = vld1_u8(pu1_pred);
                pred7 = vld1_u8(pu1_pred + 8);

                pu1_pred += i4_pred_stride;

                res0_u16x8 = vsubl_u8(src0, pred0);
                res1_u16x8 = vsubl_u8(src1, pred1);
                res2_u16x8 = vsubl_u8(src2, pred2);
                res3_u16x8 = vsubl_u8(src3, pred3);
                res4_u16x8 = vsubl_u8(src4, pred4);
                res5_u16x8 = vsubl_u8(src5, pred5);
                res6_u16x8 = vsubl_u8(src6, pred6);
                res7_u16x8 = vsubl_u8(src7, pred7);

                res0_16x8 = vreinterpretq_s16_u16(res0_u16x8);
                res1_16x8 = vreinterpretq_s16_u16(res1_u16x8);
                res2_16x8 = vreinterpretq_s16_u16(res2_u16x8);
                res3_16x8 = vreinterpretq_s16_u16(res3_u16x8);
                res4_16x8 = vreinterpretq_s16_u16(res4_u16x8);
                res5_16x8 = vreinterpretq_s16_u16(res5_u16x8);
                res6_16x8 = vreinterpretq_s16_u16(res6_u16x8);
                res7_16x8 = vreinterpretq_s16_u16(res7_u16x8);

                respred0_16x8 = vld1q_s16(pi2_res);
                respred1_16x8 = vld1q_s16(pi2_res + 8);

                pi2_res += i4_res_stride;

                respred2_16x8 = vld1q_s16(pi2_res);
                respred3_16x8 = vld1q_s16(pi2_res + 8);

                pi2_res += i4_res_stride;

                respred4_16x8 = vld1q_s16(pi2_res);
                respred5_16x8 = vld1q_s16(pi2_res + 8);

                pi2_res += i4_res_stride;

                respred6_16x8 = vld1q_s16(pi2_res);
                respred7_16x8 = vld1q_s16(pi2_res + 8);

                pi2_res += i4_res_stride;

                temp0_16x8 = vsubq_s16(res0_16x8, respred0_16x8);
                temp1_16x8 = vsubq_s16(res1_16x8, respred1_16x8);
                temp2_16x8 = vsubq_s16(res2_16x8, respred2_16x8);
                temp3_16x8 = vsubq_s16(res3_16x8, respred3_16x8);
                temp4_16x8 = vsubq_s16(res4_16x8, respred4_16x8);
                temp5_16x8 = vsubq_s16(res5_16x8, respred5_16x8);
                temp6_16x8 = vsubq_s16(res6_16x8, respred6_16x8);
                temp7_16x8 = vsubq_s16(res7_16x8, respred7_16x8);

                temp0_16x8 = vabsq_s16(temp0_16x8);
                temp1_16x8 = vabsq_s16(temp1_16x8);
                temp2_16x8 = vabsq_s16(temp2_16x8);
                temp3_16x8 = vabsq_s16(temp3_16x8);
                temp4_16x8 = vabsq_s16(temp4_16x8);
                temp5_16x8 = vabsq_s16(temp5_16x8);
                temp6_16x8 = vabsq_s16(temp6_16x8);
                temp7_16x8 = vabsq_s16(temp7_16x8);

                temp0_16x8 = vaddq_s16(temp0_16x8, temp1_16x8);
                temp1_16x8 = vaddq_s16(temp2_16x8, temp3_16x8);
                temp2_16x8 = vaddq_s16(temp4_16x8, temp5_16x8);
                temp3_16x8 = vaddq_s16(temp6_16x8, temp7_16x8);

                temp0_16x8 = vaddq_s16(temp0_16x8, temp1_16x8);
                temp1_16x8 = vaddq_s16(temp2_16x8, temp3_16x8);

                temp0_16x8 = vaddq_s16(temp0_16x8, temp1_16x8);

                temp0_32x4 = vpaddlq_s16(temp0_16x8);
                temp0_32x2 = vpadd_s32(vget_low_s32(temp0_32x4), vget_high_s32(temp0_32x4));

                u4_sad += vget_lane_s32(temp0_32x2, 0);
                u4_sad += vget_lane_s32(temp0_32x2, 1);
            }
        }
    }
    else
    {
        for(i = 0; i < u4_mb_ht; i++)
        {
            for(j = 0; j < u4_mb_wd; j++)
            {
                WORD16 i2_src = pu1_src[j + i * i4_src_stride];
                WORD16 i2_pred = pu1_pred[j + i * i4_pred_stride];
                WORD16 i2_res = pi2_res[j + i * i4_res_stride];
                u4_sad += ABS(i2_src - i2_pred - i2_res);
            }
        }
    }

    return u4_sad;
}
