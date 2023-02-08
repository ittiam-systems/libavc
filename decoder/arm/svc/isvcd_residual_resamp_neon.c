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
*  isvcd_residual_resamp_neonintr.c
*
* @brief
*  Contains routines that resample for SVC resampling
*
* @author
*  Kishore
*
*  @par List of Functions:
*  - isvcd_pred_residual_recon_4x4_neonintr()
*  - isvcd_pred_residual_recon_8x8_neonintr()
*  - isvcd_pred_residual_recon_16x16_neonintr()
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

/*!
 **************************************************************************
 * \file isvcd_residual_resamp_neonintr.c
 *
 * \brief
 *    Contains routines that resample for SVC resampling
 *
 * Detailed_description
 *
 * \date
 *
 *
 * \author : kishore
 **************************************************************************
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
/*  Function Name : isvcd_residual_luma_dyadic_neonintr                       */
/*                                                                           */
/*  Description   : this fucntion does the upsampling of luma residuals for  */
/*                  Dyadic cases                                             */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : Residual upsampling context      */
/*                  pu1_inp_data : input 8 bit data pointer                  */
/*                  i4_inp_data_stride : input buffer stride                 */
/*                  pi2_out_res : output 16 bit buffer pointer               */
/*                  i4_out_res_stride : Output buffer stride                 */
/*                  pu1_inp_bitmap : input packed sign bit data pointer      */
/*                  i4_inp_bitmap_stride : sign bit buffer stride            */
/*                  ps_ref_mb_mode : reference mb mode pointer of base layer */
/*                  ps_coord : mb co-ordinate pointer                        */
/*  Globals       : none                                                     */
/*  Processing    : it does the upsampling with fixed phase values and       */
/*                  reference layer transform size                           */
/*  Outputs       : Upsampled residuals for luma                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/

void isvcd_residual_luma_dyadic_neonintr(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                         WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                         WORD32 i4_out_res_stride, mem_element_t *ps_ref_mb_mode,
                                         UWORD16 u2_mb_x, UWORD16 u2_mb_y, WORD32 i4_ref_nnz,
                                         WORD32 i4_ref_tx_size)
{
    WORD16 *pi2_refarray_buffer;
    WORD32 i4_blk_ctr;
    residual_sampling_ctxt_t *ps_ctxt;
    UNUSED(ps_ref_mb_mode);
    UNUSED(u2_mb_x);
    UNUSED(u2_mb_y);

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    pi2_refarray_buffer = ps_ctxt->pi2_refarray_buffer;

    /* based on transform size the counter and interpolation width and */
    /* height are intialised as follows                                */
    if((i4_ref_tx_size) && (0 != i4_ref_nnz))
    {
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        WORD32 i4_i, i4_j;
        /* ----------- Horizontal Interpolation ---------------- */
        int16x8_t i2_coeff_add_16x8_r0;
        int16x8_t i2_coeff_16x8_r0_0, i2_coeff_16x8_r0_1;
        int16x8_t i2_coeff_16x8_sl_r0_0, i2_coeff_16x8_sl_r0_1;
        int16x8_t result_16x8_r0_0, result_16x8_r0_1;
        int16x8_t final_result_16x8_r0_0, final_result_16x8_r0_1;

        int16x8_t i2_coeff_add_16x8_r1;
        int16x8_t i2_coeff_16x8_r1_0, i2_coeff_16x8_r1_1;
        int16x8_t i2_coeff_16x8_sl_r1_0, i2_coeff_16x8_sl_r1_1;
        int16x8_t result_16x8_r1_0, result_16x8_r1_1;
        int16x8_t final_result_16x8_r1_0, final_result_16x8_r1_1;
        int16x8x2_t result_16x8x2_t_0;

        pi2_ref_data_byte = pi2_inp_data;

        /* ----------- Horizontal Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

        for(i4_i = 0; i4_i < BLOCK_HEIGHT; i4_i += 2)
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

            result_16x8x2_t_0 = vzipq_s16(result_16x8_r0_0, result_16x8_r0_1);
            final_result_16x8_r0_0 = result_16x8x2_t_0.val[0];
            final_result_16x8_r0_1 = result_16x8x2_t_0.val[1];

            result_16x8x2_t_0 = vzipq_s16(result_16x8_r1_0, result_16x8_r1_1);
            final_result_16x8_r1_0 = result_16x8x2_t_0.val[0];
            final_result_16x8_r1_1 = result_16x8x2_t_0.val[1];

            vst1q_s32(pi4_ref_array + 1, vmovl_s16(vget_low_s16(final_result_16x8_r0_0)));
            vst1q_s32(pi4_ref_array + 5, vmovl_s16(vget_high_s16(final_result_16x8_r0_0)));
            vst1q_s32(pi4_ref_array + 9, vmovl_s16(vget_low_s16(final_result_16x8_r0_1)));
            vst1q_s32(pi4_ref_array + 13, vmovl_s16(vget_high_s16(final_result_16x8_r0_1)));
            pi4_ref_array[0] = pi2_ref_data_byte[0] << 2;
            pi4_ref_array[15] = pi2_ref_data_byte[7] << 2;
            pi4_ref_array += 16;
            pi2_ref_data_byte += i4_inp_data_stride;

            vst1q_s32(pi4_ref_array + 1, vmovl_s16(vget_low_s16(final_result_16x8_r1_0)));
            vst1q_s32(pi4_ref_array + 5, vmovl_s16(vget_high_s16(final_result_16x8_r1_0)));
            vst1q_s32(pi4_ref_array + 9, vmovl_s16(vget_low_s16(final_result_16x8_r1_1)));
            vst1q_s32(pi4_ref_array + 13, vmovl_s16(vget_high_s16(final_result_16x8_r1_1)));

            pi4_ref_array[0] = pi2_ref_data_byte[0] << 2;
            pi4_ref_array[15] = pi2_ref_data_byte[7] << 2;
            pi4_ref_array += 16;
            /* vertical loop uopdates */
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
                pi4_ref_array_temp += MB_WIDTH;
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
            if(0 != (i4_ref_nnz & 0x1))
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

                i2_res_16x8_r0_0 = vzipq_s16(i2_res_16x8_r0_0, i2_res_16x8_r0_1).val[0];
                i2_res_16x8_r1_0 = vzipq_s16(i2_res_16x8_r1_0, i2_res_16x8_r1_1).val[0];
                i2_res_16x8_r2_0 = vzipq_s16(i2_res_16x8_r2_0, i2_res_16x8_r2_1).val[0];
                i2_res_16x8_r3_0 = vzipq_s16(i2_res_16x8_r3_0, i2_res_16x8_r3_1).val[0];

                i2_coeff1_16x8_r0_0 = vshlq_n_s16(i2_coeff1_16x8_r0_0, 1);
                i2_coeff1_16x8_r1_0 = vshlq_n_s16(i2_coeff1_16x8_r1_0, 1);
                i2_coeff1_16x8_r2_0 = vshlq_n_s16(i2_coeff1_16x8_r2_0, 1);
                i2_coeff1_16x8_r3_0 = vshlq_n_s16(i2_coeff1_16x8_r3_0, 1);

                vst1q_s16(pi2_refarray_buffer + 1, i2_res_16x8_r0_0);
                vst1q_lane_s16(pi2_refarray_buffer, i2_coeff1_16x8_r0_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 7, i2_coeff1_16x8_r0_0, 3);

                vst1q_s16(pi2_refarray_buffer + 9, i2_res_16x8_r1_0);
                vst1q_lane_s16(pi2_refarray_buffer + 8, i2_coeff1_16x8_r1_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 15, i2_coeff1_16x8_r1_0, 3);

                vst1q_s16(pi2_refarray_buffer + 17, i2_res_16x8_r2_0);
                vst1q_lane_s16(pi2_refarray_buffer + 16, i2_coeff1_16x8_r2_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 23, i2_coeff1_16x8_r2_0, 3);

                vst1q_s16(pi2_refarray_buffer + 25, i2_res_16x8_r3_0);
                vst1q_lane_s16(pi2_refarray_buffer + 24, i2_coeff1_16x8_r3_0, 0);
                vst1q_lane_s16(pi2_refarray_buffer + 31, i2_coeff1_16x8_r3_0, 3);

                {
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

                    i4_horz_add_32x4_r1_1 =
                        vaddl_s16(i4_horz_samp_16x4_r0_1, i4_horz_samp_16x4_r1_1);
                    i4_horz_add_32x4_r1_2 =
                        vaddl_s16(i4_horz_samp_16x4_r0_2, i4_horz_samp_16x4_r1_2);

                    i4_horz_add_32x4_r2_1 =
                        vaddl_s16(i4_horz_samp_16x4_r1_1, i4_horz_samp_16x4_r2_1);
                    i4_horz_add_32x4_r2_2 =
                        vaddl_s16(i4_horz_samp_16x4_r1_2, i4_horz_samp_16x4_r2_2);

                    i4_horz_add_32x4_r3_1 =
                        vaddl_s16(i4_horz_samp_16x4_r2_1, i4_horz_samp_16x4_r3_1);
                    i4_horz_add_32x4_r3_2 =
                        vaddl_s16(i4_horz_samp_16x4_r2_2, i4_horz_samp_16x4_r3_2);

                    i4_horz_samp_32x4_r0_1 = vshll_n_s16(i4_horz_samp_16x4_r0_1, 1);
                    i4_horz_samp_32x4_r0_2 = vshll_n_s16(i4_horz_samp_16x4_r0_2, 1);

                    i4_horz_samp_32x4_r1_1 = vshll_n_s16(i4_horz_samp_16x4_r1_1, 1);
                    i4_horz_samp_32x4_r1_2 = vshll_n_s16(i4_horz_samp_16x4_r1_2, 1);

                    i4_horz_samp_32x4_r2_1 = vshll_n_s16(i4_horz_samp_16x4_r2_1, 1);
                    i4_horz_samp_32x4_r2_2 = vshll_n_s16(i4_horz_samp_16x4_r2_2, 1);

                    i4_horz_samp_32x4_r3_1 = vshll_n_s16(i4_horz_samp_16x4_r3_1, 1);
                    i4_horz_samp_32x4_r3_2 = vshll_n_s16(i4_horz_samp_16x4_r3_2, 1);

                    i4_horz_res_32x4_r1_1 =
                        vaddq_s32(i4_horz_samp_32x4_r0_1, i4_horz_add_32x4_r1_1);
                    i4_horz_res_32x4_r1_2 =
                        vaddq_s32(i4_horz_samp_32x4_r0_2, i4_horz_add_32x4_r1_2);

                    i4_horz_res_32x4_r2_1 =
                        vaddq_s32(i4_horz_samp_32x4_r1_1, i4_horz_add_32x4_r1_1);
                    i4_horz_res_32x4_r2_2 =
                        vaddq_s32(i4_horz_samp_32x4_r1_2, i4_horz_add_32x4_r1_2);

                    i4_horz_res_32x4_r3_1 =
                        vaddq_s32(i4_horz_samp_32x4_r1_1, i4_horz_add_32x4_r2_1);
                    i4_horz_res_32x4_r3_2 =
                        vaddq_s32(i4_horz_samp_32x4_r1_2, i4_horz_add_32x4_r2_2);

                    i4_horz_res_32x4_r4_1 =
                        vaddq_s32(i4_horz_samp_32x4_r2_1, i4_horz_add_32x4_r2_1);
                    i4_horz_res_32x4_r4_2 =
                        vaddq_s32(i4_horz_samp_32x4_r2_2, i4_horz_add_32x4_r2_2);

                    i4_horz_res_32x4_r5_1 =
                        vaddq_s32(i4_horz_samp_32x4_r2_1, i4_horz_add_32x4_r3_1);
                    i4_horz_res_32x4_r5_2 =
                        vaddq_s32(i4_horz_samp_32x4_r2_2, i4_horz_add_32x4_r3_2);

                    i4_horz_res_32x4_r6_1 =
                        vaddq_s32(i4_horz_samp_32x4_r3_1, i4_horz_add_32x4_r3_1);
                    i4_horz_res_32x4_r6_2 =
                        vaddq_s32(i4_horz_samp_32x4_r3_2, i4_horz_add_32x4_r3_2);

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

                    pi2_out_res += BLOCK_WIDTH;
                }
            }
            else
            {
                pi2_out_res += BLOCK_WIDTH;
            }

            /* Block level loop updates */
            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= SUB_BLOCK_WIDTH;
                pi2_inp_data += (i4_inp_data_stride * SUB_BLOCK_HEIGHT);
                pi2_out_res -= MB_WIDTH;
                pi2_out_res += (i4_out_res_stride * BLOCK_HEIGHT);
                i4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += SUB_BLOCK_WIDTH;
            }

            i4_ref_nnz >>= 1;
        } /* end of loop over all the blocks */
    }

    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_residual_neonintr                       */
/*                                                                           */
/*  Description   : this fucntion does the upsampling of residuals.          */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : Residual upsampling context      */
/*                  pu1_inp_data : input 8 bit data pointer                  */
/*                  i4_inp_data_stride : input buffer stride                 */
/*                  pi2_out_res : output 16 bit buffer pointer               */
/*                  i4_out_res_stride : Output buffer stride                 */
/*                  pu1_inp_bitmap : input packed sign bit data pointer      */
/*                  i4_inp_bitmap_stride : sign bit buffer stride            */
/*                  ps_ref_mb_mode : reference mb mode pointer of base layer */
/*                  ps_coord : mb co-ordinate pointer                        */
/*  Globals       : none                                                     */
/*  Processing    : it does the upsampling with fixed phase values and       */
/*                  reference layer transform size                           */
/*  Outputs       : Upsampled residuals.                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 05 2021   Dolan          creation                              */
/*                                                                           */
/*****************************************************************************/

void isvcd_interpolate_residual_neonintr(void *pv_residual_samp_ctxt, WORD16 *pi2_out,
                                         WORD32 i4_out_stride, WORD32 i4_refarray_wd,
                                         UWORD16 u2_mb_x, UWORD16 u2_mb_y, WORD32 i4_chroma_flag)
{
    residual_sampling_ctxt_t *ps_ctxt;
    residual_samp_map_ctxt_t *ps_map_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;

    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    WORD32 i4_temp_array_ht;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD16 *pi2_ref_array;
    UWORD8 *pu1_ref_x_ptr_incr, *pu1_ref_y_ptr_incr;

    UWORD8 arr_y_ref_pos[16] = {0};
    UWORD8 arr_x_ref_pos[16] = {0};
    UWORD8 arr_x_ref_pos_low[16] = {0};
    UWORD8 arr_x_phase[16] = {0};
    UWORD8 arr_y_phase[16] = {0};
    UWORD8 *pi1_y_ref_pos;
    UWORD8 *pi1_x_ref_pos;
    UWORD8 *pi1_x_ref_pos_low;
    UWORD8 *pi1_y_phase;
    UWORD8 *pi1_x_phase;

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi2_ref_array = ps_ctxt->pi2_refarray_buffer;
    pu1_ref_x_ptr_incr = ps_ctxt->pu1_ref_x_ptr_incr;
    pu1_ref_y_ptr_incr = ps_ctxt->pu1_ref_y_ptr_incr;

    if(1 == i4_chroma_flag)
        ps_map_ctxt = &ps_lyr_ctxt->s_chroma_map_ctxt;
    else
        ps_map_ctxt = &ps_lyr_ctxt->s_luma_map_ctxt;

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

    i4_temp_array_ht = i4_mb_ht;
    i4_frm_mb_y = u2_mb_y * i4_mb_ht;
    i4_frm_mb_x = u2_mb_x * i4_mb_wd;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation                                                */
    /* --------------------------------------------------------------------- */
    if(i4_chroma_flag == 0)
    {
        int16x8_t ref_arr_16x8_r0_0, ref_arr_16x8_r0_1;
        int16x8_t ref_arr_16x8_r1_0, ref_arr_16x8_r1_1;
        uint8x16_t x_ref_pos_mask_r0_0, x_ref_rnd_mask_r0_0;
        uint8x16_t u1_incr_8x16_r0_0, x_ref_pos_mask_temp_r0_0, u1_incr_not_8x16_r0_0,
            u1_y_incr_8x16_r0_0, phs_mask_8x16_0;
        uint8x16_t u1_incr_8x16_r1_0, x_ref_pos_mask_temp_r1_0, u1_incr_not_8x16_r1_0;
        int16x8_t ref_arr_temp0_16x8_r0_0, res_16x8_r0_0, vert_res_16x8_r0_0;
        int16x8_t ref_arr_temp0_16x8_r1_0, res_16x8_r1_0, vert_res_16x8_r1_0;
        int16x8_t ref_arr_temp1_16x8_r0_0;
        int16x8_t ref_arr_temp1_16x8_r1_0;

        uint8x16_t x_ref_pos_mask_temp_r0_1, u1_incr_not_8x16_r0_1;
        uint8x16_t x_ref_pos_mask_temp_r1_1, u1_incr_not_8x16_r1_1;
        int16x8_t ref_arr_temp0_16x8_r0_1, res_16x8_r0_1, vert_res_16x8_r0_1;
        int16x8_t ref_arr_temp0_16x8_r1_1, res_16x8_r1_1, vert_res_16x8_r1_1;
        int16x8_t ref_arr_temp1_16x8_r0_1;
        int16x8_t ref_arr_temp1_16x8_r1_1;

        uint16x8_t u1_y_incr_16x8_r0_0, u1_y_incr_16x8_r0_1;

        uint8x16_t u1_incr_not_8x16_r0_0_even, u1_incr_not_8x16_r1_0_even,
            x_ref_pos_mask_temp_r0_0_even, x_ref_pos_mask_temp_r1_0_even;
        uint8x16_t u1_incr_not_8x16_r0_0_odd, u1_incr_not_8x16_r1_0_odd,
            x_ref_pos_mask_temp_r0_0_odd, x_ref_pos_mask_temp_r1_0_odd;
        uint8x16x2_t u1_incr_not_8x16_2, x_ref_pos_mask_temp;
        int16x8_t prev_res_16x8_r0_0;
        int16x8_t prev_res_16x8_r1_0;
        int16x8_t prev_res_16x8_r0_1;
        int16x8_t prev_res_16x8_r1_1;
        uint8x8x2_t u1_incr_8x8x2_t;
        uint8x8_t u1_incr_8x8_t0, u1_incr_8x8_t1;
        uint16x8_t u1_prev_y_incr_16x8_r0_0;
        uint16x8_t u1_prev_y_incr_16x8_r0_1;

        WORD32 zero_r0_r1 = 0;

        int32x4_t res_32x4_l_0, res_32x4_h_0;
        int32x4_t res_32x4_l_1, res_32x4_h_1;
        int16x8_t res_16x8_l, res_16x8_h;
        uint16x8_t phs_mask_16x8_0, phs_mask_16x8_1;
        int16x8_t const_16_16x8, phs_mask_16min_16x8_0;
        int16x8_t dup_val_1, dup_val_2, dup_val_3, dup_val_4, dup_val_5, dup_abs;
        uint8x16_t phs_mask_div8_8x16_0, mid_indx_8x16;
        int16x8_t phs_mask_16min_16x8_1;
        uint16x8_t ones = vdupq_n_u16(0xFFFF);
        uint8x16_t const_ones = vdupq_n_u8(1);
        uint8x8x2_t u1_temp_8x8x2_t;
        uint8x8_t u1_temp_8x8_t0, u1_temp_8x8_t1;

        WORD16 *pi2_ref_array_temp;
        UWORD8 *pu1_ref_x_ptr_incr_temp, *pu1_ref_y_ptr_incr_temp;
        WORD32 i4_y_phase;
        WORD32 out_stride_temp;
        WORD32 strt_indx_h = 0;

        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            arr_y_phase[i4_y] = (UWORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            arr_y_ref_pos[i4_y] = (UWORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
        }
        pi1_y_ref_pos = arr_y_ref_pos;
        pi1_y_phase = arr_y_phase;

        strt_indx_h = (ps_x_pos_phase[8 + i4_frm_mb_x].i2_ref_pos);
        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            arr_x_ref_pos[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            arr_x_phase[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
        }

        pi1_x_ref_pos = arr_x_ref_pos;
        pi1_x_phase = arr_x_phase;

        phs_mask_8x16_0 = vld1q_u8((pi1_x_phase));
        phs_mask_16x8_0 = vmovl_u8(vget_low_u8(phs_mask_8x16_0));
        phs_mask_16x8_1 = vmovl_u8(vld1_u8((pi1_x_phase + 8)));
        x_ref_pos_mask_r0_0 = vld1q_u8((pi1_x_ref_pos));
        const_16_16x8 = vdupq_n_s16(16);
        phs_mask_div8_8x16_0 = vshrq_n_u8(phs_mask_8x16_0, 3);

        phs_mask_16min_16x8_0 = vsubq_s16(const_16_16x8, vreinterpretq_s16_u16(phs_mask_16x8_0));
        phs_mask_16min_16x8_1 = vsubq_s16(const_16_16x8, vreinterpretq_s16_u16(phs_mask_16x8_1));

        x_ref_rnd_mask_r0_0 = vaddq_u8(x_ref_pos_mask_r0_0, phs_mask_div8_8x16_0);
        mid_indx_8x16 = vdupq_n_u8((strt_indx_h << 1));

        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            if((i4_y > 0) && (pi1_y_ref_pos[i4_y] == pi1_y_ref_pos[i4_y - 1]))
            {
                if(!zero_r0_r1)
                {
                    res_16x8_l = vdupq_n_s16(0);
                    res_16x8_h = vdupq_n_s16(0);

                    out_stride_temp = (i4_y * i4_out_stride);
                    vst1q_s16((pi2_out + out_stride_temp), res_16x8_l);
                    vst1q_s16((pi2_out + out_stride_temp + 8), res_16x8_h);
                    continue;
                }

                res_16x8_r0_0 = prev_res_16x8_r0_0;
                res_16x8_r1_0 = prev_res_16x8_r1_0;
                res_16x8_r0_1 = prev_res_16x8_r0_1;
                res_16x8_r1_1 = prev_res_16x8_r1_1;

                u1_y_incr_16x8_r0_0 = u1_prev_y_incr_16x8_r0_0;
                u1_y_incr_16x8_r0_1 = u1_prev_y_incr_16x8_r0_1;
            }
            else
            {
                pi2_ref_array_temp = pi2_ref_array + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                pu1_ref_x_ptr_incr_temp =
                    pu1_ref_x_ptr_incr + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                ref_arr_16x8_r0_0 = vld1q_s16((pi2_ref_array_temp));
                ref_arr_16x8_r1_0 = vld1q_s16((pi2_ref_array_temp + i4_refarray_wd));
                ref_arr_16x8_r0_1 = vld1q_s16((pi2_ref_array_temp + strt_indx_h));
                ref_arr_16x8_r1_1 = vld1q_s16((pi2_ref_array_temp + i4_refarray_wd + strt_indx_h));

                dup_val_1 = vabsq_s16(ref_arr_16x8_r0_0);
                dup_val_2 = vabsq_s16(ref_arr_16x8_r1_0);
                dup_val_3 = vabsq_s16(ref_arr_16x8_r0_1);
                dup_val_4 = vabsq_s16(ref_arr_16x8_r1_1);
                dup_val_5 = vqaddq_s16(dup_val_1, dup_val_2);
                dup_val_1 = vqaddq_s16(dup_val_3, dup_val_4);
                dup_abs = vqaddq_s16(dup_val_1, dup_val_5);
                zero_r0_r1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] ||
                             dup_abs[5] || dup_abs[6] || dup_abs[7];
                if(zero_r0_r1)
                {
                    u1_incr_8x16_r0_0 = vld1q_u8((pu1_ref_x_ptr_incr_temp));
                    u1_incr_8x16_r1_0 = vld1q_u8((pu1_ref_x_ptr_incr_temp + i4_refarray_wd));
                    u1_incr_8x8x2_t.val[0] = vget_low_u8(u1_incr_8x16_r0_0);
                    u1_incr_8x8x2_t.val[1] = vget_high_u8(u1_incr_8x16_r0_0);
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x16_r0_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                    u1_incr_8x8x2_t.val[0] = vget_low_u8(u1_incr_8x16_r1_0);
                    u1_incr_8x8x2_t.val[1] = vget_high_u8(u1_incr_8x16_r1_0);
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x16_r1_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);
                    u1_incr_not_8x16_r0_0 = vbicq_u8(phs_mask_div8_8x16_0, u1_incr_8x16_r0_0);
                    u1_incr_not_8x16_r1_0 = vbicq_u8(phs_mask_div8_8x16_0, u1_incr_8x16_r1_0);

                    u1_incr_not_8x16_r0_0 = vaddq_u8(u1_incr_not_8x16_r0_0, x_ref_pos_mask_r0_0);
                    u1_incr_not_8x16_r1_0 = vaddq_u8(u1_incr_not_8x16_r1_0, x_ref_pos_mask_r0_0);

                    x_ref_pos_mask_temp_r0_0 = vaddq_u8(u1_incr_not_8x16_r0_0, u1_incr_8x16_r0_0);
                    x_ref_pos_mask_temp_r1_0 = vaddq_u8(u1_incr_not_8x16_r1_0, u1_incr_8x16_r1_0);

                    u1_incr_not_8x16_r0_0_even = vshlq_n_u8(u1_incr_not_8x16_r0_0, 1);
                    u1_incr_not_8x16_r1_0_even = vshlq_n_u8(u1_incr_not_8x16_r1_0, 1);
                    x_ref_pos_mask_temp_r0_0_even = vshlq_n_u8(x_ref_pos_mask_temp_r0_0, 1);
                    x_ref_pos_mask_temp_r1_0_even = vshlq_n_u8(x_ref_pos_mask_temp_r1_0, 1);

                    u1_incr_not_8x16_r0_0_odd = vaddq_u8(u1_incr_not_8x16_r0_0_even, const_ones);
                    u1_incr_not_8x16_r1_0_odd = vaddq_u8(u1_incr_not_8x16_r1_0_even, const_ones);
                    x_ref_pos_mask_temp_r0_0_odd =
                        vaddq_u8(x_ref_pos_mask_temp_r0_0_even, const_ones);
                    x_ref_pos_mask_temp_r1_0_odd =
                        vaddq_u8(x_ref_pos_mask_temp_r1_0_even, const_ones);

                    u1_incr_not_8x16_2 =
                        vzipq_u8(u1_incr_not_8x16_r0_0_even, u1_incr_not_8x16_r0_0_odd);

                    u1_incr_not_8x16_r0_0 = u1_incr_not_8x16_2.val[0];
                    u1_incr_not_8x16_r0_1 = u1_incr_not_8x16_2.val[1];

                    u1_incr_not_8x16_2 =
                        vzipq_u8(u1_incr_not_8x16_r1_0_even, u1_incr_not_8x16_r1_0_odd);
                    u1_incr_not_8x16_r1_0 = u1_incr_not_8x16_2.val[0];
                    u1_incr_not_8x16_r1_1 = u1_incr_not_8x16_2.val[1];

                    x_ref_pos_mask_temp =
                        vzipq_u8(x_ref_pos_mask_temp_r0_0_even, x_ref_pos_mask_temp_r0_0_odd);
                    x_ref_pos_mask_temp_r0_0 = x_ref_pos_mask_temp.val[0];
                    x_ref_pos_mask_temp_r0_1 = x_ref_pos_mask_temp.val[1];

                    x_ref_pos_mask_temp =
                        vzipq_u8(x_ref_pos_mask_temp_r1_0_even, x_ref_pos_mask_temp_r1_0_odd);
                    x_ref_pos_mask_temp_r1_0 = x_ref_pos_mask_temp.val[0];
                    x_ref_pos_mask_temp_r1_1 = x_ref_pos_mask_temp.val[1];
                    u1_incr_not_8x16_r0_1 = vsubq_u8(u1_incr_not_8x16_r0_1, mid_indx_8x16);
                    u1_incr_not_8x16_r1_1 = vsubq_u8(u1_incr_not_8x16_r1_1, mid_indx_8x16);
                    x_ref_pos_mask_temp_r0_1 = vsubq_u8(x_ref_pos_mask_temp_r0_1, mid_indx_8x16);
                    x_ref_pos_mask_temp_r1_1 = vsubq_u8(x_ref_pos_mask_temp_r1_1, mid_indx_8x16);

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r0_0));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r0_0));
                    ref_arr_temp0_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r1_0));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r1_0));
                    ref_arr_temp0_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r0_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r0_0));
                    ref_arr_temp1_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r1_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r1_0));
                    ref_arr_temp1_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_1)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r0_1));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r0_1));
                    ref_arr_temp0_16x8_r0_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_1)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r1_1));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r1_1));
                    ref_arr_temp0_16x8_r1_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_1)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r0_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r0_1));
                    ref_arr_temp1_16x8_r0_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_1)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_1)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r1_1));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r1_1));
                    ref_arr_temp1_16x8_r1_1 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    res_16x8_r0_0 = vmulq_s16(ref_arr_temp0_16x8_r0_0, phs_mask_16min_16x8_0);
                    res_16x8_r1_0 = vmulq_s16(ref_arr_temp0_16x8_r1_0, phs_mask_16min_16x8_0);
                    res_16x8_r0_0 = vmlaq_s16(res_16x8_r0_0, ref_arr_temp1_16x8_r0_0,
                                              vreinterpretq_s16_u16(phs_mask_16x8_0));
                    res_16x8_r1_0 = vmlaq_s16(res_16x8_r1_0, ref_arr_temp1_16x8_r1_0,
                                              vreinterpretq_s16_u16(phs_mask_16x8_0));
                    res_16x8_r0_1 = vmulq_s16(ref_arr_temp0_16x8_r0_1, phs_mask_16min_16x8_1);
                    res_16x8_r1_1 = vmulq_s16(ref_arr_temp0_16x8_r1_1, phs_mask_16min_16x8_1);
                    res_16x8_r0_1 = vmlaq_s16(res_16x8_r0_1, ref_arr_temp1_16x8_r0_1,
                                              vreinterpretq_s16_u16(phs_mask_16x8_1));
                    res_16x8_r1_1 = vmlaq_s16(res_16x8_r1_1, ref_arr_temp1_16x8_r1_1,
                                              vreinterpretq_s16_u16(phs_mask_16x8_1));

                    pu1_ref_y_ptr_incr_temp =
                        pu1_ref_y_ptr_incr + (pi1_y_ref_pos[i4_y] * i4_refarray_wd);
                    u1_y_incr_8x16_r0_0 = vld1q_u8((pu1_ref_y_ptr_incr_temp));

                    u1_incr_8x8x2_t.val[0] = vget_low_u8(u1_y_incr_8x16_r0_0);
                    u1_incr_8x8x2_t.val[1] = vget_high_u8(u1_y_incr_8x16_r0_0);
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_rnd_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_rnd_mask_r0_0));
                    u1_y_incr_8x16_r0_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);
                    u1_y_incr_16x8_r0_0 = vmovl_u8(vget_low_u8(u1_y_incr_8x16_r0_0));
                    u1_y_incr_16x8_r0_1 = vmovl_u8(vget_high_u8(u1_y_incr_8x16_r0_0));
                    u1_y_incr_16x8_r0_0 = vtstq_u16(u1_y_incr_16x8_r0_0, ones);
                    u1_y_incr_16x8_r0_1 = vtstq_u16(u1_y_incr_16x8_r0_1, ones);

                    prev_res_16x8_r0_0 = res_16x8_r0_0;
                    prev_res_16x8_r1_0 = res_16x8_r1_0;
                    prev_res_16x8_r0_1 = res_16x8_r0_1;
                    prev_res_16x8_r1_1 = res_16x8_r1_1;

                    u1_prev_y_incr_16x8_r0_0 = u1_y_incr_16x8_r0_0;
                    u1_prev_y_incr_16x8_r0_1 = u1_y_incr_16x8_r0_1;
                }
            }

            if(!zero_r0_r1)
            {
                res_16x8_l = vdupq_n_s16(0);
                res_16x8_h = vdupq_n_s16(0);
            }
            else
            {
                i4_y_phase = pi1_y_phase[i4_y];

                if((i4_y_phase) >> 3)
                {
                    vert_res_16x8_r0_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r0_0, res_16x8_r1_0);
                    vert_res_16x8_r1_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r1_0, res_16x8_r1_0);
                    vert_res_16x8_r0_1 =
                        vbslq_s16(u1_y_incr_16x8_r0_1, res_16x8_r0_1, res_16x8_r1_1);
                    vert_res_16x8_r1_1 =
                        vbslq_s16(u1_y_incr_16x8_r0_1, res_16x8_r1_1, res_16x8_r1_1);
                }
                else
                {
                    vert_res_16x8_r0_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r0_0, res_16x8_r0_0);
                    vert_res_16x8_r1_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r1_0, res_16x8_r0_0);
                    vert_res_16x8_r0_1 =
                        vbslq_s16(u1_y_incr_16x8_r0_1, res_16x8_r0_1, res_16x8_r0_1);
                    vert_res_16x8_r1_1 =
                        vbslq_s16(u1_y_incr_16x8_r0_1, res_16x8_r1_1, res_16x8_r0_1);
                }

                res_32x4_l_0 = vmull_n_s16(vget_low_s16(vert_res_16x8_r0_0), 16 - i4_y_phase);
                res_32x4_l_0 =
                    vmlal_n_s16(res_32x4_l_0, vget_low_s16(vert_res_16x8_r1_0), i4_y_phase);

                res_32x4_l_1 = vmull_n_s16(vget_high_s16(vert_res_16x8_r0_0), 16 - i4_y_phase);
                res_32x4_l_1 =
                    vmlal_n_s16(res_32x4_l_1, vget_high_s16(vert_res_16x8_r1_0), i4_y_phase);
                res_32x4_h_0 = vmull_n_s16(vget_low_s16(vert_res_16x8_r0_1), 16 - i4_y_phase);
                res_32x4_h_0 =
                    vmlal_n_s16(res_32x4_h_0, vget_low_s16(vert_res_16x8_r1_1), i4_y_phase);
                res_32x4_h_1 = vmull_n_s16(vget_high_s16(vert_res_16x8_r0_1), 16 - i4_y_phase);
                res_32x4_h_1 =
                    vmlal_n_s16(res_32x4_h_1, vget_high_s16(vert_res_16x8_r1_1), i4_y_phase);

                res_32x4_l_0 = vrshrq_n_s32(res_32x4_l_0, 8);
                res_32x4_l_1 = vrshrq_n_s32(res_32x4_l_1, 8);
                res_32x4_h_0 = vrshrq_n_s32(res_32x4_h_0, 8);
                res_32x4_h_1 = vrshrq_n_s32(res_32x4_h_1, 8);

                res_16x8_l = vcombine_s16(vmovn_s32(res_32x4_l_0), vmovn_s32(res_32x4_l_1));
                res_16x8_h = vcombine_s16(vmovn_s32(res_32x4_h_0), vmovn_s32(res_32x4_h_1));
            }

            out_stride_temp = (i4_y * i4_out_stride);
            vst1q_s16((pi2_out + out_stride_temp), res_16x8_l);
            vst1q_s16((pi2_out + out_stride_temp + 8), res_16x8_h);
        }
    }
    else
    {
        int16x8_t ref_arr_16x8_r0_0;
        int16x8_t ref_arr_16x8_r1_0;
        uint8x16_t x_ref_pos_mask_r0_0, x_ref_rnd_mask_r0_0;
        uint16x8_t u1_incr_16x8_r0_0, u1_incr_mask_16x8_r0_0, phs_mask_16x8_0,
            u1_incr_not_16x8_r0_0, u1_y_incr_16x8_r0_0;
        uint16x8_t u1_incr_16x8_r1_0, u1_incr_mask_16x8_r1_0, u1_incr_not_16x8_r1_0;
        uint8x16_t u1_incr_8x16_r0_0, x_ref_pos_mask_temp_r0_0, u1_incr_mask_8x16_r0_0,
            u1_incr_not_8x16_r0_0, u1_y_incr_8x16_r0_0;
        uint8x16_t u1_incr_8x16_r1_0, x_ref_pos_mask_temp_r1_0, u1_incr_mask_8x16_r1_0,
            u1_incr_not_8x16_r1_0;
        int16x8_t ref_arr_temp0_16x8_r0_0, res_16x8_r0_0, vert_res_16x8_r0_0;
        int16x8_t ref_arr_temp0_16x8_r1_0, res_16x8_r1_0, vert_res_16x8_r1_0;
        int16x8_t ref_arr_temp1_16x8_r0_0;
        int16x8_t ref_arr_temp1_16x8_r1_0;

        int32x4_t res_32x4_l_0;
        int32x4_t res_32x4_l_1;
        int16x8_t out_16x8_0, out_16x8_1;
        uint16x8_t phs_mask_div8_16x8_0, phs_mask_div8_msb_16x8_0;
        int16x8_t const_16_16x8, phs_mask_16min_16x8_0;
        uint8x8x2_t u1_temp_8x8x2_t;
        uint8x8_t u1_temp_8x8_t0, u1_temp_8x8_t1;

        uint16x8_t chroma_mask_16x8 = vreinterpretq_u16_u32(vdupq_n_u32(0x0000ffff));
        uint16x8_t ones = vdupq_n_u16(0xFFFF);

        WORD16 *pi2_ref_array_temp;
        UWORD8 *pu1_ref_x_ptr_incr_temp, *pu1_ref_y_ptr_incr_temp;
        WORD32 i4_y_phase;
        uint8x8x2_t u1_incr_8x8x2_t;
        uint8x8_t u1_incr_8x8_t0, u1_incr_8x8_t1;
        int16x8_t prev_res_16x8_r0_0;
        int16x8_t prev_res_16x8_r1_0;
        int16x8_t dup_val_1, dup_val_2, dup_abs;
        uint16x8_t u1_prev_y_incr_16x8_r0_0;

        WORD32 out_stride_temp;
        WORD32 zero_r0_r1 = 0;
        WORD32 i4_x2 = 0;
        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            arr_y_phase[i4_y] = (WORD8) ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            arr_y_ref_pos[i4_y] = (WORD8) (ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos);
        }
        pi1_y_ref_pos = arr_y_ref_pos;
        pi1_y_phase = arr_y_phase;

        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            arr_x_ref_pos[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            arr_x_phase[i4_x] = (WORD8) ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
            i4_x2 = i4_x << 1;
            arr_x_ref_pos_low[i4_x2] = (arr_x_ref_pos[i4_x]) << 1;
            arr_x_ref_pos_low[i4_x2 + 1] = arr_x_ref_pos_low[i4_x2] + 1;
        }

        pi1_x_ref_pos_low = arr_x_ref_pos_low;
        pi1_x_phase = arr_x_phase;

        phs_mask_16x8_0 = vmovl_u8(vld1_u8((pi1_x_phase)));
        x_ref_pos_mask_r0_0 = vld1q_u8((pi1_x_ref_pos_low));
        const_16_16x8 = vdupq_n_s16(16);
        phs_mask_div8_16x8_0 = vshrq_n_u16(phs_mask_16x8_0, 3);
        phs_mask_div8_msb_16x8_0 = vsliq_n_u16(phs_mask_div8_16x8_0, phs_mask_div8_16x8_0, 8);

        phs_mask_16min_16x8_0 = vsubq_s16(const_16_16x8, vreinterpretq_s16_u16(phs_mask_16x8_0));

        x_ref_rnd_mask_r0_0 = vaddq_u8(
            x_ref_pos_mask_r0_0, vreinterpretq_u8_u16(vshlq_n_u16(phs_mask_div8_msb_16x8_0, 1)));
        for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
        {
            if((i4_y > 0) && (pi1_y_ref_pos[i4_y] == pi1_y_ref_pos[i4_y - 1]))
            {
                if(!zero_r0_r1)
                {
                    res_32x4_l_0 = vdupq_n_s32(0);
                    res_32x4_l_1 = vdupq_n_s32(0);

                    out_stride_temp = (i4_y * i4_out_stride);

                    out_16x8_0 = vld1q_s16(pi2_out + out_stride_temp);
                    out_16x8_1 = vld1q_s16(pi2_out + out_stride_temp + 8);
                    out_16x8_0 = vbslq_s16(chroma_mask_16x8, vreinterpretq_s16_s32(res_32x4_l_0),
                                           out_16x8_0);
                    out_16x8_1 = vbslq_s16(chroma_mask_16x8, vreinterpretq_s16_s32(res_32x4_l_1),
                                           out_16x8_1);
                    vst1q_s16((pi2_out + out_stride_temp), out_16x8_0);
                    vst1q_s16((pi2_out + out_stride_temp + 8), out_16x8_1);
                    continue;
                }

                res_16x8_r0_0 = prev_res_16x8_r0_0;
                res_16x8_r1_0 = prev_res_16x8_r1_0;

                u1_y_incr_16x8_r0_0 = u1_prev_y_incr_16x8_r0_0;
            }
            else
            {
                pi2_ref_array_temp = pi2_ref_array + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                pu1_ref_x_ptr_incr_temp =
                    pu1_ref_x_ptr_incr + ((pi1_y_ref_pos[i4_y]) * i4_refarray_wd);
                ref_arr_16x8_r0_0 = vld1q_s16((pi2_ref_array_temp));
                ref_arr_16x8_r1_0 = vld1q_s16((pi2_ref_array_temp + i4_refarray_wd));

                dup_val_1 = vabsq_s16(ref_arr_16x8_r0_0);
                dup_val_2 = vabsq_s16(ref_arr_16x8_r1_0);
                dup_abs = vqaddq_s16(dup_val_1, dup_val_2);
                zero_r0_r1 = dup_abs[0] || dup_abs[1] || dup_abs[2] || dup_abs[3] || dup_abs[4] ||
                             dup_abs[5] || dup_abs[6] || dup_abs[7];
                if(zero_r0_r1)
                {
                    u1_incr_16x8_r0_0 = (vmovl_u8(vld1_u8((pu1_ref_x_ptr_incr_temp))));
                    u1_incr_16x8_r1_0 =
                        (vmovl_u8(vld1_u8((pu1_ref_x_ptr_incr_temp + i4_refarray_wd))));
                    u1_incr_8x8x2_t.val[0] = vget_low_u8(vreinterpretq_u8_u16(u1_incr_16x8_r0_0));
                    u1_incr_8x8x2_t.val[1] = vget_high_u8(vreinterpretq_u8_u16(u1_incr_16x8_r0_0));
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x16_r0_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                    u1_incr_8x8x2_t.val[0] = vget_low_u8(vreinterpretq_u8_u16(u1_incr_16x8_r1_0));
                    u1_incr_8x8x2_t.val[1] = vget_high_u8(vreinterpretq_u8_u16(u1_incr_16x8_r1_0));
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_pos_mask_r0_0));
                    u1_incr_8x16_r1_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                    u1_incr_16x8_r0_0 = vreinterpretq_u16_u8(u1_incr_8x16_r0_0);
                    u1_incr_16x8_r1_0 = vreinterpretq_u16_u8(u1_incr_8x16_r1_0);
                    u1_incr_mask_16x8_r0_0 = vsliq_n_u16(u1_incr_16x8_r0_0, u1_incr_16x8_r0_0, 8);
                    u1_incr_mask_16x8_r1_0 = vsliq_n_u16(u1_incr_16x8_r1_0, u1_incr_16x8_r1_0, 8);

                    u1_incr_not_16x8_r0_0 =
                        vbicq_u16(phs_mask_div8_msb_16x8_0, u1_incr_mask_16x8_r0_0);
                    u1_incr_not_16x8_r1_0 =
                        vbicq_u16(phs_mask_div8_msb_16x8_0, u1_incr_mask_16x8_r1_0);

                    u1_incr_mask_8x16_r0_0 =
                        vreinterpretq_u8_u16(vshlq_n_u16(u1_incr_mask_16x8_r0_0, 1));
                    u1_incr_mask_8x16_r1_0 =
                        vreinterpretq_u8_u16(vshlq_n_u16(u1_incr_mask_16x8_r1_0, 1));
                    u1_incr_not_8x16_r0_0 =
                        vreinterpretq_u8_u16(vshlq_n_u16(u1_incr_not_16x8_r0_0, 1));
                    u1_incr_not_8x16_r1_0 =
                        vreinterpretq_u8_u16(vshlq_n_u16(u1_incr_not_16x8_r1_0, 1));

                    u1_incr_not_8x16_r0_0 = vaddq_u8(u1_incr_not_8x16_r0_0, x_ref_pos_mask_r0_0);
                    u1_incr_not_8x16_r1_0 = vaddq_u8(u1_incr_not_8x16_r1_0, x_ref_pos_mask_r0_0);

                    x_ref_pos_mask_temp_r0_0 =
                        vaddq_u8(u1_incr_not_8x16_r0_0, u1_incr_mask_8x16_r0_0);
                    x_ref_pos_mask_temp_r1_0 =
                        vaddq_u8(u1_incr_not_8x16_r1_0, u1_incr_mask_8x16_r1_0);

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r0_0));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r0_0));
                    ref_arr_temp0_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8_t0 = vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(u1_incr_not_8x16_r1_0));
                    u1_temp_8x8_t1 = vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(u1_incr_not_8x16_r1_0));
                    ref_arr_temp0_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r0_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r0_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r0_0));
                    ref_arr_temp1_16x8_r0_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));

                    u1_temp_8x8x2_t.val[0] =
                        vreinterpret_u8_s8(vget_low_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8x2_t.val[1] =
                        vreinterpret_u8_s8(vget_high_s8(vreinterpretq_s8_s16(ref_arr_16x8_r1_0)));
                    u1_temp_8x8_t0 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_low_u8(x_ref_pos_mask_temp_r1_0));
                    u1_temp_8x8_t1 =
                        vtbl2_u8(u1_temp_8x8x2_t, vget_high_u8(x_ref_pos_mask_temp_r1_0));
                    ref_arr_temp1_16x8_r1_0 = vreinterpretq_s16_s8(vcombine_s8(
                        vreinterpret_s8_u8(u1_temp_8x8_t0), vreinterpret_s8_u8(u1_temp_8x8_t1)));
                    res_16x8_r0_0 = vmulq_s16(ref_arr_temp0_16x8_r0_0, phs_mask_16min_16x8_0);
                    res_16x8_r1_0 = vmulq_s16(ref_arr_temp0_16x8_r1_0, phs_mask_16min_16x8_0);
                    res_16x8_r0_0 = vmlaq_s16(res_16x8_r0_0, ref_arr_temp1_16x8_r0_0,
                                              vreinterpretq_s16_u16(phs_mask_16x8_0));
                    res_16x8_r1_0 = vmlaq_s16(res_16x8_r1_0, ref_arr_temp1_16x8_r1_0,
                                              vreinterpretq_s16_u16(phs_mask_16x8_0));

                    pu1_ref_y_ptr_incr_temp =
                        pu1_ref_y_ptr_incr + (pi1_y_ref_pos[i4_y] * i4_refarray_wd);
                    u1_y_incr_16x8_r0_0 = vmovl_u8(vld1_u8((pu1_ref_y_ptr_incr_temp)));
                    u1_incr_8x8x2_t.val[0] = vget_low_u8(vreinterpretq_u8_u16(u1_y_incr_16x8_r0_0));
                    u1_incr_8x8x2_t.val[1] =
                        vget_high_u8(vreinterpretq_u8_u16(u1_y_incr_16x8_r0_0));
                    u1_incr_8x8_t0 = vtbl2_u8(u1_incr_8x8x2_t, vget_low_u8(x_ref_rnd_mask_r0_0));
                    u1_incr_8x8_t1 = vtbl2_u8(u1_incr_8x8x2_t, vget_high_u8(x_ref_rnd_mask_r0_0));
                    u1_y_incr_8x16_r0_0 = vcombine_u8(u1_incr_8x8_t0, u1_incr_8x8_t1);

                    u1_y_incr_16x8_r0_0 = vreinterpretq_u16_u8(u1_y_incr_8x16_r0_0);

                    u1_y_incr_16x8_r0_0 = vtstq_u16(u1_y_incr_16x8_r0_0, ones);

                    prev_res_16x8_r0_0 = res_16x8_r0_0;
                    prev_res_16x8_r1_0 = res_16x8_r1_0;

                    u1_prev_y_incr_16x8_r0_0 = u1_y_incr_16x8_r0_0;
                }
            }

            if(!zero_r0_r1)
            {
                res_32x4_l_0 = vdupq_n_s32(0);
                res_32x4_l_1 = vdupq_n_s32(0);
            }
            else
            {
                i4_y_phase = pi1_y_phase[i4_y];

                if((i4_y_phase) >> 3)
                {
                    vert_res_16x8_r0_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r0_0, res_16x8_r1_0);
                    vert_res_16x8_r1_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r1_0, res_16x8_r1_0);
                }
                else
                {
                    vert_res_16x8_r0_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r0_0, res_16x8_r0_0);
                    vert_res_16x8_r1_0 =
                        vbslq_s16(u1_y_incr_16x8_r0_0, res_16x8_r1_0, res_16x8_r0_0);
                }
                res_32x4_l_0 = vmull_n_s16(vget_low_s16(vert_res_16x8_r0_0), 16 - i4_y_phase);
                res_32x4_l_0 =
                    vmlal_n_s16(res_32x4_l_0, vget_low_s16(vert_res_16x8_r1_0), i4_y_phase);

                res_32x4_l_1 = vmull_n_s16(vget_high_s16(vert_res_16x8_r0_0), 16 - i4_y_phase);
                res_32x4_l_1 =
                    vmlal_n_s16(res_32x4_l_1, vget_high_s16(vert_res_16x8_r1_0), i4_y_phase);

                res_32x4_l_0 = vrshrq_n_s32(res_32x4_l_0, 8);
                res_32x4_l_1 = vrshrq_n_s32(res_32x4_l_1, 8);
            }
            out_stride_temp = (i4_y * i4_out_stride);

            out_16x8_0 = vld1q_s16(pi2_out + out_stride_temp);
            out_16x8_1 = vld1q_s16(pi2_out + out_stride_temp + 8);
            out_16x8_0 =
                vbslq_s16(chroma_mask_16x8, vreinterpretq_s16_s32(res_32x4_l_0), out_16x8_0);
            out_16x8_1 =
                vbslq_s16(chroma_mask_16x8, vreinterpretq_s16_s32(res_32x4_l_1), out_16x8_1);
            vst1q_s16((pi2_out + out_stride_temp), out_16x8_0);
            vst1q_s16((pi2_out + out_stride_temp + 8), out_16x8_1);
        }
    }
    return;
} /* End of Interpolation Function */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_reflayer_const_non_boundary_mb_neonintr    */
/*                                                                           */
/*  Description   :                                                          */
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
/*         25 11 2021   Dolan           creation                             */
/*                                                                           */
/*****************************************************************************/

void isvcd_residual_reflayer_const_non_boundary_mb_neonintr(
    WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride, WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
    WORD32 i4_refarray_ht, WORD32 i4_ref_mb_type_q0, WORD32 i4_ref_mb_type_q1,
    WORD32 i4_ref_mb_type_q2, WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
    WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag)
{
    WORD32 i4_y;
    WORD16 *pi2_ref_data_byte;
    WORD16 *pi2_ref_array_temp;

    if(i4_chroma_flag == 0)
    {
        WORD16 index_0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        int16x8_t ref_mb_type_16x8_q0, ref_mb_type_16x8_q1, ref_mb_type_16x8_q2,
            ref_mb_type_16x8_q3, mb_quard1_part_x_16x8;
        int16x8_t ref_mb_type_16x8_0, ref_mb_type_16x8_1;
        int16x8_t ref_mb_type_16x8_low_0, ref_mb_type_16x8_low_1;
        uint16x8_t mb_type_mask_16x8_0, mb_type_mask_16x8_1;
        uint16x8_t mb_type_mask_16x8_low_0, mb_type_mask_16x8_low_1;
        uint16x8_t mask_16x8_0;

        int16x8_t index_arr_0;
        int16x8_t inp_data_16x8_0, inp_data_16x8_1;
        int16x8_t res_16x8_0, res_16x8_1;
        int16x8_t one_16x8 = vdupq_n_s16(1);
        int16x8_t zero_16x8 = vdupq_n_s16(0);

        index_arr_0 = vld1q_s16(&index_0[0]);

        ref_mb_type_16x8_q0 = vdupq_n_s16(i4_ref_mb_type_q0);
        ref_mb_type_16x8_q1 = vdupq_n_s16(i4_ref_mb_type_q1);
        ref_mb_type_16x8_q2 = vdupq_n_s16(i4_ref_mb_type_q2);
        ref_mb_type_16x8_q3 = vdupq_n_s16(i4_ref_mb_type_q3);
        if((i4_mb_quard1_part_x >= i4_refarray_wd) && (i4_mb_quard1_part_y >= i4_refarray_ht))
        {
            // Quard 0
            ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
            ref_mb_type_16x8_1 = ref_mb_type_16x8_q0;
            mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
            mb_type_mask_16x8_1 = mb_type_mask_16x8_0;
        }
        else if((i4_mb_quard1_part_y >= (i4_refarray_ht - 1)) &&
                (i4_mb_quard1_part_x < i4_refarray_wd))
        {
            // Quard 0 & 1
            if(i4_mb_quard1_part_x == 8)
            {
                ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                ref_mb_type_16x8_1 = ref_mb_type_16x8_q1;
            }
            else if(i4_mb_quard1_part_x < 8)
            {
                mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x));
                mask_16x8_0 =
                    vcltq_s16(index_arr_0, mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

                ref_mb_type_16x8_0 =
                    vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);
                ref_mb_type_16x8_1 = ref_mb_type_16x8_q1;
            }
            else
            {
                mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x - 8));
                mask_16x8_0 =
                    vcleq_s16(index_arr_0, mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

                ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                ref_mb_type_16x8_1 =
                    vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);
            }

            mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
            mb_type_mask_16x8_1 = vceqq_s16(ref_mb_type_16x8_1, one_16x8);
        }
        else
        {
            if(i4_mb_quard1_part_x >= i4_refarray_wd)
            {
                ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                ref_mb_type_16x8_1 = ref_mb_type_16x8_q0;

                ref_mb_type_16x8_low_0 = ref_mb_type_16x8_q2;
                ref_mb_type_16x8_low_1 = ref_mb_type_16x8_q2;
            }
            else
            {
                // Quard 0, 1, 2, 3
                if(i4_mb_quard1_part_x == 8)
                {
                    ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                    ref_mb_type_16x8_1 = ref_mb_type_16x8_q1;

                    ref_mb_type_16x8_low_0 = ref_mb_type_16x8_q2;
                    ref_mb_type_16x8_low_1 = ref_mb_type_16x8_q3;
                }
                else if(i4_mb_quard1_part_x < 8)
                {
                    mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x));
                    mask_16x8_0 =
                        vcltq_s16(index_arr_0, mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

                    ref_mb_type_16x8_0 =
                        vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);
                    ref_mb_type_16x8_1 = ref_mb_type_16x8_q1;

                    ref_mb_type_16x8_low_0 =
                        vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q2, ref_mb_type_16x8_q3);
                    ref_mb_type_16x8_low_1 = ref_mb_type_16x8_q3;
                }
                else
                {
                    mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x - 8));
                    mask_16x8_0 =
                        vcltq_s16(index_arr_0, mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

                    ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                    ref_mb_type_16x8_1 =
                        vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);

                    ref_mb_type_16x8_low_0 = ref_mb_type_16x8_q2;
                    ref_mb_type_16x8_low_1 =
                        vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q2, ref_mb_type_16x8_q3);
                }
                mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
                mb_type_mask_16x8_1 = vceqq_s16(ref_mb_type_16x8_1, one_16x8);

                mb_type_mask_16x8_low_0 = vceqq_s16(ref_mb_type_16x8_low_0, one_16x8);
                mb_type_mask_16x8_low_1 = vceqq_s16(ref_mb_type_16x8_low_1, one_16x8);
            }
        }

        if(i4_mb_quard1_part_y < i4_refarray_ht - 1)
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = vld1q_s16((pi2_ref_data_byte));
                inp_data_16x8_1 = vld1q_s16((pi2_ref_data_byte + 8));
                if(i4_y < i4_mb_quard1_part_y)
                {
                    res_16x8_0 = vbslq_s16(mb_type_mask_16x8_0, inp_data_16x8_0, zero_16x8);
                    res_16x8_1 = vbslq_s16(mb_type_mask_16x8_1, inp_data_16x8_1, zero_16x8);
                }
                else
                {
                    res_16x8_0 = vbslq_s16(mb_type_mask_16x8_low_0, inp_data_16x8_0, zero_16x8);
                    res_16x8_1 = vbslq_s16(mb_type_mask_16x8_low_1, inp_data_16x8_1, zero_16x8);
                }
                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                vst1q_s16((pi2_ref_array_temp), res_16x8_0);
                vst1q_s16((pi2_ref_array_temp + 8), res_16x8_1);
            }
        }
        else
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8_0 = vld1q_s16((pi2_ref_data_byte));
                inp_data_16x8_1 = vld1q_s16((pi2_ref_data_byte + 8));

                res_16x8_0 = vbslq_s16(mb_type_mask_16x8_0, inp_data_16x8_0, zero_16x8);
                res_16x8_1 = vbslq_s16(mb_type_mask_16x8_1, inp_data_16x8_1, zero_16x8);

                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                vst1q_s16((pi2_ref_array_temp), res_16x8_0);
                vst1q_s16((pi2_ref_array_temp + 8), res_16x8_1);
            }
        }
    }
    else
    {
        WORD16 index_0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        int16x8_t ref_mb_type_16x8_q0, ref_mb_type_16x8_q1, ref_mb_type_16x8_q2,
            ref_mb_type_16x8_q3, mb_quard1_part_x_16x8;
        int16x8_t ref_mb_type_16x8_0;
        int16x8_t ref_mb_type_16x8_low_0;
        uint16x8_t mb_type_mask_16x8_0;
        uint16x8_t mb_type_mask_16x8_low_0;
        uint16x8_t mask_16x8_0;

        int16x8_t index_arr_0;
        int16x8x2_t inp_data_16x8x2;
        int16x8_t inp_data_16x8;
        int16x8_t res_16x8_0;
        int16x8_t one_16x8 = vdupq_n_s16(1);
        int16x8_t zero_16x8 = vdupq_n_s16(0);
        index_arr_0 = vld1q_s16(&index_0[0]);

        ref_mb_type_16x8_q0 = vdupq_n_s16(i4_ref_mb_type_q0);
        ref_mb_type_16x8_q1 = vdupq_n_s16(i4_ref_mb_type_q1);
        ref_mb_type_16x8_q2 = vdupq_n_s16(i4_ref_mb_type_q2);
        ref_mb_type_16x8_q3 = vdupq_n_s16(i4_ref_mb_type_q3);
        if((i4_mb_quard1_part_x >= i4_refarray_wd) && (i4_mb_quard1_part_y >= i4_refarray_ht))
        {
            // Quard 0
            ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
            mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
        }
        else if((i4_mb_quard1_part_y >= (i4_refarray_ht - 1)) &&
                (i4_mb_quard1_part_x < i4_refarray_wd))
        {
            // Quard 0 & 1
            mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x));
            mask_16x8_0 = vcltq_s16(index_arr_0,
                                    mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

            ref_mb_type_16x8_0 = vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);
            mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
        }
        else
        {
            if(i4_mb_quard1_part_x >= i4_refarray_wd)
            {
                ref_mb_type_16x8_0 = ref_mb_type_16x8_q0;
                ref_mb_type_16x8_low_0 = ref_mb_type_16x8_q2;
            }
            else
            {
                mb_quard1_part_x_16x8 = vdupq_n_s16((i4_mb_quard1_part_x));
                mask_16x8_0 =
                    vcltq_s16(index_arr_0, mb_quard1_part_x_16x8);  // return 1 if a<b, else 0

                ref_mb_type_16x8_0 =
                    vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q0, ref_mb_type_16x8_q1);
                ref_mb_type_16x8_low_0 =
                    vbslq_s16(mask_16x8_0, ref_mb_type_16x8_q2, ref_mb_type_16x8_q3);

                mb_type_mask_16x8_0 = vceqq_s16(ref_mb_type_16x8_0, one_16x8);
                mb_type_mask_16x8_low_0 = vceqq_s16(ref_mb_type_16x8_low_0, one_16x8);
            }
        }

        if(i4_mb_quard1_part_y < i4_refarray_ht - 1)
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8x2 = vld2q_s16((pi2_ref_data_byte));
                inp_data_16x8 = inp_data_16x8x2.val[0];

                if(i4_y < i4_mb_quard1_part_y)
                {
                    res_16x8_0 = vbslq_s16(mb_type_mask_16x8_0, inp_data_16x8, zero_16x8);
                }
                else
                {
                    res_16x8_0 = vbslq_s16(mb_type_mask_16x8_low_0, inp_data_16x8, zero_16x8);
                }
                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                vst1q_s16((pi2_ref_array_temp), res_16x8_0);
            }
        }
        else
        {
            for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_y * i4_inp_data_stride);
                inp_data_16x8x2 = vld2q_s16((pi2_ref_data_byte));
                inp_data_16x8 = inp_data_16x8x2.val[0];

                res_16x8_0 = vbslq_s16(mb_type_mask_16x8_0, inp_data_16x8, zero_16x8);

                pi2_ref_array_temp = pi2_ref_array + (i4_y * i4_refarray_wd);
                vst1q_s16((pi2_ref_array_temp), res_16x8_0);
            }
        }
    }
}
