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
******************************************************************************
* @file ih264e_downscaler_neon.c
*
* @brief
*  This file contains the ARMV8 SIMD version of the function which does
*  horizontal scaling and transpose
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - ih264e_horizontal_downscale_and_transpose_av8()
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <stdlib.h>
#include <arm_neon.h>

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvc_defs.h"
#include "isvce_defs.h"
#include "isvc_structs.h"
#include "isvce_downscaler_private_defs.h"

void isvce_horizontal_downscale_and_transpose_neon(
    downscaler_ctxt_t *ps_scaler, buffer_container_t *ps_src, buffer_container_t *ps_dst,
    FILTER_COEFF_ARRAY pai1_filters, UWORD32 u4_blk_wd, UWORD32 u4_blk_ht, UWORD8 u1_is_chroma)
{
    WORD32 i, j;
    UWORD8 u1_phase;
    UWORD8 *pu1_src_j, *pu1_dst_j;
    UWORD8 *pu1_in_pixel;
    UWORD8 *pu1_out_pixel;
    WORD8 *pi1_filter_grid;
    UWORD16 u2_full_pixel_inc;
    UWORD32 u4_num_iterations_vertical_by_16, u4_num_iterations_vertical_by_8;
    UWORD32 u4_rem_vert_loop_by_8, u4_rem_vert_loop_by_4;
    UWORD32 u4_rem_vert_loop;
    UWORD32 u4_height_finished;

    uint8x8_t reg_8x8_src_r0, reg_8x8_src_r1, reg_8x8_src_r2, reg_8x8_src_r3, reg_8x8_src_r4,
        reg_8x8_src_r5, reg_8x8_src_r6, reg_8x8_src_r7;

    uint16x8_t reg_16x8_src_r0, reg_16x8_src_r1, reg_16x8_src_r2, reg_16x8_src_r3, reg_16x8_src_r4,
        reg_16x8_src_r5, reg_16x8_src_r6, reg_16x8_src_r7;

    int16x8_t reg_16x8_mul_r0, reg_16x8_mul_r1, reg_16x8_mul_r2, reg_16x8_mul_r3, reg_16x8_mul_r4,
        reg_16x8_mul_r5, reg_16x8_mul_r6, reg_16x8_mul_r7;

    int32x4_t reg_32x4_sum_r0, reg_32x4_sum_r1, reg_32x4_sum_r2, reg_32x4_sum_r3, reg_32x4_sum_r4,
        reg_32x4_sum_r5, reg_32x4_sum_r6, reg_32x4_sum_r7;

    int32x4_t reg_32x4_sum_r01, reg_32x4_sum_r23, reg_32x4_sum_r45, reg_32x4_sum_r67,
        reg_32x4_sum_r89, reg_32x4_sum_r1011, reg_32x4_sum_r1213, reg_32x4_sum_r1415;

    uint8x8_t reg_8x8_src_r8, reg_8x8_src_r9, reg_8x8_src_r10, reg_8x8_src_r11, reg_8x8_src_r12,
        reg_8x8_src_r13, reg_8x8_src_r14, reg_8x8_src_r15;

    uint16x8_t reg_16x8_src_r8, reg_16x8_src_r9, reg_16x8_src_r10, reg_16x8_src_r11,
        reg_16x8_src_r12, reg_16x8_src_r13, reg_16x8_src_r14, reg_16x8_src_r15;

    int16x8_t reg_16x8_mul_r8, reg_16x8_mul_r9, reg_16x8_mul_r10, reg_16x8_mul_r11,
        reg_16x8_mul_r12, reg_16x8_mul_r13, reg_16x8_mul_r14, reg_16x8_mul_r15;

    int32x4_t reg_32x4_sum_r8, reg_32x4_sum_r9, reg_32x4_sum_r10, reg_32x4_sum_r11,
        reg_32x4_sum_r12, reg_32x4_sum_r13, reg_32x4_sum_r14, reg_32x4_sum_r15;

    uint8x16_t reg_8x16_src_r0, reg_8x16_src_r1, reg_8x16_src_r2, reg_8x16_src_r3, reg_8x16_src_r4,
        reg_8x16_src_r5, reg_8x16_src_r6, reg_8x16_src_r7;

    uint16x8_t reg_16x8_src_cb_r0, reg_16x8_src_cb_r1, reg_16x8_src_cb_r2, reg_16x8_src_cb_r3,
        reg_16x8_src_cb_r4, reg_16x8_src_cb_r5, reg_16x8_src_cb_r6, reg_16x8_src_cb_r7;

    uint16x8_t reg_16x8_src_cr_r0, reg_16x8_src_cr_r1, reg_16x8_src_cr_r2, reg_16x8_src_cr_r3,
        reg_16x8_src_cr_r4, reg_16x8_src_cr_r5, reg_16x8_src_cr_r6, reg_16x8_src_cr_r7;

    int16x8_t reg_16x8_mul_cb_r0, reg_16x8_mul_cb_r1, reg_16x8_mul_cb_r2, reg_16x8_mul_cb_r3,
        reg_16x8_mul_cb_r4, reg_16x8_mul_cb_r5, reg_16x8_mul_cb_r6, reg_16x8_mul_cb_r7;

    int16x8_t reg_16x8_mul_cr_r0, reg_16x8_mul_cr_r1, reg_16x8_mul_cr_r2, reg_16x8_mul_cr_r3,
        reg_16x8_mul_cr_r4, reg_16x8_mul_cr_r5, reg_16x8_mul_cr_r6, reg_16x8_mul_cr_r7;

    int32x4_t reg_32x4_sum_cb_r0, reg_32x4_sum_cb_r1, reg_32x4_sum_cb_r2, reg_32x4_sum_cb_r3,
        reg_32x4_sum_cb_r4, reg_32x4_sum_cb_r5, reg_32x4_sum_cb_r6, reg_32x4_sum_cb_r7;

    int32x4_t reg_32x4_sum_cr_r0, reg_32x4_sum_cr_r1, reg_32x4_sum_cr_r2, reg_32x4_sum_cr_r3,
        reg_32x4_sum_cr_r4, reg_32x4_sum_cr_r5, reg_32x4_sum_cr_r6, reg_32x4_sum_cr_r7;

    int32x4_t reg_32x4_sum_cb_r01, reg_32x4_sum_cb_r23, reg_32x4_sum_cb_r45, reg_32x4_sum_cb_r67;
    uint16x4_t reg_16x4_sum_cb_r01_23, reg_16x4_sum_cb_r45_67;
    uint16x8_t reg_16x8_sum_cb_r0_r7;
    uint8x8_t reg_8x8_sum_cb_r0_r7;

    int32x4_t reg_32x4_sum_cr_r01, reg_32x4_sum_cr_r23, reg_32x4_sum_cr_r45, reg_32x4_sum_cr_r67;
    uint16x4_t reg_16x4_sum_cr_r01_23, reg_16x4_sum_cr_r45_67;
    uint16x8_t reg_16x8_sum_cr_r0_r7;
    uint8x8_t reg_8x8_sum_cr_r0_r7;
    uint16x8_t reg_16x8_sum_cb_cr_r0_r3;
    uint8x8_t reg_8x8_sum_cb_cr_r0_r3;

    int32x4_t reg_32x4_sum_cb_cr_r0;
    uint16x4_t reg_16x4_sum_cb_cr_r0;

    int32x4_t reg_32x4_zero = vdupq_n_s32(0);

    uint16x4_t reg_16x4_sum_r01_23, reg_16x4_sum_r45_67;
    uint16x4_t reg_16x4_sum_r8_r11, reg_16x4_sum_r12_r15;
    uint16x8_t reg_16x8_sum_r0_r7, reg_16x8_sum_r8_r15;
    uint8x8_t reg_8x8_sum_r0_r7, reg_8x8_sum_r8_r15;
    uint8x16_t reg_8x16_sum_r0_r15;
    int8x8_t reg_8x8_filt_coeff_grid;
    int16x8_t reg_16x8_filt_coeff_grid;
    int32x4x2_t reg_32x4x2_sum_r01, reg_32x4x2_sum_r23, reg_32x4x2_sum_r45, reg_32x4x2_sum_r67;
    int32x4x2_t reg_32x4x2_sum_r89, reg_32x4x2_sum_r1011, reg_32x4x2_sum_r1213,
        reg_32x4x2_sum_r1415;
    uint8x16x2_t reg_8x16x2_src_r0, reg_8x16x2_src_r1, reg_8x16x2_src_r2, reg_8x16x2_src_r3;

    downscaler_state_t *ps_scaler_state = (downscaler_state_t *) ps_scaler->pv_scaler_state;

    UWORD32 u4_center_pixel_pos = ps_scaler_state->i4_init_offset;
    UWORD32 u4_src_vert_increments = ps_scaler_state->u4_vert_increment;
    UWORD32 u4_src_horz_increments = ps_scaler_state->u4_horz_increment;
    UWORD8 *pu1_src = (UWORD8 *) ps_src->pv_data;
    UWORD32 u4_in_stride = ps_src->i4_data_stride;
    UWORD8 *pu1_dst = (UWORD8 *) ps_dst->pv_data;
    UWORD32 u4_out_stride = ps_dst->i4_data_stride;
    UWORD32 u4_center_pixel_pos_src = u4_center_pixel_pos;

    /* Offset the input so that the input pixel to be processed
    co-incides with the centre of filter (4th coefficient)*/
    pu1_src += (1 + u1_is_chroma);

    ASSERT((1 << DOWNSCALER_Q) == u4_src_vert_increments);

    if(!u1_is_chroma)
    {
        u4_num_iterations_vertical_by_16 = u4_blk_ht >> 4;
        u4_rem_vert_loop = u4_blk_ht % 16;

        for(j = 0; j < (WORD32) u4_num_iterations_vertical_by_16; j++)
        {
            pu1_src_j = pu1_src + ((j << 4) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 4);

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);

                pi1_filter_grid = pai1_filters[u1_phase];

                /* Doing the Calculation for current Loop Count  */
                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                /******************************************************/
                /* This loop is going vertically in bottom direction */
                /* but the output pixels are stored in horizontal    */
                /* direction in transpose manner                     */
                /******************************************************/

                /* r0-r7 */
                reg_8x8_src_r0 = vld1_u8(pu1_in_pixel);
                reg_8x8_src_r1 = vld1_u8(pu1_in_pixel + u4_in_stride);
                reg_8x8_src_r2 = vld1_u8(pu1_in_pixel + 2 * u4_in_stride);
                reg_8x8_src_r3 = vld1_u8(pu1_in_pixel + 3 * u4_in_stride);
                reg_8x8_src_r4 = vld1_u8(pu1_in_pixel + 4 * u4_in_stride);
                reg_8x8_src_r5 = vld1_u8(pu1_in_pixel + 5 * u4_in_stride);
                reg_8x8_src_r6 = vld1_u8(pu1_in_pixel + 6 * u4_in_stride);
                reg_8x8_src_r7 = vld1_u8(pu1_in_pixel + 7 * u4_in_stride);

                /* r0-r7 */
                reg_16x8_src_r0 = vmovl_u8(reg_8x8_src_r0);
                reg_16x8_src_r1 = vmovl_u8(reg_8x8_src_r1);
                reg_16x8_src_r2 = vmovl_u8(reg_8x8_src_r2);
                reg_16x8_src_r3 = vmovl_u8(reg_8x8_src_r3);
                reg_16x8_src_r4 = vmovl_u8(reg_8x8_src_r4);
                reg_16x8_src_r5 = vmovl_u8(reg_8x8_src_r5);
                reg_16x8_src_r6 = vmovl_u8(reg_8x8_src_r6);
                reg_16x8_src_r7 = vmovl_u8(reg_8x8_src_r7);

                /* r8-r15 */
                reg_8x8_src_r8 = vld1_u8(pu1_in_pixel + 8 * u4_in_stride);
                reg_8x8_src_r9 = vld1_u8(pu1_in_pixel + 9 * u4_in_stride);
                reg_8x8_src_r10 = vld1_u8(pu1_in_pixel + 10 * u4_in_stride);
                reg_8x8_src_r11 = vld1_u8(pu1_in_pixel + 11 * u4_in_stride);
                reg_8x8_src_r12 = vld1_u8(pu1_in_pixel + 12 * u4_in_stride);
                reg_8x8_src_r13 = vld1_u8(pu1_in_pixel + 13 * u4_in_stride);
                reg_8x8_src_r14 = vld1_u8(pu1_in_pixel + 14 * u4_in_stride);
                reg_8x8_src_r15 = vld1_u8(pu1_in_pixel + 15 * u4_in_stride);

                reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                /*r0-r7 */
                reg_16x8_mul_r0 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r0), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r1 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r1), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r2 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r2), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r3 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r3), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r4 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r4), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r5 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r5), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r6 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r6), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r7 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r7), reg_16x8_filt_coeff_grid);

                /* r8-r15 */
                reg_16x8_src_r8 = vmovl_u8(reg_8x8_src_r8);
                reg_16x8_src_r9 = vmovl_u8(reg_8x8_src_r9);
                reg_16x8_src_r10 = vmovl_u8(reg_8x8_src_r10);
                reg_16x8_src_r11 = vmovl_u8(reg_8x8_src_r11);
                reg_16x8_src_r12 = vmovl_u8(reg_8x8_src_r12);
                reg_16x8_src_r13 = vmovl_u8(reg_8x8_src_r13);
                reg_16x8_src_r14 = vmovl_u8(reg_8x8_src_r14);
                reg_16x8_src_r15 = vmovl_u8(reg_8x8_src_r15);

                /* r0-r7 */
                reg_32x4_sum_r0 = vpaddlq_s16(reg_16x8_mul_r0);
                reg_32x4_sum_r1 = vpaddlq_s16(reg_16x8_mul_r1);
                reg_32x4_sum_r2 = vpaddlq_s16(reg_16x8_mul_r2);
                reg_32x4_sum_r3 = vpaddlq_s16(reg_16x8_mul_r3);
                reg_32x4_sum_r4 = vpaddlq_s16(reg_16x8_mul_r4);
                reg_32x4_sum_r5 = vpaddlq_s16(reg_16x8_mul_r5);
                reg_32x4_sum_r6 = vpaddlq_s16(reg_16x8_mul_r6);
                reg_32x4_sum_r7 = vpaddlq_s16(reg_16x8_mul_r7);

                /* r8-r15 */
                reg_16x8_mul_r8 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r8), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r9 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r9), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r10 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r10), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r11 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r11), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r12 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r12), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r13 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r13), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r14 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r14), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_r15 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r15), reg_16x8_filt_coeff_grid);

                /* r0-r7 */
                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r0, reg_32x4_sum_r1);
                reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_r2, reg_32x4_sum_r3);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_r4, reg_32x4_sum_r5);
                reg_32x4x2_sum_r67 = vuzpq_s32(reg_32x4_sum_r6, reg_32x4_sum_r7);

                reg_32x4_sum_r01 = vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_r23 = vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                reg_32x4_sum_r45 = vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);
                reg_32x4_sum_r67 = vaddq_s32(reg_32x4x2_sum_r67.val[0], reg_32x4x2_sum_r67.val[1]);

                /* r8-r15 */
                reg_32x4_sum_r8 = vpaddlq_s16(reg_16x8_mul_r8);
                reg_32x4_sum_r9 = vpaddlq_s16(reg_16x8_mul_r9);
                reg_32x4_sum_r10 = vpaddlq_s16(reg_16x8_mul_r10);
                reg_32x4_sum_r11 = vpaddlq_s16(reg_16x8_mul_r11);
                reg_32x4_sum_r12 = vpaddlq_s16(reg_16x8_mul_r12);
                reg_32x4_sum_r13 = vpaddlq_s16(reg_16x8_mul_r13);
                reg_32x4_sum_r14 = vpaddlq_s16(reg_16x8_mul_r14);
                reg_32x4_sum_r15 = vpaddlq_s16(reg_16x8_mul_r15);

                /* r0-r7 */
                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r01, reg_32x4_sum_r23);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_r45, reg_32x4_sum_r67);
                reg_32x4_sum_r01 = vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_r45 = vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);

                /* r8-r15 */
                reg_32x4x2_sum_r89 = vuzpq_s32(reg_32x4_sum_r8, reg_32x4_sum_r9);
                reg_32x4x2_sum_r1011 = vuzpq_s32(reg_32x4_sum_r10, reg_32x4_sum_r11);
                reg_32x4x2_sum_r1213 = vuzpq_s32(reg_32x4_sum_r12, reg_32x4_sum_r13);
                reg_32x4x2_sum_r1415 = vuzpq_s32(reg_32x4_sum_r14, reg_32x4_sum_r15);

                reg_32x4_sum_r89 = vaddq_s32(reg_32x4x2_sum_r89.val[0], reg_32x4x2_sum_r89.val[1]);
                reg_32x4_sum_r1011 =
                    vaddq_s32(reg_32x4x2_sum_r1011.val[0], reg_32x4x2_sum_r1011.val[1]);
                reg_32x4_sum_r1213 =
                    vaddq_s32(reg_32x4x2_sum_r1213.val[0], reg_32x4x2_sum_r1213.val[1]);
                reg_32x4_sum_r1415 =
                    vaddq_s32(reg_32x4x2_sum_r1415.val[0], reg_32x4x2_sum_r1415.val[1]);

                /* r0-r7 */
                reg_16x4_sum_r01_23 = vqrshrun_n_s32(reg_32x4_sum_r01, 7);
                reg_16x4_sum_r45_67 = vqrshrun_n_s32(reg_32x4_sum_r45, 7);

                /* r8-r15 */
                reg_32x4x2_sum_r89 = vuzpq_s32(reg_32x4_sum_r89, reg_32x4_sum_r1011);
                reg_32x4x2_sum_r1213 = vuzpq_s32(reg_32x4_sum_r1213, reg_32x4_sum_r1415);
                reg_32x4_sum_r89 = vaddq_s32(reg_32x4x2_sum_r89.val[0], reg_32x4x2_sum_r89.val[1]);
                reg_32x4_sum_r1213 =
                    vaddq_s32(reg_32x4x2_sum_r1213.val[0], reg_32x4x2_sum_r1213.val[1]);

                /* r0-r7 */
                reg_16x8_sum_r0_r7 = vcombine_u16(reg_16x4_sum_r01_23, reg_16x4_sum_r45_67);
                reg_8x8_sum_r0_r7 = vqmovn_u16(reg_16x8_sum_r0_r7);

                reg_16x4_sum_r8_r11 = vqrshrun_n_s32(reg_32x4_sum_r89, 7);
                reg_16x4_sum_r12_r15 = vqrshrun_n_s32(reg_32x4_sum_r1213, 7);

                reg_16x8_sum_r8_r15 = vcombine_u16(reg_16x4_sum_r8_r11, reg_16x4_sum_r12_r15);
                reg_8x8_sum_r8_r15 = vqmovn_u16(reg_16x8_sum_r8_r15);

                reg_8x16_sum_r0_r15 = vcombine_u8(reg_8x8_sum_r0_r7, reg_8x8_sum_r8_r15);

                /* r0-r7 */
                vst1q_u8(pu1_out_pixel, reg_8x16_sum_r0_r15);

                pu1_out_pixel += 16;
                pu1_in_pixel += (u4_src_vert_increments * (u4_in_stride << 4)) >> DOWNSCALER_Q;

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }

        /* Loop for the remaining height less than 16 */
        if(u4_rem_vert_loop)
        {
            u4_rem_vert_loop_by_8 = u4_rem_vert_loop >> 3;
            u4_rem_vert_loop = u4_rem_vert_loop % 8;

            u4_height_finished = (u4_num_iterations_vertical_by_16 << 4);

            pu1_src_j = pu1_src + ((u4_height_finished) *u4_in_stride);
            pu1_dst_j = pu1_dst + u4_height_finished;

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            /* 8 <= remaining height < 16 */
            if(u4_rem_vert_loop_by_8)
            {
                for(i = 0; i < (WORD32) u4_blk_wd; i++)
                {
                    u1_phase = get_filter_phase(u4_center_pixel_pos);
                    pi1_filter_grid = pai1_filters[u1_phase];

                    u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                    pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                    pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                    reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                    for(j = u4_rem_vert_loop_by_8; j > 0; j--)
                    {
                        /******************************************************/
                        /* This loop is going vertically in bottom direction */
                        /* but the output pixels are stored in horizontal    */
                        /* direction in transpose manner                     */
                        /******************************************************/

                        reg_8x8_src_r0 = vld1_u8(pu1_in_pixel);
                        reg_8x8_src_r1 = vld1_u8(pu1_in_pixel + u4_in_stride);
                        reg_8x8_src_r2 = vld1_u8(pu1_in_pixel + 2 * u4_in_stride);
                        reg_8x8_src_r3 = vld1_u8(pu1_in_pixel + 3 * u4_in_stride);
                        reg_8x8_src_r4 = vld1_u8(pu1_in_pixel + 4 * u4_in_stride);
                        reg_8x8_src_r5 = vld1_u8(pu1_in_pixel + 5 * u4_in_stride);
                        reg_8x8_src_r6 = vld1_u8(pu1_in_pixel + 6 * u4_in_stride);
                        reg_8x8_src_r7 = vld1_u8(pu1_in_pixel + 7 * u4_in_stride);

                        reg_16x8_src_r0 = vmovl_u8(reg_8x8_src_r0);
                        reg_16x8_src_r1 = vmovl_u8(reg_8x8_src_r1);
                        reg_16x8_src_r2 = vmovl_u8(reg_8x8_src_r2);
                        reg_16x8_src_r3 = vmovl_u8(reg_8x8_src_r3);
                        reg_16x8_src_r4 = vmovl_u8(reg_8x8_src_r4);
                        reg_16x8_src_r5 = vmovl_u8(reg_8x8_src_r5);
                        reg_16x8_src_r6 = vmovl_u8(reg_8x8_src_r6);
                        reg_16x8_src_r7 = vmovl_u8(reg_8x8_src_r7);
                        reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                        reg_16x8_mul_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r0),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r1 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r1),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r2 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r2),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r3 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r3),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r4 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r4),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r5 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r5),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r6 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r6),
                                                    reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_r7 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r7),
                                                    reg_16x8_filt_coeff_grid);

                        reg_32x4_sum_r0 = vpaddlq_s16(reg_16x8_mul_r0);
                        reg_32x4_sum_r1 = vpaddlq_s16(reg_16x8_mul_r1);
                        reg_32x4_sum_r2 = vpaddlq_s16(reg_16x8_mul_r2);
                        reg_32x4_sum_r3 = vpaddlq_s16(reg_16x8_mul_r3);
                        reg_32x4_sum_r4 = vpaddlq_s16(reg_16x8_mul_r4);
                        reg_32x4_sum_r5 = vpaddlq_s16(reg_16x8_mul_r5);
                        reg_32x4_sum_r6 = vpaddlq_s16(reg_16x8_mul_r6);
                        reg_32x4_sum_r7 = vpaddlq_s16(reg_16x8_mul_r7);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r0, reg_32x4_sum_r1);
                        reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_r2, reg_32x4_sum_r3);
                        reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_r4, reg_32x4_sum_r5);
                        reg_32x4x2_sum_r67 = vuzpq_s32(reg_32x4_sum_r6, reg_32x4_sum_r7);

                        reg_32x4_sum_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                        reg_32x4_sum_r23 =
                            vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                        reg_32x4_sum_r45 =
                            vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);
                        reg_32x4_sum_r67 =
                            vaddq_s32(reg_32x4x2_sum_r67.val[0], reg_32x4x2_sum_r67.val[1]);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r01, reg_32x4_sum_r23);
                        reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_r45, reg_32x4_sum_r67);
                        reg_32x4_sum_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                        reg_32x4_sum_r45 =
                            vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);

                        reg_16x4_sum_r01_23 = vqrshrun_n_s32(reg_32x4_sum_r01, 7);
                        reg_16x4_sum_r45_67 = vqrshrun_n_s32(reg_32x4_sum_r45, 7);

                        reg_16x8_sum_r0_r7 = vcombine_u16(reg_16x4_sum_r01_23, reg_16x4_sum_r45_67);
                        reg_8x8_sum_r0_r7 = vqmovn_u16(reg_16x8_sum_r0_r7);

                        vst1_u8(pu1_out_pixel, reg_8x8_sum_r0_r7);

                        pu1_out_pixel += 8;
                        pu1_in_pixel +=
                            (u4_src_vert_increments * (u4_in_stride << 3)) >> DOWNSCALER_Q;
                    }
                    /* Update the context for next Loop Count */
                    u4_center_pixel_pos += u4_src_horz_increments;
                }
            }

            /* 1 <= remaining height < 8 */
            if(u4_rem_vert_loop)
            {
                u4_height_finished =
                    ((u4_num_iterations_vertical_by_16 << 4) + (u4_rem_vert_loop_by_8 << 3));
                pu1_src_j = pu1_src + u4_height_finished * u4_in_stride;
                pu1_dst_j = pu1_dst + u4_height_finished;

                u4_center_pixel_pos = u4_center_pixel_pos_src;

                for(i = 0; i < (WORD32) u4_blk_wd; i++)
                {
                    u1_phase = get_filter_phase(u4_center_pixel_pos);
                    pi1_filter_grid = pai1_filters[u1_phase];

                    u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                    pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                    pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                    reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                    for(j = u4_rem_vert_loop; j > 0; j--)
                    {
                        /******************************************************/
                        /* This loop is going vertically in bottom direction */
                        /* but the output pixels are stored in horizontal    */
                        /* direction in transpose manner                     */
                        /******************************************************/

                        reg_8x8_src_r0 = vld1_u8(pu1_in_pixel);
                        reg_16x8_src_r0 = vmovl_u8(reg_8x8_src_r0);

                        reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                        reg_16x8_mul_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_r0),
                                                    reg_16x8_filt_coeff_grid);

                        reg_32x4_sum_r0 = vpaddlq_s16(reg_16x8_mul_r0);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r0, reg_32x4_zero);
                        reg_32x4_sum_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_r01, reg_32x4_zero);
                        reg_32x4_sum_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);

                        reg_16x4_sum_r01_23 = vqrshrun_n_s32(reg_32x4_sum_r01, 7);

                        vst1_lane_u8(pu1_out_pixel, vreinterpret_u8_u16(reg_16x4_sum_r01_23), 0);
                        pu1_out_pixel += 1;
                        pu1_in_pixel += (u4_src_vert_increments * u4_in_stride) >> DOWNSCALER_Q;
                    }
                    /* Update the context for next Loop Count */
                    u4_center_pixel_pos += u4_src_horz_increments;
                }
            }
        }
    }
    /* for chroma */
    else
    {
        u4_num_iterations_vertical_by_8 = u4_blk_ht >> 3;
        u4_rem_vert_loop = u4_blk_ht % 8;

        for(j = 0; j < (WORD32) u4_num_iterations_vertical_by_8; j++)
        {
            pu1_src_j = pu1_src + ((j << 3) * u4_in_stride);
            pu1_dst_j = pu1_dst + (j << 3);

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            for(i = 0; i < (WORD32) u4_blk_wd; i++)
            {
                u1_phase = get_filter_phase(u4_center_pixel_pos);
                pi1_filter_grid = pai1_filters[u1_phase];

                /*Doing the Calculation for current Loop Count  */
                u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                /******************************************************/
                /* This loop is going vertically in bottom direction */
                /* but the output pixels are stored in horizontal    */
                /* direction in transpose manner                     */
                /******************************************************/

                reg_8x16_src_r0 = vld1q_u8(pu1_in_pixel);
                reg_8x16_src_r1 = vld1q_u8(pu1_in_pixel + u4_in_stride);
                reg_8x16_src_r2 = vld1q_u8(pu1_in_pixel + 2 * u4_in_stride);
                reg_8x16_src_r3 = vld1q_u8(pu1_in_pixel + 3 * u4_in_stride);
                reg_8x16_src_r4 = vld1q_u8(pu1_in_pixel + 4 * u4_in_stride);
                reg_8x16_src_r5 = vld1q_u8(pu1_in_pixel + 5 * u4_in_stride);
                reg_8x16_src_r6 = vld1q_u8(pu1_in_pixel + 6 * u4_in_stride);
                reg_8x16_src_r7 = vld1q_u8(pu1_in_pixel + 7 * u4_in_stride);

                reg_8x16x2_src_r0 = vuzpq_u8(reg_8x16_src_r0, reg_8x16_src_r1);
                reg_8x16x2_src_r1 = vuzpq_u8(reg_8x16_src_r2, reg_8x16_src_r3);
                reg_8x16x2_src_r2 = vuzpq_u8(reg_8x16_src_r4, reg_8x16_src_r5);
                reg_8x16x2_src_r3 = vuzpq_u8(reg_8x16_src_r6, reg_8x16_src_r7);

                reg_16x8_src_cb_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[0]));
                reg_16x8_src_cb_r1 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r0.val[0]));
                reg_16x8_src_cb_r2 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r1.val[0]));
                reg_16x8_src_cb_r3 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r1.val[0]));
                reg_16x8_src_cb_r4 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r2.val[0]));
                reg_16x8_src_cb_r5 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r2.val[0]));
                reg_16x8_src_cb_r6 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r3.val[0]));
                reg_16x8_src_cb_r7 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r3.val[0]));

                reg_16x8_src_cr_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[1]));
                reg_16x8_src_cr_r1 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r0.val[1]));
                reg_16x8_src_cr_r2 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r1.val[1]));
                reg_16x8_src_cr_r3 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r1.val[1]));
                reg_16x8_src_cr_r4 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r2.val[1]));
                reg_16x8_src_cr_r5 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r2.val[1]));
                reg_16x8_src_cr_r6 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r3.val[1]));
                reg_16x8_src_cr_r7 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r3.val[1]));

                reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                reg_16x8_mul_cb_r0 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r0), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r1 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r1), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r2 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r2), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r3 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r3), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r4 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r4), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r5 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r5), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r6 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r6), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cb_r7 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r7), reg_16x8_filt_coeff_grid);

                reg_16x8_mul_cr_r0 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r0), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r1 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r1), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r2 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r2), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r3 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r3), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r4 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r4), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r5 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r5), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r6 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r6), reg_16x8_filt_coeff_grid);
                reg_16x8_mul_cr_r7 =
                    vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r7), reg_16x8_filt_coeff_grid);

                reg_32x4_sum_cb_r0 = vpaddlq_s16(reg_16x8_mul_cb_r0);
                reg_32x4_sum_cb_r1 = vpaddlq_s16(reg_16x8_mul_cb_r1);
                reg_32x4_sum_cb_r2 = vpaddlq_s16(reg_16x8_mul_cb_r2);
                reg_32x4_sum_cb_r3 = vpaddlq_s16(reg_16x8_mul_cb_r3);
                reg_32x4_sum_cb_r4 = vpaddlq_s16(reg_16x8_mul_cb_r4);
                reg_32x4_sum_cb_r5 = vpaddlq_s16(reg_16x8_mul_cb_r5);
                reg_32x4_sum_cb_r6 = vpaddlq_s16(reg_16x8_mul_cb_r6);
                reg_32x4_sum_cb_r7 = vpaddlq_s16(reg_16x8_mul_cb_r7);

                reg_32x4_sum_cr_r0 = vpaddlq_s16(reg_16x8_mul_cr_r0);
                reg_32x4_sum_cr_r1 = vpaddlq_s16(reg_16x8_mul_cr_r1);
                reg_32x4_sum_cr_r2 = vpaddlq_s16(reg_16x8_mul_cr_r2);
                reg_32x4_sum_cr_r3 = vpaddlq_s16(reg_16x8_mul_cr_r3);
                reg_32x4_sum_cr_r4 = vpaddlq_s16(reg_16x8_mul_cr_r4);
                reg_32x4_sum_cr_r5 = vpaddlq_s16(reg_16x8_mul_cr_r5);
                reg_32x4_sum_cr_r6 = vpaddlq_s16(reg_16x8_mul_cr_r6);
                reg_32x4_sum_cr_r7 = vpaddlq_s16(reg_16x8_mul_cr_r7);

                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_r0, reg_32x4_sum_cb_r1);
                reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_cb_r2, reg_32x4_sum_cb_r3);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_cb_r4, reg_32x4_sum_cb_r5);
                reg_32x4x2_sum_r67 = vuzpq_s32(reg_32x4_sum_cb_r6, reg_32x4_sum_cb_r7);

                reg_32x4_sum_cb_r01 =
                    vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_cb_r23 =
                    vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                reg_32x4_sum_cb_r45 =
                    vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);
                reg_32x4_sum_cb_r67 =
                    vaddq_s32(reg_32x4x2_sum_r67.val[0], reg_32x4x2_sum_r67.val[1]);

                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_r01, reg_32x4_sum_cb_r23);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_cb_r45, reg_32x4_sum_cb_r67);
                reg_32x4_sum_cb_r01 =
                    vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_cb_r45 =
                    vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);

                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cr_r0, reg_32x4_sum_cr_r1);
                reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_cr_r2, reg_32x4_sum_cr_r3);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_cr_r4, reg_32x4_sum_cr_r5);
                reg_32x4x2_sum_r67 = vuzpq_s32(reg_32x4_sum_cr_r6, reg_32x4_sum_cr_r7);

                reg_32x4_sum_cr_r01 =
                    vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_cr_r23 =
                    vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                reg_32x4_sum_cr_r45 =
                    vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);
                reg_32x4_sum_cr_r67 =
                    vaddq_s32(reg_32x4x2_sum_r67.val[0], reg_32x4x2_sum_r67.val[1]);

                reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cr_r01, reg_32x4_sum_cr_r23);
                reg_32x4x2_sum_r45 = vuzpq_s32(reg_32x4_sum_cr_r45, reg_32x4_sum_cr_r67);
                reg_32x4_sum_cr_r01 =
                    vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                reg_32x4_sum_cr_r45 =
                    vaddq_s32(reg_32x4x2_sum_r45.val[0], reg_32x4x2_sum_r45.val[1]);

                reg_16x4_sum_cb_r01_23 = vqrshrun_n_s32(reg_32x4_sum_cb_r01, 7);
                reg_16x4_sum_cb_r45_67 = vqrshrun_n_s32(reg_32x4_sum_cb_r45, 7);

                reg_16x4_sum_cr_r01_23 = vqrshrun_n_s32(reg_32x4_sum_cr_r01, 7);
                reg_16x4_sum_cr_r45_67 = vqrshrun_n_s32(reg_32x4_sum_cr_r45, 7);

                reg_16x8_sum_cb_r0_r7 =
                    vcombine_u16(reg_16x4_sum_cb_r01_23, reg_16x4_sum_cb_r45_67);
                reg_16x8_sum_cr_r0_r7 =
                    vcombine_u16(reg_16x4_sum_cr_r01_23, reg_16x4_sum_cr_r45_67);

                reg_8x8_sum_cb_r0_r7 = vqmovn_u16(reg_16x8_sum_cb_r0_r7);
                reg_8x8_sum_cr_r0_r7 = vqmovn_u16(reg_16x8_sum_cr_r0_r7);

                vst1_u8(pu1_out_pixel, reg_8x8_sum_cb_r0_r7);
                vst1_u8(pu1_out_pixel + u4_out_stride, reg_8x8_sum_cr_r0_r7);

                pu1_out_pixel += 8;

                pu1_in_pixel += (u4_src_vert_increments * (u4_in_stride << 3)) >> DOWNSCALER_Q;

                /* Update the context for next Loop Count */
                u4_center_pixel_pos += u4_src_horz_increments;
            }
        }

        /* Loop for the remaining height less than 8 */
        if(u4_rem_vert_loop)
        {
            u4_rem_vert_loop_by_4 = u4_rem_vert_loop >> 2;
            u4_rem_vert_loop = u4_rem_vert_loop % 4;
            u4_height_finished = (u4_num_iterations_vertical_by_8 << 3);
            pu1_src_j = pu1_src + ((u4_height_finished) *u4_in_stride);
            pu1_dst_j = pu1_dst + u4_height_finished;

            u4_center_pixel_pos = u4_center_pixel_pos_src;

            /* 4<= remaining height < 8 */
            if(u4_rem_vert_loop_by_4)
            {
                for(i = 0; i < (WORD32) u4_blk_wd; i++)
                {
                    u1_phase = get_filter_phase(u4_center_pixel_pos);
                    pi1_filter_grid = pai1_filters[u1_phase];

                    u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                    pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                    pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                    reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                    for(j = u4_rem_vert_loop_by_4; j > 0; j--)
                    {
                        /******************************************************/
                        /* This loop is going vertically in bottom direction */
                        /* but the output pixels are stored in horizontal    */
                        /* direction in transpose manner                     */
                        /******************************************************/

                        reg_8x16_src_r0 = vld1q_u8(pu1_in_pixel);
                        reg_8x16_src_r1 = vld1q_u8(pu1_in_pixel + u4_in_stride);
                        reg_8x16_src_r2 = vld1q_u8(pu1_in_pixel + 2 * u4_in_stride);
                        reg_8x16_src_r3 = vld1q_u8(pu1_in_pixel + 3 * u4_in_stride);

                        reg_8x16x2_src_r0 = vuzpq_u8(reg_8x16_src_r0, reg_8x16_src_r1);
                        reg_8x16x2_src_r1 = vuzpq_u8(reg_8x16_src_r2, reg_8x16_src_r3);

                        reg_16x8_src_cb_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[0]));
                        reg_16x8_src_cb_r1 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r0.val[0]));
                        reg_16x8_src_cb_r2 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r1.val[0]));
                        reg_16x8_src_cb_r3 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r1.val[0]));

                        reg_16x8_src_cr_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[1]));
                        reg_16x8_src_cr_r1 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r0.val[1]));
                        reg_16x8_src_cr_r2 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r1.val[1]));
                        reg_16x8_src_cr_r3 = vmovl_u8(vget_high_u8(reg_8x16x2_src_r1.val[1]));

                        reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                        reg_16x8_mul_cb_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r0),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cb_r1 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r1),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cb_r2 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r2),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cb_r3 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r3),
                                                       reg_16x8_filt_coeff_grid);

                        reg_16x8_mul_cr_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r0),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cr_r1 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r1),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cr_r2 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r2),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cr_r3 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r3),
                                                       reg_16x8_filt_coeff_grid);

                        reg_32x4_sum_cb_r0 = vpaddlq_s16(reg_16x8_mul_cb_r0);
                        reg_32x4_sum_cb_r1 = vpaddlq_s16(reg_16x8_mul_cb_r1);
                        reg_32x4_sum_cb_r2 = vpaddlq_s16(reg_16x8_mul_cb_r2);
                        reg_32x4_sum_cb_r3 = vpaddlq_s16(reg_16x8_mul_cb_r3);

                        reg_32x4_sum_cr_r0 = vpaddlq_s16(reg_16x8_mul_cr_r0);
                        reg_32x4_sum_cr_r1 = vpaddlq_s16(reg_16x8_mul_cr_r1);
                        reg_32x4_sum_cr_r2 = vpaddlq_s16(reg_16x8_mul_cr_r2);
                        reg_32x4_sum_cr_r3 = vpaddlq_s16(reg_16x8_mul_cr_r3);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_r0, reg_32x4_sum_cb_r1);
                        reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_cb_r2, reg_32x4_sum_cb_r3);
                        reg_32x4_sum_cb_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                        reg_32x4_sum_cb_r23 =
                            vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_r01, reg_32x4_sum_cb_r23);
                        reg_32x4_sum_cb_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cr_r0, reg_32x4_sum_cr_r1);
                        reg_32x4x2_sum_r23 = vuzpq_s32(reg_32x4_sum_cr_r2, reg_32x4_sum_cr_r3);
                        reg_32x4_sum_cr_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);
                        reg_32x4_sum_cr_r23 =
                            vaddq_s32(reg_32x4x2_sum_r23.val[0], reg_32x4x2_sum_r23.val[1]);
                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cr_r01, reg_32x4_sum_cr_r23);
                        reg_32x4_sum_cr_r01 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);

                        reg_16x4_sum_cb_r01_23 = vqrshrun_n_s32(reg_32x4_sum_cb_r01, 7);
                        reg_16x4_sum_cr_r01_23 = vqrshrun_n_s32(reg_32x4_sum_cr_r01, 7);

                        reg_16x8_sum_cb_cr_r0_r3 =
                            vcombine_u16(reg_16x4_sum_cb_r01_23, reg_16x4_sum_cr_r01_23);
                        reg_8x8_sum_cb_cr_r0_r3 = vmovn_u16(reg_16x8_sum_cb_cr_r0_r3);
                        vst1_lane_u32((uint32_t *) (pu1_out_pixel),
                                      vreinterpret_u32_u8(reg_8x8_sum_cb_cr_r0_r3), 0);
                        vst1_lane_u32((uint32_t *) (pu1_out_pixel + u4_out_stride),
                                      vreinterpret_u32_u8(reg_8x8_sum_cb_cr_r0_r3), 1);

                        pu1_out_pixel += 4;

                        pu1_in_pixel +=
                            (u4_src_vert_increments * (u4_in_stride << 2)) >> DOWNSCALER_Q;
                    }
                    /* Update the context for next Loop Count */
                    u4_center_pixel_pos += u4_src_horz_increments;
                }
            }

            /* 1<= remaining height < 4 */
            if(u4_rem_vert_loop)
            {
                u4_height_finished =
                    ((u4_num_iterations_vertical_by_8 << 3) + (u4_rem_vert_loop_by_4 << 2));
                pu1_src_j = pu1_src + u4_height_finished * u4_in_stride;
                pu1_dst_j = pu1_dst + u4_height_finished;

                u4_center_pixel_pos = u4_center_pixel_pos_src;
                for(i = 0; i < (WORD32) u4_blk_wd; i++)
                {
                    u1_phase = get_filter_phase(u4_center_pixel_pos);
                    pi1_filter_grid = pai1_filters[u1_phase];

                    u2_full_pixel_inc = u4_center_pixel_pos >> DOWNSCALER_Q;

                    pu1_in_pixel = pu1_src_j + (u2_full_pixel_inc << u1_is_chroma);

                    pu1_out_pixel = pu1_dst_j + ((i << u1_is_chroma) * u4_out_stride);

                    reg_8x8_filt_coeff_grid = vld1_s8(pi1_filter_grid);

                    for(j = u4_rem_vert_loop; j > 0; j--)
                    {
                        /******************************************************/
                        /* This loop is going vertically in bottom direction */
                        /* but the output pixels are stored in horizontal    */
                        /* direction in transpose manner                     */
                        /******************************************************/

                        reg_8x16_src_r0 = vld1q_u8(pu1_in_pixel);

                        reg_8x16x2_src_r0 = vuzpq_u8(reg_8x16_src_r0, reg_8x16_src_r0);
                        reg_16x8_src_cb_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[0]));
                        reg_16x8_src_cr_r0 = vmovl_u8(vget_low_u8(reg_8x16x2_src_r0.val[1]));

                        reg_16x8_filt_coeff_grid = vmovl_s8(reg_8x8_filt_coeff_grid);

                        reg_16x8_mul_cb_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cb_r0),
                                                       reg_16x8_filt_coeff_grid);
                        reg_16x8_mul_cr_r0 = vmulq_s16(vreinterpretq_s16_u16(reg_16x8_src_cr_r0),
                                                       reg_16x8_filt_coeff_grid);

                        reg_32x4_sum_cb_r0 = vpaddlq_s16(reg_16x8_mul_cb_r0);
                        reg_32x4_sum_cr_r0 = vpaddlq_s16(reg_16x8_mul_cr_r0);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_r0, reg_32x4_sum_cr_r0);
                        reg_32x4_sum_cb_cr_r0 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);

                        reg_32x4x2_sum_r01 = vuzpq_s32(reg_32x4_sum_cb_cr_r0, reg_32x4_zero);
                        reg_32x4_sum_cb_cr_r0 =
                            vaddq_s32(reg_32x4x2_sum_r01.val[0], reg_32x4x2_sum_r01.val[1]);

                        reg_16x4_sum_cb_cr_r0 = vqrshrun_n_s32(reg_32x4_sum_cb_cr_r0, 7);
                        vst1_lane_u8((pu1_out_pixel), vreinterpret_u8_u16(reg_16x4_sum_cb_cr_r0),
                                     0);
                        vst1_lane_u8((pu1_out_pixel + u4_out_stride),
                                     vreinterpret_u8_u16(reg_16x4_sum_cb_cr_r0), 2);

                        pu1_out_pixel += 1;

                        pu1_in_pixel += (u4_src_vert_increments * (u4_in_stride)) >> DOWNSCALER_Q;
                    }

                    /* Update the context for next Loop Count */
                    u4_center_pixel_pos += u4_src_horz_increments;
                }
            }
        }
    }
}
