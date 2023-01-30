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
* @file isvce_svc_rc_utils_neon.c
*
* @brief
*  This file contains the neom SIMD version of the function which computes
*  gradient per pixel value being used in Init Qp
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - isvce_get_gpp_neon()
*
* @remarks
*  None
*
*******************************************************************************
*/

#include <arm_neon.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "isvc_structs.h"
#include "isvce_rc_utils_private_defs.h"

/**
*******************************************************************************
*
* @brief
*   get gpp function
*
* @par Description:
*   computes gradient per pixel value for a given frame
*
* @param[in] ps_input_buf
*  pointer to yuv buffer properties
*
* @returns
*  calculated gpp value
*
* @remarks
*  none
*
*******************************************************************************
*/

DOUBLE isvce_get_gpp_neon(yuv_buf_props_t *ps_input_buf)
{
    UWORD8 *pu1_input_buf;
    UWORD32 i, j, k;
    UWORD32 u4_width, u4_height, i4_input_stride;
    DOUBLE d_gpp_y, d_gpp_u, d_gpp_v, d_gpp;

    uint8x8_t reg_8x8_src_r0, reg_8x8_src_r1, reg_8x8_src_r2, reg_8x8_src_r3, reg_8x8_src_r4,
        reg_8x8_src_r5, reg_8x8_src_r6, reg_8x8_src_r7, reg_8x8_src_r8;
    uint8x8_t reg_8x8_src_right_r0, reg_8x8_src_right_r1, reg_8x8_src_right_r2,
        reg_8x8_src_right_r3, reg_8x8_src_right_r4, reg_8x8_src_right_r5, reg_8x8_src_right_r6,
        reg_8x8_src_right_r7;
    uint16x8_t reg_16x8_abs_diff_y, reg_16x8_abs_diff_uv;
    uint64x2_t reg_64x2_gpp_y, reg_64x2_gpp_uv;

    uint8x8_t reg_8x8_shuffle = {0, 2, 4, 6, 1, 3, 5, 7};
    uint16x8_t reg_16x8_and_mask_y = {0xffff, 0xffff, 0xffff, 0xffff,
                                      0xffff, 0xffff, 0xffff, 0x0000};
    uint16x8_t reg_16x8_and_mask_uv = {0xffff, 0xffff, 0xffff, 0x0000,
                                       0xffff, 0xffff, 0xffff, 0x0000};
    uint32x4_t reg_32x4_abs_diff_hadd_y = vdupq_n_u32(0);
    uint32x4_t reg_32x4_abs_diff_hadd_uv = vdupq_n_u32(0);

    d_gpp_y = 0;
    d_gpp_u = 0;
    d_gpp_v = 0;
    d_gpp = 0;
    pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[0].pv_data;
    i4_input_stride = ps_input_buf->as_component_bufs[0].i4_data_stride;
    u4_width = ps_input_buf->u4_width;
    u4_height = ps_input_buf->u4_height;

    ASSERT((u4_width % 8) == 0);

    /***********************************************************/
    /* For Luma -                                              */
    /* This code block calculates gpp value for luma by adding */
    /* the absolute difference between the current pixel and   */
    /* it's immediate right pixel with the absolute difference */
    /* between the current pixel and it's immediate bottom     */
    /* pixel and accumulating for every pixel in the frame.    */
    /***********************************************************/
    /* -8 in the checks below since right column and bottow row being used for gradients, */
    /* and last row and column are ignored for gradient computation. */
    /* Note that input is not required to be padded */
    for(i = 0; i < u4_height - 8; i += 8)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
            reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
            reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
            reg_8x8_src_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j);
            reg_8x8_src_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j);
            reg_8x8_src_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j);
            reg_8x8_src_r8 = vld1_u8(pu1_input_buf + (i4_input_stride * 8) + j);

            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 1);
            reg_8x8_src_right_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j + 1);
            reg_8x8_src_right_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j + 1);
            reg_8x8_src_right_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j + 1);
            reg_8x8_src_right_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j + 1);
            reg_8x8_src_right_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j + 1);
            reg_8x8_src_right_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j + 1);
            reg_8x8_src_right_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j + 1);

            reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_r2);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_r3);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_r4);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r4, reg_8x8_src_r5);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r5, reg_8x8_src_r6);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r6, reg_8x8_src_r7);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r7, reg_8x8_src_r8);

            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_right_r1);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_right_r2);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_right_r3);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r4, reg_8x8_src_right_r4);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r5, reg_8x8_src_right_r5);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r6, reg_8x8_src_right_r6);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r7, reg_8x8_src_right_r7);

            reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 7   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_y                  */
        /************************************************************/
        ASSERT((u4_width - j) == 8);
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
        reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
        reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
        reg_8x8_src_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j);
        reg_8x8_src_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j);
        reg_8x8_src_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j);
        reg_8x8_src_r8 = vld1_u8(pu1_input_buf + (i4_input_stride * 8) + j);

        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 1);
        reg_8x8_src_right_r1 = vext_u8(reg_8x8_src_r1, reg_8x8_src_r1, 1);
        reg_8x8_src_right_r2 = vext_u8(reg_8x8_src_r2, reg_8x8_src_r2, 1);
        reg_8x8_src_right_r3 = vext_u8(reg_8x8_src_r3, reg_8x8_src_r3, 1);
        reg_8x8_src_right_r4 = vext_u8(reg_8x8_src_r4, reg_8x8_src_r4, 1);
        reg_8x8_src_right_r5 = vext_u8(reg_8x8_src_r5, reg_8x8_src_r5, 1);
        reg_8x8_src_right_r6 = vext_u8(reg_8x8_src_r6, reg_8x8_src_r6, 1);
        reg_8x8_src_right_r7 = vext_u8(reg_8x8_src_r7, reg_8x8_src_r7, 1);

        reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_r2);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_r3);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_r4);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r4, reg_8x8_src_r5);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r5, reg_8x8_src_r6);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r6, reg_8x8_src_r7);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r7, reg_8x8_src_r8);

        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_right_r1);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_right_r2);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_right_r3);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r4, reg_8x8_src_right_r4);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r5, reg_8x8_src_right_r5);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r6, reg_8x8_src_right_r6);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r7, reg_8x8_src_right_r7);

        reg_16x8_abs_diff_y = vandq_u16(reg_16x8_abs_diff_y, reg_16x8_and_mask_y);

        reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);

        pu1_input_buf += (i4_input_stride * 8);
    }

    /* Loop for remaining height less than 8 */
    /*    4 <= remaining_height < 8          */
    for(k = i; k < u4_height - 4; k += 4, i += 4)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
            reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
            reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 1);
            reg_8x8_src_right_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j + 1);
            reg_8x8_src_right_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j + 1);
            reg_8x8_src_right_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j + 1);

            reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_r2);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_r3);
            reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_r4);

            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_right_r1);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_right_r2);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_right_r3);

            reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 7   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_y                  */
        /************************************************************/
        ASSERT((u4_width - j) == 8);
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
        reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
        reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);

        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 1);
        reg_8x8_src_right_r1 = vext_u8(reg_8x8_src_r1, reg_8x8_src_r1, 1);
        reg_8x8_src_right_r2 = vext_u8(reg_8x8_src_r2, reg_8x8_src_r2, 1);
        reg_8x8_src_right_r3 = vext_u8(reg_8x8_src_r3, reg_8x8_src_r3, 1);

        reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_r2);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_r3);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_r4);

        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r1, reg_8x8_src_right_r1);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r2, reg_8x8_src_right_r2);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r3, reg_8x8_src_right_r3);

        reg_16x8_abs_diff_y = vandq_u16(reg_16x8_abs_diff_y, reg_16x8_and_mask_y);

        reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);

        pu1_input_buf += (i4_input_stride * 4);
    }

    /* Loop for remaining height less than 4 */
    /*    0 <= remaining_height < 4          */
    for(k = i; k < u4_height - 1; k++)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 1);

            reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_y =
                vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);

            reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 7   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_y                  */
        /************************************************************/
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 1);

        reg_16x8_abs_diff_y = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_y = vabal_u8(reg_16x8_abs_diff_y, reg_8x8_src_r0, reg_8x8_src_right_r0);

        reg_16x8_abs_diff_y = vandq_u16(reg_16x8_abs_diff_y, reg_16x8_and_mask_y);

        reg_32x4_abs_diff_hadd_y = vpadalq_u16(reg_32x4_abs_diff_hadd_y, reg_16x8_abs_diff_y);

        pu1_input_buf += i4_input_stride;
    }

    /* Pairwise add reg_32x4_abs_diff_hadd_y to get final gpp value */
    reg_64x2_gpp_y = vpaddlq_u32(reg_32x4_abs_diff_hadd_y);
    d_gpp_y = vgetq_lane_u64(reg_64x2_gpp_y, 0);
    d_gpp_y += vgetq_lane_u64(reg_64x2_gpp_y, 1);

    pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[1].pv_data;
    i4_input_stride = ps_input_buf->as_component_bufs[1].i4_data_stride;

    /***************************************************************/
    /* For Chroma -                                                */
    /* This code block first deinterleaves the Cb and Cr values,   */
    /* calculates gpp value for both Cb and Cr separately by       */
    /* adding the absolute difference between the current pixel    */
    /* and it's immediate right pixel with the absolute            */
    /* difference between the current pixel and it's immediate     */
    /* bottom pixel and accumulating for every pixel in the frame. */
    /***************************************************************/
    for(i = 0; i < (u4_height >> 1) - 8; i += 8)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
            reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
            reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
            reg_8x8_src_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j);
            reg_8x8_src_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j);
            reg_8x8_src_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j);
            reg_8x8_src_r8 = vld1_u8(pu1_input_buf + (i4_input_stride * 8) + j);

            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 2);
            reg_8x8_src_right_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j + 2);
            reg_8x8_src_right_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j + 2);
            reg_8x8_src_right_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j + 2);
            reg_8x8_src_right_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j + 2);
            reg_8x8_src_right_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j + 2);
            reg_8x8_src_right_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j + 2);
            reg_8x8_src_right_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j + 2);

            /* separating u and v */
            reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
            reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
            reg_8x8_src_r2 = vtbl1_u8(reg_8x8_src_r2, reg_8x8_shuffle);
            reg_8x8_src_r3 = vtbl1_u8(reg_8x8_src_r3, reg_8x8_shuffle);
            reg_8x8_src_r4 = vtbl1_u8(reg_8x8_src_r4, reg_8x8_shuffle);
            reg_8x8_src_r5 = vtbl1_u8(reg_8x8_src_r5, reg_8x8_shuffle);
            reg_8x8_src_r6 = vtbl1_u8(reg_8x8_src_r6, reg_8x8_shuffle);
            reg_8x8_src_r7 = vtbl1_u8(reg_8x8_src_r7, reg_8x8_shuffle);
            reg_8x8_src_r8 = vtbl1_u8(reg_8x8_src_r8, reg_8x8_shuffle);
            reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);
            reg_8x8_src_right_r1 = vtbl1_u8(reg_8x8_src_right_r1, reg_8x8_shuffle);
            reg_8x8_src_right_r2 = vtbl1_u8(reg_8x8_src_right_r2, reg_8x8_shuffle);
            reg_8x8_src_right_r3 = vtbl1_u8(reg_8x8_src_right_r3, reg_8x8_shuffle);
            reg_8x8_src_right_r4 = vtbl1_u8(reg_8x8_src_right_r4, reg_8x8_shuffle);
            reg_8x8_src_right_r5 = vtbl1_u8(reg_8x8_src_right_r5, reg_8x8_shuffle);
            reg_8x8_src_right_r6 = vtbl1_u8(reg_8x8_src_right_r6, reg_8x8_shuffle);
            reg_8x8_src_right_r7 = vtbl1_u8(reg_8x8_src_right_r7, reg_8x8_shuffle);

            reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_r2);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_r3);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_r4);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r4, reg_8x8_src_r5);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r5, reg_8x8_src_r6);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r6, reg_8x8_src_r7);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r7, reg_8x8_src_r8);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_right_r1);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_right_r2);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_right_r3);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r4, reg_8x8_src_right_r4);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r5, reg_8x8_src_right_r5);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r6, reg_8x8_src_right_r6);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r7, reg_8x8_src_right_r7);

            reg_32x4_abs_diff_hadd_uv =
                vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 6   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_uv                 */
        /************************************************************/
        ASSERT((u4_width - j) == 8);
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
        reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
        reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
        reg_8x8_src_r5 = vld1_u8(pu1_input_buf + (i4_input_stride * 5) + j);
        reg_8x8_src_r6 = vld1_u8(pu1_input_buf + (i4_input_stride * 6) + j);
        reg_8x8_src_r7 = vld1_u8(pu1_input_buf + (i4_input_stride * 7) + j);
        reg_8x8_src_r8 = vld1_u8(pu1_input_buf + (i4_input_stride * 8) + j);
        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 2);
        reg_8x8_src_right_r1 = vext_u8(reg_8x8_src_r1, reg_8x8_src_r1, 2);
        reg_8x8_src_right_r2 = vext_u8(reg_8x8_src_r2, reg_8x8_src_r2, 2);
        reg_8x8_src_right_r3 = vext_u8(reg_8x8_src_r3, reg_8x8_src_r3, 2);
        reg_8x8_src_right_r4 = vext_u8(reg_8x8_src_r4, reg_8x8_src_r4, 2);
        reg_8x8_src_right_r5 = vext_u8(reg_8x8_src_r5, reg_8x8_src_r5, 2);
        reg_8x8_src_right_r6 = vext_u8(reg_8x8_src_r6, reg_8x8_src_r6, 2);
        reg_8x8_src_right_r7 = vext_u8(reg_8x8_src_r7, reg_8x8_src_r7, 2);

        /* separating u and v */
        reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
        reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
        reg_8x8_src_r2 = vtbl1_u8(reg_8x8_src_r2, reg_8x8_shuffle);
        reg_8x8_src_r3 = vtbl1_u8(reg_8x8_src_r3, reg_8x8_shuffle);
        reg_8x8_src_r4 = vtbl1_u8(reg_8x8_src_r4, reg_8x8_shuffle);
        reg_8x8_src_r5 = vtbl1_u8(reg_8x8_src_r5, reg_8x8_shuffle);
        reg_8x8_src_r6 = vtbl1_u8(reg_8x8_src_r6, reg_8x8_shuffle);
        reg_8x8_src_r7 = vtbl1_u8(reg_8x8_src_r7, reg_8x8_shuffle);
        reg_8x8_src_r8 = vtbl1_u8(reg_8x8_src_r8, reg_8x8_shuffle);
        reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);
        reg_8x8_src_right_r1 = vtbl1_u8(reg_8x8_src_right_r1, reg_8x8_shuffle);
        reg_8x8_src_right_r2 = vtbl1_u8(reg_8x8_src_right_r2, reg_8x8_shuffle);
        reg_8x8_src_right_r3 = vtbl1_u8(reg_8x8_src_right_r3, reg_8x8_shuffle);
        reg_8x8_src_right_r4 = vtbl1_u8(reg_8x8_src_right_r4, reg_8x8_shuffle);
        reg_8x8_src_right_r5 = vtbl1_u8(reg_8x8_src_right_r5, reg_8x8_shuffle);
        reg_8x8_src_right_r6 = vtbl1_u8(reg_8x8_src_right_r6, reg_8x8_shuffle);
        reg_8x8_src_right_r7 = vtbl1_u8(reg_8x8_src_right_r7, reg_8x8_shuffle);

        reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_r2);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_r3);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_r4);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r4, reg_8x8_src_r5);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r5, reg_8x8_src_r6);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r6, reg_8x8_src_r7);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r7, reg_8x8_src_r8);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_right_r1);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_right_r2);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_right_r3);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r4, reg_8x8_src_right_r4);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r5, reg_8x8_src_right_r5);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r6, reg_8x8_src_right_r6);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r7, reg_8x8_src_right_r7);

        reg_16x8_abs_diff_uv = vandq_u16(reg_16x8_abs_diff_uv, reg_16x8_and_mask_uv);

        reg_32x4_abs_diff_hadd_uv = vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);

        pu1_input_buf += (i4_input_stride * 8);
    }

    /* Loop for remaining height less than 8 */
    /*    4 <= remaining_height < 8          */
    for(k = i; k < (u4_height >> 1) - 4; k += 4, i += 4)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
            reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
            reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 2);
            reg_8x8_src_right_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j + 2);
            reg_8x8_src_right_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j + 2);
            reg_8x8_src_right_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j + 2);

            /* separating u and v */
            reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
            reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
            reg_8x8_src_r2 = vtbl1_u8(reg_8x8_src_r2, reg_8x8_shuffle);
            reg_8x8_src_r3 = vtbl1_u8(reg_8x8_src_r3, reg_8x8_shuffle);
            reg_8x8_src_r4 = vtbl1_u8(reg_8x8_src_r4, reg_8x8_shuffle);
            reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);
            reg_8x8_src_right_r1 = vtbl1_u8(reg_8x8_src_right_r1, reg_8x8_shuffle);
            reg_8x8_src_right_r2 = vtbl1_u8(reg_8x8_src_right_r2, reg_8x8_shuffle);
            reg_8x8_src_right_r3 = vtbl1_u8(reg_8x8_src_right_r3, reg_8x8_shuffle);

            reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_r2);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_r3);
            reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_r4);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_right_r1);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_right_r2);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_right_r3);

            reg_32x4_abs_diff_hadd_uv =
                vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 6   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_uv                 */
        /************************************************************/
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_r2 = vld1_u8(pu1_input_buf + (i4_input_stride * 2) + j);
        reg_8x8_src_r3 = vld1_u8(pu1_input_buf + (i4_input_stride * 3) + j);
        reg_8x8_src_r4 = vld1_u8(pu1_input_buf + (i4_input_stride * 4) + j);
        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 2);
        reg_8x8_src_right_r1 = vext_u8(reg_8x8_src_r1, reg_8x8_src_r1, 2);
        reg_8x8_src_right_r2 = vext_u8(reg_8x8_src_r2, reg_8x8_src_r2, 2);
        reg_8x8_src_right_r3 = vext_u8(reg_8x8_src_r3, reg_8x8_src_r3, 2);

        /* separating u and v */
        reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
        reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
        reg_8x8_src_r2 = vtbl1_u8(reg_8x8_src_r2, reg_8x8_shuffle);
        reg_8x8_src_r3 = vtbl1_u8(reg_8x8_src_r3, reg_8x8_shuffle);
        reg_8x8_src_r4 = vtbl1_u8(reg_8x8_src_r4, reg_8x8_shuffle);
        reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);
        reg_8x8_src_right_r1 = vtbl1_u8(reg_8x8_src_right_r1, reg_8x8_shuffle);
        reg_8x8_src_right_r2 = vtbl1_u8(reg_8x8_src_right_r2, reg_8x8_shuffle);
        reg_8x8_src_right_r3 = vtbl1_u8(reg_8x8_src_right_r3, reg_8x8_shuffle);

        reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_r2);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_r3);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_r4);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r1, reg_8x8_src_right_r1);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r2, reg_8x8_src_right_r2);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r3, reg_8x8_src_right_r3);

        reg_16x8_abs_diff_uv = vandq_u16(reg_16x8_abs_diff_uv, reg_16x8_and_mask_uv);

        reg_32x4_abs_diff_hadd_uv = vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);

        pu1_input_buf += (i4_input_stride * 4);
    }

    /* Loop for remaining height less than 4 */
    /*    0 <= remaining_height < 4          */
    for(k = i; k < (u4_height >> 1) - 1; k++)
    {
        for(j = 0; j < u4_width - 8; j += 8)
        {
            reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
            reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
            reg_8x8_src_right_r0 = vld1_u8(pu1_input_buf + j + 2);

            /* separating u and v */
            reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
            reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
            reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);

            reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
            reg_16x8_abs_diff_uv =
                vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);

            reg_32x4_abs_diff_hadd_uv =
                vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);
        }

        /************************************************************/
        /* Remaining width -                                        */
        /* Since Last pixel is not getting processed, remaining 6   */
        /* pixels are getting processed separately by performing    */
        /* and operations with reg_16x8_and_mask_uv                 */
        /************************************************************/
        reg_8x8_src_r0 = vld1_u8(pu1_input_buf + j);
        reg_8x8_src_r1 = vld1_u8(pu1_input_buf + i4_input_stride + j);
        reg_8x8_src_right_r0 = vext_u8(reg_8x8_src_r0, reg_8x8_src_r0, 2);

        /* separating u and v */
        reg_8x8_src_r0 = vtbl1_u8(reg_8x8_src_r0, reg_8x8_shuffle);
        reg_8x8_src_r1 = vtbl1_u8(reg_8x8_src_r1, reg_8x8_shuffle);
        reg_8x8_src_right_r0 = vtbl1_u8(reg_8x8_src_right_r0, reg_8x8_shuffle);

        reg_16x8_abs_diff_uv = vabdl_u8(reg_8x8_src_r0, reg_8x8_src_r1);
        reg_16x8_abs_diff_uv = vabal_u8(reg_16x8_abs_diff_uv, reg_8x8_src_r0, reg_8x8_src_right_r0);

        reg_16x8_abs_diff_uv = vandq_u16(reg_16x8_abs_diff_uv, reg_16x8_and_mask_uv);

        reg_32x4_abs_diff_hadd_uv = vpadalq_u16(reg_32x4_abs_diff_hadd_uv, reg_16x8_abs_diff_uv);

        pu1_input_buf += i4_input_stride;
    }

    /* Pairwise add u4_abd_hadd_uv to get final gpp_u and gpp_v value */
    reg_64x2_gpp_uv = vpaddlq_u32(reg_32x4_abs_diff_hadd_uv);
    d_gpp_u = vgetq_lane_u64(reg_64x2_gpp_uv, 0);
    d_gpp_v = vgetq_lane_u64(reg_64x2_gpp_uv, 1);

    d_gpp_y /= (u4_width * u4_height);
    d_gpp_u /= ((u4_width / 2) * (u4_height / 2));
    d_gpp_v /= ((u4_width / 2) * (u4_height / 2));

    d_gpp = (DOUBLE) ((WT_LUMA_GPP * d_gpp_y) + d_gpp_u + d_gpp_v) / WT_TOTAL_GPP;

    return d_gpp;
}
