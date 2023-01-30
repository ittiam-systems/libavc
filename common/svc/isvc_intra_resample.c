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
/*!
 **************************************************************************
 * \file isvcd_resamp_svc.c
 *
 * \brief
 *    Contains routines that resample for SVC resampling
 *
 * Detailed_description
 *
 * \date
 *
 *
 * \author
 **************************************************************************
 */
#include <assert.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "isvc_macros.h"
#include "ih264_platform_macros.h"
#include "isvc_intra_resample.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"

#define NUM_SEGMENTS 16
#define NUM_INTRA_SAMP_FXNS 32
#define INTERPOL_FILTER_SIZE_LUMA 64
#define INTERPOL_FILTER_SIZE_CHROMA 32

typedef void(PF_INTRA_SAMP_PADDING)(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                    UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                    UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                    WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                    WORD32 i4_corner_pixel_available);

static const WORD8 g_ai1_interp_filter_luma[INTERPOL_FILTER_SIZE_LUMA] = {
    0,  -1, -2, -3, -3, -4, -4, -3, -3, -3, -2, -1, -1, -1, -1, -1, 32, 32, 31, 30, 28, 26,
    24, 22, 19, 16, 14, 11, 8,  6,  4,  2,  0,  2,  4,  6,  8,  11, 14, 16, 19, 22, 24, 26,
    28, 30, 31, 32, 0,  -1, -1, -1, -1, -1, -2, -3, -3, -3, -4, -4, -3, -3, -2, -1};

static const UWORD8 g_au1_interp_filter_chroma[INTERPOL_FILTER_SIZE_CHROMA] = {
    32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8,  6,  4,  2,
    0,  2,  4,  6,  8,  10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};

static const UWORD32 gu4_valid_segs_lookup[NUM_SEGMENTS] = {
    0x0F000000, 0xCF000000, 0x3F000000, 0xFF000000, 0x0F000000, 0xCF000000, 0x3F000000, 0xFF000000,
    0x0F000000, 0x8F000000, 0x6F000000, 0xEF000000, 0x1F000000, 0x9F000000, 0x7F000000, 0xFF000000};

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_copy_data                                           */
/*                                                                           */
/*  Description   : this module copies the data from source to destination   */
/*                  the amount of data to be copied is passed as input       */
/*                                                                           */
/*  Inputs        : pu1_src : pointer to the source buffer                   */
/*                  u2_src_stride : source buffer stride                     */
/*                  pu1_dst : pointer to the destination buffer              */
/*                  u2_dst_stride : destination buffer stride                */
/*                  u4_num_bytes : number of bytes to be copied              */
/*                  u4_num_lines : number of lines to be copied              */
/*  Globals       : none                                                     */
/*  Processing    : it does a memcpy from source to destination              */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : both buffers are assumed to be 2-D buffers               */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         29 04 2009   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
/** \brief performs the 2-D memory transfer  */
static void isvc_copy_data(UWORD8 *pu1_src, WORD32 i4_src_stride, UWORD8 *pu1_dst,
                           WORD32 i4_dst_stride, WORD32 i4_num_bytes, WORD32 i4_num_lines)
{
    WORD32 i4_vert_lines;
    ASSERT(NULL != pu1_src);
    ASSERT(NULL != pu1_dst);

    for(i4_vert_lines = 0; i4_vert_lines < i4_num_lines; i4_vert_lines++)
    {
        memcpy(pu1_dst, pu1_src, i4_num_bytes);
        pu1_src += i4_src_stride;
        pu1_dst += i4_dst_stride;
    }
}

static void isvc_copy_data_semiplanr(UWORD8 *pu1_src, WORD32 i4_src_stride, UWORD8 *pu1_dst1,
                                     UWORD8 *pu1_dst2, WORD32 i4_dst_stride, WORD32 i4_num_bytes,
                                     WORD32 i4_num_lines)
{
    WORD32 i4_vert_lines, u4_i;

    ASSERT(NULL != pu1_src);
    ASSERT(NULL != pu1_dst1);
    ASSERT(NULL != pu1_dst2);

    for(i4_vert_lines = 0; i4_vert_lines < i4_num_lines; i4_vert_lines++)
    {
        for(u4_i = 0; u4_i < i4_num_bytes; u4_i++)
        {
            *(pu1_dst1 + u4_i) = *(pu1_src + (2 * u4_i));
            *(pu1_dst2 + u4_i) = *(pu1_src + (2 * u4_i) + 1);
        }
        pu1_src += i4_src_stride;
        pu1_dst1 += i4_dst_stride;
        pu1_dst2 += i4_dst_stride;
    }
}

static void isvc_left_right_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                    UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                    UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                    WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                    WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(i1_yd_index);
    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + i1_xd_index;

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - i1_xd_index;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst, *pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
        pu1_src += i4_refarray_stride;
    }
}

static void isvc_left_right_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                           WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                           UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                           WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                           WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i1_yd_index);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_src_cb = pu1_dst_cb + i1_xd_index;

    pu1_dst_cr = pu1_refarray_2 + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i1_xd_index;

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_CHROMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_CHROMA);
    pu1_dst_cb = pu1_src_cb - i1_xd_index;
    pu1_dst_cr = pu1_src_cr - i1_xd_index;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst_cb, *pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;
        pu1_src_cb += i4_refarray_stride;

        memset(pu1_dst_cr, *pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
        pu1_src_cr += i4_refarray_stride;
    }
}

static void isvc_top_bot_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                 UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                 UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                 WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                 WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(i1_xd_index);
    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + (i1_yd_index * i4_refarray_stride);

    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - (i1_yd_index * i4_refarray_stride);

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memcpy(pu1_dst, pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
    }
}

static void isvc_top_bot_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                        WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                        UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                        WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                        WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i1_xd_index);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;

    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_src_cb = pu1_dst_cb + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i4_tmp;

    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_CHROMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_CHROMA);

    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_dst_cb = pu1_src_cb - i4_tmp;

    pu1_dst_cr = pu1_src_cr - i4_tmp;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memcpy(pu1_dst_cb, pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;

        memcpy(pu1_dst_cr, pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
    }
}

static void isvc_diag_reconstruction(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                     UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                     UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                     WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                     WORD32 i4_corner_pixel_available)
{
    WORD32 i4_i;
    UWORD8 *pu1_src_1, *pu1_src_2, *pu1_dst;
    UWORD8 u1_filter_delay_buf[18];
    UWORD8 u1_out_buf[16];
    WORD32 i4_width, i4_height;
    WORD32 i4_x_off, i4_y_off;
    WORD32 i4_block_size = BLK_SIZE;

    UNUSED(pu1_refarray_2);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src_1 = pu1_dst + i1_xd_index;
    pu1_src_2 = pu1_dst + (i1_yd_index * i4_refarray_stride);

    i4_width = MAX(u1_seg_wd, (((i4_mb_adjoin_x >> 3) ^ 1) * i4_block_size));
    i4_height = MAX(u1_seg_ht, (((i4_mb_adjoin_y >> 4) ^ 1) * i4_block_size));

    i4_x_off = (i4_width - u1_seg_wd);
    i4_y_off = (i4_height - u1_seg_ht);

    if(i1_xd_index > 0 && i1_yd_index > 0)
    {
        /* Quadrant 1 Processing */

        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf[0], pu1_src_2, (i4_width + 1));
        for(i4_i = i4_height; i4_i > 0; i4_i--)
        {
            u1_filter_delay_buf[i4_width + i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_width] =
                (u1_filter_delay_buf[i4_width - 1] + u1_filter_delay_buf[i4_width + 1] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_height - i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else if(i1_xd_index < 0 && i1_yd_index > 0)
    {
        /* Quadrant 2 Processing */
        /* load the pixel in the filter delay buffer  */
        for(i4_i = 0; i4_i < (i4_height + 1); i4_i++)
        {
            u1_filter_delay_buf[i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        pu1_src_2 -= i4_x_off;
        memcpy(&u1_filter_delay_buf[i4_i], pu1_src_2, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_i - 1] =
                (u1_filter_delay_buf[i4_i] + u1_filter_delay_buf[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index < 0)
    {
        /* Quadrant 3 Processing */
        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf[0], pu1_src_2, (i4_width + 1));

        pu1_src_1 -= (i4_y_off * i4_refarray_stride);
        for(i4_i = 1; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf[i4_width + i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_width] =
                (u1_filter_delay_buf[i4_width - 1] + u1_filter_delay_buf[i4_width + 1] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else
    {
        /* Quadrant 4 Processing */
        /* load the pixel in the filter delay buffer  */
        pu1_src_1 += ((u1_seg_ht - 1) * i4_refarray_stride);
        for(i4_i = 0; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf[i4_i] = *pu1_src_1;
            pu1_src_1 -= i4_refarray_stride;
        }

        pu1_src_2 -= i4_x_off;
        memcpy(&u1_filter_delay_buf[i4_i], pu1_src_2, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_i - 1] =
                (u1_filter_delay_buf[i4_i] + u1_filter_delay_buf[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
}

static void isvc_diag_reconstruction_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                            WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                            UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                            WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                            WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)
{
    WORD32 i4_i;
    UWORD8 u1_filter_delay_buf_cb[18], u1_filter_delay_buf_cr[18];
    UWORD8 u1_out_buf_cb[16], u1_out_buf_cr[16];
    WORD32 i4_width, i4_height;
    WORD32 i4_x_off, i4_y_off;
    WORD32 i4_block_size = BLK_SIZE >> 1;
    UWORD8 *pu1_src_1_cb, *pu1_src_2_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_1_cr, *pu1_src_2_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;

    pu1_src_1_cb = pu1_dst_cb + i1_xd_index;
    pu1_src_1_cr = pu1_dst_cr + i1_xd_index;

    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_src_2_cb = pu1_dst_cb + i4_tmp;
    pu1_src_2_cr = pu1_dst_cr + i4_tmp;

    i4_width = MAX(u1_seg_wd, (((i4_mb_adjoin_x >> 3) ^ 1) * i4_block_size));
    i4_height = MAX(u1_seg_ht, (((i4_mb_adjoin_y >> 4) ^ 1) * i4_block_size));

    i4_x_off = (i4_width - u1_seg_wd);
    i4_y_off = (i4_height - u1_seg_ht);

    if(i1_xd_index < 0 && i1_yd_index > 0)
    {
        /* Quadrant 1 Processing */

        /* load the pixel in the filter delay buffer  */
        for(i4_i = 0; i4_i < (i4_height + 1); i4_i++)
        {
            u1_filter_delay_buf_cb[i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        pu1_src_2_cb -= i4_x_off;
        pu1_src_2_cr -= i4_x_off;

        memcpy(&u1_filter_delay_buf_cb[i4_i], pu1_src_2_cb, i4_width);
        memcpy(&u1_filter_delay_buf_cr[i4_i], pu1_src_2_cr, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_i - 1] =
                (u1_filter_delay_buf_cb[i4_i] + u1_filter_delay_buf_cb[i4_i - 2] + 1) >> 1;

            u1_filter_delay_buf_cr[i4_i - 1] =
                (u1_filter_delay_buf_cr[i4_i] + u1_filter_delay_buf_cr[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;

            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index > 0)
    {
        /* Quadrant 2 Processing */

        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf_cb[0], pu1_src_2_cb, (i4_width + 1));
        memcpy(&u1_filter_delay_buf_cr[0], pu1_src_2_cr, (i4_width + 1));

        for(i4_i = i4_height; i4_i > 0; i4_i--)
        {
            u1_filter_delay_buf_cb[i4_width + i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_width + i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_width] =
                (u1_filter_delay_buf_cb[i4_width - 1] + u1_filter_delay_buf_cb[i4_width + 1] + 1) >>
                1;

            u1_filter_delay_buf_cr[i4_width] =
                (u1_filter_delay_buf_cr[i4_width - 1] + u1_filter_delay_buf_cr[i4_width + 1] + 1) >>
                1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_height - i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;

            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_height - i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index < 0)
    {
        /* Quadrant 3 Processing */

        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf_cb[0], pu1_src_2_cb, (i4_width + 1));
        memcpy(&u1_filter_delay_buf_cr[0], pu1_src_2_cr, (i4_width + 1));

        i4_tmp = (i4_y_off * i4_refarray_stride);
        pu1_src_1_cb -= i4_tmp;
        pu1_src_1_cr -= i4_tmp;
        for(i4_i = 1; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf_cb[i4_width + i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_width + i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_width] =
                (u1_filter_delay_buf_cb[i4_width - 1] + u1_filter_delay_buf_cb[i4_width + 1] + 1) >>
                1;

            u1_filter_delay_buf_cr[i4_width] =
                (u1_filter_delay_buf_cr[i4_width - 1] + u1_filter_delay_buf_cr[i4_width + 1] + 1) >>
                1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;

            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else
    {
        /* Quadrant 4 Processing */

        /* load the pixel in the filter delay buffer  */
        i4_tmp = ((u1_seg_ht - 1) * i4_refarray_stride);
        pu1_src_1_cb += i4_tmp;
        pu1_src_1_cr += i4_tmp;

        for(i4_i = 0; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf_cb[i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb -= i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr -= i4_refarray_stride;
        }

        pu1_src_2_cb -= i4_x_off;
        pu1_src_2_cr -= i4_x_off;

        memcpy(&u1_filter_delay_buf_cb[i4_i], pu1_src_2_cb, i4_width);
        memcpy(&u1_filter_delay_buf_cr[i4_i], pu1_src_2_cr, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_i - 1] =
                (u1_filter_delay_buf_cb[i4_i] + u1_filter_delay_buf_cb[i4_i - 2] + 1) >> 1;

            u1_filter_delay_buf_cr[i4_i - 1] =
                (u1_filter_delay_buf_cr[i4_i] + u1_filter_delay_buf_cr[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;

            memcpy(pu1_dst_cr, &u1_out_buf_cr[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
}

static void isvc_diag_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                              UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                              UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                              WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                              WORD32 i4_corner_pixel_available)

{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + i1_xd_index + (i1_yd_index * i4_refarray_stride);

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);

    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - i1_xd_index - (i1_yd_index * i4_refarray_stride);

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst, *pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
    }
}

static void isvc_diag_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                     UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                     UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                     WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                     WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;

    i4_tmp = i1_xd_index + (i1_yd_index * i4_refarray_stride);
    pu1_src_cb = pu1_dst_cb + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i4_tmp;

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);

    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);

    i4_tmp = (i1_xd_index + (i1_yd_index * i4_refarray_stride));
    pu1_dst_cb = pu1_src_cb - i4_tmp;
    pu1_dst_cr = pu1_src_cr - i4_tmp;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst_cb, *pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;

        memset(pu1_dst_cr, *pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
    }
}

static PF_INTRA_SAMP_PADDING *gpf_lookup_fxns_luma[NUM_INTRA_SAMP_FXNS] = {
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &isvc_left_right_padding,
    NULL,
    &isvc_diag_reconstruction,
    NULL,
    &isvc_left_right_padding,
    NULL,
    &isvc_diag_reconstruction,
    NULL,
    NULL,
    &isvc_top_bot_padding,
    &isvc_diag_reconstruction,
    NULL,
    NULL,
    &isvc_top_bot_padding,
    &isvc_diag_reconstruction,
    NULL,
    &isvc_left_right_padding,
    &isvc_top_bot_padding,
    &isvc_diag_reconstruction,
    &isvc_diag_padding,
    &isvc_left_right_padding,
    &isvc_top_bot_padding,
    &isvc_diag_reconstruction,
};

static PF_INTRA_SAMP_PADDING *gpf_lookup_fxns_chroma[NUM_INTRA_SAMP_FXNS] = {
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &isvc_left_right_padding_chroma,
    NULL,
    &isvc_diag_reconstruction_chroma,
    NULL,
    &isvc_left_right_padding_chroma,
    NULL,
    &isvc_diag_reconstruction_chroma,
    NULL,
    NULL,
    &isvc_top_bot_padding_chroma,
    &isvc_diag_reconstruction_chroma,
    NULL,
    NULL,
    &isvc_top_bot_padding_chroma,
    &isvc_diag_reconstruction_chroma,
    NULL,
    &isvc_left_right_padding_chroma,
    &isvc_top_bot_padding_chroma,
    &isvc_diag_reconstruction_chroma,
    &isvc_diag_padding_chroma,
    &isvc_left_right_padding_chroma,
    &isvc_top_bot_padding_chroma,
    &isvc_diag_reconstruction_chroma,
};

static void isvc_get_ref_layer_avlblty_dyadic(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                                              WORD32 i4_element_size, WORD32 i4_ref_mb_x,
                                              WORD32 i4_ref_mb_y, WORD32 *pi4_avlblty,
                                              WORD8 i1_curr_slice_id, WORD8 i1_cons_intr_samp_flag)
{
    WORD8 i1_mb_mode;

    pi1_ref_mb_modes += (i4_ref_mb_y * i4_ref_mode_stride * i4_element_size);
    pi1_ref_mb_modes += (i4_ref_mb_x * i4_element_size);
    i1_mb_mode = *pi1_ref_mb_modes;
    i1_mb_mode = (i1_mb_mode < 0) ? i1_mb_mode : SVC_EXTRACT_MB_MODE(*pi1_ref_mb_modes);

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        *pi4_avlblty = 0;
    }
    else
    {
        *pi4_avlblty = 1;
    }

    if(1 == i1_cons_intr_samp_flag)
    {
        if(1 == *pi4_avlblty)
        {
            if(i1_mb_mode != i1_curr_slice_id)
            {
                *pi4_avlblty = 0;
            }
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_diagonal_construct_dyadic                          */
/*                                                                           */
/*  Description   : This function fills the unavaible pixels in the reference*/
/*                  array with diagonally constructed samples                */
/*  Inputs        : i4_x :current position in reference array X to be filled */
/*                  i4_y :current position in reference array Y to be filled */
/*                  i4_xd_index : diagonal index in horizontal direction     */
/*                  i4_yd_index : diagonal index in vertical direction       */
/*                  pu1_refarray : popinter to reference array               */
/*                  i4_refarray_wd: width of the reference array             */
/*  Globals       : none                                                     */
/*  Processing    : Fills the sample which is unavailable with filtered      */
/*                  diagonal samples                                         */
/*  Outputs       : pixel filled                                             */
/*  Returns       : constructed pixel                                        */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         03 12 2010   Nithya          creation                             */
/*                                                                           */
/*****************************************************************************/
static UWORD8 isvc_diagonal_construct_dyadic(WORD32 i4_x, WORD32 i4_y, WORD32 i4_xd_index,
                                             WORD32 i4_yd_index, UWORD8 *pu1_refarray,
                                             WORD32 i4_refarray_wd)
{
    WORD32 i4_diff_hor_ver, i4_sgn_xy;
    WORD32 i4_xc, i4_yc;
    WORD32 i4_samp1, i4_samp2, i4_samp3;
    WORD32 i4_result;
    UWORD8 *pu1_tmp;

    i4_diff_hor_ver = ABS(i4_xd_index) - ABS(i4_yd_index);
    i4_sgn_xy = SIGN(i4_xd_index * i4_yd_index);

    if(i4_diff_hor_ver > 0)
    {
        i4_xc = i4_x - (i4_sgn_xy * i4_yd_index);
        i4_yc = i4_y - i4_yd_index;

        pu1_tmp = pu1_refarray + (i4_yc * i4_refarray_wd);

        i4_samp1 = pu1_tmp[i4_xc - 1];
        i4_samp2 = pu1_tmp[i4_xc];
        i4_samp3 = pu1_tmp[i4_xc + 1];
    }
    else if(i4_diff_hor_ver < 0)
    {
        i4_xc = i4_x - i4_xd_index;
        i4_yc = i4_y - (i4_sgn_xy * i4_xd_index);

        pu1_tmp = pu1_refarray + ((i4_yc - 1) * i4_refarray_wd);

        i4_samp1 = pu1_tmp[i4_xc];
        pu1_tmp += i4_refarray_wd;
        i4_samp2 = pu1_tmp[i4_xc];
        pu1_tmp += i4_refarray_wd;
        i4_samp3 = pu1_tmp[i4_xc];
    }
    else
    {
        WORD32 i4_ref_xd, i4_ref_yd;

        i4_ref_xd = i4_x - i4_xd_index;
        i4_ref_yd = i4_y - i4_yd_index;

        i4_xc = i4_ref_xd + SIGN(i4_xd_index);
        i4_yc = i4_ref_yd + SIGN(i4_yd_index);

        pu1_tmp = pu1_refarray + (i4_ref_yd * i4_refarray_wd);

        i4_samp1 = pu1_tmp[i4_xc];
        i4_samp2 = pu1_tmp[i4_ref_xd];
        pu1_tmp = pu1_refarray + (i4_yc * i4_refarray_wd);
        i4_samp3 = pu1_tmp[i4_ref_xd];
    }

    i4_result = (i4_samp1 + (i4_samp2 << 1) + i4_samp3 + 2) >> 2;

    pu1_tmp = pu1_refarray + (i4_y * i4_refarray_wd);
    /* Store the filled sample */
    pu1_tmp[i4_x] = i4_result;

    return i4_result;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_corner_samp_dyadic                                 */
/*                                                                           */
/*  Description   : This function fills the corner sample in the reference   */
/*                  array with diagonally constructed samples                */
/*  Inputs        : i4_x :current position in reference array X to be filled */
/*                  i4_y :current position in reference array Y to be filled */
/*                  i4_xd_index : diagonal index in horizontal direction     */
/*                  i4_yd_index : diagonal index in vertical direction       */
/*                  pu1_refarray_y : pointer to luma reference array         */
/*                  pu1_refarray_cb : pointer to Cb reference array          */
/*                  pu1_refarray_cr : pointer to Cr reference array          */
/*  Globals       : none                                                     */
/*  Processing    : Fills the sample which is unavailable with filtered      */
/*                  diagonal samples                                         */
/*  Outputs       : pixel filled                                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         03 12 2010   Nithya          creation                             */
/*                                                                           */
/*****************************************************************************/
static void isvc_corner_samp_dyadic(WORD32 i4_x, WORD32 i4_y, WORD32 i4_xD, WORD32 i4_yD,
                                    UWORD8 *pu1_refarray_y, UWORD8 *pu1_refarray_cb,
                                    UWORD8 *pu1_refarray_cr)
{
    WORD32 i4_ref_xD, i4_ref_yD;
    WORD32 i4_c_ref_xD, i4_c_ref_yD;
    WORD32 i4_xc, i4_yc;
    WORD32 i4_c_xc, i4_c_yc;
    WORD32 i4_samp1, i4_samp2;
    UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

    i4_ref_xD = i4_x - i4_xD;
    i4_ref_yD = i4_y - i4_yD;

    i4_xc = i4_ref_xD + SIGN(i4_xD);
    i4_yc = i4_ref_yD + SIGN(i4_yD);

    /* Luma */
    pu1_tmp_src = pu1_refarray_y + (i4_yc * DYADIC_REF_W_Y);
    i4_samp1 = pu1_tmp_src[i4_ref_xD];
    pu1_tmp_src = pu1_refarray_y + (i4_ref_yD * DYADIC_REF_W_Y);
    i4_samp2 = pu1_tmp_src[i4_xc];
    pu1_tmp_dst = pu1_tmp_src;

    pu1_tmp_dst[i4_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;

    /* Chroma */
    i4_c_ref_xD = i4_ref_xD >> 1;
    i4_c_ref_yD = i4_ref_yD >> 1;

    i4_c_xc = i4_c_ref_xD + SIGN(i4_xD);
    i4_c_yc = i4_c_ref_yD + SIGN(i4_yD);

    /* Cb */
    pu1_tmp_src = pu1_refarray_cb + (i4_c_yc * DYADIC_REF_W_C);
    i4_samp1 = pu1_tmp_src[i4_c_ref_xD];
    pu1_tmp_src = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
    i4_samp2 = pu1_tmp_src[i4_c_xc];
    pu1_tmp_dst = pu1_tmp_src;

    pu1_tmp_dst[i4_c_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;

    /* Cr */
    pu1_tmp_src = pu1_refarray_cr + (i4_c_yc * DYADIC_REF_W_C);
    i4_samp1 = pu1_tmp_src[i4_c_ref_xD];
    pu1_tmp_src = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
    i4_samp2 = pu1_tmp_src[i4_c_xc];
    pu1_tmp_dst = pu1_tmp_src;

    pu1_tmp_dst[i4_c_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_reflayer_construction_dyadic                       */
/*                                                                           */
/*  Description   :  This function constructs the reference array buffer     */
/*                   for dyadic cases used for intra resampling of a         */
/*                   component in an MB                                      */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_mb_mode_map : ref layer mb mode buffer desc       */
/*                  pu1_inp_luma : luma input (reference layer data)         */
/*                  pu1_inp_chroma : chroma input (reference layer data)     */
/*                  i4_inp_luma_stride : luma input buffer stride            */
/*                  i4_inp_chroma_stride : chroma input buffer stride        */
/*                  i4_top : indicates whether the core 8x8 reference block  */
/*                  is one of 0 and 1 or one of 2 and 3                      */
/*                  i4_left : indicates whether the core 8x8 reference block */
/*                  is one of 0 and 2 or one of 1 and 3                      */
/*                  ps_ref_mb_coord : coordinates of the reference MB        */
/*  Globals       : none                                                     */
/*  Processing    : it fills the reference layer data if they are falling in */
/*                  INTRA MB region. If all the pixels are not filled  it    */
/*                  calls the border extension algorithm to fill them        */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         02 12 2010   Nithya          creation               */
/*                                                                           */
/*****************************************************************************/
static void isvc_reflayer_construction_dyadic(void *pv_intra_samp_ctxt,
                                              mem_element_t *ps_ref_mb_mode_map,
                                              UWORD8 *pu1_inp_luma, UWORD8 *pu1_inp_chroma,
                                              WORD32 i4_inp_luma_stride,
                                              WORD32 i4_inp_chroma_stride, WORD32 i4_top,
                                              WORD32 i4_left, UWORD16 u2_mb_x, UWORD16 u2_mb_y)
{
    enum
    {
        TOPLEFT_MASK = 1,
        LEFT_MASK = 2,
        TOP_MASK = 4,
        TOPRIGHT_MASK = 8,
        BOTTOMLEFT_MASK = 16
    };

    WORD32 i4_x, i4_y;
    WORD32 i4_x0, i4_y0;
    WORD32 i4_xc0, i4_yc0;
    WORD32 i4_ref_xD, i4_ref_yD;
    WORD32 i4_c_ref_xD, i4_c_ref_yD;

    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD8 *pi1_ref_mb_modes;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    WORD32 i4_mbaddr_y;
    WORD32 i4_mbaddr_x;

    WORD32 i4_refarray_wd_luma, i4_refarray_wd_chroma;
    WORD32 i4_refarray_ht_luma, i4_refarray_ht_chroma;
    WORD32 i4_avlblty;
    WORD8 i1_cons_intr_samp_flag;
    WORD8 i1_slice_id;
    WORD8 i1_corner_samp_avlbl_flag;
    UWORD8 u1_ny_avlblty;

    UWORD8 *pu1_refarray_luma;
    UWORD8 *pu1_refarray_cb, *pu1_refarray_cr;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode_map->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode_map->i4_element_size;

    i1_cons_intr_samp_flag = ps_lyr_ctxt->i1_constrained_intra_rsmpl_flag;

    ASSERT(NULL != pi1_ref_mb_modes);

    pu1_refarray_luma = ps_ctxt->pu1_refarray_buffer;
    pu1_refarray_cb = ps_ctxt->pu1_refarray_cb;
    pu1_refarray_cr = ps_ctxt->pu1_refarray_cr;

    i4_mbaddr_x = u2_mb_x;
    i4_mbaddr_y = u2_mb_y;

    i4_refarray_wd_luma = 20;
    i4_refarray_ht_luma = 20;

    i4_refarray_wd_chroma = i4_refarray_wd_luma >> 1;
    i4_refarray_ht_chroma = i4_refarray_ht_luma >> 1;

    if(1 == i1_cons_intr_samp_flag)
    {
        WORD8 *pi1_ref_mb_mode_tmp;
        WORD8 i1_mb_mode;

        pi1_ref_mb_mode_tmp = pi1_ref_mb_modes;
        pi1_ref_mb_mode_tmp += (i4_mbaddr_y * i4_ref_mode_stride * i4_element_size);
        pi1_ref_mb_mode_tmp += (i4_mbaddr_x * i4_element_size);
        i1_mb_mode = *pi1_ref_mb_mode_tmp;
        i1_mb_mode = (i1_mb_mode < 0) ? i1_mb_mode : SVC_EXTRACT_MB_MODE(*pi1_ref_mb_mode_tmp);

        /* The reference layer MB should be intra */
        ASSERT(i1_mb_mode >= 0);

        i1_slice_id = i1_mb_mode;
    }
    else
    {
        i1_slice_id = -1;
    }

    {
        UWORD8 *pu1_src, *pu1_dst;
        WORD32 i4_src_stride, i4_dst_stride;

        /* Copy luma */
        i4_src_stride = i4_inp_luma_stride;
        i4_dst_stride = DYADIC_REF_W_Y;
        pu1_src = pu1_inp_luma;
        pu1_dst = pu1_refarray_luma;

        isvc_copy_data(pu1_src, i4_src_stride, pu1_dst, i4_dst_stride, i4_refarray_wd_luma,
                       i4_refarray_ht_luma);

        i4_src_stride = i4_inp_chroma_stride;
        i4_dst_stride = DYADIC_REF_W_C;
        pu1_src = pu1_inp_chroma;
        isvc_copy_data_semiplanr(pu1_src, i4_src_stride, pu1_refarray_cb, pu1_refarray_cr,
                                 i4_dst_stride, i4_refarray_wd_chroma, i4_refarray_ht_chroma);
    }

    {
        /* mb_x + left, mb_y + top */
        isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                          i4_mbaddr_x + i4_left, i4_mbaddr_y + i4_top, &i4_avlblty,
                                          i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty = i4_avlblty;

        /* mb_x + left, mb_y */
        isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                          i4_mbaddr_x + i4_left, i4_mbaddr_y, &i4_avlblty,
                                          i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 1);

        /* mb_x, mb_y + top */
        isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                          i4_mbaddr_x, i4_mbaddr_y + i4_top, &i4_avlblty,
                                          i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 2);

        /* mb_x - left, mb_y + top */
        isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                          i4_mbaddr_x - i4_left, i4_mbaddr_y + i4_top, &i4_avlblty,
                                          i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 3);

        /* mb_x + left, mb_y - top */
        isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                          i4_mbaddr_x + i4_left, i4_mbaddr_y - i4_top, &i4_avlblty,
                                          i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 4);
    }

    if((TOP_MASK | TOPLEFT_MASK | LEFT_MASK) == u1_ny_avlblty)
    {
        return;
    }

    if(!(u1_ny_avlblty & (TOP_MASK | TOPLEFT_MASK | LEFT_MASK)))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Set the 4 corner samples to (x-xD,y-yD) */
        i4_x0 = 9 + (i4_left << 3) + i4_left;
        i4_y0 = 9 + (i4_top << 3) + i4_top;

        i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
        i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;

        pu1_tmp_dst1[i4_x0] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst2[i4_x0] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];

        /* Set the corner sample of Cb and Cr to (x-xD,y-yD) */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;

        i4_c_ref_yD = i4_ref_yD >> 1;
        i4_c_ref_xD = i4_ref_xD >> 1;

        pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
        pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];
    }

    if(!(u1_ny_avlblty & (TOP_MASK | TOPLEFT_MASK)))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Copy (x0,ref_yD), (x0+1,ref_yD), ..., (x0+7,ref_yD) to */
        /* (x0,y0), (x0+1,y0), ..., (x0+7,y0) and   */
        /* (x0,y0+1), (x0+1,y0+1), ..., (x0+7,y0+1) */
        i4_x0 = 2;
        i4_y0 = 9 + (i4_top << 3) + i4_top;
        if(i4_left > 0)
        {
            i4_x0 += 8;
        }
        i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;

        for(i4_x = i4_x0; i4_x < i4_x0 + 8; i4_x++)
        {
            pu1_tmp_dst1[i4_x] = pu1_tmp_src[i4_x];
            pu1_tmp_dst2[i4_x] = pu1_tmp_src[i4_x];
        }

        /* Cb and Cr copy */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;
        i4_c_ref_yD = i4_ref_yD >> 1;

        pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);

        for(i4_x = i4_xc0; i4_x < i4_xc0 + 4; i4_x++)
        {
            pu1_tmp_dst1[i4_x] = pu1_tmp_src1[i4_x];
            pu1_tmp_dst2[i4_x] = pu1_tmp_src2[i4_x];
        }
    }

    if(!(u1_ny_avlblty & (TOPLEFT_MASK | LEFT_MASK)))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Copy (ref_xD,y0) to (x0,y0) and (x0+1,y0); */
        /* copy (ref_xD,y0+1) to (x0,y0+1) and (x0+1,y0+1); ... ;*/
        /* copy (ref_xD,y0+7) to (x0,y0+7) and (x0+1,y0+7) */
        i4_x0 = 9 + (i4_left << 3) + i4_left;
        i4_y0 = 2;
        if(i4_top > 0)
        {
            i4_y0 += 8;
        }
        i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_tmp_src;

        for(i4_y = i4_y0; i4_y < i4_y0 + 8; i4_y++)
        {
            pu1_tmp_dst1[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_src += DYADIC_REF_W_Y;
            pu1_tmp_dst1 += DYADIC_REF_W_Y;
        }

        /* Cb and Cr copy */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;
        i4_c_ref_xD = i4_ref_xD >> 1;

        pu1_tmp_src1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_tmp_src1;
        pu1_tmp_src2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_tmp_src2;

        for(i4_y = i4_yc0; i4_y < i4_yc0 + 4; i4_y++)
        {
            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];
            pu1_tmp_src1 += DYADIC_REF_W_C;
            pu1_tmp_src2 += DYADIC_REF_W_C;
            pu1_tmp_dst1 += DYADIC_REF_W_C;
            pu1_tmp_dst2 += DYADIC_REF_W_C;
        }
    }

    if(!(u1_ny_avlblty & TOP_MASK))
    {
        if(!(u1_ny_avlblty & TOPRIGHT_MASK))
        {
            UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

            i4_x0 = 9 - i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;

            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            /* Copy (x0,ref_yD) and (x0+1,ref_yD) to (x0,y0) and (x0+1,y0), and */
            /* to (x0,y0+1) and (x0+1,y0+1) */
            pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
            pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);

            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_x0];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_x0 + 1];

            pu1_tmp_dst += DYADIC_REF_W_Y;

            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_x0];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_x0 + 1];

            /* Cb copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            pu1_tmp_src = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);

            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_xc0];

            /* Cr copy */
            pu1_tmp_src = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);

            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_xc0];
        }
        else
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                              i4_mbaddr_x - i4_left, i4_mbaddr_y, &i4_avlblty,
                                              i1_slice_id, i1_cons_intr_samp_flag);
            i1_corner_samp_avlbl_flag = i4_avlblty;

            i4_x0 = 9 - i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;

            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;

            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);
            i4_ref_xD = i4_x0 - (i4_left * 7) - (i4_left >> 1);

            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;

            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Fill corner sample if not available */
            if(!i1_corner_samp_avlbl_flag)
            {
                isvc_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                        pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call diagonal construction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvc_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                   DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                           DYADIC_REF_W_C);

            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                           DYADIC_REF_W_C);
        }
    }

    if(!(u1_ny_avlblty & LEFT_MASK))
    {
        if(!(u1_ny_avlblty & BOTTOMLEFT_MASK))
        {
            UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 - i4_top;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

            /* Copy (ref_xD,y0) to (x0,y0), (x0+1,y0), and  */
            /* copy (ref_xD,y0+1) to (x0,y0+1), (x0+1,y0+1) */
            pu1_tmp_src = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst = pu1_tmp_src;

            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];

            pu1_tmp_src += DYADIC_REF_W_Y;
            pu1_tmp_dst += DYADIC_REF_W_Y;

            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];

            /* Cb copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;

            pu1_tmp_src = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_tmp_src;

            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_c_ref_xD];

            /* Cr copy */
            pu1_tmp_src = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_tmp_src;

            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_c_ref_xD];
        }
        else
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            isvc_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                              i4_mbaddr_x, i4_mbaddr_y - i4_top, &i4_avlblty,
                                              i1_slice_id, i1_cons_intr_samp_flag);
            i1_corner_samp_avlbl_flag = i4_avlblty;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 - i4_top;

            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;

            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
            i4_ref_yD = i4_y0 - (i4_top * 7) - (i4_top >> 1);

            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;

            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            if(!i1_corner_samp_avlbl_flag)
            {
                isvc_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                        pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call diagonal consrtuction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvc_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                   DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                           DYADIC_REF_W_C);

            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                           DYADIC_REF_W_C);
        }
    }

    if(u1_ny_avlblty & TOPLEFT_MASK)
    {
        if(!(u1_ny_avlblty & LEFT_MASK))
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;
            UWORD8 *pu1_tmp_dst;
            UWORD8 u1_filled_samp;

            i1_corner_samp_avlbl_flag = (u1_ny_avlblty & 4) >> 2;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 2;
            i4_ref_yD = 1;
            if(i4_top > 0)
            {
                i4_y0 += 8;
                i4_ref_yD = 18;
            }

            i4_ref_xD = i4_x0 - (i4_left) - (i4_left >> 1);

            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;

            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;

            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Fill corner sample if unavailable */
            if(!i1_corner_samp_avlbl_flag)
            {
                isvc_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                        pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call the diagonal construction for the 8 rows */
            if(i4_top == i4_left)
            {
                /* if top * left = 1 */
                /* (x0,y0) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0, i4_y0, i4_xD, i4_yD,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);

                /* (x0,y0+1), ..., (x0,y0+7) and */
                /* (x0+1,y0), ..., (x0+1,y0+6)   */
                for(i4_y = i4_y0 + 1; i4_y < i4_y0 + 8; i4_y++)
                {
                    i4_yD++;
                    u1_filled_samp = isvc_diagonal_construct_dyadic(
                        i4_x0, i4_y, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
                    pu1_tmp_dst[i4_x0 + 1] = u1_filled_samp;
                    pu1_tmp_dst += DYADIC_REF_W_Y;
                }

                /* (x0+1,y0+7) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(
                    i4_x0 + 1, i4_y0 + 7, i4_xD + 1, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
            }
            else
            {
                /* top * left = -1 */
                /* (x0+1,y0) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0 + 1, i4_y0, i4_xD + 1, i4_yD,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);

                /* (x0,y0), ..., (x0,y0+6) and   */
                /* (x0+1,y0+1), ..., (x0+1,y0+7) */
                for(i4_y = i4_y0; i4_y < i4_y0 + 7; i4_y++)
                {
                    u1_filled_samp = isvc_diagonal_construct_dyadic(
                        i4_x0, i4_y, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);

                    pu1_tmp_dst += DYADIC_REF_W_Y;
                    pu1_tmp_dst[i4_x0 + 1] = u1_filled_samp;
                    i4_yD++;
                }

                /* (x0,y0+7) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0, i4_y0 + 7, i4_xD, i4_yD,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);
            }

            /* For Cb and Cr */
            for(i4_y = i4_yc0; i4_y < i4_yc0 + 4; i4_y++)
            {
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_xc0, i4_y, i4_c_xD, i4_c_yD,
                                                                pu1_refarray_cb, DYADIC_REF_W_C);
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_xc0, i4_y, i4_c_xD, i4_c_yD,
                                                                pu1_refarray_cr, DYADIC_REF_W_C);
                i4_c_yD++;
            }
        }

        if(!(u1_ny_avlblty & TOP_MASK))
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;
            UWORD8 *pu1_tmp_dst;
            UWORD8 u1_filled_samp;

            i1_corner_samp_avlbl_flag = (u1_ny_avlblty & 2) >> 1;

            i4_y0 = 9 + (i4_top << 3) + (i4_top);
            i4_x0 = 2;
            i4_ref_xD = 1;
            if(i4_left > 0)
            {
                i4_x0 += 8;
                i4_ref_xD = 18;
            }

            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;

            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;

            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            if(!i1_corner_samp_avlbl_flag)
            {
                isvc_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                        pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call the diagonal construction for the 2 rows */
            if(i4_top == i4_left)
            {
                /* if top * left = 1 */
                /* (x0,y0) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0, i4_y0, i4_xD, i4_yD,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + ((i4_y0 + 1) * DYADIC_REF_W_Y);

                /* (x0+1,y0), ..., (x0+7,y0) and */
                /* (x0,y0+1), ..., (x0+6,y0+1)   */
                for(i4_x = i4_x0 + 1; i4_x < i4_x0 + 8; i4_x++)
                {
                    i4_xD++;
                    u1_filled_samp = isvc_diagonal_construct_dyadic(
                        i4_x, i4_y0, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
                    pu1_tmp_dst[i4_x - 1] = u1_filled_samp;
                }

                /* (x0+7,y0+1) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(
                    i4_x0 + 7, i4_y0 + 1, i4_xD, i4_yD + 1, pu1_refarray_luma, DYADIC_REF_W_Y);
            }
            else
            {
                /* top * left = -1 */
                /* (x0,y0+1) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0, i4_y0 + 1, i4_xD, i4_yD + 1,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + ((i4_y0 + 1) * DYADIC_REF_W_Y);

                /* (x0,y0), ..., (x0,y0+6) and   */
                /* (x0+1,y0+1), ..., (x0+1,y0+7) */
                for(i4_x = i4_x0; i4_x < i4_x0 + 7; i4_x++)
                {
                    u1_filled_samp = isvc_diagonal_construct_dyadic(
                        i4_x, i4_y0, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);

                    pu1_tmp_dst[i4_x + 1] = u1_filled_samp;
                    i4_xD++;
                }

                /* (x0+7,y0) */
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x0 + 7, i4_y0, i4_xD, i4_yD,
                                                                pu1_refarray_luma, DYADIC_REF_W_Y);
            }

            /* For Cb and Cr */
            for(i4_x = i4_xc0; i4_x < i4_xc0 + 4; i4_x++)
            {
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x, i4_yc0, i4_c_xD, i4_c_yD,
                                                                pu1_refarray_cb, DYADIC_REF_W_C);
                u1_filled_samp = isvc_diagonal_construct_dyadic(i4_x, i4_yc0, i4_c_xD, i4_c_yD,
                                                                pu1_refarray_cr, DYADIC_REF_W_C);
                i4_c_xD++;
            }
        }
    }

    if(!(u1_ny_avlblty & TOPLEFT_MASK))
    {
        UWORD8 *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        if(u1_ny_avlblty & LEFT_MASK)
        {
            /* (mb_x+left,mb_y) available, (mb_x,mb_y+top) unavailable */
            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            /* Copy (x0,ref_yD), (x0+1,ref_yD) to  */
            /* (x0,y0), (x0+1,y0), and (x0,y0+1), (x0+1,y0+1) */
            pu1_tmp_src1 = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
            pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;

            pu1_tmp_dst1[i4_x0] = pu1_tmp_src1[i4_x0];
            pu1_tmp_dst2[i4_x0] = pu1_tmp_src1[i4_x0];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src1[i4_x0 + 1];
            pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src1[i4_x0 + 1];

            /* Cb and Cr copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);

            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_xc0];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_xc0];
        }
        else if(u1_ny_avlblty & TOP_MASK)
        {
            /* (mb_x+left,mb_y) unavailable,
               (mb_x,mb_y+top) available */
            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

            /* Copy (ref_xD,y0) to (x0,y0) and (x0+1,y0) */
            /* copy (ref_xD,y0+1) to (x0,y0+1) and (x0+1,y0+1) */
            pu1_tmp_src1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst1 = pu1_tmp_src1;
            pu1_tmp_src2 = pu1_tmp_src1 + DYADIC_REF_W_Y;
            pu1_tmp_dst2 = pu1_tmp_src2;

            pu1_tmp_dst1[i4_x0] = pu1_tmp_src1[i4_ref_xD];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src1[i4_ref_xD];
            pu1_tmp_dst2[i4_x0] = pu1_tmp_src2[i4_ref_xD];
            pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src2[i4_ref_xD];

            /* Copy Cb and Cr */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;

            pu1_tmp_src1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst1 = pu1_tmp_src1;
            pu1_tmp_src2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst2 = pu1_tmp_src2;

            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];
        }
        else if(u1_ny_avlblty & (TOP_MASK | LEFT_MASK))
        {
            /* (mb_x+left,mb_y) available,
               (mb_x,mb_y+top) available */
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_x0 = 9 + (i4_left << 3) + i4_left;

            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;

            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;

            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;

            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Call diagonal construction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvc_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                   DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                           DYADIC_REF_W_C);

            isvc_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                           DYADIC_REF_W_C);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_get_ref_layer_mbtype                               */
/*                                                                           */
/*  Description   : This function is used to find the mb type of the         */
/*                    corresponding MB in the reference layer                */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra samp context                  */
/*                  pi1_ref_mb_modes : ref mb modes buffer pointer           */
/*                  i4_ref_mode_stride : mb mode buffer stride               */
/*                  i4_x_ref : reference location X                          */
/*                  i4_y_ref : reference location Y                          */
/*                  pi4_mb_type : pointer to store the mb type               */
/*                  i4_chroma_flag : chroma flag                             */
/*  Globals       : none                                                     */
/*  Processing    : it derives the bit corresponding to reference MB and     */
/*                  stores the mbtype as INTRA if the bit is set             */
/*  Outputs       : none                                                     */
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
static WORD8 isvc_get_ref_layer_mbtype(WORD8 *pi1_ref_mb_modes, WORD32 *pi4_mb_type,
                                       WORD8 i1_curr_slice_id, WORD8 i1_cons_intr_samp_flag)
{
    WORD8 i1_intra_slice_id;
    WORD8 i1_mb_mode;

    i1_mb_mode = *pi1_ref_mb_modes;
    i1_mb_mode = (i1_mb_mode < 0) ? i1_mb_mode : SVC_EXTRACT_MB_MODE(*pi1_ref_mb_modes);

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        *pi4_mb_type = SVC_INTER_MB;
        i1_intra_slice_id = -1;
    }
    else
    {
        *pi4_mb_type = SVC_INTRA_MB;
        i1_intra_slice_id = i1_mb_mode;

        if(1 == i1_cons_intr_samp_flag)
        {
            if(i1_mb_mode != i1_curr_slice_id)
            {
                *pi4_mb_type = SVC_INTER_MB;
            }
        }
    }
    return i1_intra_slice_id;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_fill_non_ava_pixel                                  */
/*                                                                           */
/*  Description   :  This function does the core pixel level processing      */
/*                    while filling the non available pixel                  */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  i4_refarray_wd : width of the reference array            */
/*                  i4_refarray_ht : height of the reference array           */
/*                  ps_mb_coord  : current mb coord structure                */
/*                  i4_chroma_flag : chroam processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : based on the map buffer values the non available pixels  */
/*                   are filled using border extension algorithm             */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2009   vijayakumar          creation                        */
/*         07 03 2011   A.D.Almeida          Optimized the filling pixels    */
/*                                                                           */
/*****************************************************************************/
static void isvc_fill_non_avail_pixel(intra_samp_lyr_ctxt *ps_lyr_ctxt, UWORD8 *pu1_refarray_1,
                                      UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                      WORD32 i4_chroma_flag, UWORD8 u1_avail_map[4][4])
{
    WORD32 i4_x, i4_y;
    WORD32 i4_corner_pixel_available;

    seg_lookup_desc_t *ps_segments_x;
    seg_lookup_desc_t *ps_segments_y;
    seg_description_t *ps_seg_desc_x, *ps_seg_desc_y;
    seg_description_t *ps_seg_x_tmp, *ps_seg_y_tmp;
    UWORD8 u1_num_sgmts_x, u1_num_sgmts_y;

    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refmb_wd;
    WORD32 i4_refmb_ht;
    WORD32 i4_xr_index, i4_yr_index;
    WORD32 i4_j, i4_i;
    WORD32 i4_cur_x;
    UWORD32 u4_lookup_4bit, u4_lookup_5bit, u4_4thbit;
    WORD32 i4_pad_size;

    WORD32 i4_x_min;
    WORD32 i4_y_min;
    WORD32 i4_x_start_pos, i4_y_start_pos;

    UWORD8 *pu1_ref_idx_x, *pu1_ref_idx_y;

    PF_INTRA_SAMP_PADDING *pf_intra_samp_padding;
    PF_INTRA_SAMP_PADDING **pf_intra_samp_lookup;

    i4_x_offset = ps_lyr_ctxt->ps_offsets->i4_abscissa;
    i4_y_offset = ps_lyr_ctxt->ps_offsets->i4_ordinate;

    i4_refmb_wd = (MB_SIZE >> i4_chroma_flag) - 1;
    i4_refmb_ht = (MB_SIZE >> i4_chroma_flag) - 1;

    if(0 == i4_chroma_flag)
    {
        pf_intra_samp_lookup = gpf_lookup_fxns_luma;
    }
    else
    {
        pf_intra_samp_lookup = gpf_lookup_fxns_chroma;
    }

    i4_x_min = ps_lyr_ctxt->i2_x_min_pos;
    i4_y_min = ps_lyr_ctxt->i2_y_min_pos;

    i4_pad_size = 2 >> i4_chroma_flag;
    i4_x_start_pos = (i4_x_min - i4_pad_size);
    i4_y_start_pos = (i4_y_min - i4_pad_size);

    i4_xr_index = (i4_x_start_pos + i4_x_offset) & i4_refmb_wd;
    i4_yr_index = (i4_y_start_pos + i4_y_offset) & i4_refmb_ht;

    ps_segments_x = (ps_lyr_ctxt->as_seg_lookup_horz + i4_xr_index);
    ps_segments_y = (ps_lyr_ctxt->as_seg_lookup_vert + i4_yr_index);

    u1_num_sgmts_x = ps_segments_x->u1_num_segments;
    u1_num_sgmts_y = ps_segments_y->u1_num_segments;

    ps_seg_desc_x = ps_segments_x->s_segments;
    ps_seg_desc_y = ps_segments_y->s_segments;

    pu1_ref_idx_x = &(ps_lyr_ctxt->au1_refarray_x_idx[0]);
    pu1_ref_idx_y = &(ps_lyr_ctxt->au1_refarray_y_idx[0]);

    i4_cur_x = pu1_ref_idx_x[i4_x_start_pos];

    u4_4thbit = ps_segments_x->u4_start_pos;

    for(i4_j = 0; i4_j < u1_num_sgmts_y; i4_j++)
    {
        UWORD8 i4_idx_a, i4_idx_b;
        UWORD8 u1_seg_ht, u1_seg_wd;
        UWORD8 u1_mb_adjoin_x, u1_mb_adjoin_y;
        WORD8 i1_nearst_mb_bdry_x, i1_nearst_mb_bdry_y;
        UWORD32 u4_num_valid_segs;
        WORD32 i4_idx_a_plus_ny, i4_idx_b_plus_nx, i4_index;
        WORD8 i1_yd_index, i1_xd_index;

        ps_seg_y_tmp = &ps_seg_desc_y[i4_j];

        i4_y = i4_y_start_pos + ps_seg_y_tmp->u1_seg_off;
        u1_seg_ht = ps_seg_y_tmp->u1_seg_dim;
        i1_yd_index = ps_seg_y_tmp->i1_dist_idx;
        i1_nearst_mb_bdry_y = ps_seg_y_tmp->i1_nearst_mb_bdry;
        u1_mb_adjoin_y = ps_seg_y_tmp->u1_mb_adjoin;

        i4_idx_a = pu1_ref_idx_y[i4_y];

        i4_idx_a_plus_ny = (i4_idx_a + i1_nearst_mb_bdry_y);

        /* Pack the availabilities of the next three horizontal MBs in 3bit
           format and 4th bit indicating if the start position is greater
           than the mb_width/2
         */
        u4_lookup_4bit = u4_4thbit | u1_avail_map[i4_idx_a][i4_cur_x + 2] << 2 |
                         u1_avail_map[i4_idx_a][i4_cur_x + 1] << 1 |
                         u1_avail_map[i4_idx_a][i4_cur_x];

        u4_num_valid_segs = gu4_valid_segs_lookup[u4_lookup_4bit];

        i4_i = CLZ(~u4_num_valid_segs);
        u4_num_valid_segs <<= (i4_i + 1);

        for(; i4_i < u1_num_sgmts_x; i4_i++)
        {
            ps_seg_x_tmp = &ps_seg_desc_x[i4_i];

            i4_x = i4_x_start_pos + ps_seg_x_tmp->u1_seg_off;
            i4_idx_b = pu1_ref_idx_x[i4_x];

            u1_seg_wd = ps_seg_x_tmp->u1_seg_dim;
            i1_xd_index = ps_seg_x_tmp->i1_dist_idx;
            i1_nearst_mb_bdry_x = ps_seg_x_tmp->i1_nearst_mb_bdry;
            u1_mb_adjoin_x = ps_seg_x_tmp->u1_mb_adjoin;

            i4_idx_b_plus_nx = (i4_idx_b + i1_nearst_mb_bdry_x);

            /* Find the avalability of (x,y-Yd),(x-Xd,y),(x-Xd,y-Yd) and pack
                it to 3 bits.
             */
            u4_lookup_5bit = u1_avail_map[i4_idx_a_plus_ny][i4_idx_b_plus_nx] << 2 |
                             u1_avail_map[i4_idx_a_plus_ny][i4_idx_b] << 1 |
                             u1_avail_map[i4_idx_a][i4_idx_b_plus_nx] | u1_mb_adjoin_x |
                             u1_mb_adjoin_y;

            i4_corner_pixel_available = u1_avail_map[i4_idx_a_plus_ny][i4_idx_b_plus_nx];

            /* Use a function pointer table based on lookup to compute
                Left,Top,Bottom,Right,Diagonal padding.
             */
            pf_intra_samp_padding = pf_intra_samp_lookup[u4_lookup_5bit];

            if(pf_intra_samp_padding != NULL)
            {
                pf_intra_samp_padding(i4_x, i4_y, i1_xd_index, i1_yd_index, u1_seg_wd, u1_seg_ht,
                                      pu1_refarray_1, pu1_refarray_2, i4_refarray_stride,
                                      u1_mb_adjoin_x, u1_mb_adjoin_y, i4_corner_pixel_available);
            }

            /* increment to the next unavailable segment */
            i4_index = CLZ(~u4_num_valid_segs);
            u4_num_valid_segs <<= (i4_index + 1);
            i4_i += i4_index;
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_intra_resamp_generate_segment_lookup                */
/*                                                                           */
/*  Description   : This function generates segment lookup used to derive    */
/*                  segments which have to be be intra resampled             */
/*                                                                           */
/*  Inputs        : pv_lookup_table : look up table                          */
/*                  i4_dimension    : dimension of the block which is used in*/
/*                  resampling process.                                      */
/*                  i4_mb_size  : size of the mb                             */
/*  Globals       : None                                                     */
/*  Processing    : This function generates segment lookup used to derive    */
/*                  segments which have to be be intra resampled             */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues       : None                                                      */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         03 03 2011   A.D.Almeida     Creation                             */
/*                                                                           */
void isvc_intra_resamp_generate_segment_lookup(seg_lookup_desc_t *ps_seg_lookup_table,
                                               WORD32 i4_dimension, WORD32 i4_mb_size,
                                               WORD32 i4_shift_val)
{
    WORD32 i4_x;
    WORD32 i4_position, i4_dist_prev_mb, i4_dist_next_mb;
    UWORD8 u1_seg_dim;
    UWORD8 u1_num_sgmts;
    WORD32 i4_block_size = i4_mb_size >> 1;
    UWORD8 u1_offset = 0;
    seg_lookup_desc_t *ps_segments;
    seg_description_t *ps_seg_desc;

    memset(ps_seg_lookup_table, 0, i4_mb_size * sizeof(seg_lookup_desc_t));

    for(i4_x = 0; i4_x < i4_mb_size; i4_x++)
    {
        ps_segments = &ps_seg_lookup_table[i4_x];
        ps_seg_desc = ps_segments->s_segments;
        i4_position = i4_x;

        if(i4_x >= i4_block_size)
        {
            /* set the fourth bit so that later it can be directly OR ed */
            ps_segments->u4_start_pos = 8;
        }
        else
        {
            ps_segments->u4_start_pos = 0;
        }

        u1_num_sgmts = 0;
        u1_offset = 0;

        while(i4_position < (i4_x + i4_dimension))
        {
            /* check and fill the nearest mb boundry flag */
            if((i4_position & (i4_mb_size - 1)) < i4_block_size)
            {
                ps_seg_desc->i1_nearst_mb_bdry = -1;
            }
            else
            {
                ps_seg_desc->i1_nearst_mb_bdry = 1;
            }

            /* find the distance from the previous MB for start of segment*/
            i4_dist_prev_mb = (i4_position & (i4_mb_size - 1));

            ps_seg_desc->i1_dist_idx =
                ((i4_dist_prev_mb >= i4_mb_size >> 1) ? (i4_mb_size - i4_dist_prev_mb)
                                                      : -(i4_dist_prev_mb + 1));

            /* find the size of the segment */
            u1_seg_dim = (i4_block_size - (i4_position & (i4_block_size - 1)));
            i4_position += u1_seg_dim;
            if(i4_position > (i4_x + i4_dimension))
            {
                i4_position = (i4_x + i4_dimension);
                u1_seg_dim = (i4_position & (i4_block_size - 1));
            }

            /* find the distance from the next MB for end of segment */
            i4_dist_next_mb = (i4_position & (i4_mb_size - 1));

            ps_seg_desc->u1_seg_dim = u1_seg_dim;
            ps_seg_desc->u1_seg_off = u1_offset;

            /* check if the segment has a adjoining MB edge */
            if(i4_dist_prev_mb == 0)
            {
                if(0 == u1_num_sgmts)
                {
                    ps_seg_desc->u1_mb_adjoin = 0;
                }
                else
                {
                    ps_seg_desc->u1_mb_adjoin = 1 << i4_shift_val;
                }
            }
            else if(i4_dist_next_mb == 0)
            {
                if(i4_position == (i4_x + i4_dimension))
                {
                    ps_seg_desc->u1_mb_adjoin = 0;
                }
                else
                {
                    ps_seg_desc->u1_mb_adjoin = 1 << i4_shift_val;
                }
            }
            else
            {
                ps_seg_desc->u1_mb_adjoin = 0;
            }

            u1_offset += u1_seg_dim;
            u1_num_sgmts++;
            ps_seg_desc++;
        }

        /* fill the number of segments for this position */
        ps_segments->u1_num_segments = u1_num_sgmts;
    }
}

static void isvc_reflayer_construction(void *pv_intra_samp_ctxt, UWORD8 *pu1_inp_1,
                                       WORD32 i4_inp_stride, WORD32 i4_refarray_stride,
                                       mem_element_t *ps_ref_mb_mode_map, WORD32 i4_chroma_flag)
{
    WORD32 i4_x, i4_y;

    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD8 *pi1_ref_mb_modes, *pi1_ref_mb_modes_bkp_1;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    WORD32 i4_dummy;
    WORD32 i4_mb_ht, i4_mb_wd;

    /* 4x4 mb grid buffer to store the mb availablity */
    UWORD8 u1_map_buf[BLK_SIZE][BLK_SIZE];
    WORD32 i4_ref_wd;
    WORD32 i4_ref_ht;
    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refarray_wd;
    WORD32 i4_refarray_ht;
    WORD32 i4_mb_type;
    WORD8 i1_cons_intr_samp_flag;
    WORD8 i1_slice_id = 0;
    WORD32 i4_mb_wd_sft, i4_mb_ht_sft;

    WORD32 i4_unfill_check;
    UWORD8 *pu1_refarray_1, *pu1_refarray_2;

    UNUSED(i4_dummy);
    memset(&u1_map_buf[0][0], 0, 16);

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode_map->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode_map->i4_element_size;

    /* get the condtrained intra sampling flag */
    i1_cons_intr_samp_flag = ps_lyr_ctxt->i1_constrained_intra_rsmpl_flag;

    ASSERT(NULL != pi1_ref_mb_modes);

    {
        WORD32 i4_base_width = ps_lyr_ctxt->i4_ref_width;
        WORD32 i4_base_height = ps_lyr_ctxt->i4_ref_height;

        i4_ref_wd = i4_base_width >> i4_chroma_flag;
        i4_ref_ht = i4_base_height >> i4_chroma_flag;

        i4_mb_wd_sft = (MB_WIDTH_SHIFT - i4_chroma_flag);
        i4_mb_ht_sft = (MB_HEIGHT_SHIFT - i4_chroma_flag);
    }

    i4_x_offset = ps_lyr_ctxt->ps_offsets->i4_abscissa;
    i4_y_offset = ps_lyr_ctxt->ps_offsets->i4_ordinate;
    i4_refarray_wd = ps_lyr_ctxt->ps_ref_array_dims->i4_abscissa;
    i4_refarray_ht = ps_lyr_ctxt->ps_ref_array_dims->i4_ordinate;

    i4_mb_wd = (MB_SIZE >> i4_chroma_flag);
    i4_mb_ht = (MB_SIZE >> i4_chroma_flag);

    if(1 == i1_cons_intr_samp_flag)
    {
        WORD32 i4_x_min, i4_x_max;
        WORD32 i4_y_min, i4_y_max;

        i4_x_min = ps_lyr_ctxt->i2_x_min_pos;
        i4_x_max = ps_lyr_ctxt->i2_x_max_pos;
        i4_y_min = ps_lyr_ctxt->i2_y_min_pos;
        i4_y_max = ps_lyr_ctxt->i2_y_max_pos;

        i4_mb_type = SVC_INTER_MB;
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_mb_x, i4_mb_y;

            /* derive local varaibles */
            i4_y_ref = (i4_y_min + 1) + i4_y_offset;
            i4_x_ref = (i4_x_min + 1) + i4_x_offset;

            i4_mb_x = (i4_x_ref >> i4_mb_wd_sft);
            i4_mb_y = (i4_y_ref >> i4_mb_ht_sft);

            pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;

            /* get the location of the byte which has the current mb mode */
            pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_element_size);
            pi1_ref_mb_modes += (i4_mb_x * i4_element_size);
        }

        for(i4_y = (i4_y_min + 1); i4_y <= (i4_y_max - 1);)
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_distleftX, i4_rangeX;
            WORD32 i4_disttopY, i4_rangeY;

            i4_y_ref = (i4_y + i4_y_offset);
            i4_disttopY = (i4_y_ref) & (i4_mb_ht - 1);
            i4_rangeY = (i4_mb_ht - i4_disttopY);

            pi1_ref_mb_modes_bkp_1 = pi1_ref_mb_modes;

            for(i4_x = (i4_x_min + 1); i4_x <= (i4_x_max - 1);)
            {
                i4_x_ref = (i4_x + i4_x_offset);
                i4_distleftX = (i4_x_ref) & (i4_mb_wd - 1);
                i4_rangeX = (i4_mb_wd - i4_distleftX);

                /* get the referecne layer mb type */
                i1_slice_id =
                    isvc_get_ref_layer_mbtype(pi1_ref_mb_modes_bkp_1, &i4_mb_type, i1_slice_id, 0);
                if(SVC_INTRA_MB == i4_mb_type)
                {
                    break;
                }
                i4_x += i4_rangeX;
                pi1_ref_mb_modes_bkp_1 += i4_element_size;
            }

            if(SVC_INTRA_MB == i4_mb_type)
            {
                break;
            }

            i4_y += i4_rangeY;
            pi1_ref_mb_modes += (i4_ref_mode_stride * i4_element_size);
        }
    }
    else
    {
        i1_slice_id = -1;
    }

    i4_unfill_check = 0;

    /* --------------------------------------------------------------------- */
    /* Copying the data from recon buffer to refSample Array.
     */
    /* NOTE: The copying of the data from recon buffer to refSample Array    */
    /*       can be optimized by bring in data at N-MB level,thus taking	 */
    /*		 advantage of the overlapping data which now gets copied every
     * MB*/
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_x_ref_start, i4_x_ref_end;
        WORD32 i4_y_ref_start, i4_y_ref_end;
        WORD32 i4_rangeW, i4_rangeH;
        WORD32 i4_offset;
        UWORD8 *pu1_src, *pu1_dst;
        UWORD8 *pu1_dst1, *pu1_dst2;

        /* Copy (refW x refH) dimension into reference sample array */
        i4_x_ref_start = MAX(0, MIN((i4_ref_wd - 1), i4_x_offset));
        i4_x_ref_end = MAX(0, MIN((i4_ref_wd - 1), (i4_refarray_wd - 1) + i4_x_offset));
        i4_y_ref_start = MAX(0, MIN((i4_ref_ht - 1), i4_y_offset));
        i4_y_ref_end = MAX(0, MIN((i4_ref_ht - 1), (i4_refarray_ht - 1) + i4_y_offset));

        /* find the actual data to be copied */
        i4_rangeW = (i4_x_ref_end - i4_x_ref_start + 1);
        i4_rangeH = (i4_y_ref_end - i4_y_ref_start + 1);

        /* get the reconbuffer pointer and ref sample array pointer */
        i4_offset =
            (i4_x_ref_start - i4_x_offset) + ((i4_y_ref_start - i4_y_offset) * i4_refarray_stride);

        if(0 == i4_chroma_flag)
        {
            pu1_refarray_1 = ps_ctxt->pu1_refarray_buffer;
            pu1_refarray_2 = NULL;

            pu1_src = pu1_inp_1;
            pu1_dst = pu1_refarray_1 + i4_offset;

            /* Copy luma data into refsample array */
            isvc_copy_data(pu1_src, i4_inp_stride, pu1_dst, i4_refarray_stride, i4_rangeW,
                           i4_rangeH);
        }
        else
        {
            pu1_refarray_1 = ps_ctxt->pu1_refarray_buffer;
            pu1_refarray_2 = ps_ctxt->pu1_refarray_cb;

            pu1_src = pu1_inp_1;
            pu1_dst1 = pu1_refarray_1 + i4_offset;

            pu1_dst2 = pu1_refarray_2 + i4_offset;

            isvc_copy_data_semiplanr(pu1_src, i4_inp_stride, pu1_dst1, pu1_dst2, i4_refarray_stride,
                                     i4_rangeW, i4_rangeH);
        }
    }
    {
        WORD32 i4_i, i4_j;
        UWORD8 *pu1_ref_idx_x, *pu1_ref_idx_y;

        WORD32 i4_x_ref;
        WORD32 i4_y_ref;
        WORD32 i4_mb_x, i4_mb_y;

        i4_y_ref = i4_y_offset;
        i4_x_ref = i4_x_offset;

        i4_mb_x = (i4_x_ref >> i4_mb_wd_sft);
        i4_mb_y = (i4_y_ref >> i4_mb_ht_sft);

        pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;

        pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_element_size);
        pi1_ref_mb_modes += (i4_mb_x * i4_element_size);

        pu1_ref_idx_x = &(ps_lyr_ctxt->au1_refarray_x_idx[0]);
        pu1_ref_idx_y = &(ps_lyr_ctxt->au1_refarray_y_idx[0]);

        i4_j = 0;
        for(i4_y = 0; i4_y < i4_refarray_ht;)
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_distleftX, i4_rangeX;
            WORD32 i4_disttopY, i4_rangeY;

            i4_y_ref = i4_y + i4_y_offset;
            i4_disttopY = (i4_y_ref) & (i4_mb_ht - 1);
            i4_rangeY = (i4_mb_ht - i4_disttopY);

            memset(pu1_ref_idx_y, i4_j, i4_rangeY);
            pu1_ref_idx_y += i4_rangeY;

            i4_i = 0;
            pi1_ref_mb_modes_bkp_1 = pi1_ref_mb_modes;
            for(i4_x = 0; i4_x < i4_refarray_wd;)
            {
                i4_x_ref = i4_x + i4_x_offset;
                i4_distleftX = (i4_x_ref) & (i4_mb_wd - 1);
                i4_rangeX = (i4_mb_wd - i4_distleftX);

                if(0 == i4_j)
                {
                    memset(pu1_ref_idx_x, i4_i, i4_rangeX);
                    pu1_ref_idx_x += i4_rangeX;
                }

                isvc_get_ref_layer_mbtype(pi1_ref_mb_modes_bkp_1, &i4_mb_type, i1_slice_id,
                                          i1_cons_intr_samp_flag);

                if(SVC_INTRA_MB == i4_mb_type)
                {
                    u1_map_buf[i4_j][i4_i] = 1;
                    i4_dummy = 1;
                }
                else
                {
                    i4_unfill_check = 1;
                }

                i4_x = i4_x + i4_rangeX;
                i4_i++;
                pi1_ref_mb_modes_bkp_1 += i4_element_size;
            }
            i4_j++;
            i4_y = i4_y + i4_rangeY;
            pi1_ref_mb_modes += (i4_ref_mode_stride * i4_element_size);
        }
        ASSERT(1 == i4_dummy);
    }

    if(i4_unfill_check == 1)
    {
        isvc_fill_non_avail_pixel(ps_lyr_ctxt, pu1_refarray_1, pu1_refarray_2, i4_refarray_stride,
                                  i4_chroma_flag, u1_map_buf);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_interpolate_base_luma_dyadic                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  intra resampling for dyadic scaling ratios               */
/*  Inputs        : pu1_inp_buf : ptr to the 12x12 reference sample buffer   */
/*                  pi2_tmp_filt_buf : ptr to the 12x16 buffer to hold the   */
/*                  vertically interpolated data                             */
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
/*         03 12 2010   Nithya          creation                             */
/*                                                                           */
/*****************************************************************************/
void isvc_interpolate_base_luma_dyadic(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                       UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1, i4_samp_2, i4_samp_3;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp, *pu1_out;
    WORD16 *pi2_tmp;

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
    for(i4_x = 0; i4_x < 12; i4_x++)
    {
        /* y = 0, y_phase = 12 */
        i4_samp_0 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_2 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_3 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        /* since y_phase 12 for y = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_3;
        i4_rslt_1 += i4_samp_1 * i4_coeff_2;
        i4_rslt_1 += i4_samp_2 * i4_coeff_1;
        i4_rslt_1 += i4_samp_3 * i4_coeff_0;

        pi2_tmp[i4_x] = i4_rslt_1;
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 15; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = i4_samp_2;
            i4_samp_2 = i4_samp_3;
            i4_samp_3 = pu1_inp[i4_x];

            /* y_phase is 4 for odd values of y */
            /* and 12 for even values of y    */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;
            i4_rslt_1 += i4_samp_2 * i4_coeff_2;
            i4_rslt_1 += i4_samp_3 * i4_coeff_3;

            i4_rslt_2 = i4_samp_0 * i4_coeff_3;
            i4_rslt_2 += i4_samp_1 * i4_coeff_2;
            i4_rslt_2 += i4_samp_2 * i4_coeff_1;
            i4_rslt_2 += i4_samp_3 * i4_coeff_0;

            /* Storing the results */
            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;

            /* Incrementing the pointers */
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;
        }

        /* y = 15, y_phase = 4 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = i4_samp_2;
        i4_samp_2 = i4_samp_3;
        i4_samp_3 = pu1_inp[i4_x];

        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;
        i4_rslt_1 += i4_samp_2 * i4_coeff_2;
        i4_rslt_1 += i4_samp_3 * i4_coeff_3;

        pi2_tmp[i4_x] = i4_rslt_1;
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;
    }

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 16; i4_y++)
    {
        /* x = 0, x_phase = 12 */
        i4_samp_0 = *pi2_tmp++;
        i4_samp_1 = *pi2_tmp++;
        i4_samp_2 = *pi2_tmp++;
        i4_samp_3 = *pi2_tmp++;

        /* since x_phase 12 for x = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_3;
        i4_rslt_1 += i4_samp_1 * i4_coeff_2;
        i4_rslt_1 += i4_samp_2 * i4_coeff_1;
        i4_rslt_1 += i4_samp_3 * i4_coeff_0;
        i4_rslt_1 += 512;

        i4_rslt_1 >>= 10;
        pu1_out[0] = CLIPUCHAR(i4_rslt_1);

        for(i4_x = 1; i4_x < 15; i4_x += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = i4_samp_2;
            i4_samp_2 = i4_samp_3;
            i4_samp_3 = *pi2_tmp++;

            /* x_phase is 4 for odd values of x */
            /* and 12 for even values of x    */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;
            i4_rslt_1 += i4_samp_2 * i4_coeff_2;
            i4_rslt_1 += i4_samp_3 * i4_coeff_3;
            i4_rslt_1 += 512;

            i4_rslt_2 = i4_samp_0 * i4_coeff_3;
            i4_rslt_2 += i4_samp_1 * i4_coeff_2;
            i4_rslt_2 += i4_samp_2 * i4_coeff_1;
            i4_rslt_2 += i4_samp_3 * i4_coeff_0;
            i4_rslt_2 += 512;

            i4_rslt_1 >>= 10;
            i4_rslt_2 >>= 10;

            pu1_out[i4_x] = CLIPUCHAR(i4_rslt_1);
            pu1_out[i4_x + 1] = CLIPUCHAR(i4_rslt_2);
        }

        /* x = 15 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = i4_samp_2;
        i4_samp_2 = i4_samp_3;
        i4_samp_3 = *pi2_tmp++;

        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;
        i4_rslt_1 += i4_samp_2 * i4_coeff_2;
        i4_rslt_1 += i4_samp_3 * i4_coeff_3;
        i4_rslt_1 += 512;

        i4_rslt_1 >>= 10;
        pu1_out[i4_x] = CLIPUCHAR(i4_rslt_1);
        pu1_out += i4_out_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_vert_interpol_chroma_dyadic                        */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios  for */
/*                  chroma for the following ref_lyr_chroma_phase_y_plus1 and*/
/*                  chroma_phase_y_plus1:                                    */
/*                  ref_lyr    cur_lyr                                       */
/*                    0      0                                               */
/*                    1      0                                               */
/*                    1      1                                               */
/*                    1      2                                               */
/*                    2      1                                               */
/*                    2      2                                               */
/*  Inputs         : pu1_inp_buf : ptr to the 6x6 reference sample buffer    */
/*                   pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the    */
/*                   vertically interpolated data                            */
/*                   i4_phase_0 : y phase for even values of y               */
/*                   i4_phase_1 : y phase for odd values of y                */
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
/*         06 12 2010   Nithya          creation                             */
/*                                                                           */
/*****************************************************************************/
void isvc_vert_interpol_chroma_dyadic(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                      WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 16 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 16 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;

    /* Vertical interpolation */
    for(i4_x = 0; i4_x < 6; i4_x++)
    {
        /* y = 0, y_phase = phase_0 */
        i4_samp_0 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        /* since y_phase = phase_0 for y = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;

        pi2_tmp[i4_x] = i4_rslt_1;
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 7; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = pu1_inp[i4_x];

            /* y_phase is phase_1 for odd values of y */
            /* and phase_0 for even values of y      */
            i4_rslt_1 = i4_samp_0 * i4_coeff_2;
            i4_rslt_1 += i4_samp_1 * i4_coeff_3;

            i4_rslt_2 = i4_samp_0 * i4_coeff_0;
            i4_rslt_2 += i4_samp_1 * i4_coeff_1;

            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;
        }

        /* y = 7, y_phase = phase_1 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = pu1_inp[i4_x];

        i4_rslt_1 = i4_samp_0 * i4_coeff_2;
        i4_rslt_1 += i4_samp_1 * i4_coeff_3;

        pi2_tmp[i4_x] = i4_rslt_1;

        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvc_horz_interpol_chroma_dyadic                        */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  horizontal intra resampling for dyadic scaling ratios for*/
/*                  chroma with following ref_lyr_chroma_phase_x_plus1_flag  */
/*                  and chroma_phase_x_plus1_flag:                           */
/*                  ref_lyr    cur_lyr                                       */
/*                    0      0                                               */
/*                    1      0                                               */
/*                    1      1                                               */
/*  Inputs        : pi2_tmp_filt_buf : ptr to the 6x8 buffer containing the  */
/*                  vertically interpolated data                             */
/*                  pu1_out_buf : pointer to the output buffer               */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_phase_0 : x phase for even values of x                */
/*                  i4_phase_1 : x phase for odd values of x                 */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : resampled samples                                        */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 12 2010   Nithya          creation                   */
/*                                                                           */
/*****************************************************************************/
void isvc_horz_interpol_chroma_dyadic(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                      WORD32 i4_out_stride, WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_dst_stride;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 16 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 16 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_dst_stride = i4_out_stride;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y++)
    {
        /* x = 0, x_phase = phase_0 */
        i4_samp_0 = *pi2_tmp++;
        i4_samp_1 = *pi2_tmp++;

        /* since x_phase = phase_0 for x = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;

        /* Round to 8-bit value */
        i4_rslt_1 += 128;
        i4_rslt_1 >>= 8;

        pu1_out[0] = i4_rslt_1;

        for(i4_x = 1; i4_x < 7; i4_x += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = *pi2_tmp++;

            /* x_phase is phase_1 for odd values of x */
            /* and phase_0 for even values of x      */
            i4_rslt_1 = i4_samp_0 * i4_coeff_2;
            i4_rslt_1 += i4_samp_1 * i4_coeff_3;

            i4_rslt_2 = i4_samp_0 * i4_coeff_0;
            i4_rslt_2 += i4_samp_1 * i4_coeff_1;

            /* Rounding to 8-bit values */
            i4_rslt_1 += 128;
            i4_rslt_1 >>= 8;
            i4_rslt_2 += 128;
            i4_rslt_2 >>= 8;

            pu1_out[2 * i4_x] = i4_rslt_1;
            pu1_out[2 * (i4_x + 1)] = i4_rslt_2;
        }

        /* y = 7, y_phase = phase_1 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = *pi2_tmp++;

        /* since x_phase = phase_1 for x = 7 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_2;
        i4_rslt_1 += i4_samp_1 * i4_coeff_3;

        /* Round to 8-bit value */
        i4_rslt_1 += 128;
        i4_rslt_1 >>= 8;

        pu1_out[2 * 7] = i4_rslt_1;
        pu1_out += i4_dst_stride;
    }
}

static void isvc_interpolate_intra_base(void *pv_intra_samp_ctxt, UWORD8 *pu1_out,
                                        WORD32 i4_out_stride, WORD32 i4_refarray_wd,
                                        WORD32 i4_chroma_flag, WORD32 i4_refarray_flag)
{
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD32 i4_x, i4_y;
    UWORD8 *pu1_refarray;
    coordinates_t *ps_phase;

    WORD32 i4_temp_array_ht;
    WORD32 *pi4_interp_buff;
    WORD32 *pi4_interp_buff_temp;

    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;

    WORD32 i4_x_min, i4_x_max;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    if(0 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_buffer;
    }
    else
    {
        pu1_refarray = ps_ctxt->pu1_refarray_cb;
    }

    i4_mb_wd = MB_SIZE >> i4_chroma_flag;
    i4_mb_ht = MB_SIZE >> i4_chroma_flag;

    i4_x_min = ps_lyr_ctxt->i2_x_min_pos;
    i4_x_max = ps_lyr_ctxt->i2_x_max_pos;

    ps_phase = ps_lyr_ctxt->ps_phase;

    i4_temp_array_ht = i4_mb_ht;
    pi4_interp_buff = ps_ctxt->pi4_temp_interpolation_buffer;
    pi4_interp_buff_temp = pi4_interp_buff;

    for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
    {
        for(i4_x = (i4_x_min - 1); i4_x <= (i4_x_max + 2); i4_x++)
        {
            WORD32 i4_y_ref = ps_lyr_ctxt->pi4_ref_array_positions_y[i4_y];
            WORD32 i4_y_phase =
                ps_phase[(ps_lyr_ctxt->ps_mb_pos->i4_ordinate * i4_mb_ht + i4_y) % 3].i4_ordinate;
            UWORD8 *pu1_refarray_temp = pu1_refarray + i4_x + (i4_y_ref * i4_refarray_wd);

            if(0 == i4_chroma_flag)
            {
                *(pi4_interp_buff + i4_x) =
                    (g_ai1_interp_filter_luma[i4_y_phase]) *
                        (*(pu1_refarray_temp - i4_refarray_wd)) +

                    (g_ai1_interp_filter_luma[16 + i4_y_phase]) * (*(pu1_refarray_temp)) +

                    (g_ai1_interp_filter_luma[32 + i4_y_phase]) *
                        (*(pu1_refarray_temp + i4_refarray_wd)) +

                    (g_ai1_interp_filter_luma[48 + i4_y_phase]) *
                        (*(pu1_refarray_temp + (2 * i4_refarray_wd)));
            }
            else
            {
                *(pi4_interp_buff + i4_x) =
                    (g_au1_interp_filter_chroma[i4_y_phase]) * (*(pu1_refarray_temp)) +
                    (g_au1_interp_filter_chroma[16 + i4_y_phase]) *
                        (*(pu1_refarray_temp + i4_refarray_wd));
            }
        }

        pi4_interp_buff = pi4_interp_buff + i4_refarray_wd;
    }

    pi4_interp_buff = pi4_interp_buff_temp;

    for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            WORD32 i4_x_ref = ps_lyr_ctxt->pi4_ref_array_positions_y[i4_x];
            WORD32 i4_x_phase =
                ps_phase[(ps_lyr_ctxt->ps_mb_pos->i4_abscissa * MAX_REF_ARR_WD_HT + i4_x) % 3]
                    .i4_ordinate;

            pi4_interp_buff_temp = pi4_interp_buff + i4_x_ref;

            if(0 == i4_chroma_flag)
            {
                *(pu1_out + i4_x + (i4_y * i4_out_stride)) = CLIPUCHAR(
                    ((g_ai1_interp_filter_luma[i4_x_phase]) * (*(pi4_interp_buff_temp - 1)) +
                     (g_ai1_interp_filter_luma[16 + i4_x_phase]) * (*(pi4_interp_buff_temp)) +
                     (g_ai1_interp_filter_luma[32 + i4_x_phase]) * (*(pi4_interp_buff_temp + 1)) +
                     (g_ai1_interp_filter_luma[48 + i4_x_phase]) * (*(pi4_interp_buff_temp + 2)) +
                     512) >>
                    10);
            }
            else
            {
                *(pu1_out + (2 * i4_x) + (i4_y * i4_out_stride)) = CLIPUCHAR(
                    ((g_au1_interp_filter_chroma[i4_x_phase]) * (*(pi4_interp_buff_temp)) +
                     (g_au1_interp_filter_chroma[16 + i4_x_phase]) * (*(pi4_interp_buff_temp + 1)) +
                     512) >>
                    10);
            }
        }

        pi4_interp_buff = pi4_interp_buff + i4_refarray_wd;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvc_intra_samp_mb_dyadic                              */
/*                                                                           */
/*  Description   : MB level function which performs the intra resampling    */
/*                  of data of an MB (luma and chroma inclusive) for dyadic  */
/*                  scaling ratios                                           */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_mb_mode_map : ref layer mb mode map buff desc     */
/*                  ps_curr_luma : current layer out luma buffer desc        */
/*                  ps_curr_chroma : current layer out chroma buffer desc    */
/*                  x,y : current mb coorinate                               */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interpolation function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         07 12 2010   Nithya          creation                             */
/*                                                                           */
/*****************************************************************************/
void isvc_intra_samp_mb_dyadic(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                               mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode_map,
                               mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                               UWORD16 u2_mb_x, UWORD16 u2_mb_y,
                               WORD32 i4_scaled_ref_layer_left_offset,
                               WORD32 i4_scaled_ref_layer_top_offset)
{
    UWORD8 *pu1_inp_luma, *pu1_inp_chroma;
    UWORD8 *pu1_out_luma, *pu1_out_chroma;
    UWORD8 *pu1_out_cb, *pu1_out_cr;
    UWORD8 *pu1_refarray_luma, *pu1_refarray_cb, *pu1_refarray_cr;
    WORD16 *pi2_tmp_filt_buf;
    WORD32 i4_inp_luma_stride, i4_inp_chroma_stride;
    WORD32 i4_out_luma_stride, i4_out_chroma_stride;
    UWORD16 u2_mb_x_ref, u2_mb_y_ref;
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD32 i4_scaled_mb_x, i4_scaled_mb_y;
    WORD32 i4_top, i4_left;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    i4_scaled_mb_x = u2_mb_x - (i4_scaled_ref_layer_left_offset >> 4);
    i4_scaled_mb_y = u2_mb_y - (i4_scaled_ref_layer_top_offset >> 4);

    if(i4_scaled_mb_x & 0x1)
    {
        i4_left = 1;
    }
    else
    {
        i4_left = -1;
    }
    if(i4_scaled_mb_y & 0x1)
    {
        i4_top = 1;
    }
    else
    {
        i4_top = -1;
    }

    u2_mb_x_ref = (i4_scaled_mb_x >> 1);
    u2_mb_y_ref = (i4_scaled_mb_y >> 1);

    pu1_inp_luma = (UWORD8 *) ps_ref_luma->pv_buffer;
    pu1_inp_chroma = (UWORD8 *) ps_ref_chroma->pv_buffer;

    i4_inp_luma_stride = ps_ref_luma->i4_num_element_stride;
    i4_inp_chroma_stride = ps_ref_chroma->i4_num_element_stride;

    /* ------- Constructing refSampleArray ----------------------- */
    isvc_reflayer_construction_dyadic(pv_intra_samp_ctxt, ps_ref_mb_mode_map, pu1_inp_luma,
                                      pu1_inp_chroma, i4_inp_luma_stride, i4_inp_chroma_stride,
                                      i4_top, i4_left, u2_mb_x_ref, u2_mb_y_ref);

    /* --------------------------------------------------------------------- */
    /* LUMA INTERPOLATION                                     */
    /* --------------------------------------------------------------------- */
    pu1_refarray_luma = ps_ctxt->pu1_refarray_buffer;
    if(1 == i4_top)
    {
        pu1_refarray_luma += (DYADIC_REF_W_Y << 3);
    }
    if(1 == i4_left)
    {
        pu1_refarray_luma += 8;
    }
    pu1_out_luma = (UWORD8 *) ps_curr_luma->pv_buffer;
    i4_out_luma_stride = ps_curr_luma->i4_num_element_stride;
    pi2_tmp_filt_buf = (WORD16 *) ps_ctxt->pi4_temp_interpolation_buffer;

    ps_lyr_ctxt->pf_interpolate_luma(pu1_refarray_luma, pi2_tmp_filt_buf, pu1_out_luma,
                                     i4_out_luma_stride);

    /* --------------------------------------------------------------------- */
    /* CHROMA INTERPOLATION                                   */
    /* --------------------------------------------------------------------- */
    pu1_out_chroma = (UWORD8 *) ps_curr_chroma->pv_buffer;
    i4_out_chroma_stride = ps_curr_chroma->i4_num_element_stride;

    /* CB */
    pu1_out_cb = pu1_out_chroma;
    pu1_refarray_cb = ps_ctxt->pu1_refarray_cb;

    if(1 == i4_top)
    {
        pu1_refarray_cb += (DYADIC_REF_W_C << 2);
    }
    if(1 == i4_left)
    {
        pu1_refarray_cb += 4;
    }

    /* Vertical interpolation */
    ps_lyr_ctxt->pf_vert_interpol_chroma(pu1_refarray_cb, pi2_tmp_filt_buf,
                                         ps_lyr_ctxt->i4_y_phase_0, ps_lyr_ctxt->i4_y_phase_1);

    /* Horizontal interpolation */
    ps_lyr_ctxt->pf_horz_interpol_chroma(pi2_tmp_filt_buf, pu1_out_cb, i4_out_chroma_stride,
                                         ps_lyr_ctxt->i4_x_phase_0, ps_lyr_ctxt->i4_x_phase_1);

    /* CR */
    pu1_out_cr = pu1_out_chroma + 1;
    pu1_refarray_cr = ps_ctxt->pu1_refarray_cr;

    if(1 == i4_top)
    {
        pu1_refarray_cr += (DYADIC_REF_W_C << 2);
    }
    if(1 == i4_left)
    {
        pu1_refarray_cr += 4;
    }

    /* Vertical interpolation */
    ps_lyr_ctxt->pf_vert_interpol_chroma(pu1_refarray_cr, pi2_tmp_filt_buf,
                                         ps_lyr_ctxt->i4_y_phase_0, ps_lyr_ctxt->i4_y_phase_1);

    /* Horizontal interpolation */
    ps_lyr_ctxt->pf_horz_interpol_chroma(pi2_tmp_filt_buf, pu1_out_cr, i4_out_chroma_stride,
                                         ps_lyr_ctxt->i4_x_phase_0, ps_lyr_ctxt->i4_x_phase_1);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvc_intra_samp_mb                                     */
/*                                                                           */
/*  Description   : MB level function which performs the intra resampling    */
/*                  of data of an MB (luma and chroma inclusive)             */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_mb_mode_map : ref layer mb mode map buff desc     */
/*                  ps_curr_luma : current layer out luma buffer desc        */
/*                  ps_curr_chroma : current layer out chroma buffer desc    */
/*                  x,y : current mb coorinate                       */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interpolation function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         07 12 2010   Nithya          creation               */
/*                                                                           */
/*****************************************************************************/
void isvc_intra_samp_mb(void *pv_intra_samp_ctxt_luma, void *pv_intra_samp_ctxt_chroma,
                        mem_element_t *ps_ref_luma, mem_element_t *ps_ref_chroma,
                        mem_element_t *ps_ref_mb_mode_map, mem_element_t *ps_curr_luma,
                        mem_element_t *ps_curr_chroma)
{
    UWORD8 *pu1_inp_luma, *pu1_inp_chroma;
    UWORD8 *pu1_out_luma, *pu1_out_chroma;
    UWORD8 *pu1_out_cb, *pu1_out_cr;
    WORD32 i4_inp_luma_stride, i4_inp_chroma_stride;
    WORD32 i4_out_luma_stride, i4_out_chroma_stride;
    WORD32 i4_chroma_flag, i4_refarray_stride;

    intra_sampling_ctxt_t *ps_ctxt_luma;
    intra_sampling_ctxt_t *ps_ctxt_chroma;

    ps_ctxt_luma = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt_luma;
    ps_ctxt_chroma = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt_chroma;

    i4_refarray_stride = ps_ctxt_luma->i4_refarray_stride;

    pu1_inp_luma = (UWORD8 *) ps_ref_luma->pv_buffer;
    pu1_inp_chroma = (UWORD8 *) ps_ref_chroma->pv_buffer;

    i4_inp_luma_stride = ps_ref_luma->i4_num_element_stride;
    i4_inp_chroma_stride = ps_ref_chroma->i4_num_element_stride;

    pu1_out_luma = (UWORD8 *) ps_curr_luma->pv_buffer;
    i4_out_luma_stride = ps_curr_luma->i4_num_element_stride;

    i4_chroma_flag = 0;

    /* ------- Constructing refSampleArray ----------------------- */
    isvc_reflayer_construction(pv_intra_samp_ctxt_luma, pu1_inp_luma, i4_inp_luma_stride,
                               i4_refarray_stride, ps_ref_mb_mode_map, i4_chroma_flag);

    /* ---- Interpolation process for Intra_Base prediction	 ------ */
    isvc_interpolate_intra_base(pv_intra_samp_ctxt_luma, pu1_out_luma, i4_out_luma_stride,
                                i4_refarray_stride, i4_chroma_flag, 0);

    pu1_out_chroma = (UWORD8 *) ps_curr_chroma->pv_buffer;
    i4_out_chroma_stride = ps_curr_chroma->i4_num_element_stride;

    pu1_out_cb = pu1_out_chroma;
    pu1_out_cr = pu1_out_cb + 1;

    i4_refarray_stride = ps_ctxt_chroma->i4_refarray_stride;

    i4_chroma_flag = 1;

    /* ------- Constructing refSampleArray ----------------------- */
    isvc_reflayer_construction(pv_intra_samp_ctxt_chroma, pu1_inp_chroma, i4_inp_chroma_stride,
                               i4_refarray_stride, ps_ref_mb_mode_map, i4_chroma_flag);

    /* ---- Cb Interpolation process for Intra_Base prediction	 ------ */
    isvc_interpolate_intra_base(pv_intra_samp_ctxt_chroma, pu1_out_cb, i4_out_chroma_stride,
                                i4_refarray_stride, i4_chroma_flag, 0);

    /* ---- Cr Interpolation process for Intra_Base prediction	 ------ */
    isvc_interpolate_intra_base(pv_intra_samp_ctxt_chroma, pu1_out_cr, i4_out_chroma_stride,
                                i4_refarray_stride, i4_chroma_flag, 1);
}
