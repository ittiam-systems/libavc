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
 *  isvc_mem_fns_av8.c
 *
 * @brief
 *  armv8 variants of
 * functions used for memory operations
 *
 * *******************************************************************************
 */
#include <arm_neon.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "isvc_mem_fns.h"

void isvc_memset_2d_neon(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 u1_val, WORD32 i4_blk_wd,
                         WORD32 i4_blk_ht)
{
    if(i4_blk_wd == 4)
    {
        vst1_lane_u32((UWORD32 *) pu1_dst, vreinterpret_u32_u8(vdup_n_u8(u1_val)), 0);
        pu1_dst += i4_dst_stride;

        vst1_lane_u32((UWORD32 *) pu1_dst, vreinterpret_u32_u8(vdup_n_u8(u1_val)), 0);
        pu1_dst += i4_dst_stride;

        vst1_lane_u32((UWORD32 *) pu1_dst, vreinterpret_u32_u8(vdup_n_u8(u1_val)), 0);
        pu1_dst += i4_dst_stride;

        vst1_lane_u32((UWORD32 *) pu1_dst, vreinterpret_u32_u8(vdup_n_u8(u1_val)), 0);
    }
    else if(i4_blk_wd == 8)
    {
        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
        pu1_dst += i4_dst_stride;

        vst1_u8(pu1_dst, vdup_n_u8(u1_val));
    }
    else if((i4_blk_wd % 16 == 0) && (i4_blk_ht % 16 == 0))
    {
        WORD32 i, j;
        UWORD8 *pu1_dst_col_ptr, *pu1_dst_row_ptr;
        WORD32 i4_width_by_16 = i4_blk_wd / 16;
        WORD32 i4_height_by_16 = i4_blk_ht / 16;

        for(i = 0; i < i4_height_by_16; i++)
        {
            pu1_dst_row_ptr = pu1_dst + i * 16 * i4_dst_stride;
            for(j = 0; j < i4_width_by_16; j++)
            {
                pu1_dst_col_ptr = pu1_dst_row_ptr + (j << 4);

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
                pu1_dst_col_ptr += i4_dst_stride;

                vst1q_u8(&pu1_dst_col_ptr[0], vdupq_n_u8(u1_val));
            }
        }
    }
    else
    {
        WORD32 i;

        for(i = 0; i < i4_blk_ht; i++)
        {
            memset(pu1_dst, u1_val, i4_blk_wd);
            pu1_dst += i4_dst_stride;
        }
    }
}
