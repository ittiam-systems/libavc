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
 *  isvc_mem_fns.c
 *
 * @brief
 *  Functions used for memory operations
 *
 * @author
 *  Ittiam
 *
 * @par List of Functions:
 *  isvc_memcpy()
 *  isvc_memcpy_mul_8()
 *  isvc_memset()
 *  isvc_memset_mul_8()
 *  isvc_memset_16bit()
 *  isvc_memset_16bit_mul_8()
 *  isvc_memory_alloc()
 *  isvc_memory_free()
 *
 * @remarks
 *  None
 *
 ******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/
/* System include files */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* User include files */
#include "ih264_typedefs.h"
#include "isvc_mem_fns.h"

/**
********************************************************************************
*  @brief  copies a 2d blk from one location to another
*
*  @param[out] pu1_dst : dst pointer
*
*  @param[in] i4_dst_stride: stride of destination
*
*  @param[in] pu1_src : src ptr
*
*  @param[in] i4_src_stride: stride of src
*
*  @param[in] i4_blk_wd : blk width
*
*  @param[in] i4_blk_ht : blk height
*
*  @return void
********************************************************************************
*/

void isvc_copy_2d(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 *pu1_src, WORD32 i4_src_stride,
                  WORD32 i4_blk_wd, WORD32 i4_blk_ht)
{
    WORD32 i;

    for(i = 0; i < i4_blk_ht; i++)
    {
        memmove(pu1_dst, pu1_src, i4_blk_wd * sizeof(pu1_dst[0]));

        pu1_dst += i4_dst_stride;
        pu1_src += i4_src_stride;
    }
}

/**
********************************************************************************
*  @brief  memsets a 2d blk
*
*  @param[out] pu1_dst : dst pointer
*
*  @param[in] i4_dst_stride: stride of destination
*
*  @param[in] i4_blk_wd : blk width
*
*  @param[in] i4_blk_ht : blk height
*
*  @return void
********************************************************************************
*/
void isvc_memset_2d(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 u1_val, WORD32 i4_blk_wd,
                    WORD32 i4_blk_ht)
{
    WORD32 i;

    for(i = 0; i < i4_blk_ht; i++)
    {
        memset(pu1_dst, u1_val, i4_blk_wd);

        pu1_dst += i4_dst_stride;
    }
}

/**
 *******************************************************************************
 *
 * @brief
 * Checks if any pixel in a block is non-zero
 *
 * @param[in] pu1_data
 *  UWORD8 pointer to the block to be checked
 *
 * @param[in] i4_data_strd
 *  Stride of data buffer
 *
 * @param[in] u4_wd
 *  Width of the block
 *
 * @param[in] u4_ht
 *  Height of the block
 *
 *******************************************************************************
 */
UWORD8 isvc_is_nonzero_blk(UWORD8 *pu1_data, WORD32 i4_data_strd, UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 i, j;

    for(i = 0; i < u4_ht; i++)
    {
        for(j = 0; j < u4_wd; j++)
        {
            if(pu1_data[j + i * i4_data_strd])
            {
                return 1;
            }
        }
    }

    return 0;
}
