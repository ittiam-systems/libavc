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
 * Function for copying to an interleaved destination
 *
 * @par Description:
 *    Copies the array of width 'wd' and height 'ht' from the  location pointed
 *    by 'src' to the location pointed by 'dst'
 *
 * @param[in] pu1_src
 *  UWORD8 pointer to the source
 *
 * @param[out] pu1_dst
 *  UWORD8 pointer to the destination
 *
 * @param[in] src_strd
 *  integer source stride
 *
 * @param[in] dst_strd
 *  integer destination stride
 *
 * @param[in] ht
 *  integer height of the array
 *
 * @param[in] wd
 *  integer width of the array
 *
 * @returns
 *
 * @remarks
 *  The alternate elements of src will be copied to alternate locations in dsr
 *  Other locations are not touched
 *
 *******************************************************************************
 */
void isvc_interleaved_copy(UWORD8 *pu1_src, UWORD8 *pu1_dst, WORD32 src_strd, WORD32 dst_strd,
                           WORD32 ht, WORD32 wd)
{
    WORD32 row, col;
    wd *= 2;

    for(row = 0; row < ht; row++)
    {
        for(col = 0; col < wd; col += 2)
        {
            pu1_dst[col] = pu1_src[col];
        }

        pu1_src += src_strd;
        pu1_dst += dst_strd;
    }
}

/**
 *******************************************************************************
 *
 * @brief
 * Function for copying to an interleaved destination
 *
 * @par Description:
 *    Copies the array of width 'wd' and height 'ht' from the  location pointed
 *    by 'src' to the location pointed by 'dst'
 *
 * @param[in] pu1_src
 *  UWORD8 pointer to the source
 *
 * @param[out] pu1_dst
 *  UWORD8 pointer to the destination
 *
 * @param[in] src_strd
 *  integer source stride
 *
 * @param[in] dst_strd
 *  integer destination stride
 *
 * @param[in] ht
 *  integer height of the array
 *
 * @param[in] wd
 *  integer width of the array
 *
 * @returns
 *
 * @remarks
 *  The alternate elements of src will be copied to alternate locations in dsr
 *  Other locations are not touched
 *
 *******************************************************************************
 */
void isvc_16bit_interleaved_copy(WORD16 *pi2_src, WORD16 *pi2_dst, WORD32 src_strd, WORD32 dst_strd,
                                 WORD32 ht, WORD32 wd)
{
    WORD32 row, col;
    wd *= 2;

    for(row = 0; row < ht; row++)
    {
        for(col = 0; col < wd; col += 2)
        {
            pi2_dst[col] = pi2_src[col];
        }

        pi2_src += src_strd;
        pi2_dst += dst_strd;
    }
}

/**
 *******************************************************************************
 *
 * @brief
 * Function for memsetting to an interleaved destination
 *
 * @par Description:
 *    Memsets the array of width 'wd' and height 'ht' pointed by 'src'
 *
 * @param[in] pu1_src
 *  UWORD8 pointer to the source
 *
 * @param[in] src_strd
 *  integer source stride
 *
 * @param[in] value
 *  Value to set
 *
 * @param[in] ht
 *  integer height of the array
 *
 * @param[in] wd
 *  integer width of the array
 *
 * @returns
 *
 * @remarks
 *  The alternate elements of src will be copied to alternate locations in dsr
 *  Other locations are not touched
 *
 *******************************************************************************
 */
void isvc_16bit_interleaved_memset(WORD16 *pi2_src, WORD32 i4_src_strd, WORD16 i2_value,
                                   WORD32 i4_wd, WORD32 i4_ht)
{
    WORD32 row, col;

    i4_wd *= 2;

    for(row = 0; row < i4_ht; row++)
    {
        for(col = 0; col < i4_wd; col += 2)
        {
            pi2_src[col] = i2_value;
        }

        pi2_src += i4_src_strd;
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
