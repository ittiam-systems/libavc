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
* @file
*  isvce_core_coding.h
*
* @brief
*  This file contains extern declarations of core coding routines
*
* @author
*  ittiam
*
* @remarks
*  none
******************************************************************************
*/

#ifndef _ISVCE_CORE_CODING_H_
#define _ISVCE_CORE_CODING_H_

#include "isvce_structs.h"

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/**
******************************************************************************
*  @brief      Enable/Disable Hadamard transform of DC Coeff's
******************************************************************************
*/
#define DISABLE_DC_TRANSFORM 0
#define ENABLE_DC_TRANSFORM 1

/**
*******************************************************************************
*  @brief bit masks for DC and AC control flags
*******************************************************************************
*/

#define DC_COEFF_CNT_LUMA_MB 16
#define NUM_4X4_BLKS_LUMA_MB_ROW 4
#define NUM_LUMA4x4_BLOCKS_IN_MB 16
#define NUM_CHROMA4x4_BLOCKS_IN_MB 8

#define SIZE_4X4_BLK_HRZ TRANS_SIZE_4
#define SIZE_4X4_BLK_VERT TRANS_SIZE_4

#define CNTRL_FLAG_DC_MASK_LUMA 0x0000FFFF
#define CNTRL_FLAG_AC_MASK_LUMA 0xFFFF0000

#define CNTRL_FLAG_AC_MASK_CHROMA_U 0xF0000000
#define CNTRL_FLAG_DC_MASK_CHROMA_U 0x0000F000

#define CNTRL_FLAG_AC_MASK_CHROMA_V 0x0F000000
#define CNTRL_FLAG_DC_MASK_CHROMA_V 0x00000F00

#define CNTRL_FLAG_AC_MASK_CHROMA (CNTRL_FLAG_AC_MASK_CHROMA_U | CNTRL_FLAG_AC_MASK_CHROMA_V)
#define CNTRL_FLAG_DC_MASK_CHROMA (CNTRL_FLAG_DC_MASK_CHROMA_U | CNTRL_FLAG_DC_MASK_CHROMA_V)

#define CNTRL_FLAG_DCBLK_MASK_CHROMA 0x0000C000

/**
*******************************************************************************
*  @brief macros for transforms
*******************************************************************************
*/
#define DEQUEUE_BLKID_FROM_CONTROL(u4_cntrl, blk_lin_id) \
    {                                                    \
        blk_lin_id = CLZ(u4_cntrl);                      \
        u4_cntrl &= (0x7FFFFFFF >> blk_lin_id);          \
    };

#define IND2SUB_LUMA_MB(u4_blk_id, i4_offset_x, i4_offset_y) \
    {                                                        \
        i4_offset_x = (u4_blk_id % 4) << 2;                  \
        i4_offset_y = (u4_blk_id / 4) << 2;                  \
    }

#define IS_V_BLK(u4_blk_id) ((u4_blk_id) > 3)

#define IND2SUB_CHROMA_MB(u4_blk_id, i4_offset_x, i4_offset_y)        \
    {                                                                 \
        i4_offset_x = ((u4_blk_id & 0x1) << 3) + IS_V_BLK(u4_blk_id); \
        i4_offset_y = (u4_blk_id & 0x2) << 1;                         \
    }

/* Typedefs */

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

extern FT_CORE_CODING isvce_code_luma_intra_macroblock_16x16;

extern FT_CORE_CODING isvce_code_luma_intra_macroblock_4x4;

extern FT_CORE_CODING isvce_code_luma_intra_macroblock_4x4_rdopt_on;

extern FT_CORE_CODING isvce_code_chroma_intra_macroblock_8x8;

extern FT_CORE_CODING isvce_code_luma_inter_macroblock_16x16;

extern FT_CORE_CODING isvce_code_chroma_inter_macroblock_8x8;

#endif
