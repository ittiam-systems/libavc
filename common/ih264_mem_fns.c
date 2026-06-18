/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
 *  ih264_mem_fns.c
 *
 * @brief
 *  Functions used for memory operations
 *
 * @author
 *  ittiam
 *
 * @par List of Functions:
 *  ih264_memset_16bit
 *  ih264_memset_16bit_mul_8
 *
 * @remarks
 *  none
 *
 ******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System Include Files */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* User Include Files */
#include "ih264_typedefs.h"

#include "ih264_debug.h"
#include "ih264_mem_fns.h"

/**
*******************************************************************************
*
* @brief memset of 16bit data of a 8,16 or 32 bytes
*
* @par Description
*  Does memset of 16bit data for 8,16 or 32 number of bytes
*
* @param[in] pu2_dst
*  pointer to the destination
*
* @param[in] value
*  value used for memset
*
* @param[in] num_words
*  number of words to set
*
* @returns none
*
* @remarks none
*
*******************************************************************************
*/
void ih264_memset_16bit(UWORD16 *pu2_dst, UWORD16 value, UWORD32 num_words)
{
    UWORD32 i;

    for(i = 0; i < num_words; i++)
    {
        *pu2_dst++ = value;
    }
}

void ih264_memset_16bit_mul_8(UWORD16 *pu2_dst,
                              UWORD16 value,
                              UWORD32 num_words)
{
    UWORD32 i;

    ASSERT(num_words % 8 == 0);
    for(i = 0; i < num_words; i++)
    {
        *pu2_dst++ = value;
    }
}

