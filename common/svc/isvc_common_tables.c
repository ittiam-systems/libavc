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
*  isvc_common_tables.c
*
* @brief
*  Contains common global tables
*
* @author
*  Harish M
*
* @par List of Functions:
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* User include files */
#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_macros.h"
#include "isvc_structs.h"
#include "ih264_common_tables.h"
#include "isvc_common_tables.h"

/*****************************************************************************/
/* Extern global definitions                                                 */
/*****************************************************************************/

/**
 ******************************************************************************
 * @brief  while encoding, basing on the input configuration parameters, the
 * the level of the bitstream is computed basing on the table below.
 * input  : table_idx
 * output : level_idc or cpb size
 * @remarks Table A-1 – level table limits
 ******************************************************************************
 */
const level_tables_t gas_isvc_lvl_tbl[16] = {
    {IH264_LEVEL_10, 1485, 99, 396, 64, 175, 64},
    {IH264_LEVEL_1B, 1485, 99, 396, 128, 350, 64},
    {IH264_LEVEL_11, 3000, 396, 900, 192, 500, 128},
    {IH264_LEVEL_12, 6000, 396, 2376, 384, 1000, 128},
    {IH264_LEVEL_13, 11880, 396, 2376, 768, 2000, 128},
    {IH264_LEVEL_20, 11880, 396, 2376, 2000, 2000, 128},
    {IH264_LEVEL_21, 19800, 792, 4752, 4000, 4000, 256},
    {IH264_LEVEL_22, 20250, 1620, 8100, 4000, 4000, 256},
    {IH264_LEVEL_30, 40500, 1620, 8100, 10000, 10000, 256},
    {IH264_LEVEL_31, 108000, 3600, 18000, 14000, 14000, 512},
    {IH264_LEVEL_32, 216000, 5120, 20480, 20000, 20000, 512},
    {IH264_LEVEL_40, 245760, 8192, 32768, 20000, 25000, 512},
    {IH264_LEVEL_41, 245760, 8192, 32768, 50000, 62500, 512},
    {IH264_LEVEL_42, 522240, 8704, 34816, 50000, 62500, 512},
    {IH264_LEVEL_50, 589824, 22080, 110400, 135000, 135000, 512},
    {IH264_LEVEL_51, 983040, 36864, 184320, 240000, 240000, 512},
};
