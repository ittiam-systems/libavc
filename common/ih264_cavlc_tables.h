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
******************************************************************************
* @file
*  ih264_cavlc_tables.h
*
* @brief
*  This file contains enumerations, macros and extern declarations of H264
*  cavlc tables
*
* @author
*  ittiam
*
* @remarks
*  none
******************************************************************************
*/

#ifndef _IH264_CAVLC_TABLES_H_
#define _IH264_CAVLC_TABLES_H_

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/
/**
******************************************************************************
 *  @brief  maximum zeros left
******************************************************************************
 */
#define MAX_ZERO_LEFT 6

/*****************************************************************************/
/* global declarations                                                       */
/*****************************************************************************/

extern const UWORD8 gu1_cbp_map_tables[48][2];
extern const UWORD8 gu1_code_coeff_token_table[3][4][16];
extern const UWORD8 gu1_size_coeff_token_table[3][4][16];
extern const UWORD8 gu1_code_coeff_token_table_chroma[4][4];
extern const UWORD8 gu1_size_coeff_token_table_chroma[4][4];
extern const UWORD8 gu1_threshold_vlc_level[6];
extern const UWORD8 gu1_size_zero_table[135];
extern const UWORD8 gu1_code_zero_table[135];
extern const UWORD8 gu1_size_zero_table_chroma[9];
extern const UWORD8 gu1_code_zero_table_chroma[9];
extern const UWORD8 gu1_index_zero_table[15];
extern const UWORD8 gu1_size_run_table[42];
extern const UWORD8 gu1_code_run_table[42];
extern const UWORD8 gu1_index_run_table[7];

#endif /* _IH264_CAVLC_TABLES_H_ */
