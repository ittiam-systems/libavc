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
* @file isvc_cabac_tables.h
*
* @brief
*  This file contains enumerations, macros and extern declarations of H264
*  cabac tables
*
* @author
*  Ittiam
*
* @remarks
*  none
******************************************************************************
*/

#ifndef _ISVC_CABAC_TABLES_H_
#define _ISVC_CABAC_TABLES_H_

#include "ih264_cabac_tables.h"
/**
******************************************************************************
*  @brief  max range of cabac contexts in H264 (0-459)
******************************************************************************
*/
#define NUM_SVC_CABAC_CTXTS 467

extern const UWORD32 (*gau4_isvc_cabac_table)[4];

/*****************************************************************************/
/* Cabac tables for context initialization depending upon type of Slice,     */
/* cabac init Idc value and Qp.                                              */
/*****************************************************************************/
extern const UWORD8 gau1_isvc_cabac_ctxt_init_table[NUM_CAB_INIT_IDC_PLUS_ONE][QP_RANGE]
                                                   [NUM_SVC_CABAC_CTXTS];

#endif
