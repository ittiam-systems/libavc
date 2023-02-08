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
 *  isvcd_tables.h
 *
 * @brief
 *  Declaration of all tables used by h264 decoder
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_TABLES_H_
#define _ISVCD_TABLES_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_cabac.h"
#include "isvcd_cabac.h"

/*****************************************************************************/
/* Cabac tables for context initialization depending upon type of Slice,     */
/* cabac init Idc value and Qp.                                              */
/*****************************************************************************/
extern const UWORD8 gau1_isvcd_cabac_ctxt_init_table[NUM_CAB_INIT_IDC_PLUS_ONE][QP_RANGE]
                                                    [NUM_CABAC_CTXTS_SVC];

#endif /*TABLES_H*/
