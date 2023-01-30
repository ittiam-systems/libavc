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
*  isvce_deblk.h
*
* @brief
*  This file contains extern declarations of deblocking routines
*
* @author
*  ittiam
*
* @remarks
*  none
******************************************************************************
*/

#ifndef _ISVCE_DEBLK_H_
#define _ISVCE_DEBLK_H_

#include "ih264_typedefs.h"
#include "isvce_structs.h"

#define CSBP_LEFT_BLOCK_MASK 0x1111
#define CSBP_RIGHT_BLOCK_MASK 0x8888

#define NUM_EDGES_IN_MB 4

extern void isvce_compute_bs(isvce_process_ctxt_t *ps_proc, UWORD8 u1_inter_layer_deblk_flag);

extern void isvce_deblock_mb(isvce_process_ctxt_t *ps_proc, isvce_deblk_ctxt_t *ps_deblk,
                             UWORD8 u1_inter_layer_deblk_flag);

#endif
