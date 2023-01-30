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
*  isvce_mc.h
*
* @brief
*  This file contains declarations of routines that perform motion compensation
*  of luma and chroma macroblocks.
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/
#ifndef _ISVCE_MC_H_
#define _ISVCE_MC_H_

/**
******************************************************************************
*
* @brief
*  performs motion compensation for a luma mb for the given mv.
*
* @par Description
*  This routine performs motion compensation of an inter mb. When the inter
*  mb mode is P16x16, there is no need to copy 16x16 unit from reference buffer
*  to pred buffer. In this case the function returns pointer and stride of the
*  ref. buffer and this info is used in place of pred buffer else where.
*  In other cases, the pred buffer is populated via copy / filtering + copy
*  (q pel cases) and returned.
*
* @param[in] ps_proc
*  pointer to current proc ctxt
*
* @return  none
*
* @remarks Assumes half pel buffers for the entire frame are populated.
*
******************************************************************************
*/
extern void isvce_motion_comp_luma(isvce_process_ctxt_t *ps_proc, buffer_container_t *ps_pred);

/**
******************************************************************************
*
* @brief
*  performs motion compensation for chroma mb
*
* @par   Description
*  Copies a MB of data from the reference buffer (Full pel, half pel or q pel)
*  according to the motion vectors given
*
* @param[in] ps_proc
*  pointer to current proc ctxt
*
* @return  none
*
* @remarks Assumes half pel and quarter pel buffers for the entire frame are
*  populated.
******************************************************************************
*/
extern void isvce_motion_comp_chroma(isvce_process_ctxt_t *ps_proc, buffer_container_t *ps_pred);

#endif
