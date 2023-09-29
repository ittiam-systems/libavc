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
*  ih264e_mc.h
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

#ifndef _IH264E_MC_H_
#define _IH264E_MC_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_motion_comp_luma(process_ctxt_t *ps_proc, UWORD8 **pu1_pseudo_pred,
                             WORD32 *pi4_pseudo_pred_strd);

void ih264e_motion_comp_chroma(process_ctxt_t *ps_proc);

#endif /* _IH264E_MC_H_ */
