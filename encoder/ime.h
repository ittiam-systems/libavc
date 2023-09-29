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
 *  ime.h
 *
 * @brief
 *  Contains declarations of me functions
 *
 * @author
 *  Ittiam
 *
 * @remarks
 *
 *******************************************************************************
 */

#ifndef _IME_H_
#define _IME_H_

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/**
******************************************************************************
 *  @brief      Number of iterations before exiting during diamond search
******************************************************************************
 */
#define NUM_LAYERS 16

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

extern void ime_diamond_search_16x16(me_ctxt_t *ps_me_ctxt, WORD32 i4_reflist);

extern void ime_evaluate_init_srchposn_16x16(me_ctxt_t *ps_me_ctxt,
                                             WORD32 i4_reflist);

extern void ime_full_pel_motion_estimation_16x16(me_ctxt_t *ps_me_ctxt,
                                                 WORD32 i4_ref_list);

extern void ime_sub_pel_motion_estimation_16x16(me_ctxt_t *ps_me_ctxt,
                                                WORD32 i4_reflist);

extern void ime_compute_skip_cost(me_ctxt_t *ps_me_ctxt,
                                  ime_mv_t *ps_skip_mv,
                                  mb_part_ctxt *ps_smb_part_info,
                                  UWORD32 u4_use_stat_sad,
                                  WORD32 i4_reflist,
                                  WORD32 is_slice_type_b);


#endif /* _IME_H_ */
