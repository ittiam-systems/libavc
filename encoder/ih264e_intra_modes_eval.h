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
*  ih264e_intra_modes_eval.h
*
* @brief
*  This file contains declarations of routines that perform rate distortion
*  analysis on a macroblock if coded as intra.
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_INTRA_MODES_EVAL_H_
#define _IH264E_INTRA_MODES_EVAL_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_derive_nghbr_avbl_of_mbs(process_ctxt_t *ps_proc_ctxt);

UWORD8 ih264e_derive_ngbr_avbl_of_mb_partitions(block_neighbors_t *s_ngbr_avbl,
                                                WORD8 i1_pel_pos_x,
                                                WORD8 i1_pel_pos_y);

void ih264e_evaluate_intra16x16_modes_for_least_cost_rdoptoff(process_ctxt_t *ps_proc);

void ih264e_evaluate_intra8x8_modes_for_least_cost_rdoptoff(process_ctxt_t *ps_proc);

void ih264e_evaluate_intra4x4_modes_for_least_cost_rdoptoff(process_ctxt_t *ps_proc);

void ih264e_evaluate_intra4x4_modes_for_least_cost_rdopton(process_ctxt_t *ps_proc);

void ih264e_evaluate_chroma_intra8x8_modes_for_least_cost_rdoptoff(process_ctxt_t *ps_proc);

void ih264e_evaluate_intra8x8_modes_for_least_cost_rdoptoff(process_ctxt_t *ps_proc_ctxt);

typedef void ih264e_evaluate_intra_modes_ft(UWORD8 *pu1_src,
                                            UWORD8 *pu1_ngbr_pels_i16,
                                            UWORD8 *pu1_dst,
                                            UWORD32 src_strd,
                                            UWORD32 dst_strd,
                                            WORD32 u4_n_avblty,
                                            UWORD32 *u4_intra_mode,
                                            WORD32 *pu4_sadmin,
                                            UWORD32 u4_valid_intra_modes);

typedef void ih264e_evaluate_intra_4x4_modes_ft(UWORD8 *pu1_src,
                                                UWORD8 *pu1_ngbr_pels,
                                                UWORD8 *pu1_dst,
                                                UWORD32 src_strd,
                                                UWORD32 dst_strd,
                                                WORD32 u4_n_avblty,
                                                UWORD32 *u4_intra_mode,
                                                WORD32 *pu4_sadmin,
                                                UWORD32 u4_valid_intra_modes,
                                                UWORD32 u4_lambda,
                                                UWORD32 u4_predictd_mode);

/* C Declarations */
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra16x16_modes;
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra_chroma_modes;
ih264e_evaluate_intra_4x4_modes_ft ih264e_evaluate_intra_4x4_modes;

/* A9 Declarations */
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra16x16_modes_a9q;
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra_chroma_modes_a9q;
ih264e_evaluate_intra_4x4_modes_ft ih264e_evaluate_intra_4x4_modes_a9q;

/* AV8 Declarations */
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra16x16_modes_av8;
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra_chroma_modes_av8;
ih264e_evaluate_intra_4x4_modes_ft ih264e_evaluate_intra_4x4_modes_av8;

/* SSSE3 Declarations */
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra16x16_modes_ssse3;
ih264e_evaluate_intra_modes_ft ih264e_evaluate_intra_chroma_modes_ssse3;
ih264e_evaluate_intra_4x4_modes_ft ih264e_evaluate_intra_4x4_modes_ssse3;

#endif /* _IH264E_INTRA_MODES_EVAL_H_ */
