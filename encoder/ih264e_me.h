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
 *  ih264e_me.h
 *
 * @brief
 *  Contains declarations of h264 me
 *
 * @author
 *  ittiam
 *
 * @remarks
 *
 *******************************************************************************
 */

#ifndef _IH264E_ME_H_
#define _IH264E_ME_H_


/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/**
******************************************************************************
 *  @brief     Skip Bias value for P slice
******************************************************************************
 */
#define SKIP_BIAS_P 0

/**
******************************************************************************
 *  @brief     Skip Bias value for B slice
******************************************************************************
 */
#define SKIP_BIAS_B 0


/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/

/**
 ******************************************************************************
 *  @brief      compute median of 3 elements (a, b, c) and store the output
 *  in to result. This is used for mv prediction
 ******************************************************************************
 */
#define MEDIAN(a, b, c, result) if (a > b){\
                                    if (b > c)\
                                        result = b;\
                                    else {\
                                        if (a > c)\
                                            result = c;\
                                        else \
                                            result = a;\
                                    }\
                                }\
                                else {\
                                    if (c > b)\
                                        result = b;\
                                    else {\
                                        if (c > a)\
                                            result = c;\
                                        else \
                                            result = a;\
                                    }\
                                }

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_init_mv_bits(me_ctxt_t *ps_me);

ih264e_skip_params_ft  ih264e_find_pskip_params;

ih264e_skip_params_ft  ih264e_find_pskip_params_me;

ih264e_skip_params_ft  ih264e_find_bskip_params;

ih264e_skip_params_ft  ih264e_find_bskip_params_me;

void ih264e_get_mv_predictor(enc_pu_t *ps_left_mb_pu, enc_pu_t *ps_top_row_pu,
                             enc_pu_mv_t *ps_pred_mv, WORD32 i4_ref_list);

ih264e_compute_me_ft  ih264e_compute_me_multi_reflist;

ih264e_compute_me_ft  ih264e_compute_me_single_reflist;

void ih264e_init_me(process_ctxt_t *ps_proc);

void ih264e_compute_me_nmb(process_ctxt_t *ps_proc, UWORD32 u4_nmb_count);

void ih264e_mv_pred(process_ctxt_t *ps_proc, WORD32 i4_reflist);

void ih264e_mv_pred_me(process_ctxt_t *ps_proc, WORD32 i4_ref_list);

#endif /* _IH264E_ME_H_ */
