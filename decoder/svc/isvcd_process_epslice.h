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
 *  isvcd_process_epslice.h
 *
 * @brief
 *  Contains declarations of routines that decode an EP slice type
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PROCESS_EPSLICE_H_
#define _ISVCD_PROCESS_EPSLICE_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"

WORD32 isvcd_parse_epslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_pmb_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                             UWORD32 u4_num_mbsNby2);

WORD32 isvcd_parse_pmb_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                             UWORD32 u4_num_mbsNby2);
WORD32 isvcd_process_ibl_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                            UWORD32 u4_mb_num, UWORD8 u1_inter_intra_mode);

WORD32 isvcd_process_residual_resample_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_mb_info_t *ps_cur_mb_info);

WORD32 isvcd_process_inter_mb_no_rsd_pred_non_target(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                     dec_mb_info_t *ps_cur_mb_info,
                                                     UWORD8 u1_inference_mode);

WORD32 isvcd_process_inter_mb_rsd_pred_non_target(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_mb_info_t *ps_cur_mb_info,
                                                  UWORD8 u1_inference_mode,
                                                  UWORD16 *pu2_res_luma_csbp);

WORD32 isvcd_process_inter_mb_rsd_pred_target_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mb_num,
                                                  UWORD8 u1_inference_mode,
                                                  UWORD16 *pu2_res_luma_csbp);

WORD32 isvcd_mv_pred_ref_tfr_nby2_epmb(dec_struct_t *ps_dec, UWORD32 u4_num_mbs,
                                       UWORD32 u4_num_mbsNby2);

WORD32 isvcd_decode_recon_tfr_nmb_non_base_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                               UWORD32 u4_mb_idx, UWORD32 u4_num_mbs,
                                               UWORD32 u4_num_mbs_next, UWORD32 u4_tfr_n_mb,
                                               UWORD32 u4_end_of_row);
WORD32 isvcd_decode_recon_tfr_nmb_base_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD32 u4_mb_idx,
                                           UWORD32 u4_num_mbs, UWORD32 u4_num_mbs_next,
                                           UWORD32 u4_tfr_n_mb, UWORD32 u4_end_of_row);
void isvcd_retrive_infer_mode_mv(svc_dec_lyr_struct_t *ps_svc_lyr_dec, mv_pred_t *ps_mvpred,
                                 UWORD8 u1_lx, UWORD8 u1_sub_mb_num);

void isvcd_update_intra_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                            dec_mb_info_t *ps_cur_mb_info);

void isvcd_update_ipcm_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info);

void isvcd_update_ibl_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_mb_info_t *ps_cur_mb_info);
void isvcd_update_inter_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                            dec_mb_info_t *ps_cur_mb_info,
                                            UWORD8 u1_inference_mode);
void isvcd_ii_pred_mb(void *pv_svc_dec, dec_mb_info_t *ps_cur_mb_info);

WORD32 isvcd_mark_err_slice_skip(svc_dec_lyr_struct_t *ps_svc_lyr_dec, WORD32 num_mb_skip,
                                 UWORD8 u1_is_idr_slice, UWORD16 u2_frame_num,
                                 pocstruct_t *ps_cur_poc, WORD32 prev_slice_err);

void isvcd_one_to_one(svc_dec_lyr_struct_t *ps_svc_lyr_dec, struct pic_buffer_t *ps_col_pic,
                      directmv_t *ps_direct, UWORD8 u1_wd_x, WORD32 u2_sub_mb_ofst,
                      dec_mb_info_t *ps_cur_mb_info);

WORD32 isvcd_process_ii_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                           dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num);
#endif /* _ISVCD_PROCESS_EPSLICE_H_ */