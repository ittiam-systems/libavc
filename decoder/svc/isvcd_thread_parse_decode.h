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
 *  isvcd_thread_parse_decode.h
 *
 * @brief
 *  Contains routines that for multi-thread decoder
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_THREAD_PARSE_DECPDE_H_
#define _ISVCD_THREAD_PARSE_DECPDE_H_
WORD32 isvcd_decode_recon_tfr_nmb_thread(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD8 u1_num_mbs,
                                         UWORD8 u1_num_mbs_next, UWORD8 u1_end_of_row);
void isvcd_decode_picture_thread(svc_dec_lyr_struct_t *ps_svc_lyr_dec);
WORD32 isvcd_decode_slice_thread(svc_dec_lyr_struct_t *ps_svc_lyr_dec);
void ih264d_compute_bs_non_mbaff_thread(dec_struct_t *ps_dec, dec_mb_info_t *ps_cur_mb_info,
                                        UWORD32 u4_mb_num);

#endif /* _ISVCD_THREAD_PARSE_DECPDE_H_ */
