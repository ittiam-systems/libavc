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
 *  isvcd_process_ebslice.h
 *
 * @brief
 *  Contains declarations of routines that decode an EB slice type
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PROCESS_EBSLICE_H_
#define _ISVCD_PROCESS_EBSLICE_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"

WORD32 isvcd_parse_ebslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_bslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_mv_pred_ref_tfr_nby2_ebmb(dec_struct_t *ps_svc_lyr_dec, UWORD32 u4_mb_idx,
                                       UWORD32 u4_num_mbs);

WORD32 isvcd_parse_bmb_non_direct_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                        dec_mb_info_t *ps_cur_mb_info,
                                        dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                                        UWORD32 u4_num_mbsNby2);

WORD32 isvcd_parse_bmb_non_direct_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                        dec_mb_info_t *ps_cur_mb_info,
                                        dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                                        UWORD32 u4_num_mbsNby2);

WORD32 isvcd_decode_spatial_direct(dec_struct_t *ps_svc_lyr_dec, UWORD8 u1_wd_x,
                                   dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mb_num);
#endif /* _ISVCD_PROCESS_EBSLICE_H_ */