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
 *  isvcd_parse_slice.h
 *
 * @brief
 *  Contains routines that decode a Enhancement slice type
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PARSE_SLICE_H_
#define _ISVCD_PARSE_SLICE_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_tables.h"

WORD32 isvcd_parse_islice_data_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                     dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_islice_data_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                     dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_bmb_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                             UWORD32 u4_num_mbsNby2);

WORD32 isvcd_parse_bmb_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                             UWORD32 u4_num_mbsNby2);

WORD32 isvcd_parse_imb_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num,
                             UWORD8 u1_mb_type);

WORD32 isvcd_parse_eislice_data_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                      dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_eislice_data_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                      dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_imb_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_type);

WORD32 isvcd_parse_inter_slice_data_cavlc_enh_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_slice_params_t *ps_slice,
                                                  UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_inter_slice_data_cabac_enh_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_slice_params_t *ps_slice,
                                                  UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_inter_slice_data_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_slice_params_t *ps_slice,
                                          UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_inter_slice_data_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_slice_params_t *ps_slice,
                                          UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_eislice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_islice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_pslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice);

WORD32 isvcd_parse_decode_slice(UWORD8 u1_is_idr_slice, UWORD8 u1_nal_ref_idc,
                                svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_parse_slice_header(svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_set_default_slice_header_ext(svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_start_of_pic(svc_dec_lyr_struct_t *ps_svc_lyr_dec, WORD32 i4_poc,
                          pocstruct_t *ps_temp_poc, UWORD16 u2_frame_num, dec_pic_params_t *ps_pps);

WORD32 isvcd_parse_decode_slice_ext_nal(UWORD8 u1_is_idr_slice, UWORD8 u1_nal_ref_idc,
                                        svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_parse_interlayer_resamp_func_init(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                               UWORD16 u2_first_mb_in_slice);

#endif /* _ISVCD_PARSE_SLICE_H_ */
