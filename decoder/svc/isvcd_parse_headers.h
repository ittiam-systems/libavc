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
 *  isvcd_parse_headers.h
 *
 * @brief
 *  Contains declarations high level syntax[above slice] parsing routines
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PARSE_HEADERS_H_
#define _ISVCD_PARSE_HEADERS_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "isvcd_structs.h"

WORD32 isvcd_parse_subset_sps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm);

WORD32 isvcd_parse_nal_unit(svc_dec_lyr_struct_t *dec_hdl, UWORD8 u1_nal_ref_idc);

WORD32 isvcd_dec_ref_base_pic_marking(
    dec_ref_base_pic_marking_params_t *ps_ref_base_pic_marking_svc_ext,
    dec_bit_stream_t *ps_bitstrm);

WORD32 isvcd_parse_pps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm);

WORD32 isvcd_parse_sps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm);

#endif /*_ISVCD_PARSE_HEADERS_H_ */