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
 *  isvcd_utils.h
 *
 * @brief
 *  Contains declaration of routines that handle of start and end of
 *    pic processing
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_UTILS_H_
#define _ISVCD_UTILS_H_

#include "ih264d_defs.h"
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"

WORD16 isvcd_allocate_dynamic_bufs(svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD16 isvcd_free_dynamic_bufs(svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_init_pic(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_frame_num, WORD32 i4_poc,
                      dec_pic_params_t *ps_pps);

WORD32 isvcd_decode_gaps_in_frame_num(dec_struct_t *ps_dec, UWORD16 u2_frame_num);

WORD32 isvcd_decode_pic_order_cnt(
    UWORD8 u1_is_idr_slice, UWORD32 u2_frame_num, pocstruct_t *ps_prev_poc, pocstruct_t *ps_cur_poc,
    dec_slice_params_t *ps_cur_slice, /*!< Pointer to current slice Params*/
    dec_pic_params_t *ps_pps, UWORD8 u1_nal_ref_idc, UWORD8 u1_bottom_field_flag,
    UWORD8 u1_field_pic_flag, WORD32 *pi4_poc, dec_struct_t *ps_dec);

void isvcd_init_dpb_ref_bufs(dec_struct_t *ps_dec);

#endif /*_ISVCD_UTILS_H_*/