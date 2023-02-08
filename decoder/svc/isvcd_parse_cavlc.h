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
 *  isvcd_parse_cavlc.h
 *
 * @brief
 *  Declaration of UVLC and CAVLC functions
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PARSE_CAVLC_H_
#define _ISVCD_PARSE_CAVLC_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "isvcd_structs.h"
#include "ih264d_cabac.h"

void isvcd_parse_bmb_ref_index_cavlc_range1(UWORD32 u4_num_part, dec_bit_stream_t *ps_bitstrm,
                                            WORD8 *pi1_ref_idx,
                                            UWORD32 u4_num_ref_idx_active_minus1,
                                            UWORD8 *pu1_motion_prediction_flag);

WORD32 isvcd_parse_bmb_ref_index_cavlc(UWORD32 u4_num_part, dec_bit_stream_t *ps_bitstrm,
                                       WORD8 *pi1_ref_idx, UWORD32 u4_num_ref_idx_active_minus1,
                                       UWORD8 *pu1_motion_prediction_flag);

WORD32 isvcd_parse_pmb_ref_index_cavlc(UWORD32 u4_num_part, dec_bit_stream_t *ps_bitstrm,
                                       WORD8 *pi1_ref_idx, UWORD32 u4_num_ref_idx_active_minus1,
                                       UWORD8 *pu1_motion_prediction_flag);

void isvcd_parse_pmb_ref_index_cavlc_range1(UWORD32 u4_num_part, dec_bit_stream_t *ps_bitstrm,
                                            WORD8 *pi1_ref_idx,
                                            UWORD32 u4_num_ref_idx_active_minus1,
                                            UWORD8 *pu1_motion_prediction_flag);

#endif /*_ISVCD_PARSE_CAVLC_H_ */