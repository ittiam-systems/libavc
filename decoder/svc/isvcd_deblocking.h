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
 *  isvcd_deblocking.h
 *
 * @brief
 *  Declarations of deblocking functions
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_DEBLOCKING_H_
#define _ISVCD_DEBLOCKING_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"

void isvcd_compute_bs_non_mbaff(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                                const UWORD16 u2_mbxn_mb);

void isvcd_compute_bs_non_mbaff_target_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info, const UWORD16 u2_mbxn_mb);

void isvcd_compute_bs_non_mbaff_medial_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info, const UWORD16 u2_mbxn_mb);

void isvcd_compute_bs_non_mbaff_target_lyr_no_inter_layer(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                          dec_mb_info_t *ps_cur_mb_info,
                                                          const UWORD16 u2_mbxn_mb);

#endif /*_ISVCD_DEBLOCKING_H_ */