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
 *  isvcd_mb_utils.h
 *
 * @brief
 *  Contains declarations of the utility functions needed to decode SVCD MB
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_MB_UTILS_H_
#define _ISVCD_MB_UTILS_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"

UWORD32 isvcd_get_mb_info_cabac_nonmbaff(dec_struct_t *ps_dec, const UWORD16 u2_cur_mb_address,
                                         dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mbskip);

UWORD32 isvcd_get_mb_info_cavlc_nonmbaff(dec_struct_t *ps_dec, const UWORD16 u2_cur_mb_address,
                                         dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mbskip_run);

#endif /* _ISVCD_MB_UTILS_H_ */