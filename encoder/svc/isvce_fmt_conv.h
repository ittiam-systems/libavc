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
*  ih264e_fmt_conv.h
*
* @brief
*  The file contains extern declarations of color space conversion routines
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_FMT_CONV_H_
#define _ISVCE_FMT_CONV_H_

#include "ih264e_fmt_conv.h"
#include "isvce_structs.h"

IH264E_ERROR_T isvce_fmt_conv(isvce_codec_t *ps_codec, svc_au_buf_t *ps_pic, UWORD8 *pu1_y_dst,
                              UWORD8 *pu1_u_dst, UWORD8 *pu1_v_dst, UWORD32 u4_dst_y_strd,
                              UWORD32 u4_dst_uv_strd, WORD32 cur_row, WORD32 num_rows);

#endif
