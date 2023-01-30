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
*  isvce_mode_stat_visualiser.h
*
* @brief
*  Contains function declarations for function declared in
*  isvce_mode_stat_visualiser.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_MODE_STAT_VISUALISER_H_
#define _ISVCE_MODE_STAT_VISUALISER_H_
#if ENABLE_MODE_STAT_VISUALISER

#include <stdio.h>

#include "ih264_typedefs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"

typedef struct mode_stat_visualiser_t
{
    FILE *ps_output_file;

    yuv_buf_props_t s_frame_buf;

} mode_stat_visualiser_t;

extern UWORD32 isvce_get_msv_ctxt_size(UWORD32 u4_wd, UWORD32 u4_ht);

extern void isvce_msv_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_msv_ctxt_delete(mode_stat_visualiser_t *ps_mode_stat_visualiser);

extern void isvce_msv_get_input_frame(mode_stat_visualiser_t *ps_mode_stat_visualiser,
                                      isvce_inp_buf_t *ps_inp_buf);

extern void isvce_msv_dump_visualisation(mode_stat_visualiser_t *ps_mode_stat_visualiser);

extern void isvce_msv_set_mode(mode_stat_visualiser_t *ps_mode_stat_visualiser,
                               isvce_mb_info_t *ps_mb_info, coordinates_t *ps_mb_pos);
#endif

#endif
