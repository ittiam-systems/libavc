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
 *  isvcd_vui.h
 *
 * @brief
 *  This file contains routines to parse SEI NAL's
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_VUI_H_
#define _ISVCD_VUI_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"

#define MAX_VUI_EXT_NUM_ENTRIES 1024

typedef struct
{
    UWORD32 u4_vui_ext_num_entries_minus1;
    UWORD8 u1_vui_ext_dependency_id[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_quality_id[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_temporal_id[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_timing_info_present_flag[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD32 u4_vui_ext_num_units_in_tick[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD32 u4_vui_ext_time_scale[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_fixed_frame_rate_flag[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_nal_hrd_params_present_flag[MAX_VUI_EXT_NUM_ENTRIES];
    hrd_t s_nal_hrd[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_vcl_hrd_params_present_flag[MAX_VUI_EXT_NUM_ENTRIES];
    hrd_t s_vcl_hrd[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_low_delay_hrd_flag[MAX_VUI_EXT_NUM_ENTRIES];
    UWORD8 u1_vui_ext_pic_struct_present_flag[MAX_VUI_EXT_NUM_ENTRIES];
} svc_vui_ext_t;

WORD32 isvcd_parse_vui_ext_parametres(svc_vui_ext_t *ps_svc_vui_ext, dec_bit_stream_t *ps_bitstrm);

#endif /*_ISVCD_VUI_H_ */
