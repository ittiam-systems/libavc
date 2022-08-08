/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

#ifndef _IMVC_STRUCTS_H_
#define _IMVC_STRUCTS_H_

#include "ih264_typedefs.h"
#include "imvc_defs.h"

typedef struct nalu_mvc_ext_t
{
    UWORD8 u1_non_idr_flag;

    UWORD8 u1_priority_id;

    UWORD16 u2_view_id;

    UWORD8 u1_temporal_id;

    UWORD8 u1_anchor_pic_flag;

    UWORD8 u1_inter_view_flag;

} nalu_mvc_ext_t;

typedef struct mvc_ivp_ref_data_t
{
    UWORD8 u1_num_refs;

    UWORD16 au2_ref_view_ids[MAX_NUM_IVP_REFS];

} mvc_ivp_ref_data_t;

typedef struct mvc_op_data_t
{
    UWORD8 u1_temporal_id;

    UWORD16 u2_num_ops;

    UWORD16 u2_num_target_views;

    UWORD16 au2_target_view_ids[MAX_NUM_VIEWS];

    /* Counter for num target views and views each target is dependent on */
    UWORD16 u2_num_views;

} mvc_op_data_t;

typedef struct mvc_level_info_t
{
    UWORD32 u4_level_idc;

    mvc_op_data_t as_mvc_op_data[MAX_NUM_OPERATING_POINTS];

} mvc_level_info_t;

typedef struct sps_mvc_ext_t
{
    UWORD16 u2_num_views;

    UWORD16 au2_view_ids[MAX_NUM_VIEWS];

    /* 0 => L0; 1 => L1 */
    mvc_ivp_ref_data_t as_anchor_ref_data[2][MAX_NUM_VIEWS];

    /* 0 => L0; 1 => L1 */
    mvc_ivp_ref_data_t as_non_anchor_ref_data[2][MAX_NUM_VIEWS];

    UWORD8 u1_num_level_values_signalled;

    mvc_level_info_t as_mvc_level_info[MAX_NUM_LEVEL_VALUES_SIGNALLED];

} sps_mvc_ext_t;

typedef struct mvc_vui_ext_t
{
    UWORD16 u2_vui_mvc_num_ops;

    UWORD8 u1_vui_mvc_temporal_id[MAX_NUM_OPERATING_POINTS];

    UWORD16 u2_vui_mvc_num_target_output_views[MAX_NUM_OPERATING_POINTS];

    UWORD16 u2_vui_mvc_view_id[MAX_NUM_OPERATING_POINTS][MAX_NUM_VIEWS];

    UWORD8 u1_vui_mvc_timing_info_present_flag[MAX_NUM_OPERATING_POINTS];

    UWORD32 u4_vui_mvc_num_units_in_tick[MAX_NUM_OPERATING_POINTS];

    UWORD32 u4_vui_mvc_time_scale[MAX_NUM_OPERATING_POINTS];

    UWORD8 u1_vui_mvc_fixed_frame_rate_flag[MAX_NUM_OPERATING_POINTS];

    UWORD8 u1_vui_mvc_nal_hrd_parameters_present_flag[MAX_NUM_OPERATING_POINTS];

    UWORD8 u1_vui_mvc_vcl_hrd_parameters_present_flag[MAX_NUM_OPERATING_POINTS];

    UWORD8 u1_vui_mvc_low_delay_hrd_flag[MAX_NUM_OPERATING_POINTS];

    UWORD8 u1_vui_mvc_pic_struct_present_flag[MAX_NUM_OPERATING_POINTS];

} mvc_vui_ext_t;

typedef struct buffer_container_t
{
    void *pv_data;

    WORD32 i4_data_stride;

} buffer_container_t;

typedef struct yuv_buf_props_t
{
    buffer_container_t as_component_bufs[NUM_COMPONENTS];

    UWORD8 u1_bit_depth;

    UWORD16 u2_width;

    UWORD16 u2_height;

} yuv_buf_props_t;

typedef struct iv_mvc_yuv_buf_t
{
    yuv_buf_props_t as_view_buf_props[MAX_NUM_VIEWS];

} iv_mvc_yuv_buf_t;

typedef struct coordinates_t
{
    WORD32 i4_abscissa;

    WORD32 i4_ordinate;
} coordinates_t;

typedef struct offsets_t
{
    UWORD16 u2_top_offset;

    UWORD16 u2_bottom_offset;

    UWORD16 u2_left_offset;

    UWORD16 u2_right_offset;
} offsets_t;

#endif
