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
*  isvc_structs.h
*
* @brief
*  Contains struct definition used for SVC
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVC_STRUCTS_H_
#define _ISVC_STRUCTS_H_

#include "ih264_typedefs.h"
#include "iv2.h"
#include "ih264_defs.h"
#include "ih264_structs.h"
#include "isvc_defs.h"

typedef struct buffer_container_t
{
    void *pv_data;

    WORD32 i4_data_stride;

} buffer_container_t;

typedef struct yuv_buf_props_t
{
    buffer_container_t as_component_bufs[NUM_COMPONENTS];

    IV_COLOR_FORMAT_T e_color_format;

    UWORD32 u4_width;

    UWORD32 u4_height;

    UWORD8 u1_bit_depth;
} yuv_buf_props_t;

typedef struct nal_unit_header_t
{
    UWORD8 u1_nal_ref_idc;

    UWORD8 u1_nal_unit_type;
} nal_unit_header_t;

typedef struct coordinates_t
{
    WORD32 i4_abscissa;

    WORD32 i4_ordinate;
} coordinates_t;

typedef struct svc_au_buf_t
{
    /* Array of structs that contain properties of the buffers used for storing */
    yuv_buf_props_t *ps_layer_yuv_buf_props;

    /* Temporal ID */
    WORD8 i1_temporal_id;

    /* Num Spatial Layers */
    UWORD8 u1_num_spatial_layers;

    /* Resolution ration b/w spatial layers */
    DOUBLE d_spatial_res_ratio;

    /* absolute value of POC */
    WORD32 i4_abs_poc;

    /* POC % MaxPicOrderCntLSB */
    WORD32 i4_poc_lsb;

    /* Lower 32 bits of time stamp */
    UWORD32 u4_timestamp_low;

    /* Higher 32 bits of time stamp */
    UWORD32 u4_timestamp_high;

    /* Is Pic used as refPic for future frames? */
    WORD32 i4_used_as_ref;

    /* frame_num in the slice header */
    WORD32 i4_frame_num;

    /*
     *  0: Top Field
     *  1: Bottom Field
     */
    WORD8 i1_field_type;

    /* buffer ID from frame buffer manager */
    WORD32 i4_buf_id;

} svc_au_buf_t;

typedef struct svc_nalu_ext_t
{
    nal_unit_header_t s_nalu_header;

    /* idr_flag */
    UWORD8 u1_idr_flag;

    /* priority_id (Range = [0, 63]) */
    UWORD8 u1_priority_id;

    /* no_inter_layer_pred_flag */
    UWORD8 u1_no_inter_layer_pred_flag;

    /* dependency_id (Range = [0, 7]) */
    UWORD8 u1_dependency_id;

    /* quality_id (Range = [0, 15]) */
    UWORD8 u1_quality_id;

    /* temporal_id (Range = [0, 7]) */
    UWORD8 u1_temporal_id;

    /* use_ref_base_pic_flag */
    UWORD8 u1_use_ref_base_pic_flag;

    /* discardable_flag */
    UWORD8 u1_discardable_flag;

    /* output_flag */
    UWORD8 u1_output_flag;

    /* reserved_three_2bits */
    UWORD8 u1_reserved_three_2bits;

} svc_nalu_ext_t;

typedef struct svc_vui_ext_t
{
    /* specifies the maximum layers in the SVC bitstream */
    UWORD32 u4_vui_ext_num_entries_minus1;

    /* specifies the dependency ID for each layer */
    UWORD8 u1_vui_ext_dependency_id[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the quality ID for each layer */
    UWORD8 u1_vui_ext_quality_id[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the temporal ID for each layer */
    UWORD8 u1_vui_ext_temporal_id[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the timing_info_present_flag value of the i-th sub-bitstream */
    UWORD8 u1_vui_ext_timing_info_present_flag[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the num_units_in_tick value of the i-th sub-bitstream */
    UWORD32 u4_vui_ext_num_units_in_tick[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the time_scale value of the i-th sub-bitstream */
    UWORD32 u4_vui_ext_time_scale[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the fixed_frame_rate_flag value of the i-th sub-bitstream */
    UWORD8 u1_vui_ext_fixed_frame_rate_flag[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the nal_hrd_parameters_present_flag value of the i-th */
    UWORD8 u1_vui_ext_nal_hrd_params_present_flag[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the vcl_hrd_parameters_present_flag value of the i-th */
    UWORD8 u1_vui_ext_vcl_hrd_params_present_flag[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the low_delay_hrd_flag value of the i-th sub-bitstream */
    UWORD8 u1_vui_ext_low_delay_hrd_flag[MAX_VUI_EXT_NUM_ENTRIES];

    /* specifies the pic_struct_present_flag value of the i-th sub-bitstream */
    UWORD8 u1_vui_ext_pic_struct_present_flag[MAX_VUI_EXT_NUM_ENTRIES];

} svc_vui_ext_t;

typedef struct sps_svc_ext_t
{
    /* inter_layer_deblocking_filter_control_present_flag */
    UWORD8 u1_inter_layer_deblocking_filter_control_present_flag;

    /* extended_spatial_scalability_idc */
    UWORD8 u1_extended_spatial_scalability_idc;

    /* chroma_phase_x_plus1_flag */
    UWORD8 u1_chroma_phase_x_plus1;

    /* chroma_phase_y_plus1 */
    UWORD8 u1_chroma_phase_y_plus1;

    /* seq_ref_layer_chroma_phase_x_plus1_flag */
    UWORD8 u1_seq_ref_layer_chroma_phase_x_plus1_flag;

    /* seq_ref_layer_chroma_phase_y_plus1 */
    UWORD8 u1_seq_ref_layer_chroma_phase_y_plus1;

    /* seq_scaled_ref_layer_left_offset */
    WORD32 i4_seq_scaled_ref_layer_left_offset;

    /* seq_scaled_ref_layer_top_offset */
    WORD32 i4_seq_scaled_ref_layer_top_offset;

    /* seq_scaled_ref_layer_right_offset */
    WORD32 i4_seq_scaled_ref_layer_right_offset;

    /* seq_scaled_ref_layer_bottom_offset */
    WORD32 i4_seq_scaled_ref_layer_bottom_offset;

    /* seq_tcoeff_level_prediction_flag */
    WORD8 i1_seq_tcoeff_level_prediction_flag;

    /* adaptive_tcoeff_level_prediction_flag */
    WORD8 i1_adaptive_tcoeff_level_prediction_flag;

    /* slice_header_restriction_flag */
    WORD8 i1_slice_header_restriction_flag;

} sps_svc_ext_t;

typedef struct subset_sps_t
{
    /* SPS structure */
    sps_t s_sps;

    /* Structure containing flags specific to SVC SPS */
    sps_svc_ext_t s_sps_svc_ext;

    /* svc_vui_parameters_present_flag */
    WORD8 i1_svc_vui_parameters_present_flag;

    svc_vui_ext_t s_svc_vui;

    /* additional_extension2_data_flag */
    WORD8 i1_additional_extension2_flag;

} subset_sps_t;

typedef struct svc_slice_header_t
{
    /* ref_layer_dq_id */
    UWORD32 u4_ref_layer_dq_id;

    /* disable_inter_layer_deblocking_filter_idc */
    UWORD32 u4_disable_inter_layer_deblocking_filter_idc;

    /* inter_layer_slice_alpha_c0_offset_div2 */
    WORD32 i4_inter_layer_slice_alpha_c0_offset_div2;

    /* inter_layer_slice_beta_offset_div2 */
    WORD32 i4_inter_layer_slice_beta_offset_div2;

    /* constrained_intra_resampling_flag */
    WORD8 i1_constrained_intra_resampling_flag;

    /* ref_layer_chroma_phase_x_plus1_flag */
    WORD8 i1_ref_layer_chroma_phase_x_plus1_flag;

    /* ref_layer_chroma_phase_y_plus1 */
    WORD8 i1_ref_layer_chroma_phase_y_plus1;

    /* scaled_ref_layer_left_offset */
    WORD32 i4_scaled_ref_layer_left;

    /* scaled_ref_layer_top_offset */
    WORD32 i4_scaled_ref_layer_top;

    /* scaled_ref_layer_right_offset */
    WORD32 i4_scaled_ref_layer_right;

    /* scaled_ref_layer_bottom_offset */
    WORD32 i4_scaled_ref_layer_bottom;

    /* slice_skip_flag */
    WORD8 i1_slice_skip_flag;

    /* num_mbs_in_slice_minus1 */
    UWORD32 u4_num_mbs_in_slice_minus1;

    /* adaptive_base_mode_flag */
    WORD8 i1_adaptive_base_mode_flag;

    /* default_base_mode_flag */
    WORD8 i1_default_base_mode_flag;

    /* adaptive_motion_prediction_flag */
    WORD8 i1_adaptive_motion_prediction_flag;

    /* default_motion_prediction_flag */
    WORD8 i1_default_motion_prediction_flag;

    /* adaptive_residual_prediction_flag */
    WORD8 i1_adaptive_residual_prediction_flag;

    /* default_residual_prediction_flag */
    WORD8 i1_default_residual_prediction_flag;

    /* tcoeff_level_prediction_flag */
    WORD8 i1_tcoeff_level_prediction_flag;

    /* scan_idx_start */
    UWORD32 u4_scan_idx_start;

    /* scan_idx_end */
    UWORD32 u4_scan_idx_end;

    WORD32 i4_store_ref_base_pic_flag;

    slice_header_t s_slice_header;
} svc_slice_header_t;

#endif
