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
/**
*******************************************************************************
* @file
*  imvcd.h
*
* @brief
*  This file contains all the necessary structure and  enumeration
* definitions needed for the Application  Program Interface(API) of the
* Ittiam MVC Decoder
*
*******************************************************************************
*/

#ifndef _IMVCD_H_
#define _IMVCD_H_
#include <stdbool.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"

/* Extern functions */
extern IV_API_CALL_STATUS_T imvcd_api_function(iv_obj_t *ps_dec_hdl, void *pv_ip, void *pv_op);

/* Typedefs */
typedef enum IMVCD_CTL_SUB_CMDS
{
    IMVCD_CTL_SET_NUM_CORES  = IVD_CMD_CTL_CODEC_SUBCMD_START,
    IMVCD_CTL_SET_PROCESSOR  = IVD_CMD_CTL_CODEC_SUBCMD_START + 1,
    IMVCD_CTL_GET_VUI_PARAMS = IVD_CMD_CTL_CODEC_SUBCMD_START + 2,
    IMVCD_CTL_DEGRADE        = IVD_CMD_CTL_CODEC_SUBCMD_START + 3,

} IMVCD_CTL_SUB_CMDS;

typedef struct imvcd_create_ip_t
{
    ivd_create_ip_t s_ivd_ip;

} imvcd_create_ip_t;

typedef struct imvcd_create_op_t
{
    ivd_create_op_t s_ivd_op;

} imvcd_create_op_t;

typedef struct imvcd_delete_ip_t
{
    ivd_delete_ip_t s_ivd_ip;

} imvcd_delete_ip_t;

typedef struct imvcd_delete_op_t
{
    ivd_delete_op_t s_ivd_op;

} imvcd_delete_op_t;

typedef struct imvcd_video_decode_ip_t
{
    ivd_video_decode_ip_t s_ivd_ip;

} imvcd_video_decode_ip_t;

typedef struct imvcd_video_decode_op_t
{
    ivd_video_decode_op_t s_ivd_op;

    iv_yuv_buf_t *ps_view_disp_bufs;

} imvcd_video_decode_op_t;

typedef struct imvcd_set_config_ip_t
{
    ivd_ctl_set_config_ip_t s_ivd_ip;
} imvcd_set_config_ip_t;

typedef struct imvcd_set_config_op_t
{
    ivd_ctl_set_config_op_t s_ivd_op;
} imvcd_set_config_op_t;

typedef struct imvcd_set_num_cores_ip_t
{
    UWORD32 u4_size;

    IVD_API_COMMAND_TYPE_T e_cmd;

    IVD_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

    UWORD32 u4_num_cores;
} imvcd_set_num_cores_ip_t;

typedef struct imvcd_set_num_cores_op_t
{
    UWORD32 u4_size;

    UWORD32 u4_error_code;
} imvcd_set_num_cores_op_t;

typedef struct imvcd_set_arch_ip_t
{
    UWORD32 u4_size;

    IVD_API_COMMAND_TYPE_T e_cmd;

    IVD_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

    IVD_ARCH_T e_arch;

    IVD_SOC_T e_soc;

} imvcd_set_arch_ip_t;

typedef struct imvcd_set_arch_op_t
{
    UWORD32 u4_size;

    UWORD32 u4_error_code;
} imvcd_set_arch_op_t;

typedef struct imvcd_set_degrade_mode_ip_t
{
    UWORD32 u4_size;

    IVD_API_COMMAND_TYPE_T e_cmd;

    IVD_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

    /**
     * Pictures that are are degraded
     * 0 : No degrade
     * 1 : Only on non-reference frames
     * 2 : Use interval specified by u4_nondegrade_interval
     * 3 : All non-key frames
     * 4 : All frames
     */
    WORD32 i4_degrade_pics;

    /**
     * Interval for pictures which are completely decoded without any degradation
     */
    WORD32 i4_nondegrade_interval;

    /**
     * bit position (lsb is zero): Type of degradation
     * 1 : Disable deblocking
     * 2 : Faster inter prediction filters
     * 3 : Fastest inter prediction filters
     */
    WORD32 i4_degrade_type;

} imvcd_set_degrade_mode_ip_t;

typedef struct imvcd_set_degrade_mode_op_t
{
    UWORD32 u4_size;

    UWORD32 u4_error_code;
} imvcd_set_degrade_mode_op_t;

typedef struct imvcd_flush_dec_ip_t
{
    ivd_ctl_flush_ip_t s_ivd_ip;
} imvcd_flush_dec_ip_t;

typedef struct imvcd_flush_dec_op_t
{
    ivd_ctl_flush_op_t s_ivd_op;
} imvcd_flush_dec_op_t;

typedef struct imvcd_get_buf_info_ip_t
{
    ivd_ctl_getbufinfo_ip_t s_ivd_ip;
} imvcd_get_buf_info_ip_t;

typedef struct ivd_mvc_buf_info_t
{
    UWORD16 u2_num_views;
} ivd_mvc_buf_info_t;

typedef struct imvcd_get_buf_info_op_t
{
    ivd_ctl_getbufinfo_op_t s_ivd_op;

    ivd_mvc_buf_info_t s_mvc_buf_info;

} imvcd_get_buf_info_op_t;

typedef struct imvcd_get_vui_ip_t
{
    UWORD32 u4_size;

    IVD_API_COMMAND_TYPE_T e_cmd;

    IVD_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

} imvcd_get_vui_ip_t;

typedef struct imvcd_get_vui_op_t
{
    bool b_is_vui_available;

    UWORD32 u4_error_code;

    UWORD8 u1_aspect_ratio_idc;

    UWORD16 u2_sar_width;

    UWORD16 u2_sar_height;

    UWORD8 u1_overscan_appropriate_flag;

    UWORD8 u1_video_format;

    UWORD8 u1_video_full_range_flag;

    UWORD8 u1_colour_primaries;

    UWORD8 u1_tfr_chars;

    UWORD8 u1_matrix_coeffs;

    UWORD8 u1_cr_top_field;

    UWORD8 u1_cr_bottom_field;

    UWORD32 u4_num_units_in_tick;

    UWORD32 u4_time_scale;

    UWORD8 u1_fixed_frame_rate_flag;

    UWORD8 u1_nal_hrd_params_present;

    UWORD8 u1_vcl_hrd_params_present;

    UWORD8 u1_low_delay_hrd_flag;

    UWORD8 u1_pic_struct_present_flag;

    UWORD8 u1_bitstream_restriction_flag;

    UWORD8 u1_mv_over_pic_boundaries_flag;

    UWORD32 u4_max_bytes_per_pic_denom;

    UWORD32 u4_max_bits_per_mb_denom;

    UWORD32 u4_log2_max_mv_length_horz;

    UWORD32 u4_log2_max_mv_length_vert;

    UWORD32 u4_num_reorder_frames;

    UWORD32 u4_max_dec_frame_buffering;

} imvcd_get_vui_op_t;

#endif
