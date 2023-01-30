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
/*****************************************************************************/
/*                                                                           */
/*  File Name         : isvce.h                                              */
/*                                                                           */
/*  Description       : This file contains all the necessary structure and   */
/*                      enumeration definitions needed for the Application   */
/*                      Program Interface(API) of the Ittiam SVC Encoder     */
/*                                                                           */
/*  List of Functions : isvce_api_function                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _ISVCE_H_
#define _ISVCE_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>

#include "iv2.h"
#include "ive2.h"

    /*****************************************************************************/
    /* Enums                                                                     */
    /*****************************************************************************/
    typedef enum ISVCE_API_COMMAND_TYPE_T
    {
        ISVCE_CMD_VIDEO_NA = 0x7FFFFFFF,
        ISVCE_CMD_GET_NUM_MEM_REC = 0x0,
        ISVCE_CMD_FILL_NUM_MEM_REC = 0x1,
        ISVCE_CMD_RETRIEVE_MEMREC = 0x2,
        ISVCE_CMD_INIT = 0x3,
        ISVCE_CMD_EXTENSIONS = 0x100,
        ISVCE_CMD_VIDEO_CTL,
        ISVCE_CMD_VIDEO_ENCODE
    } ISVCE_API_COMMAND_TYPE_T;

    typedef enum ISVCE_CONTROL_API_COMMAND_TYPE_T
    {
        ISVCE_CMD_CT_NA = 0x7FFFFFFF,
        ISVCE_CMD_CTL_SETDEFAULT = 0x0,
        ISVCE_CMD_CTL_SET_DIMENSIONS = 0x1,
        ISVCE_CMD_CTL_SET_FRAMERATE = 0x2,
        ISVCE_CMD_CTL_SET_BITRATE = 0x3,
        ISVCE_CMD_CTL_SET_FRAMETYPE = 0x4,
        ISVCE_CMD_CTL_SET_QP = 0x5,
        ISVCE_CMD_CTL_SET_ENC_MODE = 0x6,
        ISVCE_CMD_CTL_SET_VBV_PARAMS = 0x7,
        ISVCE_CMD_CTL_SET_AIR_PARAMS = 0x8,
        ISVCE_CMD_CTL_SET_ME_PARAMS = 0X9,
        ISVCE_CMD_CTL_SET_GOP_PARAMS = 0XA,
        ISVCE_CMD_CTL_SET_PROFILE_PARAMS = 0XB,
        ISVCE_CMD_CTL_SET_DEBLOCK_PARAMS = 0XC,
        ISVCE_CMD_CTL_SET_IPE_PARAMS = 0XD,
        ISVCE_CMD_CTL_SET_VUI_PARAMS = 0XE,
        ISVCE_CMD_CTL_SET_NUM_CORES = 0x30,
        ISVCE_CMD_CTL_RESET = 0xA0,
        ISVCE_CMD_CTL_FLUSH = 0xB0,
        ISVCE_CMD_CTL_GETBUFINFO = 0xC0,
        ISVCE_CMD_CTL_GETVERSION = 0xC1,
        ISVCE_CMD_CTL_SET_SEI_MDCV_PARAMS = 0xD0,
        ISVCE_CMD_CTL_SET_SEI_CLL_PARAMS = 0xD1,
        ISVCE_CMD_CTL_SET_SEI_AVE_PARAMS = 0xD2,
        ISVCE_CMD_CTL_SET_SEI_CCV_PARAMS = 0xD3,
        ISVCE_CMD_CTL_GET_ENC_FRAME_DIMENSIONS = 0xE1
    } ISVCE_CONTROL_API_COMMAND_TYPE_T;

    /*****************************************************************************/
    /* Extended Structures                                                       */
    /*****************************************************************************/

    /*****************************************************************************/
    /*  Get Number of Memory Records                                             */
    /*****************************************************************************/
    typedef struct svc_inp_params_t
    {
        /**
         * Num Temporal Layers
         */
        UWORD8 u1_num_temporal_layers;

        /**
         * Num Spatial Layers
         */
        UWORD8 u1_num_spatial_layers;

        /**
         * Resolution ration b/w spatial layers
         */
        DOUBLE d_spatial_res_ratio;

    } svc_inp_params_t;

    typedef struct isvce_num_mem_rec_ip_t
    {
        iv_num_mem_rec_ip_t s_ive_ip;
    } isvce_num_mem_rec_ip_t;

    typedef struct isvce_num_mem_rec_op_t
    {
        iv_num_mem_rec_op_t s_ive_op;
    } isvce_num_mem_rec_op_t;

    /*****************************************************************************/
    /*  Fill Memory Records                                                      */
    /*****************************************************************************/

    typedef struct isvce_fill_mem_rec_ip_t
    {
        iv_fill_mem_rec_ip_t s_ive_ip;

        svc_inp_params_t s_svc_inp_params;

        UWORD32 u4_wd;

        UWORD32 u4_ht;

    } isvce_fill_mem_rec_ip_t;

    typedef struct isvce_fill_mem_rec_op_t
    {
        iv_fill_mem_rec_op_t s_ive_op;
    } isvce_fill_mem_rec_op_t;

    /*****************************************************************************/
    /*  Retrieve Memory Records                                                  */
    /*****************************************************************************/

    typedef struct isvce_retrieve_mem_rec_ip_t
    {
        iv_retrieve_mem_rec_ip_t s_ive_ip;
    } isvce_retrieve_mem_rec_ip_t;

    typedef struct isvce_retrieve_mem_rec_op_t
    {
        iv_retrieve_mem_rec_op_t s_ive_op;
    } isvce_retrieve_mem_rec_op_t;

    /*****************************************************************************/
    /*   Initialize encoder                                                      */
    /*****************************************************************************/

    typedef struct isvce_init_ip_t
    {
        ive_init_ip_t s_ive_ip;

        svc_inp_params_t s_svc_inp_params;

        UWORD32 *pu4_max_bitrate;

        UWORD32 u4_wd;

        UWORD32 u4_ht;

        bool b_use_default_vui;

        bool b_nalu_info_export_enable;

    } isvce_init_ip_t;

    typedef struct isvce_init_op_t
    {
        ive_init_op_t s_ive_op;
    } isvce_init_op_t;

    /*****************************************************************************/
    /*   Video control  Flush                                                    */
    /*****************************************************************************/

    typedef struct isvce_ctl_flush_ip_t
    {
        ive_ctl_flush_ip_t s_ive_ip;
    } isvce_ctl_flush_ip_t;

    typedef struct isvce_ctl_flush_op_t
    {
        ive_ctl_flush_op_t s_ive_op;
    } isvce_ctl_flush_op_t;

    /*****************************************************************************/
    /*   Video control reset                                                     */
    /*****************************************************************************/

    typedef struct isvce_ctl_reset_ip_t
    {
        ive_ctl_reset_ip_t s_ive_ip;
    } isvce_ctl_reset_ip_t;

    typedef struct isvce_ctl_reset_op_t
    {
        ive_ctl_reset_op_t s_ive_op;
    } isvce_ctl_reset_op_t;

    /*****************************************************************************/
    /*   Video control:Get Buf Info                                              */
    /*****************************************************************************/

    typedef struct isvce_ctl_getbufinfo_ip_t
    {
        ive_ctl_getbufinfo_ip_t s_ive_ip;
    } isvce_ctl_getbufinfo_ip_t;

    typedef struct isvce_ctl_getbufinfo_op_t
    {
        ive_ctl_getbufinfo_op_t s_ive_op;

        UWORD32 au4_min_rec_buf_size[IVE_MAX_IO_BUFFER_COMPONENTS];

        UWORD32 u4_rec_comp_cnt;

        UWORD32 u4_min_rec_bufs;

        UWORD32 u4_min_nalu_info_bufs;

        UWORD32 u4_min_nalu_info_buf_size;
    } isvce_ctl_getbufinfo_op_t;

    /*****************************************************************************/
    /*   Video control:Get Version Info                                          */
    /*****************************************************************************/

    typedef struct isvce_ctl_getversioninfo_ip_t
    {
        ive_ctl_getversioninfo_ip_t s_ive_ip;
    } isvce_ctl_getversioninfo_ip_t;

    typedef struct isvce_ctl_getversioninfo_op_t
    {
        ive_ctl_getversioninfo_op_t s_ive_op;
    } isvce_ctl_getversioninfo_op_t;

    /*****************************************************************************/
    /*   Video control:Set default params                                       */
    /*****************************************************************************/

    typedef struct isvce_ctl_setdefault_ip_t
    {
        ive_ctl_setdefault_ip_t s_ive_ip;
    } isvce_ctl_setdefault_ip_t;

    typedef struct isvce_ctl_setdefault_op_t
    {
        ive_ctl_setdefault_op_t s_ive_op;
    } isvce_ctl_setdefault_op_t;

    /*****************************************************************************/
    /*   Video control  Set IPE params                                           */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_ipe_params_ip_t
    {
        ive_ctl_set_ipe_params_ip_t s_ive_ip;
    } isvce_ctl_set_ipe_params_ip_t;

    typedef struct isvce_ctl_set_ipe_params_op_t
    {
        ive_ctl_set_ipe_params_op_t s_ive_op;
    } isvce_ctl_set_ipe_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set Frame dimensions                                     */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_dimensions_ip_t
    {
        ive_ctl_set_dimensions_ip_t s_ive_ip;
    } isvce_ctl_set_dimensions_ip_t;

    typedef struct isvce_ctl_set_dimensions_op_t
    {
        ive_ctl_set_dimensions_op_t s_ive_op;
    } isvce_ctl_set_dimensions_op_t;

    /*   Video control - Get Enc Frame dimensions */
    typedef struct isvce_ctl_get_enc_dimensions_ip_t
    {
        UWORD32 u4_inp_frame_wd;

        UWORD32 u4_inp_frame_ht;
    } isvce_ctl_get_enc_dimensions_ip_t;

    typedef struct isvce_ctl_get_enc_dimensions_op_t
    {
        UWORD32 u4_error_code;

        UWORD32 u4_enc_frame_wd;

        UWORD32 u4_enc_frame_ht;

    } isvce_ctl_get_enc_dimensions_op_t;

    /*****************************************************************************/
    /*   Video control  Set Frame rates                                          */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_frame_rate_ip_t
    {
        ive_ctl_set_frame_rate_ip_t s_ive_ip;
    } isvce_ctl_set_frame_rate_ip_t;

    typedef struct isvce_ctl_set_frame_rate_op_t
    {
        ive_ctl_set_frame_rate_op_t s_ive_op;
    } isvce_ctl_set_frame_rate_op_t;

    /*****************************************************************************/
    /*   Video control  Set Bitrate                                              */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_bitrate_ip_t
    {
        ive_ctl_set_bitrate_ip_t s_ive_ip;

        UWORD32 *pu4_target_bitrate;
    } isvce_ctl_set_bitrate_ip_t;

    typedef struct isvce_ctl_set_bitrate_op_t
    {
        ive_ctl_set_bitrate_op_t s_ive_op;
    } isvce_ctl_set_bitrate_op_t;

    /*****************************************************************************/
    /*   Video control  Set Frame type                                           */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_frame_type_ip_t
    {
        ive_ctl_set_frame_type_ip_t s_ive_ip;
    } isvce_ctl_set_frame_type_ip_t;

    typedef struct isvce_ctl_set_frame_type_op_t
    {
        ive_ctl_set_frame_type_op_t s_ive_op;
    } isvce_ctl_set_frame_type_op_t;

    /*****************************************************************************/
    /*   Video control  Set Encode mode                                          */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_enc_mode_ip_t
    {
        ive_ctl_set_enc_mode_ip_t s_ive_ip;
    } isvce_ctl_set_enc_mode_ip_t;

    typedef struct isvce_ctl_set_enc_mode_op_t
    {
        ive_ctl_set_enc_mode_op_t s_ive_op;
    } isvce_ctl_set_enc_mode_op_t;

    /*****************************************************************************/
    /*   Video control  Set QP                                                   */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_qp_ip_t
    {
        ive_ctl_set_qp_ip_t s_ive_ip;

        UWORD32 *pu4_i_qp;

        UWORD32 *pu4_i_qp_max;

        UWORD32 *pu4_i_qp_min;

        UWORD32 *pu4_p_qp;

        UWORD32 *pu4_p_qp_max;

        UWORD32 *pu4_p_qp_min;

        UWORD32 *pu4_b_qp;

        UWORD32 *pu4_b_qp_max;

        UWORD32 *pu4_b_qp_min;

    } isvce_ctl_set_qp_ip_t;

    typedef struct isvce_ctl_set_qp_op_t
    {
        ive_ctl_set_qp_op_t s_ive_op;
    } isvce_ctl_set_qp_op_t;

    /*****************************************************************************/
    /*   Video control  Set AIR params                                           */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_air_params_ip_t
    {
        ive_ctl_set_air_params_ip_t s_ive_ip;
    } isvce_ctl_set_air_params_ip_t;

    typedef struct isvce_ctl_set_air_params_op_t
    {
        ive_ctl_set_air_params_op_t s_ive_op;
    } isvce_ctl_set_air_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set VBV params                                           */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_vbv_params_ip_t
    {
        ive_ctl_set_vbv_params_ip_t s_ive_ip;

        UWORD32 *pu4_vbv_buffer_delay;
    } isvce_ctl_set_vbv_params_ip_t;

    typedef struct isvce_ctl_set_vbv_params_op_t
    {
        ive_ctl_set_vbv_params_op_t s_ive_op;
    } isvce_ctl_set_vbv_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set Processor Details                                    */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_num_cores_ip_t
    {
        ive_ctl_set_num_cores_ip_t s_ive_ip;
    } isvce_ctl_set_num_cores_ip_t;

    typedef struct isvce_ctl_set_num_cores_op_t
    {
        ive_ctl_set_num_cores_op_t s_ive_op;
    } isvce_ctl_set_num_cores_op_t;

    /*****************************************************************************/
    /*   Video control  Set Motion estimation params                             */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_me_params_ip_t
    {
        ive_ctl_set_me_params_ip_t s_ive_ip;
    } isvce_ctl_set_me_params_ip_t;

    typedef struct isvce_ctl_set_me_params_op_t
    {
        ive_ctl_set_me_params_op_t s_ive_op;
    } isvce_ctl_set_me_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set GOP params                                           */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_gop_params_ip_t
    {
        ive_ctl_set_gop_params_ip_t s_ive_ip;
    } isvce_ctl_set_gop_params_ip_t;

    typedef struct isvce_ctl_set_gop_params_op_t
    {
        ive_ctl_set_gop_params_op_t s_ive_op;
    } isvce_ctl_set_gop_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set Deblock params                                       */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_deblock_params_ip_t
    {
        ive_ctl_set_deblock_params_ip_t s_ive_ip;
    } isvce_ctl_set_deblock_params_ip_t;

    typedef struct isvce_ctl_set_deblock_params_op_t
    {
        ive_ctl_set_deblock_params_op_t s_ive_op;
    } isvce_ctl_set_deblock_params_op_t;

    /*****************************************************************************/
    /*   Video control  Set Profile params                                       */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_profile_params_ip_t
    {
        ive_ctl_set_profile_params_ip_t s_ive_ip;
    } isvce_ctl_set_profile_params_ip_t;

    typedef struct isvce_ctl_set_profile_params_op_t
    {
        ive_ctl_set_profile_params_op_t s_ive_op;
    } isvce_ctl_set_profile_params_op_t;

    /*****************************************************************************/
    /*   Synchronous video encode call                                           */
    /*****************************************************************************/
    typedef struct isvce_nalu_info_buf_t
    {
        /* For each NALU, following info will be copied as a csv string - */
        /* 'type,length,SId,TID,isIDR,isFirstSliceInLayer,isLastSliceInLayer' */
        UWORD8 *pu1_buf;

        UWORD32 u4_num_bytes;

        UWORD32 u4_buf_size;
    } isvce_nalu_info_buf_t;

    typedef struct isvce_video_encode_ip_t
    {
        ive_video_encode_ip_t s_ive_ip;

        isvce_nalu_info_buf_t *ps_nalu_info_buf;

    } isvce_video_encode_ip_t;

    typedef struct isvce_video_encode_op_t
    {
        ive_video_encode_op_t s_ive_op;

        bool b_is_nalu_info_present;

        isvce_nalu_info_buf_t *ps_nalu_info_buf;

    } isvce_video_encode_op_t;

    /*****************************************************************************/
    /*   Video usability information                                             */
    /*****************************************************************************/
    typedef struct isvce_vui_ip_t
    {
        /** size of the structure  */
        UWORD32 u4_size;

        /** Command type : ISVCE_CMD_VIDEO_CTL  */
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        /** Sub command type : ISVCE_CMD_CTL_SET_GOP_PARAMS */
        ISVCE_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

        /** indicates the presence of aspect_ratio */
        UWORD8 u1_aspect_ratio_info_present_flag;

        /** specifies the aspect ratio of the luma samples */
        UWORD8 u1_aspect_ratio_idc;

        /** width of the luma samples. user dependent */
        UWORD16 u2_sar_width;

        /** Height of the luma samples. user dependent */
        UWORD16 u2_sar_height;

        /** if 1, specifies that the overscan_appropriate_flag is present
         * if 0, the preferred display method for the video signal is unspecified */
        UWORD8 u1_overscan_info_present_flag;

        /** if 1,indicates that the cropped decoded pictures output
         * are suitable for display using overscan */
        UWORD8 u1_overscan_appropriate_flag;

        /** if 1 specifies that video_format, video_full_range_flag and
         * colour_description_present_flag are present */
        UWORD8 u1_video_signal_type_present_flag;

        /** pal, secam, ntsc, ...  */
        UWORD8 u1_video_format;

        /** indicates the black level and range of the luma and chroma signals */
        UWORD8 u1_video_full_range_flag;

        /** if 1,specifies that colour_primaries, transfer_characteristics
         * and matrix_coefficients are present */
        UWORD8 u1_colour_description_present_flag;

        /** indicates the chromaticity coordinates of the source primaries  */
        UWORD8 u1_colour_primaries;

        /** indicates the opto-electronic transfer characteristic of the source picture */
        UWORD8 u1_transfer_characteristics;

        /** the matrix coefficients used in deriving luma and chroma signals
         * from the green, blue, and red primaries */
        UWORD8 u1_matrix_coefficients;

        /** if 1, specifies that chroma_sample_loc_type_top_field and
         * chroma_sample_loc_type_bottom_field are present */
        UWORD8 u1_chroma_loc_info_present_flag;

        /** location of chroma samples */
        UWORD8 u1_chroma_sample_loc_type_top_field;

        UWORD8 u1_chroma_sample_loc_type_bottom_field;

        /**  Indicates the presence of the num_units_in_ticks, time_scale flag */
        UWORD8 u1_vui_timing_info_present_flag;

        /**  Number of units that correspond to one increment of the
         *   clock. Indicates the  resolution */
        UWORD32 u4_vui_num_units_in_tick;

        /**  The number of time units that pass in one second */
        UWORD32 u4_vui_time_scale;

        /** Flag indicating that time difference between two frames is a constant */
        UWORD8 u1_fixed_frame_rate_flag;

        /** Indicates the presence of NAL HRD parameters */
        UWORD8 u1_nal_hrd_parameters_present_flag;

        /** Indicates the presence of VCL HRD parameters */
        UWORD8 u1_vcl_hrd_parameters_present_flag;

        /** Specifies the HRD operational mode */
        UWORD8 u1_low_delay_hrd_flag;

        /** Indicates presence of SEI messages which include pic_struct syntax element */
        UWORD8 u1_pic_struct_present_flag;

        /** 1, specifies that the following cvs bitstream restriction parameters are present */
        UWORD8 u1_bitstream_restriction_flag;

        /** if 0, indicates that no pel outside the pic boundaries and
         * no sub-pels derived using pels outside the pic boundaries is used for inter prediction */
        UWORD8 u1_motion_vectors_over_pic_boundaries_flag;

        /** Indicates a number of bytes not exceeded by the sum of the sizes of the VCL NAL units
         * associated with any coded picture */
        UWORD8 u1_max_bytes_per_pic_denom;

        /** Indicates an upper bound for the number of bits of coding_unit() data */
        UWORD8 u1_max_bits_per_mb_denom;

        /** Indicate the maximum absolute value of a decoded horizontal MV component
         * in quarter-pel luma units */
        UWORD8 u1_log2_max_mv_length_horizontal;

        /** Indicate the maximum absolute value of a decoded vertical MV component
         * in quarter-pel luma units */
        UWORD8 u1_log2_max_mv_length_vertical;

        /** Max number of frames that are not synchronized in display and decode order */
        UWORD8 u1_num_reorder_frames;

        /** specifies required size of the HRD DPB in units of frame buffers */
        UWORD8 u1_max_dec_frame_buffering;

    } isvce_vui_ip_t;

    typedef struct isvce_vui_op_t
    {
        /** size of the structure                                           */
        UWORD32 u4_size;

        /** Return error code                                               */
        UWORD32 u4_error_code;
    } isvce_vui_op_t;

    /*****************************************************************************/
    /*    Video control  Set SEI MDCV params                                     */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_sei_mdcv_params_ip_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Command type : ISVCE_CMD_VIDEO_CTL                                  */
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        /** Sub command type : ISVCE_CMD_CTL_SET_SEI_MDCV_PARAMS                */
        ISVCE_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

        /** mastering display color volume info present flag                  */
        UWORD8 u1_sei_mdcv_params_present_flag;

        /** Array to store the display_primaries_x values                     */
        UWORD16 au2_display_primaries_x[3];

        /** Array to store the display_primaries_y values                     */
        UWORD16 au2_display_primaries_y[3];

        /** Variable to store the white point x value                         */
        UWORD16 u2_white_point_x;

        /** Variable to store the white point y value                         */
        UWORD16 u2_white_point_y;

        /** Variable to store the max display mastering luminance value       */
        UWORD32 u4_max_display_mastering_luminance;

        /** Variable to store the min display mastering luminance value       */
        UWORD32 u4_min_display_mastering_luminance;

        /** Lower 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_low;

        /** Upper 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_high;

    } isvce_ctl_set_sei_mdcv_params_ip_t;

    typedef struct isvce_ctl_set_sei_mdcv_params_op_t
    {
        /** size of the structure                                           */
        UWORD32 u4_size;

        /** Return error code                                               */
        UWORD32 u4_error_code;

    } isvce_ctl_set_sei_mdcv_params_op_t;

    /*****************************************************************************/
    /*    Video control  Set SEI CLL params                                      */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_sei_cll_params_ip_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Command type : ISVCE_CMD_VIDEO_CTL                                  */
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        /** Sub command type : ISVCE_CMD_CTL_SET_SEI_CLL_PARAMS                 */
        ISVCE_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

        /** content light level info present flag                             */
        UWORD8 u1_sei_cll_params_present_flag;

        /** The maximum pixel intensity of all samples                        */
        UWORD16 u2_max_content_light_level;

        /** The average pixel intensity of all samples                        */
        UWORD16 u2_max_pic_average_light_level;

        /** Lower 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_low;

        /** Upper 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_high;

    } isvce_ctl_set_sei_cll_params_ip_t;

    typedef struct isvce_ctl_set_sei_cll_params_op_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Return error code                                                 */
        UWORD32 u4_error_code;

    } isvce_ctl_set_sei_cll_params_op_t;

    /*****************************************************************************/
    /*    Video control  Set SEI AVE params                                      */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_sei_ave_params_ip_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Command type : ISVCE_CMD_VIDEO_CTL                                  */
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        /** Sub command type : ISVCE_CMD_CTL_SET_SEI_AVE_PARAMS                 */
        ISVCE_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

        /** ambient viewing environment info present flag                     */
        UWORD8 u1_sei_ave_params_present_flag;

        /** specifies the environmental illluminance of the ambient viewing
         * environment                                                        */
        UWORD32 u4_ambient_illuminance;

        /** specify the normalized x chromaticity coordinates of the
         * environmental ambient light in the nominal viewing environment     */
        UWORD16 u2_ambient_light_x;

        /** specify the normalized y chromaticity coordinates of the
         * environmental ambient light in the nominal viewing environment     */
        UWORD16 u2_ambient_light_y;

        /** Lower 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_low;

        /** Upper 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_high;

    } isvce_ctl_set_sei_ave_params_ip_t;

    typedef struct isvce_ctl_set_sei_ave_params_op_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Return error code                                                 */
        UWORD32 u4_error_code;

    } isvce_ctl_set_sei_ave_params_op_t;

    /*****************************************************************************/
    /*    Video control  Set SEI CCV params                                      */
    /*****************************************************************************/
    typedef struct isvce_ctl_set_sei_ccv_params_ip_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Command type : ISVCE_CMD_VIDEO_CTL                                  */
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        /** Sub command type : ISVCE_CMD_CTL_SET_SEI_CCV_PARAMS                 */
        ISVCE_CONTROL_API_COMMAND_TYPE_T e_sub_cmd;

        /** content color volume info present flag                            */
        UWORD8 u1_sei_ccv_params_present_flag;

        /** Flag used to control persistence of CCV SEI messages              */
        UWORD8 u1_ccv_cancel_flag;

        /** specifies the persistence of the CCV SEI message for the
         * current layer                                                      */
        UWORD8 u1_ccv_persistence_flag;

        /** specifies the presence of syntax elements ccv_primaries_x
         * and ccv_primaries_y                                                */
        UWORD8 u1_ccv_primaries_present_flag;

        /** specifies that the syntax element ccv_min_luminance_value
         * is present                                                         */
        UWORD8 u1_ccv_min_luminance_value_present_flag;

        /** specifies that the syntax element ccv_max_luminance_value
         *  is present                                                        */
        UWORD8 u1_ccv_max_luminance_value_present_flag;

        /** specifies that the syntax element ccv_avg_luminance_value
         *  is present                                                        */
        UWORD8 u1_ccv_avg_luminance_value_present_flag;

        /** shall be equal to 0 in bitstreams conforming to this version.
         * Other values for reserved_zero_2bits are reserved for future use   */
        UWORD8 u1_ccv_reserved_zero_2bits;

        /** specify the normalized x chromaticity coordinates of the colour
         * primary component c of the nominal content colour volume           */
        WORD32 ai4_ccv_primaries_x[3];

        /** specify the normalized y chromaticity coordinates of the colour
         * primary component c of the nominal content colour volume           */
        WORD32 ai4_ccv_primaries_y[3];

        /** specifies the normalized minimum luminance value                  */
        UWORD32 u4_ccv_min_luminance_value;

        /** specifies the normalized maximum luminance value                  */
        UWORD32 u4_ccv_max_luminance_value;

        /** specifies the normalized average luminance value                  */
        UWORD32 u4_ccv_avg_luminance_value;

        /** Lower 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_low;

        /** Upper 32bits of time stamp corresponding to input buffer,
         * from which this command takes effect                               */
        UWORD32 u4_timestamp_high;

    } isvce_ctl_set_sei_ccv_params_ip_t;

    typedef struct isvce_ctl_set_sei_ccv_params_op_t
    {
        /** size of the structure                                             */
        UWORD32 u4_size;

        /** Return error code                                                 */
        UWORD32 u4_error_code;

    } isvce_ctl_set_sei_ccv_params_op_t;

    /* The enum values should not have greater than 8 bits as this is assigned to WORD8 */
    typedef enum IV_MB_TYPE_T
    {
        INTRA16x16 = 0,
        INTRA4x4,
        INTER16x16
    } IV_MB_TYPE_T;

    /*****************************************************************************/
    /*   Pic info structures                                                     */
    /*****************************************************************************/
    typedef struct isvce_pic_info1_t
    {
        /** Qp  */
        UWORD32 u4_qp;

        /** Pic Type */
        IV_PICTURE_CODING_TYPE_T e_frame_type;

    } isvce_pic_info1_t;

    /*****************************************************************************/
    /*   MB info structures                                                     */
    /*****************************************************************************/
    typedef struct isvce_mv_t
    {
        /** MV X    */
        WORD16 i2_mv_x;

        /** MV Y    */
        WORD16 i2_mv_y;
    } isvce_mv_t;

    typedef struct isvce_mb_info1_t
    {
        /** Intra / Inter    */
        WORD8 i1_mb_type;

        union
        {
            isvce_mv_t as_mv[1];

            /** Intra mode */
            WORD8 ai1_intra_mode[1];
        };
    } isvce_mb_info1_t;

    typedef struct isvce_mb_info2_t
    {
        /** Intra / Inter    */
        WORD8 i1_mb_type;

        /** SAD     */
        UWORD16 u2_sad;

        union
        {
            isvce_mv_t as_mv[1];

            /** Intra mode */
            WORD8 ai1_intra_mode[1];
        };

    } isvce_mb_info2_t;

    typedef struct isvce_mb_info3_t
    {
        /** Intra / Inter    */
        WORD8 i1_mb_type;

        union
        {
            isvce_mv_t as_mv[4];

            /** Intra mode */
            WORD8 ai1_intra_mode[16];
        };

    } isvce_mb_info3_t;

    typedef struct isvce_mb_info4_t
    {
        /** Intra / Inter    */
        WORD8 i1_mb_type;

        /** Intra Mode      */
        WORD8 i1_intra_mode;

        /** SAD     */
        UWORD16 u2_sad;

        union
        {
            isvce_mv_t as_mv[16];

            /** Intra mode */
            WORD8 ai1_intra_mode[16];
        };

    } isvce_mb_info4_t;

    /* Add any new structures to the following union. It is used to calculate the
     * max size needed for allocation of memory */
    typedef struct isvce_api_mb_info_t
    {
        union
        {
            isvce_mb_info1_t s_mb_info1;
            isvce_mb_info2_t s_mb_info2;
            isvce_mb_info3_t s_mb_info3;
            isvce_mb_info4_t s_mb_info4;
        };
    } isvce_api_mb_info_t;

    typedef struct isvce_pic_info2_t
    {
        /** Qp  */
        UWORD32 u4_qp;

        /** Pic Type */
        IV_PICTURE_CODING_TYPE_T e_frame_type;

        /** Disable deblock level (0: Enable completely, 3: Disable completely */
        UWORD32 u4_disable_deblock_level;

    } isvce_pic_info2_t;

    typedef struct isvce_api_cmds_t
    {
        ISVCE_API_COMMAND_TYPE_T e_cmd;

        ISVCE_CONTROL_API_COMMAND_TYPE_T e_ctl_cmd;
    } isvce_api_cmds_t;

    extern IV_STATUS_T isvce_api_function(iv_obj_t *ps_handle, void *pv_api_ip, void *pv_api_op,
                                          isvce_api_cmds_t *ps_iv_api_cmds);

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif

#endif
