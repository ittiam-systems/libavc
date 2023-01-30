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
*  isvce_structs.h
*
* @brief
*  Contains struct definition used for SVC encoding
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_STRUCTS_H_
#define _ISVCE_STRUCTS_H_

#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264_defs.h"
#include "ih264_deblk_edge_filters.h"
#include "isvc_inter_pred_filters.h"
#include "ithread.h"
#include "isvc_defs.h"
#include "isvc_mem_fns.h"
#include "isvc_cabac_tables.h"
#include "isvc_trans_quant_itrans_iquant.h"

/* Dependencies of ime_structs.h */
#include "ime_defs.h"
#include "ime_distortion_metrics.h"

/* Dependencies of ih264e_cabac_structs.h */
#include "ih264_cabac_tables.h"

/* Dependencies of ih264e_structs.h */
#include "ih264e_error.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "ih264e_bitstream.h"
#include "ih264e_cabac_structs.h"
#include "ih264e_defs.h"
#include "ime_structs.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"

#include "ih264e_structs.h"
#include "isvce_cabac_structs.h"
#include "isvce_defs.h"
#include "isvce_downscaler.h"
#include "isvce_interface_structs.h"
#include "isvce_nalu_stat_aggregator.h"
#include "isvce_pred_structs.h"
#include "isvce_rc_utils.h"

#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"

typedef struct svc_params_t
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

} svc_params_t;

typedef struct svc_layer_data_t
{
    /**
     * Array of structs that contain mode_info per MB for every MB per layer
     */
    isvce_mb_info_t *ps_mb_info;

    UWORD32 *pu4_num_pus_in_mb;

} svc_layer_data_t;

typedef struct svc_au_data_t
{
    /**
     * Array of structs that contain layer-wise data used for svc prediction
     */
    svc_layer_data_t *ps_svc_layer_data;

    /**
     * Absolute POC for the current MV Bank
     */
    WORD32 i4_abs_poc;

    /**
     * Buffer Id
     */
    WORD32 i4_buf_id;

} svc_au_data_t;

typedef struct isvce_inp_buf_t
{
    /* App's buffer */
    isvce_raw_inp_buf_t s_inp_props;

    /* A copy of SVC parameters */
    svc_params_t s_svc_params;

    /**
     * Array of structs that contain properties of the buffers used for storing
     * layer-wise YUV data
     */
    yuv_buf_props_t as_layer_yuv_buf_props[MAX_NUM_SPATIAL_LAYERS];

} isvce_inp_buf_t;

typedef struct mb_intra_modes_t
{
    UWORD8 au1_intra_modes[MAX_PU_IN_MB];
} mb_intra_modes_t;

typedef struct nbr_info_t
{
    isvce_mb_info_t *ps_top_row_mb_info;

    isvce_mb_info_t *ps_left_mb_info;

    mb_intra_modes_t *ps_top_mb_intra_modes;

    mb_intra_modes_t *ps_left_mb_intra_modes;

} nbr_info_t;

typedef struct svc_nbr_info_t
{
    /**
     * Array of structs that contain properties of the buffers used for storing
     * layer-wise neighbour info
     */
    nbr_info_t *ps_layer_nbr_info;
} svc_nbr_info_t;

typedef struct layer_resampler_props_t
{
    UWORD32 u4_shift_x;

    UWORD32 u4_shift_y;

    UWORD32 u4_scale_x;

    UWORD32 u4_scale_y;

    WORD32 i4_offset_x;

    WORD32 i4_offset_y;

    WORD32 i4_add_x;

    WORD32 i4_add_y;

    WORD32 i4_delta_x;

    WORD32 i4_delta_y;

    WORD32 i4_refphase_x;

    WORD32 i4_refphase_y;

    WORD32 i4_phase_x;

    WORD32 i4_phase_y;

    UWORD32 u4_sub_wd;

    UWORD32 u4_sub_ht;

    UWORD32 u4_mb_wd;

    UWORD32 u4_mb_ht;

} layer_resampler_props_t;

typedef struct svc_ilp_data_t
{
    /* Pointer to current AU buf */
    svc_au_data_t *ps_svc_au_data;

    /* Array of bufs corresponding to numSpatialLayers */
    layer_resampler_props_t *aps_layer_resampler_props[NUM_SP_COMPONENTS];

    /* Array of bufs corresponding to numSpatialLayers */
    yuv_buf_props_t *ps_intra_recon_bufs;

    /* Array of bufs corresponding to numSpatialLayers */
    yuv_buf_props_t *ps_residual_bufs;
} svc_ilp_data_t;

typedef struct ilp_mv_t
{
    isvce_enc_pu_mv_t as_mv[ENC_MAX_PU_IN_MB][NUM_PRED_DIRS];

    MBTYPES_T e_mb_type;

    PRED_MODE_T ae_pred_mode[ENC_MAX_PU_IN_MB];
} ilp_mv_t;

typedef struct ilp_me_cands_t
{
    isvce_enc_pu_mv_t as_mv[MAX_PU_IN_MB + MAX_ILP_MV_IN_NBR_RGN][NUM_PRED_DIRS];

    MBTYPES_T e_mb_type[MAX_PU_IN_MB + MAX_ILP_MV_IN_NBR_RGN];

    PRED_MODE_T ae_pred_mode[MAX_PU_IN_MB + MAX_ILP_MV_IN_NBR_RGN];

    UWORD32 u4_num_ilp_mvs;

    UWORD32 u4_num_ilp_mvs_incl_nbrs;
} ilp_me_cands_t;

typedef struct isvce_cfg_params_t
{
    /** maximum width for which codec should request memory requirements    */
    UWORD32 u4_max_wd;

    /** maximum height for which codec should request memory requirements   */
    UWORD32 u4_max_ht;

    /** Maximum number of reference frames                                  */
    UWORD32 u4_max_ref_cnt;

    /** Maximum number of reorder frames                                    */
    UWORD32 u4_max_reorder_cnt;

    /** Maximum level supported                                             */
    UWORD32 u4_max_level;

    /** Input color format                                                  */
    IV_COLOR_FORMAT_T e_inp_color_fmt;

    /** Flag to enable/disable - To be used only for debugging/testing      */
    UWORD32 u4_enable_recon;

    /** Recon color format                                                  */
    IV_COLOR_FORMAT_T e_recon_color_fmt;

    /** Encoder Speed preset - Value between 0 (slowest) and 100 (fastest)  */
    IVE_SPEED_CONFIG u4_enc_speed_preset;

    /** Rate control mode                                                   */
    IVE_RC_MODE_T e_rc_mode;

    /** Maximum frame rate to be supported                                  */
    UWORD32 u4_max_framerate;

    /** Maximum bitrate to be supported                                     */
    UWORD32 au4_max_bitrate[MAX_NUM_SPATIAL_LAYERS];

    /** Maximum number of consecutive  B frames                             */
    UWORD32 u4_num_bframes;

    /** Content type Interlaced/Progressive                                 */
    IV_CONTENT_TYPE_T e_content_type;

    /** Maximum search range to be used in X direction                      */
    UWORD32 u4_max_srch_rng_x;

    /** Maximum search range to be used in Y direction                      */
    UWORD32 u4_max_srch_rng_y;

    /** Slice Mode                                                          */
    IVE_SLICE_MODE_T e_slice_mode;

    /** Slice parameter                                                     */
    UWORD32 u4_slice_param;

    /** Processor architecture                                          */
    IV_ARCH_T e_arch;

    /** SOC details                                                     */
    IV_SOC_T e_soc;

    /** Input width to be sent in bitstream                                */
    UWORD32 u4_disp_wd;

    /** Input height to be sent in bitstream                               */
    UWORD32 u4_disp_ht;

    /** Input width                                                     */
    UWORD32 u4_wd;

    /** Input height                                                    */
    UWORD32 u4_ht;

    /** Input stride                                                    */
    UWORD32 u4_strd;

    /** Source frame rate                                               */
    UWORD32 u4_src_frame_rate;

    /** Target frame rate                                               */
    UWORD32 u4_tgt_frame_rate;

    /** Target bitrate in kilobits per second                           */
    UWORD32 au4_target_bitrate[MAX_NUM_SPATIAL_LAYERS];

    /** Force current frame type                                        */
    IV_PICTURE_CODING_TYPE_T e_frame_type;

    /** Encoder mode                                                    */
    IVE_ENC_MODE_T e_enc_mode;

    /** Set initial Qp for I pictures                                   */
    UWORD32 au4_i_qp[MAX_NUM_SPATIAL_LAYERS];

    /** Set initial Qp for P pictures                                   */
    UWORD32 au4_p_qp[MAX_NUM_SPATIAL_LAYERS];

    /** Set initial Qp for B pictures                                   */
    UWORD32 au4_b_qp[MAX_NUM_SPATIAL_LAYERS];

    /** Set minimum Qp for I pictures                                   */
    UWORD32 au4_i_qp_min[MAX_NUM_SPATIAL_LAYERS];

    /** Set maximum Qp for I pictures                                   */
    UWORD32 au4_i_qp_max[MAX_NUM_SPATIAL_LAYERS];

    /** Set minimum Qp for P pictures                                   */
    UWORD32 au4_p_qp_min[MAX_NUM_SPATIAL_LAYERS];

    /** Set maximum Qp for P pictures                                   */
    UWORD32 au4_p_qp_max[MAX_NUM_SPATIAL_LAYERS];

    /** Set minimum Qp for B pictures                                   */
    UWORD32 au4_b_qp_min[MAX_NUM_SPATIAL_LAYERS];

    /** Set maximum Qp for B pictures                                   */
    UWORD32 au4_b_qp_max[MAX_NUM_SPATIAL_LAYERS];

    /** Adaptive intra refresh mode                                     */
    IVE_AIR_MODE_T e_air_mode;

    /** Adaptive intra refresh period in frames                         */
    UWORD32 u4_air_refresh_period;

    /** VBV buffer delay                                                */
    UWORD32 au4_vbv_buffer_delay[MAX_NUM_SPATIAL_LAYERS];

    /** Number of cores to be used                                      */
    UWORD32 u4_num_cores;

    /** ME speed preset - Value between 0 (slowest) and 100 (fastest)      */
    UWORD32 u4_me_speed_preset;

    /** Flag to enable/disable half pel motion estimation               */
    UWORD32 u4_enable_hpel;

    /** Flag to enable/disable quarter pel motion estimation            */
    UWORD32 u4_enable_qpel;

    /** Flag to enable/disable intra 4x4 analysis                       */
    UWORD32 u4_enable_intra_4x4;

    /** Flag to enable/disable intra 8x8 analysis                       */
    UWORD32 u4_enable_intra_8x8;

    /** Flag to enable/disable intra 16x16 analysis                     */
    UWORD32 u4_enable_intra_16x16;

    /** Flag to enable/disable fast SAD approximation                   */
    UWORD32 u4_enable_fast_sad;

    /*flag to enable/disable alternate reference frames                 */
    UWORD32 u4_enable_alt_ref;

    /*Flag to enable/disable computation of SATDQ in ME*/
    UWORD32 u4_enable_satqd;

    /*Minimum SAD to search for*/
    WORD32 i4_min_sad;

    /** Maximum search range in X direction for farthest reference      */
    UWORD32 u4_srch_rng_x;

    /** Maximum search range in Y direction for farthest reference      */
    UWORD32 u4_srch_rng_y;

    /** I frame interval                                                */
    UWORD32 u4_i_frm_interval;

    /** IDR frame interval                                              */
    UWORD32 u4_idr_frm_interval;

    /** Disable deblock level (0: Enable completely, 3: Disable completely */
    UWORD32 u4_disable_deblock_level;

    /** Profile                                                         */
    IV_PROFILE_T e_profile;

    /** Lower 32bits of time stamp corresponding to input buffer,
     * from which this command takes effect                             */
    UWORD32 u4_timestamp_low;

    /** Upper 32bits of time stamp corresponding to input buffer,
     * from which this command takes effect                             */
    UWORD32 u4_timestamp_high;

    /** Flag to say if the current config parameter set is valid
     * Will be zero to start with and will be set to 1, when configured
     * Once encoder uses the parameter set, this will be set to zero */
    UWORD32 u4_is_valid;

    /** Command associated with this config param set */
    ISVCE_CONTROL_API_COMMAND_TYPE_T e_cmd;

    /** Input width in mbs                                                    */
    UWORD32 i4_wd_mbs;

    /** Input height in mbs                                                   */
    UWORD32 i4_ht_mbs;

    /** entropy coding mode flag                                              */
    UWORD32 u4_entropy_coding_mode;

    /** enable weighted prediction                                            */
    UWORD32 u4_weighted_prediction;

    /** Pic info type */
    UWORD32 u4_pic_info_type;
    /**
     * MB info type
     */
    UWORD32 u4_isvce_mb_info_type;

    /** VUI structure                                                         */
    vui_t s_vui;

    /** SEI structure                                                         */
    sei_params_t s_sei;

    /** Flag to enable/disable VUI from header                          */
    UWORD32 u4_disable_vui;

    /** SVC params                                                            */
    svc_params_t s_svc_params;

    bool b_nalu_info_export_enable;

} isvce_cfg_params_t;

typedef struct mb_qp_ctxt_t
{
    UWORD8 u1_cur_mb_qp;

} mb_qp_ctxt_t;

typedef struct isvce_entropy_ctxt_t
{
    /**
     * Pointer to the cabac context
     */
    isvce_cabac_ctxt_t *ps_cabac;

    mb_qp_ctxt_t *ps_mb_qp_ctxt;

    /**
     * start of frame / start of slice flag
     */
    WORD32 i4_sof;

    /**
     * end of frame / end of slice flag
     */
    WORD32 i4_eof;

    /**
     * generate header upon request
     */
    WORD32 i4_gen_header;

    WORD32 i4_gen_subset_sps;

    /**
     * Pointer to base of sequence parameter set structure array
     */
    sps_t *ps_sps_base;

    /**
     * Pointer to base of Picture parameter set structure array
     */
    pps_t *ps_pps_base;

    /**
     * Current slice idx
     */
    WORD32 i4_cur_slice_idx;

    /**
     * Points to the array of slice indices which is used to identify the
     * independent slice to which each MB in a frame belongs.
     */
    UWORD8 *pu1_slice_idx;

    /**
     * Pointer to base of svc_nalu_ext structure array
     */
    svc_nalu_ext_t *ps_svc_nalu_ext_base;

    /**
     * Pointer to base of subset sequence parameter set structure array
     */
    subset_sps_t *ps_subset_sps_base;

    /**
     * Pointer to base of slice header structure array
     */
    slice_header_t *ps_slice_hdr_base;

    /**
     * Pointer to base of SVC slice header structure array
     */
    svc_slice_header_t *ps_svc_slice_hdr_base;

    /**
     * entropy status
     */
    UWORD8 *pu1_entropy_map;

    /**
     * MB's x position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_x;

    /**
     * MB's y position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_y;

    /**
     * MB start address
     */
    WORD32 i4_mb_cnt;

    /**
     * MB start address
     */
    WORD32 i4_mb_start_add;

    /**
     * MB end address
     */
    WORD32 i4_mb_end_add;

    /**
     * Input width in mbs
     */
    WORD32 i4_wd_mbs;

    /**
     * Input height in mbs
     */
    WORD32 i4_ht_mbs;

    /**
     * Bitstream structure
     */
    bitstrm_t *ps_bitstrm;

#if ENABLE_RE_ENC_AS_SKIP
    bitstrm_t *ps_bitstrm_after_slice_hdr;
#endif

    /**
     *  transform_8x8_mode_flag
     */
    WORD8 i1_transform_8x8_mode_flag;

    /**
     *  entropy_coding_mode_flag
     */
    WORD8 u1_entropy_coding_mode_flag;

    /**
     * Pointer to the top row nnz for luma
     */
    UWORD8 (*pu1_top_nnz_luma)[4];

    /**
     * left nnz for luma
     */
    UWORD32 u4_left_nnz_luma;

    /**
     * Pointer to zero runs before for the mb
     */
    UWORD8 au1_zero_run[16];

    /**
     * Pointer to the top row nnz for chroma
     */
    UWORD8 (*pu1_top_nnz_cbcr)[4];

    /**
     * left nnz for chroma
     */
    UWORD8 u4_left_nnz_cbcr;

    /**
     * Pointer frame level mb subblock coeff data
     */
    void *pv_pic_mb_coeff_data;

    /**
     * Pointer to mb subblock coeff data and number of subblocks and scan idx
     * Incremented each time a coded subblock is processed
     */
    void *pv_mb_coeff_data;

    /**
     * Pointer frame level mb header data
     */
    void *pv_pic_mb_header_data;

    /**
     * Pointer to mb header data and
     * incremented each time a coded mb is encoded
     */
    void *pv_mb_header_data;

    /**
     * Error code during parse stage
     */
    IH264E_ERROR_T i4_error_code;

    /**
     * Void pointer to job context
     */
    void *pv_proc_jobq, *pv_entropy_jobq;

    /**
     * Flag to signal end of frame
     */
    WORD32 i4_end_of_frame;

    /**
     * Abs POC count of the frame
     */
    WORD32 i4_abs_pic_order_cnt;

    /**
     * mb skip run
     */
    WORD32 *pi4_mb_skip_run;

    /**
     * Flag to signal end of sequence
     */
    UWORD32 u4_is_last;

    /**
     * Lower 32bits of time-stamp corresponding to the buffer being encoded
     */
    UWORD32 u4_timestamp_low;

    /**
     * Upper 32bits of time-stamp corresponding to the buffer being encoded
     */
    UWORD32 u4_timestamp_high;

    /**
     * Current Picture count - used for synchronization
     */
    WORD32 i4_pic_cnt;

    /**
     * Number of bits consumed by header for I and P mb types
     */
    UWORD32 u4_header_bits[MAX_MB_TYPE];

    /**
     * Number of bits consumed by residue for I and P mb types
     */
    UWORD32 u4_residue_bits[MAX_MB_TYPE];

    UWORD8 u1_spatial_layer_id;

} isvce_entropy_ctxt_t;

/**
 ******************************************************************************
 *  @brief      Rate control related variables
 ******************************************************************************
 */
typedef struct isvce_rate_control_ctxt_t
{
    void *apps_rate_control_api[MAX_NUM_SPATIAL_LAYERS];

    void *pps_frame_time;

    void *pps_time_stamp;

    void *pps_pd_frm_rate;

    /**
     * frame rate pull down
     */
    WORD32 pre_encode_skip[MAX_CTXT_SETS];

    /**
     * skip frame (cbr)
     */
    WORD32 post_encode_skip[MAX_CTXT_SETS];

    /**
     * rate control type
     */
    rc_type_e e_rc_type;

    /**
     * pic type
     */
    picture_type_e e_pic_type;

    /**
     * rc utils context
     */
    svc_rc_utils_ctxt_t s_rc_utils;

    /**
     * intra cnt in previous frame
     */
    WORD32 ai4_num_intra_in_prev_frame[MAX_NUM_SPATIAL_LAYERS];

    /**
     * avg activity of prev frame
     */
    WORD32 ai4_avg_activity[MAX_NUM_SPATIAL_LAYERS];

} isvce_rate_control_ctxt_t;

typedef struct
{
    /**
     * mb type and mode
     */
    UWORD8 u1_mb_type_mode;

    /**
     * CBP
     */
    UWORD8 u1_cbp;

    /**
     * MB qp delta
     */
    UWORD8 u1_mb_qp;

    /**
     * Element to align structure to 2 byte boundary
     */
    UWORD8 u1_pad;

    UWORD8 u1_base_mode_flag;

    UWORD8 u1_residual_prediction_flag;

} isvce_mb_hdr_common_t;

/**
******************************************************************************
*  @brief      macro block info for I4x4 MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

    /**
     * Sub block modes, 2 modes per byte
     */
    UWORD8 au1_sub_blk_modes[8];
} isvce_mb_hdr_i4x4_t;

/**
******************************************************************************
*  @brief      macro block info for I8x8 MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

    /**
     * Sub block modes, 2 modes per byte
     */
    UWORD8 au1_sub_blk_modes[2];
} isvce_mb_hdr_i8x8_t;

/**
******************************************************************************
*  @brief      macro block info for I16x16 MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

} isvce_mb_hdr_i16x16_t;

/**
******************************************************************************
*  @brief      macro block info for P16x16 MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

    /**
     * MV
     */
    WORD16 ai2_mvd[2];

    UWORD8 u1_mvp_idx;
} isvce_mb_hdr_p16x16_t;

/**
******************************************************************************
*  @brief      macro block info for PSKIP MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

} isvce_mb_hdr_pskip_t;

/**
******************************************************************************
*  @brief      macro block info for B16x16 MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

    /**
     * MV
     */
    WORD16 ai2_mvd[NUM_PRED_DIRS][2];

    UWORD8 au1_mvp_idx[NUM_PRED_DIRS];
} isvce_mb_hdr_b16x16_t;

/**
******************************************************************************
*  @brief      macro block info for BDIRECT MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

} isvce_mb_hdr_bdirect_t;

/**
******************************************************************************
*  @brief      macro block info for PSKIP MB
******************************************************************************
*/
typedef struct
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

} isvce_mb_hdr_bskip_t;

/**
******************************************************************************
*  @brief      macro block info for IBL MB
******************************************************************************
*/
typedef struct isvce_mb_hdr_base_mode_t
{
    /**
     * Common MB header params
     */
    isvce_mb_hdr_common_t common;

} isvce_mb_hdr_base_mode_t;

/**
******************************************************************************
*  @brief      Union of mb_hdr structures for size calculation
*  and to access first few common elements
******************************************************************************
*/

typedef union isvce_mb_hdr_t
{
    isvce_mb_hdr_i4x4_t mb_hdr_i4x4;
    isvce_mb_hdr_i8x8_t mb_hdr_i8x8;
    isvce_mb_hdr_i16x16_t mb_hdr_i16x16;
    isvce_mb_hdr_p16x16_t mb_hdr_p16x16;
    isvce_mb_hdr_pskip_t mb_hdr_pskip;
    isvce_mb_hdr_b16x16_t mb_hdr_b16x16;
    isvce_mb_hdr_bdirect_t mb_hdr_bdirect;
    isvce_mb_hdr_bskip_t mb_hdr_bskip;
    isvce_mb_hdr_base_mode_t mb_hdr_base_mode;
} isvce_mb_hdr_t;

typedef struct isvce_bs_ctxt_t
{
    /**
     * MB's x position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_x;

    /**
     * MB's y position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_y;

    /**
     * MB's x position within a Slice in raster scan in MB units
     */
    WORD32 i4_mb_slice_x;

    /**
     * MB's y position within a Slice in raster scan in MB units
     */
    WORD32 i4_mb_slice_y;

    /**
     * Vertical strength, Two bits per edge.
     * Stored in format. BS[15] | BS[14] | .. |BS[0]
     */
    UWORD32 *pu4_pic_vert_bs;

    UWORD32 *pu4_intra_base_vert_bs;

    /**
     * Boundary strength, Two bits per edge.
     * Stored in format. BS[15] | BS[14] | .. |BS[0]
     */
    UWORD32 *pu4_pic_horz_bs;

    UWORD32 *pu4_intra_base_horz_bs;

    /**
     *  Qp array stored for each mb
     */
    UWORD8 *pu1_pic_qp;

} isvce_bs_ctxt_t;

typedef struct isvce_deblk_ctxt_t
{
    /**
     * MB's x position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_x;

    /**
     * MB's y position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_y;

    /**
     * structure that contains BS and QP frame level arrays
     */
    isvce_bs_ctxt_t s_bs_ctxt;

    /*
     * Recon Buffers
     */
    yuv_buf_props_t s_rec_pic_buf_props;

    /**
     *  Points to the array of slice indices which is used to identify the slice
     *  to which each MB in a frame belongs.
     */
    UWORD8 *pu1_slice_idx;

} isvce_deblk_ctxt_t;

/**
**************************************************************************
*   @brief   isvce_me_ctxt_t
*
*   Structure encapsulating the parameters used in the motion estimation
*   context
**************************************************************************
*/
typedef struct isvce_me_ctxt_t
{
    /**
     * Ref pointer to current MB luma for each ref list
     */
    UWORD8 *apu1_ref_buf_luma[MAX_NUM_REFLIST];

    /**
     * Src pointer to current MB luma
     */
    UWORD8 *pu1_src_buf_luma;

    /**
     * source stride
     * (strides for luma and chroma are the same)
     */
    WORD32 i4_src_strd;

    /**
     * recon stride
     * (strides for luma and chroma are the same)
     */
    WORD32 ai4_rec_strd[MAX_NUM_REFLIST];

    /**
     * Offset for half pel x plane from the pic buf
     */
    UWORD32 u4_half_x_offset;

    /**
     * Offset for half pel y plane from half x plane
     */
    UWORD32 u4_half_y_offset;

    /**
     * Offset for half pel xy plane from half y plane
     */
    UWORD32 u4_half_xy_offset;

    /**
     *  Search range in the X, Y axis in terms of pixels
     */
    WORD32 ai2_srch_boundaries[2];

    /**
     *  Search range in the north direction in terms of pixels
     */
    WORD32 i4_srch_range_n;

    /**
     *  Search range in the south direction in terms of pixels
     */
    WORD32 i4_srch_range_s;

    /**
     *  Search range in the east direction in terms of pixels
     */
    WORD32 i4_srch_range_e;

    /**
     *  Search range in the west direction in terms of pixels
     */
    WORD32 i4_srch_range_w;

    /**
     * left mb motion vector
     */
    ime_mv_t s_left_mv;

    /**
     * top left mb motion vector
     */
    ime_mv_t s_top_left_mv;

    /*
     * ilp MVs for ME candidates *
     */
    ilp_me_cands_t *ps_ilp_me_cands;

    /**
     * Number of valid candidates for the Initial search position
     */
    UWORD32 u4_num_candidates[MAX_NUM_REFLIST + 1];

    /**
     * Motion vector predictors derived from neighboring
     * blocks for each of the six block partitions
     */
    ime_mv_t as_mv_init_search[MAX_NUM_REFLIST + 1][MAX_FPEL_SEARCH_CANDIDATES];

    /**
     * mv bits
     */
    UWORD8 *pu1_mv_bits;

    /**
     * lambda (lagrange multiplier for cost computation)
     */
    UWORD32 u4_lambda_motion;

    /**
     * enabled fast sad computation
     */
    UWORD32 u4_enable_fast_sad;

    /*
     * Enable SKIP block prediction based on SATQD
     */
    UWORD32 u4_enable_stat_sad;

    /*
     * Minimum distortion to search for
     * */
    WORD32 i4_min_sad;

    /*
     * Signal that minimum sad has been reached in ME
     * */
    UWORD32 u4_min_sad_reached;

    /**
     * Flag to enable/disbale half pel motion estimation
     */
    UWORD32 u4_enable_hpel;

    /**
     * Diamond search Iteration Max Cnt
     */
    UWORD32 u4_num_layers;

    /**
     * encoder me speed
     */
    UWORD32 u4_me_speed_preset;

    UWORD32 u4_left_is_intra;

    UWORD32 u4_left_is_skip;

    /* skip_type can be PREDL0, PREDL1 or  BIPRED */
    WORD32 i4_skip_type;

    /* Biasing given for skip prediction */
    WORD32 i4_skip_bias[2];

    /**
     * Structure to store the MB partition info
     * We need 1(L0)+1(L1)+1(bi)
     */
    mb_part_ctxt as_mb_part[MAX_NUM_REFLIST + 1];
    /*
     * Threshold to compare the sad with
     */
    UWORD16 *pu2_sad_thrsh;

    /**
     * fn ptrs for compute sad routines
     */
    ime_compute_sad_ft *pf_ime_compute_sad_16x16[2];
    ime_compute_sad_ft *pf_ime_compute_sad_16x8;
    ime_compute_sad4_diamond *pf_ime_compute_sad4_diamond;
    ime_compute_sad3_diamond *pf_ime_compute_sad3_diamond;
    ime_compute_sad2_diamond *pf_ime_compute_sad2_diamond;
    ime_sub_pel_compute_sad_16x16_ft *pf_ime_sub_pel_compute_sad_16x16;

    /*
     * Function poitners for SATQD
     */
    ime_compute_sad_stat *pf_ime_compute_sad_stat_luma_16x16;

    /**
     * Qp
     */
    UWORD8 u1_mb_qp;

    /*
     * Buffers for holding subpel and bipred temp buffers
     */
    UWORD8 *apu1_subpel_buffs[SUBPEL_BUFF_CNT];

    WORD32 u4_subpel_buf_strd;

    /*
     * Buffers to store the best halfpel plane*
     */
    UWORD8 *pu1_hpel_buf;

} isvce_me_ctxt_t;

typedef struct isvce_mb_info_nmb_t
{
    UWORD32 u4_mb_type;
    UWORD32 u4_min_sad;
    UWORD32 u4_min_sad_reached;
    WORD32 i4_mb_cost;
    WORD32 i4_mb_distortion;

    isvce_enc_pu_mv_t as_skip_mv[4];

    isvce_enc_pu_mv_t as_pred_mv[2];

    block_neighbors_t s_ngbr_avbl;

    /*
     * Buffer to hold best subpel buffer in each MB of NMB
     */
    UWORD8 *pu1_best_sub_pel_buf;

    /*
     * Stride for subpel buffer
     */
    UWORD32 u4_bst_spel_buf_strd;

} isvce_mb_info_nmb_t;

typedef struct isvce_process_ctxt_t
{
    svc_params_t s_svc_params;

    /* Resolves circular dependency with svc_ilp_mv_ctxt_t */
    void *ps_svc_ilp_mv_ctxt;

    /* Resolves circular dependency with svc_res_pred_ctxt_t */
    void *ps_res_pred_ctxt;

    /* Resolves circular dependency with svc_intra_pred_ctxt_t */
    void *ps_intra_pred_ctxt;

    /* Resolves circular dependency with svc_sub_pic_rc_ctxt_t */
    void *ps_sub_pic_rc_ctxt;

    yuv_buf_props_t *ps_mb_pred_buf;

    yuv_buf_props_t *ps_mb_res_buf;

    ilp_mv_t *ps_ilp_mv;

    /**
     * entropy context
     */
    isvce_entropy_ctxt_t s_entropy;

    /**
     * me context
     */
    isvce_me_ctxt_t s_me_ctxt;

    /* Resolves circular dependency with isvce_codec_t */
    void *ps_codec;

    /**
     * N mb process contest
     */
    n_mb_process_ctxt_t s_n_mb_ctxt;

    /*
     * Src Buffers
     */
    yuv_buf_props_t s_src_buf_props;

    /*
     * Recon Buffers
     */
    yuv_buf_props_t s_rec_buf_props;

    /*
     * Reference Frame Buffers
     */
    yuv_buf_props_t as_ref_buf_props[MAX_REF_PIC_CNT];

    /*
     * Src Buffers
     */
    yuv_buf_props_t s_src_pic_buf_props;

    /*
     * Recon Buffers
     */
    yuv_buf_props_t s_rec_pic_buf_props;

    /*
     * Reference Frame Buffers
     */
    yuv_buf_props_t as_ref_pic_buf_props[MAX_REF_PIC_CNT];

    /**
     * Pointer to ME NMB info
     */
    isvce_mb_info_nmb_t *ps_nmb_info;

    isvce_mb_info_nmb_t *ps_cur_mb;

    /**
     * Offset for half pel x plane from the pic buf
     */
    UWORD32 u4_half_x_offset;

    /**
     * Offset for half pel y plane from half x plane
     */
    UWORD32 u4_half_y_offset;

    /**
     * Offset for half pel xy plane from half y plane
     */
    UWORD32 u4_half_xy_offset;

    /**
     * pred buffer pointer (temp buffer 1)
     */
    UWORD8 *pu1_pred_mb;

    /**
     * pred buffer pointer (prediction buffer for intra 16x16
     */
    UWORD8 *pu1_pred_mb_intra_16x16;

    /**
     * pred buffer pointer (prediction buffer for intra 16x16_plane
     */
    UWORD8 *pu1_pred_mb_intra_16x16_plane;

    /**
     * pred buffer pointer (prediction buffer for intra chroma
     */
    UWORD8 *pu1_pred_mb_intra_chroma;

    /**
     * pred buffer pointer (prediction buffer for intra chroma plane
     */
    UWORD8 *pu1_pred_mb_intra_chroma_plane;

    /**
     * temp. reference buffer ptr for intra 4x4 when rdopt is on
     */
    UWORD8 *pu1_ref_mb_intra_4x4;

    /**
     * prediction buffer stride
     */
    WORD32 i4_pred_strd;

    /**
     * transform buffer pointer (temp buffer 2)
     */
    WORD16 *pi2_res_buf;

    /**
     * temp. transform buffer ptr for intra 4x4 when rdopt is on
     */
    WORD16 *pi2_res_buf_intra_4x4;

    /**
     * transform buffer stride
     */
    WORD32 i4_res_strd;

    /**
     * scratch buffer for inverse transform (temp buffer 3)
     */
    void *pv_scratch_buff;

    /**
     * frame num
     */
    WORD32 i4_frame_num;

    /**
     * start address of frame / sub-frame
     */
    WORD32 i4_frame_strt_add;

    /**
     *  IDR pic
     */
    UWORD32 u4_is_idr;

    /**
     *  idr_pic_id
     */
    UWORD32 u4_idr_pic_id;

    /**
     * Input width in mbs
     */
    WORD32 i4_wd_mbs;

    /**
     * Input height in mbs
     */
    WORD32 i4_ht_mbs;

    /**
     *  slice_type
     */
    WORD32 i4_slice_type;

    /**
     * Current slice idx
     */
    WORD32 i4_cur_slice_idx;

    /**
     * MB's x position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_x;

    /**
     * MB's y position within a picture in raster scan in MB units
     */
    WORD32 i4_mb_y;

    /**
     * MB's x position within a Slice in raster scan in MB units
     */
    WORD32 i4_mb_slice_x;

    /**
     * MB's y position within a Slice in raster scan in MB units
     */
    WORD32 i4_mb_slice_y;

    /**
     * mb neighbor availability pointer
     */
    block_neighbors_t *ps_ngbr_avbl;

    /**
     * lambda (lagrange multiplier for cost computation)
     */
    UWORD32 u4_lambda;

    /**
     * mb distortion
     */
    WORD32 i4_mb_distortion;

    /**
     * mb cost
     */
    WORD32 i4_mb_cost;

    /********************************************************************/
    /* i4_ngbr_avbl_mb_16 - ngbr avbl of curr mb                        */
    /* i4_ngbr_avbl_sb_8 - ngbr avbl of all 8x8 sub blocks of curr mb   */
    /* i4_ngbr_avbl_sb_4 - ngbr avbl of all 4x4 sub blocks of curr mb   */
    /* i4_ngbr_avbl_mb_c - chroma ngbr avbl of curr mb                  */
    /********************************************************************/
    WORD32 i4_ngbr_avbl_16x16_mb;
    WORD32 ai4_neighbor_avail_8x8_subblks[4];
    UWORD8 au1_ngbr_avbl_4x4_subblks[16];
    WORD32 i4_chroma_neighbor_avail_8x8_mb;

    /**
     * array to store the mode of mb sub blocks
     */
    UWORD8 au1_intra_luma_mb_4x4_modes[16];

    /**
     * array to store the predicted mode of mb sub blks
     */
    UWORD8 au1_predicted_intra_luma_mb_4x4_modes[16];

    /**
     * macro block intra 16x16 mode
     */
    UWORD8 u1_l_i16_mode;

    /**
     * array to store the mode of the macro block intra 8x8 4 modes
     */
    UWORD8 au1_intra_luma_mb_8x8_modes[4];

    /**
     * intra chroma mb mode
     */
    UWORD8 u1_c_i8_mode;

    /********************************************************************/
    /* array to store pixels from the neighborhood for intra prediction */
    /* i16 - 16 left pels + 1 top left pel + 16 top pels = 33 pels      */
    /* i8 - 8 lpels + 1 tlpels + 8 tpels + 8 tr pels = 25 pels          */
    /* i4 - 4 lpels + 1 tlpels + 4 tpels + 4 tr pels = 13 pels          */
    /* ic - 8 left pels + 1 top left pel + 8 top pels )*2               */
    /********************************************************************/
    UWORD8 au1_ngbr_pels[34];

    /**
     * array for 8x8 intra pels filtering (temp buff 4)
     */
    UWORD8 au1_neighbor_pels_i8x8_unfiltered[25];

    /**
     * Number of sub partitons in the inter pred MB
     */
    UWORD32 u4_num_sub_partitions;

    /**
     *  Pointer to hold num PUs each MB in a picture
     */
    UWORD32 *pu4_mb_pu_cnt;

    /**
     * Pointer to the array of structures having motion vectors, size
     *  and position of sub partitions
     */
    isvce_mb_info_t *ps_mb_info;

    /**
     * Pointer to the pu of current co-located MB in list 1
     */
    isvce_mb_info_t *ps_col_mb;

    /**
     * predicted motion vector
     */
    isvce_enc_pu_mv_t *ps_skip_mv;

    /**
     * predicted motion vector
     */
    isvce_enc_pu_mv_t *ps_pred_mv;

    /**
     * top row mb syntax information base
     * In normal working scenarios, for a given context set,
     * the mb syntax info pointer is identical across all process threads.
     * But when the hard bound on slices are enabled, in multi core, frame
     * is partitioned in to sections equal to set number of cores and each
     * partition is run independently. In this scenario, a ctxt set will alone
     * appear to run multiple frames at a time. For this to occur, the common
     * pointers across the proc ctxt should disappear.
     *
     * This is done by allocating MAX_PROCESS_THREADS memory and distributing
     * across individual ctxts when byte bnd per slice is enabled.
     */
    svc_nbr_info_t s_nbr_info_base;

    nbr_info_t s_nbr_info;

    /**
     * mb neighbor availability pointer
     */
    block_neighbors_t s_ngbr_avbl;

    /**
     * coded block pattern
     */
    UWORD32 u4_cbp;

    /**
     *  number of non zero coeffs
     */
    UWORD32 au4_nnz[5];

    UWORD8 au1_chroma_nnz[2 * (NUM_4x4_IN_8x8 + 1)];

    /**
     *  number of non zero coeffs for intra 4x4 when rdopt is on
     */
    UWORD32 au4_nnz_intra_4x4[4];

    /**
     * frame qp & mb qp
     */
    UWORD8 u1_frame_qp;

    UWORD8 u1_mb_qp;

    /**
     * quantization parameters for luma & chroma planes
     */
    quant_params_t *ps_qp_params[3];

    /**
     * Pointer frame level mb subblock coeff data
     */
    void *pv_pic_mb_coeff_data;

    /**
     * Pointer to mb subblock coeff data and number of subblocks and scan idx
     * Incremented each time a coded subblock is processed
     */
    void *pv_mb_coeff_data;

    /**
     * Pointer frame level mb header data
     */
    void *pv_pic_mb_header_data;

    /**
     * Pointer to mb header data and
     * incremented each time a coded mb is encoded
     */
    void *pv_mb_header_data;

    /**
     * Signal that pic_init is called first time
     */
    WORD32 i4_first_pic_init;

    /**
     * Current MV Bank's buffer ID
     */
    WORD32 i4_cur_mv_bank_buf_id;

    /**
     * Void pointer to job context
     */
    void *pv_proc_jobq, *pv_entropy_jobq;

    /**
     * Number of MBs to be processed in the current Job
     */
    WORD32 i4_mb_cnt;

    /**
     * ID for the current context - Used for debugging
     */
    WORD32 i4_id;

    /**
     * Pointer to current picture buffer structure
     */
    svc_au_buf_t *ps_cur_pic;

    /**
     * Pointer to current picture's mv buffer structure
     */
    svc_au_data_t *ps_cur_mv_buf;

    /**
     * Flag to indicate if ps_proc was initialized at least once in a frame.
     * This is needed to handle cases where a core starts to handle format
     * conversion jobs directly
     */
    WORD32 i4_init_done;

    /**
     * Process status: one byte per MB
     */
    UWORD8 *pu1_proc_map;

    /**
     * Deblk status: one byte per MB
     */
    UWORD8 *pu1_deblk_map;

    /**
     * Process status: one byte per MB
     */
    UWORD8 *pu1_me_map;

    /*
     * Intra refresh mask.
     * Indicates if an Mb is coded in intra mode within the current AIR interval
     * NOTE Refreshes after each AIR period
     * NOTE The map is shared between process
     */
    UWORD8 *pu1_is_intra_coded;

    /**
     * Disable deblock level (0: Enable completely, 3: Disable completely
     */
    UWORD32 u4_disable_deblock_level;

    /**
     * Pointer to the structure that contains deblock context
     */
    isvce_deblk_ctxt_t s_deblk_ctxt;

    /**
     * Points to the array of slice indices which is used to identify the
     * independent slice to which each MB in a frame belongs.
     */
    UWORD8 *pu1_slice_idx;

    /**
     * Pointer to base of svc_nalu_ext structure array
     */
    svc_nalu_ext_t *ps_svc_nalu_ext_base;

    /**
     * Pointer to base of subset sequence parameter set structure array
     */
    subset_sps_t *ps_subset_sps_base;

    /**
     * Pointer to base of slice header structure array
     */
    slice_header_t *ps_slice_hdr_base;

    /**
     * Pointer to base of SVC slice header structure array
     */
    svc_slice_header_t *ps_svc_slice_hdr_base;

    /**
     * Number of mb's to process in one loop
     */
    WORD32 i4_nmb_ntrpy;

    /**
     * Number of mb's to process in one loop
     */
    UWORD32 u4_nmb_me;

    /**
     * Structure for current input buffer
     */
    isvce_inp_buf_t s_inp_buf;

    /**
     * api call cnt
     */
    WORD32 i4_encode_api_call_cnt;

    /**
     * Current Picture count - used for synchronization
     */
    WORD32 i4_pic_cnt;

    /**
     * Intermediate buffer for interpred leaf level functions
     */
    WORD32 ai16_pred1[HP_BUFF_WD * HP_BUFF_HT];

    /**
     * Reference picture for the current picture
     * TODO: Only 2 reference assumed currently
     */
    svc_au_buf_t *aps_ref_pic[MAX_REF_PIC_CNT];

    /**
     * Reference MV buff for the current picture
     */
    svc_au_data_t *aps_mv_buf[MAX_REF_PIC_CNT];

    /**
     * frame info used by RC
     */
    frame_info_t s_frame_info;

    /*
     * NOTE NOT PERSISTANT INSIDE FUNCTIONS
     * Min sad for current MB
     * will be populated initially
     * Once a sad less than eq to u4_min_sad is reached, the value will be copied
     * to the cariable
     */
    UWORD32 u4_min_sad;

    /*
     * indicates weather we have rached minimum sa or not
     */
    UWORD32 u4_min_sad_reached;

    /**
     * Current error code
     */
    WORD32 i4_error_code;

    /*
     * Enables or disables computation of recon
     */
    UWORD32 u4_compute_recon;

    /*
     * Temporary buffers to be used for subpel computation
     */
    UWORD8 *apu1_subpel_buffs[SUBPEL_BUFF_CNT];

    /*
     * Buffer holding best sub pel values
     */
    UWORD8 *pu1_best_subpel_buf;

    /*
     * Stride for buffer holding best sub pel
     */
    UWORD32 u4_bst_spel_buf_strd;

    /*
     * SVC spatial layer ID
     */
    UWORD8 u1_spatial_layer_id;
} isvce_process_ctxt_t;

typedef UWORD8 FT_CORE_CODING(isvce_process_ctxt_t *ps_proc);

typedef WORD32 FT_FIND_SKIP_PARAMS(isvce_process_ctxt_t *, WORD32);

typedef void FT_ME_ALGORITHM(isvce_process_ctxt_t *);

typedef struct enc_loop_fxns_t
{
    /**
     * luma core coding function pointer
     */
    FT_CORE_CODING *apf_luma_energy_compaction[MAX_MBTYPES];

    /**
     * chroma core coding function pointer
     */
    FT_CORE_CODING *apf_chroma_energy_compaction[2];

    /**
     * forward transform for intra blk of mb type 16x16
     */
    FT_LUMA_16X16_RESI_TRANS_DCTRANS_QUANT
    *pf_resi_trans_dctrans_quant_16x16;

    /**
     * inverse transform for intra blk of mb type 16x16
     */
    FT_LUMA_16X16_IDCTRANS_IQUANT_ITRANS_RECON
    *pf_idctrans_iquant_itrans_recon_16x16;

    /**
     * forward transform for 4x4 blk luma
     */
    FT_RESI_TRANS_QUANT *apf_resi_trans_quant_4x4[NUM_RESI_TRANS_QUANT_VARIANTS];

    /**
     * forward transform for 4x4 blk luma
     */
    FT_RESI_TRANS_QUANT
    *apf_resi_trans_quant_chroma_4x4[NUM_RESI_TRANS_QUANT_VARIANTS];

    /*
     * hadamard transform and quant for a 4x4 block
     */
    FT_HADAMARD_QUANT *pf_hadamard_quant_4x4;

    /*
     *  hadamard transform and quant for a 4x4 block
     */
    FT_HADAMARD_QUANT *pf_hadamard_quant_2x2_uv;

    /**
     * inverse transform for 4x4 blk
     */
    FT_IQ_IT_RECON *apf_iquant_itrans_recon_4x4[NUM_IQ_IT_RECON_VARIANTS];

    /**
     * inverse transform for chroma 4x4 blk
     */
    FT_IQ_IT_RECON *apf_iquant_itrans_recon_chroma_4x4[NUM_IQ_IT_RECON_VARIANTS];

    /**
     * inverse transform for 4x4 blk with only single dc coeff
     */
    FT_IQ_IT_RECON *apf_iquant_itrans_recon_4x4_dc[NUM_IQ_IT_RECON_VARIANTS];

    /**
     * inverse transform for chroma 4x4 blk with only single dc coeff
     */
    FT_IQ_IT_RECON
    *apf_iquant_itrans_recon_chroma_4x4_dc[NUM_IQ_IT_RECON_VARIANTS];

    /*
     * Inverse hadamard transform and iquant for a 4x4 block
     */
    FT_IHADAMARD_SCALING *pf_ihadamard_scaling_4x4;

    /*
     * Inverse hadamard transform and iquant for a 4x4 block
     */
    FT_IHADAMARD_SCALING *pf_ihadamard_scaling_2x2_uv;

    /**
     * forward transform for 8x8 blk
     */
    FT_RESI_TRANS_QUANT *apf_resi_trans_quant_8x8[NUM_RESI_TRANS_QUANT_VARIANTS];

    /**
     * inverse transform for 8x8 blk
     */
    FT_IQ_IT_RECON *apf_iquant_itrans_recon_8x8[NUM_IQ_IT_RECON_VARIANTS];

    FT_IQ_IT_RECON *pf_zcbf_iquant_itrans_recon_4x4;

    FT_IQ_IT_RECON *pf_chroma_zcbf_iquant_itrans_recon_4x4;

} enc_loop_fxns_t;

typedef struct inter_pred_fxns_t
{
    FT_INTER_PRED_LUMA *pf_inter_pred_luma_copy;

    FT_INTER_PRED_LUMA *pf_inter_pred_luma_horz;

    FT_INTER_PRED_LUMA *pf_inter_pred_luma_vert;

    FT_INTER_PRED_LUMA_BILINEAR *pf_inter_pred_luma_bilinear;

    FT_INTER_PRED_CHROMA *pf_inter_pred_chroma;
} inter_pred_fxns_t;

typedef struct mem_fxns_t
{
    FT_MEMCPY *pf_mem_cpy;

    FT_MEMSET *pf_mem_set;

    FT_MEMCPY *pf_mem_cpy_mul8;

    FT_MEMSET *pf_mem_set_mul8;

    FT_COPY_2D *pf_copy_2d;

    FT_MEMSET_2D *pf_memset_2d;

    FT_16BIT_INTERLEAVED_COPY *pf_16bit_interleaved_copy;

    FT_16BIT_INTERLEAVED_MEMSET *pf_16bit_interleaved_memset;

    FT_NONZERO_CHECKER *pf_nonzero_checker;

} mem_fxns_t;

typedef struct isa_dependent_fxns_t
{
    enc_loop_fxns_t s_enc_loop_fxns;

    inter_pred_fxns_t s_inter_pred_fxns;

    mem_fxns_t s_mem_fxns;
} isa_dependent_fxns_t;

/**
 * Reference set containing pointers to MV buf and pic buf
 */
typedef struct
{
    /** Picture count */
    WORD32 i4_pic_cnt;

    /** POC */
    WORD32 i4_poc;

    /** picture buffer */
    svc_au_buf_t *ps_pic_buf;

    /** mv buffer */
    svc_au_data_t *ps_svc_au_data;

} isvce_ref_set_t;

typedef struct isvce_codec_t
{
    /**
     * downscaler context
     */
    downscaler_ctxt_t s_scaler;

    svc_ilp_data_t s_svc_ilp_data;

    nalu_descriptors_t as_nalu_descriptors[MAX_NUM_SPATIAL_LAYERS];

    isa_dependent_fxns_t s_isa_dependent_fxns;

#if ENABLE_MODE_STAT_VISUALISER
    /* Resolves circular dependency with mode_stat_visualiser_t */
    void *ps_mode_stat_visualiser;
#endif

    /** enable constrained intra prediction */
    UWORD32 au4_constrained_intra_pred[MAX_NUM_SPATIAL_LAYERS];

    /**
     * Id of current pic (input order)
     */
    WORD32 i4_poc;

    /**
     * Number of encode frame API calls made
     * This variable must only be used for context selection [Read only]
     */
    WORD32 i4_encode_api_call_cnt;

    /**
     * Number of pictures encoded
     */
    WORD32 i4_pic_cnt;

    /**
     * Number of threads created
     */
    WORD32 i4_proc_thread_cnt;

    /**
     * Mutex used to keep the control calls thread-safe
     */
    void *pv_ctl_mutex;

    /**
     * Current active config parameters
     */
    isvce_cfg_params_t s_cfg;

    /**
     * Array containing the config parameter sets
     */
    isvce_cfg_params_t as_cfg[MAX_ACTIVE_CONFIG_PARAMS];

    /**
     * Color format used by encoder internally
     */
    IV_COLOR_FORMAT_T e_codec_color_format;

    /**
     * recon stride
     * (strides for luma and chroma are the same)
     */
    WORD32 i4_rec_strd;

    /**
     * Flag to enable/disable deblocking of a frame
     */
    WORD32 u4_disable_deblock_level;

    /**
     * Number of continuous frames where deblocking was disabled
     */
    WORD32 u4_disable_deblock_level_cnt;

    /**
     * frame type
     */
    PIC_TYPE_T pic_type;

    /**
     * frame qp
     */
    UWORD32 au4_frame_qp[MAX_NUM_SPATIAL_LAYERS];

    /**
     * Enable inital QP calculation based on BPP and GPP
     */
    UWORD8 u1_enable_init_qp;

    /**
     * frame num
     */
    WORD32 i4_frame_num;

    /**
     *  slice_type
     */
    WORD32 i4_slice_type;

    /*
     * Force current frame to specific type
     */
    IV_PICTURE_CODING_TYPE_T force_curr_frame_type;

    /**
     *  IDR pic
     */
    UWORD32 u4_is_idr;

    /**
     *  idr_pic_id
     */
    WORD32 i4_idr_pic_id;

    /**
     * Flush mode
     */
    WORD32 i4_flush_mode;

    /**
     * Encode header mode
     */
    WORD32 i4_header_mode;

    /**
     * Flag to indicate if header has already
     * been generated when i4_api_call_cnt 0
     */
    UWORD32 u4_header_generated;

    /**
     * Encode generate header
     */
    WORD32 i4_gen_header;

    /**
     * To signal successful completion of init
     */
    WORD32 i4_init_done;

    /**
     * To signal that at least one picture was decoded
     */
    WORD32 i4_first_pic_done;

    /**
     * Reset flag - Codec is reset if this flag is set
     */
    WORD32 i4_reset_flag;

    /**
     * Current error code
     */
    WORD32 i4_error_code;

    /**
     * threshold residue
     */
    WORD32 u4_thres_resi;

    /**
     * disable intra inter gating
     */
    UWORD32 u4_inter_gate;

    /**
     * Holds mem records passed during init.
     * This will be used to return the mem records during retrieve call
     */
    iv_mem_rec_t *ps_mem_rec_backup;

    /**
     * Flag to determine if the entropy thread is active
     */
    volatile UWORD32 au4_entropy_thread_active[MAX_CTXT_SETS];

    /**
     * Mutex used to keep the entropy calls thread-safe
     */
    void *pv_entropy_mutex;

    /**
     * Job queue buffer base
     */
    void *pv_proc_jobq_buf, *pv_entropy_jobq_buf;

    /**
     * Job Queue mem tab size
     */
    WORD32 i4_proc_jobq_buf_size, i4_entropy_jobq_buf_size;

    /**
     * Memory for svc_au_data buffer manager
     */
    void *pv_svc_au_data_store_mgr_base;

    /**
     * svc_au_data buffer manager
     */
    void *pv_svc_au_data_store_mgr;

    /**
     * Pointer to svc_au_data structure array
     */
    svc_au_data_t *ps_svc_au_data;

    /**
     * Base address for svc_au_data
     */
    svc_au_data_t *ps_svc_au_data_base;

    /**
     * svc_au_data size
     */
    WORD32 i4_svc_au_data_size;

    /**
     * Memory for Picture buffer manager for reference pictures
     */
    void *pv_ref_buf_mgr_base;

    /**
     * Picture buffer manager for reference pictures
     */
    void *pv_ref_buf_mgr;

    /**
     * Number of reference buffers added to the buffer manager
     */
    WORD32 i4_ref_buf_cnt;

    /**
     * Pointer to Pic Buf structure array
     */
    svc_au_buf_t *ps_pic_buf;

    /**
     * Base address for Picture buffer
     */
    svc_au_buf_t *ps_pic_buf_base;

    /**
     * Total pic buffer size allocated
     */
    WORD32 i4_total_pic_buf_size;

    /**
     * Memory for Buffer manager for output buffers
     */
    void *pv_out_buf_mgr_base;

    /**
     * Buffer manager for output buffers
     */
    void *pv_out_buf_mgr;

    /**
     * Current output buffer's buffer ID
     */
    WORD32 i4_out_buf_id;

    /**
     * Number of output buffers added to the buffer manager
     */
    WORD32 i4_out_buf_cnt;

    /**
     * Memory for Picture buffer manager for input buffers
     */
    void *pv_inp_buf_mgr_base;

    /**
     * Picture buffer manager for input buffers
     */
    void *pv_inp_buf_mgr;

    /**
     * Current input buffer's buffer ID
     */
    WORD32 i4_inp_buf_id;

    /**
     * Number of input buffers added to the buffer manager
     */
    WORD32 i4_inp_buf_cnt;

    /**
     * Pointer to dpb manager structure
     */
    void *pv_dpb_mgr;

    /**
     * Pointer to base of Sequence parameter set structure array
     */
    sps_t *ps_sps_base;

    /**
     * Pointer to base of Picture parameter set structure array
     */
    pps_t *ps_pps_base;

    /**
     * Pointer to base of svc_nalu_ext structure array
     */
    svc_nalu_ext_t *ps_svc_nalu_ext_base;

    /**
     * Pointer to base of subset sequence parameter set structure array
     */
    subset_sps_t *ps_subset_sps_base;

    /**
     * Pointer to base of slice header structure array
     */
    slice_header_t *ps_slice_hdr_base;

    /**
     * Pointer to base of SVC slice header structure array
     */
    svc_slice_header_t *ps_svc_slice_hdr_base;

    /**
     * packed residue coeff data size for 1 row of mbs
     */
    UWORD32 u4_size_coeff_data;

    /**
     * packed header data size for 1 row of mbs
     */
    UWORD32 u4_size_header_data;

    /**
     * Processing context - One for each processing thread
     * Create two sets, each set used for alternate frames
     */
    isvce_process_ctxt_t as_process[MAX_PROCESS_CTXT];

    /**
     * Thread handle for each of the processing threads
     */
    void *apv_proc_thread_handle[MAX_PROCESS_THREADS];

    /**
     * Thread created flag for each of the processing threads
     */
    WORD32 ai4_process_thread_created[MAX_PROCESS_THREADS];

    /**
     * Void pointer to process job context
     */
    void *pv_proc_jobq, *pv_entropy_jobq;

    /**
     * Number of MBs processed together for better instruction cache handling
     */
    WORD32 i4_proc_nmb;

    /**
     * Previous POC lsb
     */
    WORD32 i4_prev_poc_lsb;

    /**
     * Previous POC msb
     */
    WORD32 i4_prev_poc_msb;

    /**
     * Max POC lsb that has arrived till now
     */
    WORD32 i4_max_prev_poc_lsb;

    /**
     * Context for format conversion
     */
    fmt_conv_t s_fmt_conv;

    /**
     * Absolute pic order count
     */
    WORD32 i4_abs_pic_order_cnt;

    /**
     *  Pic order count of lsb
     */
    WORD32 i4_pic_order_cnt_lsb;

    /**
     * Array giving current picture being processed in each context set
     */
    WORD32 ai4_pic_cnt[MAX_CTXT_SETS];

    /*
     * Min sad to search for
     */
    UWORD32 u4_min_sad;

    /**
     * Reference picture set
     */
    isvce_ref_set_t as_ref_set[MAX_DPB_SIZE + MAX_CTXT_SETS];

    /*
     * Air pic cnt
     * Contains the number of pictures that have been encoded with air
     * This value is moudulo air refresh period
     */
    WORD32 i4_air_pic_cnt;

    /*
     * Intra refresh map
     * Stores the frames at which intra refresh should occur for a MB
     */
    UWORD16 *pu2_intr_rfrsh_map;

    /*
     * Indicates if the current frame is used as a reference frame
     */
    UWORD32 u4_is_curr_frm_ref;

    /*
     * Indicates if there can be non reference frames in the stream
     */
    WORD32 i4_non_ref_frames_in_stream;

    /*
     * Memory for color space conversion for luma plane
     */
    UWORD8 *pu1_y_csc_buf_base;

    /*
     * Memory for color space conversion foe chroma plane
     */
    UWORD8 *pu1_uv_csc_buf_base;

    /**
     * Function pointers for intra pred leaf level functions luma
     */
    pf_intra_pred apf_intra_pred_16_l[MAX_I16x16];
    pf_intra_pred apf_intra_pred_8_l[MAX_I8x8];
    pf_intra_pred apf_intra_pred_4_l[MAX_I4x4];

    /**
     * Function pointers for intra pred leaf level functions chroma
     */
    pf_intra_pred apf_intra_pred_c[MAX_CH_I8x8];

    /**
     * deblock vertical luma edge with blocking strength 4
     */
    ih264_deblk_edge_bs4_ft *pf_deblk_luma_vert_bs4;

    /**
     * deblock vertical chroma edge with blocking strength 4
     */
    ih264_deblk_chroma_edge_bs4_ft *pf_deblk_chroma_vert_bs4;

    /**
     * deblock vertical luma edge with blocking strength less than 4
     */
    ih264_deblk_edge_bslt4_ft *pf_deblk_luma_vert_bslt4;

    /**
     * deblock vertical chroma edge with blocking strength less than 4
     */
    ih264_deblk_chroma_edge_bslt4_ft *pf_deblk_chroma_vert_bslt4;

    /**
     * deblock horizontal luma edge with blocking strength 4
     */
    ih264_deblk_edge_bs4_ft *pf_deblk_luma_horz_bs4;

    /**
     * deblock horizontal chroma edge with blocking strength 4
     */
    ih264_deblk_chroma_edge_bs4_ft *pf_deblk_chroma_horz_bs4;

    /**
     * deblock horizontal luma edge with blocking strength less than 4
     */
    ih264_deblk_edge_bslt4_ft *pf_deblk_luma_horz_bslt4;

    /**
     * deblock horizontal chroma edge with blocking strength less than 4
     */
    ih264_deblk_chroma_edge_bslt4_ft *pf_deblk_chroma_horz_bslt4;

    /**
     * functions for padding
     */
    pf_pad pf_pad_top;
    pf_pad pf_pad_bottom;
    pf_pad pf_pad_left_luma;
    pf_pad pf_pad_left_chroma;
    pf_pad pf_pad_right_luma;
    pf_pad pf_pad_right_chroma;

    /**
     * fn ptrs for compute sad routines
     */
    ime_compute_sad_ft *apf_compute_sad_16x16[2];
    ime_compute_sad_ft *pf_compute_sad_16x8;

    /**
     * Function pointer for computing ME
     * 1 for PSLICE and 1 for BSLICE
     */
    FT_ME_ALGORITHM *apf_compute_me[2];

    /**
     * Function pointers for computing SKIP parameters
     */
    FT_FIND_SKIP_PARAMS *apf_find_skip_params_me[2];

    /**
     * intra mode eval -encoder level function
     */
    pf_evaluate_intra_modes pf_ih264e_evaluate_intra16x16_modes;
    pf_evaluate_intra_modes pf_ih264e_evaluate_intra_chroma_modes;
    pf_evaluate_intra_4x4_modes pf_ih264e_evaluate_intra_4x4_modes;

    /* Half pel generation function - encoder level
     *
     */
    pf_sixtapfilter_horz pf_ih264e_sixtapfilter_horz;
    pf_sixtap_filter_2dvh_vert pf_ih264e_sixtap_filter_2dvh_vert;

    /**
     * color space conversion from YUV 420P to YUV 420Sp
     */
    pf_fmt_conv_420p_to_420sp pf_ih264e_conv_420p_to_420sp;

    /**
     * color space conversion from YUV 420P to YUV 420Sp
     */
    pf_fmt_conv_422ile_to_420sp pf_ih264e_fmt_conv_422i_to_420sp;

    /**
     * write mb layer for a given slice I, P, B
     */
    IH264E_ERROR_T (*pf_write_mb_syntax_layer[2][3])(isvce_entropy_ctxt_t *ps_ent_ctxt);

    /**
     * Output buffer
     */
    isvce_out_buf_t as_out_buf[MAX_CTXT_SETS];

    /**
     * recon buffer
     */
    isvce_rec_buf_t as_rec_buf[MAX_CTXT_SETS];

    /**
     * rate control context
     */
    isvce_rate_control_ctxt_t s_rate_control;

    /**
     * input buffer queue
     */
    isvce_inp_buf_t as_inp_list[SVC_MAX_NUM_INP_FRAMES];

    /**
     * Flag to indicate if any IDR requests are pending
     */
    WORD32 i4_pending_idr_flag;

    /**
     *Flag to indicate if we have recived the last input frame
     */
    WORD32 i4_last_inp_buff_received;

    /*
     * Max num reference frames to be signaled in SPS
     */
    WORD32 i4_max_num_reference_frames;

    /**
     * backup sei params for comparison
     */
    sei_params_t s_sei;
} isvce_codec_t;

#endif
