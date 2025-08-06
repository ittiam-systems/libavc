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
 *  isvcd_structs.h
 *
 * @brief
 *  Contains structures required for decoder
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_STRUCTS_H_
#define _ISVCD_STRUCTS_H_

#include "isvcd_defs.h"
#include "isvcd_cabac.h"
#include "ih264d_structs.h"
#include "isvcd_iquant_itrans_residual_recon.h"
#include "isvcd_iquant_itrans_residual.h"
#include "isvcd_iquant_itrans.h"
#include "isvcd_pred_residual_recon.h"

#include "isvcd_nal.h"
#include "isvcd_nal_structs.h"
#include "isvcd_nal_parse.h"
#include "isvcd_nal_parse_structs.h"

#include "isvcd_intra_resamp.h"
#include "isvcd_ii_pred.h"
#include "isvcd_residual_resamp.h"

#include "isvcd_vui.h"

#define NUM_MB_PARTS 4
#define NUM_SUB_MB_PARTS 4
#define NUM_INTRA_SUB_BLOCKS 16
#define MAX_NUM_MB_PART NUM_MB_PARTS *NUM_SUB_MB_PARTS

#define ANNEX_B 0     /*!< Annex B stream*/
#define NON_ANNEX_B 1 /*!< Non Annex B RFC stream */

#define BUFFER_ALIGN_4 4

#define MAX_VCL_NAL_BUFF_SIZE (1024 * 1024 * 2)
#define MAX_NON_VCL_NAL_BUFF_SIZE (1024 * 1024)
#define MAX_SCLD_REF_LAYER_OFFSET 32768
#define MIN_SCLD_REF_LAYER_OFFSET -32768
#define MAX_SVC_NAL_UNIT_TYPE 31
/*! Nal unit svc extennsion parameters */

struct _SvcDecLyrStruct;

typedef struct
{
    UWORD8 u1_nal_ref_idc;              /** NAL ref idc of the Slice NAL unit */
    UWORD8 u1_svc_ext_flag;             /** svc nal extension */
    UWORD8 u1_idr_flag;                 /** IDR picture when dependency_id = maximum value of
                                           dependency_id */
    UWORD8 u1_priority_id;              /** priority identifier for the NAL unit */
    UWORD8 u1_no_inter_layer_pred_flag; /** Usage of the inter-layer prediction */
    UWORD8 u1_dependency_id;            /** dependency identifier for the NAL unit */
    UWORD8 u1_quality_id;               /** quality identifier for the NAL unit */
    UWORD8 u1_temporal_id;              /** temporal identifier for the NAL unit */
    UWORD8 u1_use_ref_base_pic_flag;    /** specifies reference pictures for inter
                                           prediction */
    UWORD8 u1_discardable_flag;         /** current NAL unit is not used in dependant
                                           decoding */
    UWORD8 u1_output_flag;              /** decoded picture output and removal process */

} dec_nal_unit_svc_ext_params_t;

/* Structure to contain information about reference base pic marking */
typedef struct
{
    UWORD8 u1_adaptive_ref_base_pic_marking_mode_flag;
    UWORD32 u4_memory_management_base_control_operation;
    UWORD32 u4_difference_of_base_pic_nums_minus1;
    UWORD32 u4_long_term_base_pic_num;

} dec_ref_base_pic_marking_params_t;

/*! Sequence level parameters svc extension */

typedef struct
{
    /** Presence of deblocking filter for inter-layer prediction in the slice
     * header */
    UWORD8 u1_inter_layer_deblocking_filter_control_present_flag;
    UWORD8 u1_extended_spatial_scalability_idc; /** Geometrical parameters for the
                                                   resampling processes */
    /** horizontal phase shift of the chroma components in units of half luma
     * samples*/
    UWORD8 u1_chroma_phase_x_plus1_flag;
    /** vertical phase shift of the chroma components in units of half luma
     * samples */
    UWORD8 u1_chroma_phase_y_plus1;
    /* horizontal phase shift of chroma in units of half luma samples used for
     * inter-layer prediction */
    UWORD8 u1_seq_ref_layer_chroma_phase_x_plus1_flag;
    /* vertical phase shift of chroma in units of half luma samples used for
     * inter-layer prediction */
    UWORD8 u1_seq_ref_layer_chroma_phase_y_plus1;
    WORD32 i4_seq_scaled_ref_layer_left_offset;      /** horizontal left offset */
    WORD32 i4_seq_scaled_ref_layer_top_offset;       /** vertical top offset */
    WORD32 i4_seq_scaled_ref_layer_right_offset;     /**horizontal right offset */
    WORD32 i4_seq_scaled_ref_layer_bottom_offset;    /** vertical bottom offset */
    UWORD8 u1_seq_tcoeff_level_prediction_flag;      /** presence of the syntax element
                                                        adaptive_tcoeff_level_prediction_flag
                                                      */
    UWORD8 u1_adaptive_tcoeff_level_prediction_flag; /** presence of
                                                        tcoeff_level_prediction_flag
                                                      */
    /**  specifies presence of syntax elements in slice headers that refer to the
     * subset sequence parameter set */
    UWORD8 u1_slice_header_restriction_flag;
    UWORD8 u1_svc_vui_parameters_present_flag;
    svc_vui_ext_t *ps_svc_vui_ext;
} dec_subset_seq_params_t;

typedef struct
{
    WORD32 i4_left_offset; /*!< Scaled horizontal offset of the top left
                           corner luma sample of reference layer from
                           the top left corner luma sample of the current
                           layer. (In the units of num MBs) */
    WORD32 i4_rt_offset;   /*!< Scaled horizontal offset of the bottom right
                           corner luma sample of reference layer from
                           the top left corner luma sample of the
                           current layer. (In the units of num MBs) */
    WORD32 i4_top_offset;  /*!< Scaled vertical offset of the top left
                           corner luma sample of reference layer from
                           the top left corner luma sample of the current
                           layer. (In the units of num MBs) */
    WORD32 i4_bot_offset;  /*!< Scaled vertical offset of the bottom right
                           corner luma sample of reference layer from
                           the top left corner luma sample of the current
                           layer. (In the units of num MBs) */
} dec_svc_crop_wnd_offset_t;

typedef struct
{
    UWORD32 u4_ref_layer_dq_id;
    UWORD32 u4_disable_inter_layer_deblk_filter_idc;
    WORD32 i4_inter_layer_slice_alpha_c0_offset_div2;
    WORD32 i4_inter_layer_slice_beta_offset_div2;
    UWORD8 u1_constrained_intra_resampling_flag;
    UWORD8 u1_ref_layer_chroma_phase_x_plus1_flag;
    UWORD8 u1_ref_layer_chroma_phase_y_plus1;

    WORD32 i4_scaled_ref_layer_left_offset;
    WORD32 i4_scaled_ref_layer_right_offset;
    WORD32 i4_scaled_ref_layer_top_offset;
    WORD32 i4_scaled_ref_layer_bottom_offset;
    UWORD8 u1_slice_skip_flag;
    UWORD32 u4_num_mbs_in_slice_minus1;
    UWORD8 u1_adaptive_base_mode_flag;
    UWORD8 u1_default_base_mode_flag;
    UWORD8 u1_adaptive_motion_prediction_flag;
    UWORD8 u1_default_motion_prediction_flag;
    UWORD8 u1_adaptive_residual_prediction_flag;
    UWORD8 u1_default_residual_prediction_flag;
    UWORD8 u1_tcoeff_level_prediction_flag;
    UWORD8 u1_scan_idx_start;
    UWORD8 u1_scan_idx_end;
    UWORD8 u1_base_pred_weight_table_flag;
    dec_ref_base_pic_marking_params_t s_ref_base_pic_marking_svc_ext;
    UWORD8 u1_store_ref_base_pic_flag; /* specifies when dependency_id is equal to
                                          the max value of the VCL NAL units of
                                          the coded picture */

} dec_slice_svc_ext_params_t;

/* Prefix NAL unit svc extension parameters*/

typedef struct
{
    dec_nal_unit_svc_ext_params_t s_nal_svc_ext;
    dec_ref_base_pic_marking_params_t s_ref_base_pic_marking_svc_ext;
    UWORD8 u1_store_ref_base_pic_flag;                  /* specifies when dependency_id is equal to
                                                           the max value of the VCL NAL units of
                                                           the coded picture */
    UWORD8 u1_additional_prefix_nal_unit_ext_flag;      /* To indicate whether additional nal unit
                                                           extension data flag syntax elements */
    UWORD8 u1_additional_prefix_nal_unit_ext_data_flag; /*Used for FUTURE USE */

} dec_prefix_nal_unit_svc_ext_params_t;

typedef enum
{
    LIST_0 = 0,
    LIST_1 = 1,
    NUM_REF_LISTS
} REF_LIST_T;

typedef enum
{
    INTRA_16x16 = 0,
    INTRA_NXN
} INTRA_MB_PRED_MODE_T;

/* Transform type */
typedef enum
{
    T_4X4 = 0,
    T_8X8,
    T_PCM
} TRANSFORM_TYPE_T;

typedef struct
{
    vcl_node_t *ps_top_node; /*!< VCL node corresponding to top most
                             layer in the access unit.
                             */
    vcl_node_t *ps_bot_node; /*!< VCL node corresponding to
                             bottom most layer to be decoded in the
                             access unit. This parameter keeps updated
                             during course of decoding from base layer
                             till top most layer.
                             */

    WORD32 i4_num_res_lyrs;  /*!< Number of layers with
                             different resolutions. Layers with spatial
                             resoluton change flag equal to 0 are
                             considered to be of same resolution as
                             reference layer's resolution.
                             */
    /* following 2 parameter will be updated only if      */
    /* the picture boundary is detected due to difference */
    /* in slice header syntax                             */

    UWORD16 u2_frm_num_next;  /*!< frame number of the next
                              slice after picture boundary detection
                              */
    WORD8 i1_nal_ref_id_next; /*!< nal ref id of the next slice after picture boundary detection
                              range [-1,3];  -1 says the picture boundary is detected
                              by DQID of the layers
                              */
} vcl_nal_t;

typedef struct
{
    WORD32 i4_num_non_vcl_nals;              /*!< Total number of non vcl nals that are
                                             extracted from the bitstream.
                                             */

    non_vcl_buf_hdr_t *ps_first_non_vcl_nal; /*!< This shall point to first NON VCL NAL
                                             that is extracted from the input bitstream. This
                                             shall be set to NULL if there are no VCL NALs
                                             present in the bitstream in a access unit */

} non_vcl_nal_t;

typedef struct
{
    WORD16 i2_mv_x; /*!< motion vectors in horizontal direction
                    QPEL units
                    */
    WORD16 i2_mv_y; /*!< motion vectors in vertical direction QPEL
                    units
                    */
} mot_vec_t;

typedef struct
{
    mot_vec_t as_mv[NUM_SUB_MB_PARTS];
    WORD32 i4_ref_idx;
} mb_part_mv_t;

typedef struct
{
    /*
    PRED_16X16,
    PRED_16X8,
    PRED_8X16,
    PRED_8X8,
    */
    UWORD8 u1_part_type;

    UWORD8 u1_mv_cnt;
    /*
    0-1 bits   : 1st partition
    2-3 bits   : 2nd partition
    4-5 bits   : 3rd partition
    6-7 bits   : 4th partition

    Value
    00     : B_DIRECT
    01     : L0
    10     : L1
    11     : BiPred
    */
    UWORD8 u1_pred_mode; /* Pred mode shall have valid value till
                         mode motion prediction only. Hence this value shall
                         not be used for the motion compensation
                         */
    UWORD8 au1_mot_pred_flags[2];
    UWORD8 au1_sub_mb_part_ht_wd[NUM_MB_PARTS * NUM_SUB_MB_PARTS];
    UWORD8 au1_sub_mb_num[NUM_MB_PARTS * NUM_SUB_MB_PARTS >> 1];
    WORD8 ai1_sub_mb_part_type[NUM_MB_PARTS];
    mb_part_mv_t as_mb_part_dmv[2][NUM_MB_PARTS];
} inter_mb_prms_t; /* We need to improve on commenting */

typedef struct
{
    /*
    0   : I_16x16
    1   : I_NXN
    */
    WORD32 i4_pred_mode;
    WORD32 i4_chroma_intra_pred_mode;
    WORD8 ai1_rem_intra_pred_mode[NUM_INTRA_SUB_BLOCKS];
} intra_mb_prms_t; /* We need to improve on commenting */

typedef union
{
    inter_mb_prms_t s_inter;
    intra_mb_prms_t s_intra;
} mb_prms_ext_t;

typedef struct
{
    UWORD8 u1_chroma_nnz; /*! NNZs of Chroma. Here each bit corresonds
                          to a NNZs of 4x4 sub block. Lower 4 bits are
                          used for Cb and upper are used for Cr */
    UWORD16 u2_luma_nnz;  /*! NNZs of Luma. Here each bit corresonds
                          to a NNZs of 4x4 sub block in raster scan
                          order. */
    WORD8 i1_mb_mode;     /*! MB mode of an MB */

    WORD8 i1_tx_size;     /*! transform size of an MB */

    WORD8 i1_slice_id;
} inter_lyr_mb_prms_t;

/* the following 2 structures are used store certain parameters  */
/* across units for each layer or each dependency layer          */
typedef struct
{
    WORD32 i4_updated_sts; /*!< flag to indicate whether the
                           params have been updated
                           */
    WORD32 i4_ref_dq_id;   /*!< place to hold the ref_dqid of previous
                           access unit
                           */
    WORD32 i4_nal_ref_id;  /*!< place to hold the nal_ref_id of previous
                           access unit
                           */
    UWORD16 u2_frm_num;    /*!< place to hold the frame number of
                           previous access unit  will be used to
                           handle Errors in "frame_num"
                           syntax elements
                           */
} prev_au_prms_t;

typedef struct
{
    WORD32 i4_updated_sts; /*!< flag to indicate whether the
                           params have been updated
                           */

    UWORD8 u1_pps_id;      /*!< PPS ID of an access unit for a particular
                           layer. will be used in concealment of
                           Errors in next access unit
                           */

    UWORD8 u1_sps_id;      /*!< SPS ID of an access unit for a particular
                           layer. will be used in concealment of
                           Errors in next access unit
                           */

} prev_au_sps_pps_t;

typedef enum
{
    /*CABAC SVC related flags*/
    CABAC_BASE_MODE_FLAG = 460,
    CABAC_MOT_PRED_FLAG0 = 463,
    CABAC_MOT_PRED_FLAG1 = 464,
    CABAC_RES_PRED_FLAG = 465

} svc_cabac_table_num_t;

typedef struct
{
    dec_seq_params_t *ps_seq;
    dec_subset_seq_params_t s_sps_svc_ext;

    /* sequence associated frame paramateres*/
    WORD32 i4_reorder_depth;
    UWORD16 u2_disp_height;
    UWORD16 u2_disp_width;
    UWORD16 u2_pic_wd;
    UWORD16 u2_pic_ht;
    UWORD16 u2_frm_wd_y;
    UWORD16 u2_frm_ht_y;
    UWORD16 u2_frm_wd_uv;
    UWORD16 u2_frm_ht_uv;
    UWORD8 u1_pad_len_y_v;
    UWORD8 u1_pad_len_cr_v;
    UWORD16 u2_crop_offset_y;
    UWORD16 u2_crop_offset_uv;
} dec_svc_seq_params_t;

typedef struct
{
    /*svc related flags*/
    UWORD8 u1_base_mode_flag;
    UWORD8 u1_residual_prediction_flag;
    UWORD8 u1_crop_window_flag;
    UWORD8 au1_motion_pred_flag[2];
} dec_svc_mb_info_t;

typedef struct _SvcDecLyrStruct
{
    dec_struct_t s_dec;

    /*Pred + Res = Target when csbp is zero*/
    ih264_pred_residual_recon_ft *pf_pred_residual_recon_luma_4x4;

    ih264_pred_residual_recon_ft *pf_pred_residual_recon_luma_8x8;

    ih264_pred_residual_recon_ft *pf_pred_residual_recon_luma_16x16;

    ih264_pred_residual_recon_chroma_ft *pf_pred_residual_recon_chroma_4x4;

    ih264_pred_residual_recon_chroma_ft *pf_pred_residual_recon_chroma_8x8;

    /* IT + Res + Recon*/
    ih264_iquant_itrans_residual_recon_ft *pf_iquant_itrans_residual_recon_luma_4x4;

    ih264_iquant_itrans_residual_recon_ft *pf_iquant_itrans_residual_recon_luma_4x4_dc;

    ih264_iquant_itrans_residual_recon_ft *pf_iquant_itrans_residual_recon_luma_8x8;

    ih264_iquant_itrans_residual_recon_ft *pf_iquant_itrans_residual_recon_luma_8x8_dc;

    ih264_iquant_itrans_residual_recon_chroma_ft *pf_iquant_itrans_residual_recon_chroma_4x4;

    ih264_iquant_itrans_residual_recon_chroma_ft *pf_iquant_itrans_residual_recon_chroma_4x4_dc;

    /* Res nnz*/
    ih264_residual_ft *pf_residual_luma_4x4;
    ih264_residual_ft *pf_residual_luma_8x8;
    ih264_residual_ft *pf_residual_luma_16x16;

    ih264_residual_chroma_ft *pf_residual_chroma_cb_cr_8x8;

    /*IT + residual */
    ih264_iquant_itrans_residual_ft *pf_iquant_itrans_residual_luma_4x4;

    ih264_iquant_itrans_residual_ft *pf_iquant_itrans_residual_luma_4x4_dc;

    ih264_iquant_itrans_residual_ft *pf_iquant_itrans_residual_luma_8x8;

    ih264_iquant_itrans_residual_ft *pf_iquant_itrans_residual_luma_8x8_dc;

    ih264_iquant_itrans_residual_chroma_ft *pf_iquant_itrans_residual_chroma_4x4;

    ih264_iquant_itrans_residual_chroma_ft *pf_iquant_itrans_residual_chroma_4x4_dc;

    /* IT */
    ih264_iquant_itrans_ft *pf_iquant_itrans_luma_4x4;

    ih264_iquant_itrans_ft *pf_iquant_itrans_luma_4x4_dc;

    ih264_iquant_itrans_ft *pf_iquant_itrans_luma_8x8;

    ih264_iquant_itrans_ft *pf_iquant_itrans_luma_8x8_dc;

    ih264_iquant_itrans_chroma_ft *pf_iquant_itrans_chroma_4x4;

    ih264_iquant_itrans_chroma_ft *pf_iquant_itrans_chroma_4x4_dc;

    /**
     *SVC extension parsing strcture place holders
     */
    dec_nal_unit_svc_ext_params_t *ps_nal_svc_ext;
    dec_prefix_nal_unit_svc_ext_params_t s_pre_nal_unit_svc_ext;
    dec_svc_crop_wnd_offset_t *ps_crop_wnd_offset;
    UWORD8 *apu1_crop_wnd_flag[MAX_DEP_LYRS_IN_RES];

    /**
     *contexts for the CABAC related parsing
     */
    bin_ctxt_model_t *ps_base_mode_flag;
    bin_ctxt_model_t *ps_motion_prediction_flag_l0;
    bin_ctxt_model_t *ps_motion_prediction_flag_l1;
    bin_ctxt_model_t *ps_residual_prediction_flag;

    /**
     * Function pointers to read Params common to CAVLC and CABAC
     */
    WORD32(*pf_parse_inter_mb_svc_ext)
    (struct _SvcDecLyrStruct *ps_dec, dec_mb_info_t *ps_cur_mb_info,
     dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD32 u4_mb_num, UWORD32 u4_num_mbsNby2);

    WORD32(*pf_parse_inter_slice_svc_ext)
    (struct _SvcDecLyrStruct *ps_dec, dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

    /**
     * Function pointers to parse inter slice data
     */

    WORD32(*pf_parse_svc_inter_slice)
    (struct _SvcDecLyrStruct *ps_dec, dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice);

    /* inter layer precition buffers */

    /* 4x4 level */
    mv_pred_t *ps_il_pred_mv_bank_buf_base;

    /* 16x16 level */
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms_base;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms_frm_start;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms_cur_mb;
    UWORD16 u2_inter_lyr_mb_prms_stride;
    UWORD32 u4_inter_lyr_mb_prms_size; /* in Bytes */

    /* full frame size : -255 -255 */
    WORD16 *pi2_il_residual_resample_luma_base;
    WORD16 *pi2_il_residual_resample_chroma_base;
    WORD16 *pi2_il_residual_resample_mb_luma_frm_start;
    WORD16 *pi2_il_residual_resample_mb_chroma_frm_start;

    UWORD16 u2_residual_resample_luma_stride;
    UWORD16 u2_residual_resample_chroma_stride;
    UWORD32 u4_residual_resample_luma_size;   /* in Bytes */
    UWORD32 u4_residual_resample_chroma_size; /* in Bytes */

    mv_pred_t *ps_il_pred_mv_bank_buf_cur_mb;

    UWORD8 *pu1_crop_wnd_flag;
    /*
     * Layer info flag - Base layer; Intermediate Enhancement Layers; Target
     * Enhacement Layer.
     */
    UWORD8 u1_layer_identifier;
    /* layer id of the current layer */
    UWORD8 u1_layer_id;
    /* flag to indicate if spatial layers are dyadic */
    UWORD8 u1_dyadic_flag;
    /* flag to indicate if current layer is base layer */
    UWORD8 u1_base_res_flag;
    /* reference layer for inter layer prediction, no quality layers */
    UWORD8 u1_ref_layer_id;

    UWORD8 u1_restricted_res_change_flag;

    res_prms_t s_res_prms;

    void *pv_ref_lyr_offset;
    void *pv_mode_mv_sample_ctxt;
    void *pv_ii_pred_ctxt;
    void *pv_residual_sample_ctxt;

    void *pv_intra_sample_ctxt;
    /*!< projected locations buffer pointer exported by Intra Upsampling module
    for luma this buffer contains the projected offsets and window width in
    reference layer for each MB (in horizontal direction) of current resolution
    layer.*/
    ref_mb_map_t *ps_intsam_luma_map_horz;

    /*!< projected locations buffer pointer exported by Intra Upsampling module
      for chroma this buffer contains the projected offsets and window width in
      reference layer
      for each MB (in horizontal direction) of current resolution layer.*/
    ref_mb_map_t *ps_intsam_chroma_map_horz;

    /*!< projected locations  buffer pointer exported by Intra Upsampling module
       for luma this buffer contains the projected offsets and window width in
       reference layer for each MB (in vertical direction) of current resolution
       layer. */
    ref_mb_map_t *ps_intsam_luma_map_vert;

    /*!<  projected locations buffer pointer exported by Intra Upsampling module
       for chroma this buffer contains the projected offsets and window width in
       reference layer for each MB (in vertical direction) of current resolution
       layer.  */
    ref_mb_map_t *ps_intsam_chroma_map_vert;

    /*!< projected locations buffer pointer exported by Residual Upsampling module
    for luma. this buffer contains the projected offsets and window width in
    reference layer for each MB (in horizontal direction) of current resolution
    layer. */
    ref_mb_map_t *ps_ressam_luma_map_horz;

    /*!< projected locations buffer pointer exported by Residual Upsampling module
        for chroma. this buffer contains the projected offsets and window width in
       reference layer
        for each MB (in horizontal direction) of current resolution layer. */
    ref_mb_map_t *ps_ressam_chroma_map_horz;

    /*!< projected locationscbuffer pointercexported by Residual Upsampling
    modulec for chroma. this buffer contains the projected offsets and window
    width in reference layer for each MB (in vertical direction) ofv    current
    resolution layer. */
    ref_mb_map_t *ps_ressam_luma_map_vert;

    /*!< projected locationscbuffer pointerccexported by Residual Upsampling
    module for chroma.cthis buffer contains the projected offsets and window width
    in reference layer for each MB (in vertical direction) of current resolution
    layer.*/
    ref_mb_map_t *ps_ressam_chroma_map_vert;

    /* pointer to decoder layer referered by current layer */
    void *ps_dec_svc_ref_layer;
    /* pointer to master context */
    void *ps_svcd_ctxt;

    UWORD8 u1_inter_lyr_disable_dblk_filter_idc;
    WORD8 i1_inter_lyr_slice_alpha_c0_offset;
    WORD8 i1_inter_lyr_slice_beta_offset;

    UWORD8 *pu1_ii_resamp_buffer_luma;
    UWORD8 *pu1_ii_resamp_buffer_chroma;

    dec_slice_svc_ext_params_t s_svc_slice_params;
    dec_svc_seq_params_t *ps_subset_sps;
    dec_svc_seq_params_t *ps_cur_subset_sps;
    void *pv_scratch_subset_sps;

    /* Variables Required for N MB design */
    dec_svc_mb_info_t *ps_svc_nmb_info;

    dec_svc_mb_info_t *ps_svc_frm_mb_info;

    void (*pf_svc_compute_bs)(struct _SvcDecLyrStruct *ps_svc_lyr_dec,
                              struct _DecMbInfo *ps_cur_mb_info, const UWORD16 u2_mbxn_mb);

    UWORD16 *pu2_frm_res_luma_csbp;
    WORD32 i4_frm_res_luma_csbp_stride;

    UWORD8 *pu1_svc_base_mode_flag;
    WORD32 i4_frm_svc_base_mode_cabac_stride;
    WORD32 i4_frm_svc_base_mode_cabac_size;
    UWORD32 u4_pps_id_for_layer;
    UWORD8 u1_error_in_cur_frame;
    UWORD8 u1_res_init_done;
    WORD32 pic_width;
    WORD32 pic_height;
} svc_dec_lyr_struct_t;

typedef struct
{
    /* common parameters for all layers in SVC */
    UWORD32 u4_num_cores;
    IVD_ARCH_T e_processor_arch;
    IVD_SOC_T e_processor_soc;
    UWORD8 u1_target_layer_id;
    UWORD8 u1_cur_layer_id;

    /* dcode context for all layers in SVC */
    svc_dec_lyr_struct_t *ps_svc_dec_lyr;

    dec_pic_params_t *ps_pps;
    dec_seq_params_t *ps_sps;
    dec_svc_seq_params_t *ps_subset_sps;
    struct _sei *ps_sei;
    struct _sei *ps_sei_parse;

    /* attributes related to set tgt layer api func */
    WORD32 u1_tgt_dep_id;
    WORD32 u1_tgt_quality_id;
    WORD32 u1_tgt_temp_id;
    WORD32 u1_tgt_priority_id;

    ref_lyr_scaled_offset_t as_ref_lyr_offsets[MAX_NUM_RES_LYRS];

    void *pv_ref_lyr_offset;
    void *pv_mode_mv_sample_ctxt;
    void *pv_ii_pred_ctxt;
    void *pv_residual_sample_ctxt;
    void *pv_intra_sample_ctxt;

    void *pv_nal_parse_ctxt;
    non_vcl_nal_t s_non_vcl_nal; /*!< NON VCL nal structure */
    vcl_nal_t s_vcl_nal;         /*!< VCL nal structure */

    /*!< array to store the Did of bottom most layer in each resolution */
    WORD32 ai4_dq_id_map[MAX_NUM_RES_LYRS];
    WORD32 i4_error_code;
    void *pv_vcl_nal_buff;
    void *pv_non_vcl_nal_buff;

    /*!< array of structure to store the reference layer DQID,
    poc syntax and frame num, for each depedency id
    present in an access unit this will be used as reference
    for the next access unit */
    prev_au_prms_t as_au_prms_dep[MAX_DEPENDENCY_LYRS];
    /*!< array to store the pps id for each layer in a resolution */
    prev_au_sps_pps_t as_pps_sps_prev_au[MAX_TOTAL_LYRS];

    WORD32 i4_eos_flag;
    UWORD8 u1_prev_num_res_layers;
    UWORD32 u4_num_sps_ctr;
    UWORD32 u4_num_pps_ctr;
    UWORD8 u1_parse_nal_unit_error;
    UWORD8 u1_exit_till_next_IDR;
    UWORD8 u1_pre_parse_in_flush;
    WORD32 pic_width;
    WORD32 pic_height;
} svc_dec_ctxt_t;

#endif /*_ISVCD_STRUCTS_H_*/