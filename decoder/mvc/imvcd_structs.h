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
#ifndef _IMVCD_STRUCTS_H_
#define _IMVCD_STRUCTS_H_

#include <stdbool.h>

#include "ih264_typedefs.h"
#include "imvcd.h"
#include "ih264_error.h"
#include "ih264_buf_mgr.h"
#include "ih264_disp_mgr.h"
#include "ih264d_sei.h"
#include "ih264d_structs.h"
#include "imvc_defs.h"
#include "imvc_structs.h"
#include "imvcd_defs.h"

/* structs */
typedef struct mvc_au_mv_pred_t
{
    mv_pred_t *aps_mvs[MAX_NUM_VIEWS];

    /* colZeroFlag | // 0th bit
       field_flag  | // 1st bit
       XX          | // 2:3 bit don't cares
       subMbMode   | // 4:5 bit
       MbMode      | // 6:7 bit */
    UWORD8 *apu1_mode_descriptors[MAX_NUM_VIEWS];

} mvc_au_mv_pred_t;

typedef struct ivp_data_t
{
    bool b_is_ivp_ref;

    /* Due to the structuring of dpb_mgr_t, */
    /* mvc_pic_buffer_t used for referencing ought to contain */
    /* all data in indices corresponding to view_id doing the referencing. */
    /* This struct, and this variable in particular, is used for identifying */
    /* the reference view's view_id */

    UWORD16 u2_ref_view_id;
} ivp_data_t;

typedef struct mvc_au_buffer_t
{
    /** pic_buffer for all views */
    yuv_buf_props_t as_view_buffers[MAX_NUM_VIEWS];

    /** display offsets for all views */
    offsets_t as_disp_offsets[MAX_NUM_VIEWS];

    ivp_data_t s_ivp_data;

    /** SEI data */
    sei s_sei_pic;

    /** AU MV Data */
    mvc_au_mv_pred_t *ps_au_mv_data;

    /* It will contain information about types of slices */
    UWORD32 au4_pack_slc_typ[MAX_NUM_VIEWS];

    /** Width of the display luma frame in pixels */
    UWORD16 u2_disp_width;

    /** Height of the display luma frame in pixels */
    UWORD16 u2_disp_height;

    /** Time at which frame has to be displayed */
    UWORD32 u4_time_stamp;

    /** (1: short 0: long) term ref pic */
    bool b_is_short_term_ref;

    /** frame / field / complementary field pair */
    UWORD8 u1_pic_type;

    /** Idx into the picBufAPI array */
    WORD32 i4_pic_buf_id;

    WORD32 i4_mv_buf_id;

    WORD32 i4_poc;

    WORD32 i4_frame_num;

    /* Derived based on '8.2.4.1' */
    WORD32 i4_pic_num;

    /** minPOC */
    WORD32 i4_avg_poc;

    /*Same as u1_pic_type..u1_pic_type gets overwritten whereas this doesnot get
    overwritten
    ...stores the pictype of frame/complementary field pair/ mbaff */
    UWORD8 u1_picturetype;

    UWORD8 u1_long_term_frm_idx;

    UWORD8 u1_long_term_pic_num;

    /* Refer to SEI table D-1 */
    UWORD8 u1_pic_struct;

} mvc_au_buffer_t;

typedef struct mvc_au_buf_mgr_t
{
    void *pv_mem;

    buf_mgr_t *ps_buf_mgr_ctxt;

    void *pv_au_buf_base;

    mvc_au_buffer_t *aps_buf_id_to_au_buf_map[MAX_DISP_BUFS_NEW];

    UWORD8 au1_au_buf_id_to_mv_buf_id_map[MAX_DISP_BUFS_NEW];

    UWORD8 au1_au_buf_ref_flag[MAX_DISP_BUFS_NEW];

} mvc_au_buf_mgr_t;

typedef struct mvc_au_mv_pred_buf_mgr_t
{
    void *pv_mem;

    buf_mgr_t *ps_buf_mgr_ctxt;

    void *pv_au_mv_pred_buf_base;

    mvc_au_mv_pred_t *aps_buf_id_to_mv_pred_buf_map[MAX_DISP_BUFS_NEW];

} mvc_au_mv_pred_buf_mgr_t;

typedef struct subset_sps_t
{
    dec_seq_params_t s_sps_data;

    sps_mvc_ext_t s_sps_mvc_ext;

    mvc_vui_ext_t s_mvc_vui_ext;

    offsets_t s_disp_offsets;

    UWORD8 u1_mvc_vui_parameters_present_flag;

} subset_sps_t;

typedef struct ref_pic_list_mod_data_t
{
    UWORD8 au1_num_active_refs[2];

    UWORD8 au1_ref_pic_list_modification_flag_lx[2];

    UWORD8 au1_modification_of_pic_nums_idc[2][MVC_MAX_REF_PICS + 1];

    WORD32 ai4_abs_diff_pic_num_minus1[2][MVC_MAX_REF_PICS + 1];

    WORD32 ai4_long_term_pic_num[2][MVC_MAX_REF_PICS + 1];

    WORD32 ai4_abs_diff_view_idx_minus1[2][MVC_MAX_REF_PICS + 1];
} ref_pic_list_mod_data_t;

typedef struct mvc_dec_ctxt_t
{
    dec_struct_t s_view_dec_ctxt;

    iv_mvc_yuv_buf_t s_out_buffer;

    /* Resolves circular dependency with mvc_dpb_manager_t */
    void *ps_dpb_mgr;

    subset_sps_t as_subset_sps[MAX_NUM_SEQ_PARAMS];

    /* Indexed via viewOrderID */
    nalu_mvc_ext_t as_nalu_mvc_ext[MAX_NUM_VIEWS];

    /* Indexed via viewOrderID */
    dec_slice_params_t as_slices[MAX_NUM_VIEWS];

    ref_pic_list_mod_data_t as_ref_pic_list_mod_data[MAX_NUM_VIEWS];

    subset_sps_t *aps_pps_id_to_subset_sps_map[MAX_NUM_PIC_PARAMS];

    disp_mgr_t s_mvc_disp_buf_mgr;

    mvc_au_buf_mgr_t s_mvc_au_buf_mgr;

    mvc_au_mv_pred_buf_mgr_t s_mvc_au_mv_pred_buf_mgr;

    mvc_au_buffer_t *ps_cur_au;

    AVC_EXT_NALU_ID_T ae_nalu_id[MAX_NUM_VIEWS];

    UWORD8 au1_nal_ref_idc[MAX_NUM_VIEWS];

    UWORD32 u4_num_aus_decoded;

    UWORD16 u2_num_views;

    UWORD16 u2_num_views_decoded;

    UWORD8 u1_num_sps;

    UWORD8 u1_num_subset_sps;

    UWORD8 u1_num_pps;

    bool b_header_only_decode;

    bool b_flush_enabled;

} mvc_dec_ctxt_t;

#endif
