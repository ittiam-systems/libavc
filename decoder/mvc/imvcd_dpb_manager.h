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

/*****************************************************************************/
/*                                                                           */
/*  File Name         : imvcd_nalu_parser.h                                  */
/*                                                                           */
/*  Description       : Functions for MVC NALU parsing                       */
/*                                                                           */
/*****************************************************************************/

#ifndef _IMVCD_DPB_MANAGER_H_
#define _IMVCD_DPB_MANAGER_H_
#include <stdbool.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "ih264_error.h"
#include "ih264_buf_mgr.h"
#include "ih264_disp_mgr.h"
#include "ih264d_dpb_manager.h"
#include "imvcd_defs.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

#define NUM_MMCO_CMD_IDS ((RESET_ALL_PICTURES) + 1)

typedef struct mvc_dpb_info_t
{
    mvc_au_buffer_t *ps_au_buf;

    struct mvc_dpb_info_t *ps_prev_short;

    struct mvc_dpb_info_t *ps_prev_long;

    struct field_t s_top_field;

    struct field_t s_bot_field;

    bool b_used_as_ref;

} mvc_dpb_info_t;

typedef struct display_buf_info_t
{
    WORD32 i4_poc;

    WORD32 i4_poc_buf_id;

    WORD32 i4_frame_num;
} display_buf_info_t;

typedef struct dpb_ivp_ctxt_t
{
    sps_mvc_ext_t *ps_sps_mvc_ext;

    nalu_mvc_ext_t *ps_nalu_mvc_exts;

    UWORD8 au1_au_buf_ids[MVC_MAX_REF_PICS];

    UWORD8 au1_mv_buf_ids[MVC_MAX_REF_PICS];

    UWORD32 u4_num_ivp_refs;
} dpb_ivp_ctxt_t;

typedef struct mvc_dpb_manager_t
{
    /** DPB in default index order */
    mvc_au_buffer_t *aps_def_dpb[MVC_MAX_REF_PICS];

    /** DPB in reordered index order, 0-fwd,1-bwd */
    mvc_au_buffer_t *aps_mod_dpb[2][MVC_MAX_REF_PICS];

    /** DPB in reordered index order, 0-fwd,1-bwd */
    mvc_au_buffer_t as_init_dpb[2][MVC_MAX_REF_PICS];

    /** Replicates view level data in 'aps_mod_dpb' */
    pic_buffer_t *aps_view_mod_dpb[2][MVC_MAX_REF_PICS];

    /** Replicates view level data in 'aps_init_dpb' */
    pic_buffer_t as_view_init_dpb[2][MVC_MAX_REF_PICS];

    mvc_dpb_info_t as_dpb_info[MVC_MAX_REF_PICS];

    display_buf_info_t as_display_buf_info[MAX_FRAMES];

    dpb_ivp_ctxt_t s_dpb_ivp_ctxt;

    mvc_dpb_info_t *ps_dpb_st_head;

    mvc_dpb_info_t *ps_dpb_lt_head;

    mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr;

    mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr;

    disp_mgr_t *ps_disp_buf_mgr;

    WORD32 ai4_gaps_start_frm_num[MAX_FRAMES];

    WORD32 ai4_gaps_end_frm_num[MAX_FRAMES];

    WORD8 ai1_gaps_per_seq[MAX_FRAMES];

    UWORD8 au1_num_active_st_refs[2];

    UWORD8 au1_num_active_lt_refs[2];

    WORD32 i4_max_pic_num;

    WORD32 i4_display_delay;

    WORD32 i4_cur_display_seq;

    UWORD16 u2_num_views;

    UWORD8 u1_num_st_ref_bufs;

    UWORD8 u1_num_lt_ref_bufs;

    UWORD8 u1_max_lt_frame_idx;

    UWORD8 u1_num_gaps;

    WORD8 i1_poc_buf_id_entries;

    WORD8 i1_gaps_deleted;

    UWORD8 u1_mmco_error_in_seq;

} mvc_dpb_manager_t;

/* Function declarations */
extern void imvcd_init_dpb_mgr(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                               mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr,
                               disp_mgr_t *ps_disp_buf_mgr);

extern WORD32 imvcd_dpb_assign_display_seq(mvc_dpb_manager_t *ps_dpb_mgr);

extern WORD32 imvcd_dpb_insert_pic_in_display_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                                   WORD32 i4_display_poc, UWORD32 u4_frame_num,
                                                   WORD32 i4_buf_id);

extern WORD32 imvcd_dpb_do_mmco_for_gaps(mvc_dpb_manager_t *ps_dpb_mgr, UWORD8 u1_num_ref_frames);

extern void imvcd_dpb_delete_nonref_nondisplay_pics(mvc_dpb_manager_t *ps_dpb_mgr);

extern void imvcd_reset_dpb(mvc_dpb_manager_t *ps_dpb_mgr);

extern void imvcd_dpb_release_display_bufs(mvc_dpb_manager_t *ps_dpb_mgr);

extern void imvcd_assign_pic_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_max_frame_num,
                                 WORD32 i4_cur_frame_num,
                                 bool b_are_gaps_in_frame_num_value_allowed);

extern UWORD8 imvcd_dpb_st_lt_deduplicator(mvc_dpb_manager_t *ps_dpb_mgr);

extern WORD32 imvcd_init_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                      nalu_mvc_ext_t *ps_cur_nalu_mvc_ext,
                                      mvc_au_buffer_t *ps_cur_au, UWORD16 u2_view_order_id);

extern WORD32 imvcd_dpb_reorder_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                             nalu_mvc_ext_t *ps_cur_nalu_mvc_ext,
                                             mvc_au_buffer_t *ps_cur_au,
                                             ref_pic_list_mod_data_t *ps_ref_pic_list_mod_data,
                                             UWORD16 u2_view_order_id);

extern WORD32 imvcd_dpb_insert_st_node(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buffer_t *ps_au_buf);

extern WORD32 imvcd_dpb_delete_st_node_or_make_lt(mvc_dpb_manager_t *ps_dpb_mgr,
                                                  WORD32 i4_frame_num, UWORD32 u4_lt_idx);

extern WORD32 imvcd_dpb_do_mmco(dpb_commands_t *ps_dpb_cmds, mvc_dpb_manager_t *ps_dpb_mgr,
                                mvc_au_buffer_t *ps_cur_au, UWORD8 u1_max_num_ref_frames,
                                UWORD8 u1_curr_pic_in_err);

extern WORD32 imvcd_dpb_update_default_index_list(mvc_dpb_manager_t *ps_dpb_mgr);

extern void imvcd_dpb_set_display_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_display_num);

extern void imvcd_dpb_set_max_pic_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_max_pic_num);

extern void imvcd_dpb_set_num_views(mvc_dpb_manager_t *ps_dpb_mgr, UWORD16 u2_num_views);

extern void imvcd_dpb_set_display_delay(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_display_delay);

extern void imvcd_dpb_init_au_bufs(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buffer_t *ps_cur_au);

extern void imvcd_dpb_init_view_bufs(mvc_dpb_manager_t *ps_dpb_mgr, UWORD16 u2_view_order_id,
                                     UWORD16 u2_view_id);

extern void imvcd_dpb_init_ivp_ctxt(mvc_dpb_manager_t *ps_dpb_mgr, sps_mvc_ext_t *ps_sps_mvc_ext,
                                    nalu_mvc_ext_t *ps_nalu_mvc_exts);

extern void imvcd_dpb_reset_ivp_ctxt(mvc_dpb_manager_t *ps_dpb_mgr);

extern pic_buffer_t **imvcd_dpb_get_view_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                                      UWORD16 u2_view_order_id, UWORD16 u2_view_id,
                                                      UWORD8 u1_pred_dir);

extern bool imvcd_dpb_is_diff_poc_valid(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_curr_poc);

#endif
