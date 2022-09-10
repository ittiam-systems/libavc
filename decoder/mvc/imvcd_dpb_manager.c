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

#include "ih264_typedefs.h"
#include "ih264d_error_handler.h"
#include "imvcd_dpb_manager.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

void imvcd_dpb_set_display_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_display_num)
{
    ps_dpb_mgr->i4_cur_display_seq = i4_display_num;
}

void imvcd_dpb_set_max_pic_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_max_pic_num)
{
    ps_dpb_mgr->i4_max_pic_num = i4_max_pic_num;
}

void imvcd_dpb_set_num_views(mvc_dpb_manager_t *ps_dpb_mgr, UWORD16 u2_num_views)
{
    ps_dpb_mgr->u2_num_views = u2_num_views;
}

void imvcd_dpb_set_display_delay(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_display_delay)
{
    ps_dpb_mgr->i4_display_delay = i4_display_delay;
}

void imvcd_dpb_init_au_bufs(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buffer_t *ps_cur_au)
{
    WORD32 i;

    for(i = 0; i < 2; i++)
    {
        ps_dpb_mgr->as_init_dpb[i][0] = ps_cur_au[0];
    }
}

void imvcd_dpb_init_view_bufs(mvc_dpb_manager_t *ps_dpb_mgr, UWORD16 u2_view_order_id,
                              UWORD16 u2_view_id)
{
    WORD32 i;

    for(i = 0; i < 2; i++)
    {
        imvcd_convert_au_buf_to_view_buf(&ps_dpb_mgr->as_init_dpb[i][0],
                                         &ps_dpb_mgr->as_view_init_dpb[i][0], u2_view_order_id,
                                         u2_view_id);
    }
}

void imvcd_dpb_init_ivp_ctxt(mvc_dpb_manager_t *ps_dpb_mgr, sps_mvc_ext_t *ps_sps_mvc_ext,
                             nalu_mvc_ext_t *ps_nalu_mvc_exts)
{
    ps_dpb_mgr->s_dpb_ivp_ctxt.ps_nalu_mvc_exts = ps_nalu_mvc_exts;
    ps_dpb_mgr->s_dpb_ivp_ctxt.ps_sps_mvc_ext = ps_sps_mvc_ext;

    ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs = 0;
}

void imvcd_dpb_reset_ivp_ctxt(mvc_dpb_manager_t *ps_dpb_mgr)
{
    UWORD32 i;

    for(i = 0; i < ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs; i++)
    {
        ih264_buf_mgr_release(ps_dpb_mgr->ps_mvc_au_buf_mgr->ps_buf_mgr_ctxt,
                              ps_dpb_mgr->s_dpb_ivp_ctxt.au1_au_buf_ids[i],
                              BUF_MGR_REF | BUF_MGR_IO);

        ih264_buf_mgr_release(ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr->ps_buf_mgr_ctxt,
                              ps_dpb_mgr->s_dpb_ivp_ctxt.au1_mv_buf_ids[i],
                              BUF_MGR_REF | BUF_MGR_IO);
    }

    ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs = 0;
}

pic_buffer_t **imvcd_dpb_get_view_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                               UWORD16 u2_view_order_id, UWORD16 u2_view_id,
                                               UWORD8 u1_pred_dir)
{
    WORD32 i;

    UWORD8 u1_num_ref_bufs = ps_dpb_mgr->au1_num_active_st_refs[u1_pred_dir] +
                             ps_dpb_mgr->au1_num_active_lt_refs[u1_pred_dir];

    for(i = 0; i < u1_num_ref_bufs; i++)
    {
        imvcd_convert_au_buf_to_view_buf(ps_dpb_mgr->aps_mod_dpb[u1_pred_dir][i],
                                         &ps_dpb_mgr->as_view_init_dpb[u1_pred_dir][i],
                                         u2_view_order_id, u2_view_id);

        ps_dpb_mgr->aps_view_mod_dpb[u1_pred_dir][i] =
            &ps_dpb_mgr->as_view_init_dpb[u1_pred_dir][i];
    }

    return ps_dpb_mgr->aps_view_mod_dpb[u1_pred_dir];
}

void imvcd_init_dpb_mgr(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                        mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr,
                        disp_mgr_t *ps_disp_buf_mgr)
{
    WORD32 i, j, k, l;

    mvc_dpb_info_t *ps_dpb_info = ps_dpb_mgr->as_dpb_info;

    for(i = 0; i < 2; i++)
    {
        mvc_au_buffer_t *ps_init_dpb = ps_dpb_mgr->as_init_dpb[i];
        pic_buffer_t *ps_view_init_dpb = ps_dpb_mgr->as_view_init_dpb[i];

        for(j = 0; j < MVC_MAX_REF_PICS; j++)
        {
            for(k = 0; k < MAX_NUM_VIEWS; k++)
            {
                for(l = 0; l < NUM_COMPONENTS; l++)
                {
                    ps_init_dpb->as_view_buffers[k].as_component_bufs[l].pv_data = NULL;
                }
            }

            ps_view_init_dpb->pu1_buf1 = NULL;
            ps_view_init_dpb->pu1_buf2 = NULL;
            ps_view_init_dpb->pu1_buf3 = NULL;

            ps_dpb_mgr->aps_mod_dpb[i][j] = ps_init_dpb;
            ps_dpb_mgr->aps_view_mod_dpb[i][j] = ps_view_init_dpb;

            ps_init_dpb++;
            ps_view_init_dpb++;
        }
    }

    for(i = 0; i < MVC_MAX_REF_PICS; i++)
    {
        ps_dpb_info[i].b_used_as_ref = false;
        ps_dpb_info[i].ps_prev_short = NULL;
        ps_dpb_info[i].ps_prev_long = NULL;
        ps_dpb_info[i].ps_au_buf = NULL;
        ps_dpb_info[i].s_top_field.u1_reference_info = UNUSED_FOR_REF;
        ps_dpb_info[i].s_bot_field.u1_reference_info = UNUSED_FOR_REF;
        ps_dpb_info[i].s_top_field.u1_long_term_frame_idx = MVC_MAX_REF_PICS + 1;
        ps_dpb_info[i].s_bot_field.u1_long_term_frame_idx = MVC_MAX_REF_PICS + 1;
    }

    ps_dpb_mgr->u1_num_st_ref_bufs = ps_dpb_mgr->u1_num_lt_ref_bufs = 0;
    ps_dpb_mgr->ps_dpb_st_head = NULL;
    ps_dpb_mgr->ps_dpb_lt_head = NULL;
    ps_dpb_mgr->i1_gaps_deleted = 0;
    ps_dpb_mgr->i1_poc_buf_id_entries = 0;
    ps_dpb_mgr->u1_mmco_error_in_seq = 0;
    ps_dpb_mgr->u1_num_gaps = 0;
    ps_dpb_mgr->i4_display_delay = 0;
    ps_dpb_mgr->i4_cur_display_seq = 0;

    for(i = 0; i < MAX_FRAMES; i++)
    {
        ps_dpb_mgr->ai4_gaps_start_frm_num[i] = INVALID_FRAME_NUM;
        ps_dpb_mgr->ai4_gaps_end_frm_num[i] = 0;
        ps_dpb_mgr->ai1_gaps_per_seq[i] = 0;
        ps_dpb_mgr->as_display_buf_info[i].i4_poc_buf_id = -1;
        ps_dpb_mgr->as_display_buf_info[i].i4_poc = INT32_MAX;
        ps_dpb_mgr->as_display_buf_info[i].i4_frame_num = 0;
    }

    ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs = 0;
    ps_dpb_mgr->s_dpb_ivp_ctxt.ps_nalu_mvc_exts = NULL;
    ps_dpb_mgr->s_dpb_ivp_ctxt.ps_sps_mvc_ext = NULL;

    ps_dpb_mgr->ps_mvc_au_buf_mgr = ps_mvc_au_buf_mgr;
    ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr = ps_mvc_au_mv_pred_buf_mgr;
    ps_dpb_mgr->ps_disp_buf_mgr = ps_disp_buf_mgr;
}

WORD32 imvcd_dpb_assign_display_seq(mvc_dpb_manager_t *ps_dpb_mgr)
{
    WORD32 i;

    display_buf_info_t *ps_display_buf_info = ps_dpb_mgr->as_display_buf_info;

    WORD32 i4_min_poc = INT32_MAX;
    WORD32 i4_min_poc_buf_id = -1;
    WORD32 i4_min_index = -1;

    if(ps_dpb_mgr->i1_poc_buf_id_entries >= ps_dpb_mgr->i4_display_delay)
    {
        for(i = 0; i < MAX_FRAMES; i++)
        {
            if((-1 != ps_display_buf_info[i].i4_poc_buf_id) &&
               (DO_NOT_DISP != ps_display_buf_info[i].i4_poc_buf_id))
            {
                /* Checking for <= is necessary to handle cases where there is one
                   valid buffer with poc set to 0x7FFFFFFF. */
                if(ps_display_buf_info[i].i4_poc <= i4_min_poc)
                {
                    i4_min_poc = ps_display_buf_info[i].i4_poc;
                    i4_min_poc_buf_id = ps_display_buf_info[i].i4_poc_buf_id;
                    i4_min_index = i;
                }
            }
        }

        if((i4_min_index != -1) && (DO_NOT_DISP != i4_min_poc_buf_id))
        {
            ps_dpb_mgr->i4_cur_display_seq++;

            ih264_disp_mgr_add(
                ps_dpb_mgr->ps_disp_buf_mgr, i4_min_poc_buf_id, ps_dpb_mgr->i4_cur_display_seq,
                ps_dpb_mgr->ps_mvc_au_buf_mgr->aps_buf_id_to_au_buf_map[i4_min_poc_buf_id]);

            ps_display_buf_info[i4_min_index].i4_poc_buf_id = -1;
            ps_display_buf_info[i4_min_index].i4_poc = 0x7fffffff;

            ps_dpb_mgr->i1_poc_buf_id_entries--;
        }
        else if(DO_NOT_DISP == i4_min_poc_buf_id)
        {
            return ERROR_GAPS_IN_FRM_NUM;
        }
    }

    return OK;
}

WORD32 imvcd_dpb_insert_pic_in_display_list(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_display_poc,
                                            UWORD32 u4_frame_num, WORD32 i4_buf_id)
{
    WORD32 i;

    display_buf_info_t *ps_display_buf_info = ps_dpb_mgr->as_display_buf_info;

    for(i = 0; i < MAX_FRAMES; i++)
    {
        /* Find an empty slot */
        if(ps_display_buf_info[i].i4_poc_buf_id == -1)
        {
            if(GAP_FRAME_NUM == ps_display_buf_info[i].i4_frame_num)
            {
                ps_dpb_mgr->i1_gaps_deleted--;
            }
            else
            {
                ps_dpb_mgr->i1_poc_buf_id_entries++;
            }

            ps_display_buf_info[i].i4_poc_buf_id = i4_buf_id;
            ps_display_buf_info[i].i4_poc = i4_display_poc;
            ps_display_buf_info[i].i4_frame_num = u4_frame_num;

            break;
        }
    }

    if(MAX_FRAMES == i)
    {
        return ERROR_GAPS_IN_FRM_NUM;
    }

    return OK;
}

static WORD32 imvcd_dpb_delete_gap_frm_sliding(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_pic_num,
                                               UWORD8 *pu1_del_node)
{
    WORD32 i, j, j_min;
    WORD8 i1_gap_idx;
    WORD32 *pi4_gaps_start_frm_num, *pi4_gaps_end_frm_num, i4_gap_frame_num;
    WORD32 i4_start_frm_num, i4_end_frm_num;
    WORD32 i4_max_pic_num;
    WORD32 i4_frm_num, i4_gap_frm_num_min;

    /* find the least frame num from gaps and current DPB node    */
    /* Delete the least one                                       */
    *pu1_del_node = 1;

    if(0 == ps_dpb_mgr->u1_num_gaps)
    {
        return OK;
    }

    pi4_gaps_start_frm_num = ps_dpb_mgr->ai4_gaps_start_frm_num;
    pi4_gaps_end_frm_num = ps_dpb_mgr->ai4_gaps_end_frm_num;
    i4_gap_frame_num = INVALID_FRAME_NUM;
    i4_max_pic_num = ps_dpb_mgr->i4_max_pic_num;

    i1_gap_idx = -1;

    if(INVALID_FRAME_NUM != i4_pic_num)
    {
        i4_gap_frame_num = i4_pic_num;

        for(i = 0; i < MAX_FRAMES; i++)
        {
            i4_start_frm_num = pi4_gaps_start_frm_num[i];

            if(INVALID_FRAME_NUM != i4_start_frm_num)
            {
                i4_end_frm_num = pi4_gaps_end_frm_num[i];

                if(i4_end_frm_num < i4_max_pic_num)
                {
                    if(i4_start_frm_num <= i4_gap_frame_num)
                    {
                        i4_gap_frame_num = i4_start_frm_num;
                        i1_gap_idx = i;
                    }
                }
                else
                {
                    if(((i4_start_frm_num <= i4_gap_frame_num) &&
                        (i4_gap_frame_num <= i4_max_pic_num)) ||
                       ((i4_start_frm_num >= i4_gap_frame_num) &&
                        ((i4_gap_frame_num + i4_max_pic_num) >= i4_end_frm_num)))
                    {
                        i4_gap_frame_num = i4_start_frm_num;
                        i1_gap_idx = i;
                    }
                }
            }
        }
    }
    else
    {
        /* no valid short term buffers, delete one gap from the least start */
        /* of gap sequence                                                  */
        i4_gap_frame_num = pi4_gaps_start_frm_num[0];
        i1_gap_idx = 0;

        for(i = 1; i < MAX_FRAMES; i++)
        {
            if(INVALID_FRAME_NUM != pi4_gaps_start_frm_num[i])
            {
                if(pi4_gaps_start_frm_num[i] < i4_gap_frame_num)
                {
                    i4_gap_frame_num = pi4_gaps_start_frm_num[i];
                    i1_gap_idx = i;
                }
            }
        }
        if(INVALID_FRAME_NUM == i4_gap_frame_num)
        {
            return ERROR_DBP_MANAGER_T;
        }
    }

    if(-1 != i1_gap_idx)
    {
        /* find least frame_num in the poc_map, which is in this range */
        i4_start_frm_num = pi4_gaps_start_frm_num[i1_gap_idx];

        if(i4_start_frm_num < 0)
        {
            i4_start_frm_num += i4_max_pic_num;
        }

        i4_end_frm_num = pi4_gaps_end_frm_num[i1_gap_idx];

        if(i4_end_frm_num < 0)
        {
            i4_end_frm_num += i4_max_pic_num;
        }

        i4_gap_frm_num_min = INT32_MIN;
        j_min = MAX_FRAMES;

        for(j = 0; j < MAX_FRAMES; j++)
        {
            i4_frm_num = ps_dpb_mgr->as_display_buf_info[j].i4_frame_num;

            if((i4_start_frm_num <= i4_frm_num) && (i4_end_frm_num >= i4_frm_num))
            {
                if(i4_frm_num < i4_gap_frm_num_min)
                {
                    j_min = j;
                    i4_gap_frm_num_min = i4_frm_num;
                }
            }
        }

        if(j_min != MAX_FRAMES)
        {
            ps_dpb_mgr->as_display_buf_info[j_min].i4_poc_buf_id = -1;
            ps_dpb_mgr->as_display_buf_info[j_min].i4_poc = 0x7fffffff;
            ps_dpb_mgr->as_display_buf_info[j_min].i4_frame_num = GAP_FRAME_NUM;

            ps_dpb_mgr->i1_gaps_deleted++;
            ps_dpb_mgr->ai1_gaps_per_seq[i1_gap_idx]--;
            ps_dpb_mgr->u1_num_gaps--;
            *pu1_del_node = 0;

            if(0 == ps_dpb_mgr->ai1_gaps_per_seq[i1_gap_idx])
            {
                ps_dpb_mgr->ai4_gaps_start_frm_num[i1_gap_idx] = INVALID_FRAME_NUM;
                ps_dpb_mgr->ai4_gaps_end_frm_num[i1_gap_idx] = 0;
            }
        }
    }

    return OK;
}

WORD32 imvcd_dpb_do_mmco_for_gaps(mvc_dpb_manager_t *ps_dpb_mgr, UWORD8 u1_num_ref_frames)
{
    mvc_dpb_info_t *ps_next_dpb;

    WORD32 i;
    WORD32 i4_error_code;
    UWORD8 u1_num_gaps;
    UWORD8 u1_num_st_ref_bufs, u1_num_lt_ref_bufs, u1_del_node;

    WORD32 i4_frame_gaps = 1;

    // Sliding window - implements 8.2.5.3, flush out buffers
    u1_num_st_ref_bufs = ps_dpb_mgr->u1_num_st_ref_bufs;
    u1_num_lt_ref_bufs = ps_dpb_mgr->u1_num_lt_ref_bufs;

    while(1)
    {
        u1_num_gaps = ps_dpb_mgr->u1_num_gaps;

        if((u1_num_st_ref_bufs + u1_num_lt_ref_bufs + u1_num_gaps + i4_frame_gaps) >
           u1_num_ref_frames)
        {
            if(0 == (u1_num_st_ref_bufs + u1_num_gaps))
            {
                i4_frame_gaps = 0;

                ps_dpb_mgr->u1_num_gaps = (u1_num_ref_frames - u1_num_lt_ref_bufs);
            }
            else
            {
                u1_del_node = 1;
                ps_next_dpb = ps_dpb_mgr->ps_dpb_st_head;

                if(u1_num_st_ref_bufs > 1)
                {
                    for(i = 1; i < (u1_num_st_ref_bufs - 1); i++)
                    {
                        if(ps_next_dpb == NULL)
                        {
                            return ERROR_DBP_MANAGER_T;
                        }

                        ps_next_dpb = ps_next_dpb->ps_prev_short;
                    }

                    if(ps_next_dpb->ps_prev_short->ps_prev_short != NULL)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    if(u1_num_gaps)
                    {
                        i4_error_code = imvcd_dpb_delete_gap_frm_sliding(
                            ps_dpb_mgr, ps_next_dpb->ps_prev_short->ps_au_buf->i4_pic_num,
                            &u1_del_node);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }

                    if(u1_del_node)
                    {
                        u1_num_st_ref_bufs--;
                        ps_next_dpb->ps_prev_short->b_used_as_ref = false;
                        ps_next_dpb->ps_prev_short->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                        ps_next_dpb->ps_prev_short->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                        imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                            ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                            ps_next_dpb->ps_prev_short->ps_au_buf->i4_pic_buf_id);

                        ps_next_dpb->ps_prev_short->ps_au_buf = NULL;
                        ps_next_dpb->ps_prev_short = NULL;
                    }
                }
                else
                {
                    if(u1_num_st_ref_bufs)
                    {
                        if(u1_num_gaps)
                        {
                            i4_error_code = imvcd_dpb_delete_gap_frm_sliding(
                                ps_dpb_mgr, ps_next_dpb->ps_au_buf->i4_pic_num, &u1_del_node);

                            if(i4_error_code != OK)
                            {
                                return i4_error_code;
                            }
                        }

                        if(u1_del_node)
                        {
                            u1_num_st_ref_bufs--;
                            ps_next_dpb->b_used_as_ref = false;
                            ps_next_dpb->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                            ps_next_dpb->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                            imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                                ps_next_dpb->ps_au_buf->i4_pic_buf_id);

                            ps_next_dpb->ps_au_buf = NULL;
                            ps_next_dpb = NULL;
                            ps_dpb_mgr->ps_dpb_st_head = NULL;
                            ps_dpb_mgr->u1_num_st_ref_bufs = u1_num_st_ref_bufs;
                        }
                    }
                    else
                    {
                        i4_error_code = imvcd_dpb_delete_gap_frm_sliding(
                            ps_dpb_mgr, INVALID_FRAME_NUM, &u1_del_node);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }

                        if(u1_del_node)
                        {
                            return ERROR_DBP_MANAGER_T;
                        }
                    }
                }
            }
        }
        else
        {
            ps_dpb_mgr->u1_num_gaps += i4_frame_gaps;

            break;
        }
    }

    ps_dpb_mgr->u1_num_st_ref_bufs = u1_num_st_ref_bufs;

    return OK;
}

void imvcd_dpb_delete_nonref_nondisplay_pics(mvc_dpb_manager_t *ps_dpb_mgr)
{
    WORD32 i;

    display_buf_info_t *ps_display_buf_info = ps_dpb_mgr->as_display_buf_info;

    /* remove all gaps marked as unused for ref */
    for(i = 0; (i < MAX_FRAMES) && ps_dpb_mgr->i1_gaps_deleted; i++)
    {
        if(GAP_FRAME_NUM == ps_display_buf_info[i].i4_frame_num)
        {
            ps_dpb_mgr->i1_gaps_deleted--;
            ps_dpb_mgr->i1_poc_buf_id_entries--;
            ps_display_buf_info[i].i4_poc_buf_id = -1;
            ps_display_buf_info[i].i4_poc = 0x7fffffff;
            ps_display_buf_info[i].i4_frame_num = 0;
        }
    }
}

void imvcd_reset_dpb(mvc_dpb_manager_t *ps_dpb_mgr)
{
    WORD32 i;

    mvc_dpb_info_t *ps_dpb_info = ps_dpb_mgr->as_dpb_info;

    for(i = 0; i < MVC_MAX_REF_PICS; i++)
    {
        if(ps_dpb_info[i].b_used_as_ref)
        {
            imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                ps_dpb_info[i].ps_au_buf->i4_pic_buf_id);

            ps_dpb_info[i].b_used_as_ref = false;
            ps_dpb_info[i].ps_prev_short = NULL;
            ps_dpb_info[i].ps_prev_long = NULL;
            ps_dpb_info[i].ps_au_buf = NULL;
            ps_dpb_info[i].s_top_field.u1_reference_info = UNUSED_FOR_REF;
            ps_dpb_info[i].s_bot_field.u1_reference_info = UNUSED_FOR_REF;
            ps_dpb_info[i].s_top_field.u1_long_term_frame_idx = MVC_MAX_REF_PICS + 1;
            ps_dpb_info[i].s_bot_field.u1_long_term_frame_idx = MVC_MAX_REF_PICS + 1;
        }
    }

    ps_dpb_mgr->u1_num_st_ref_bufs = ps_dpb_mgr->u1_num_lt_ref_bufs = 0;
    ps_dpb_mgr->ps_dpb_st_head = NULL;
    ps_dpb_mgr->ps_dpb_lt_head = NULL;
    ps_dpb_mgr->u1_mmco_error_in_seq = 0;

    /* release all gaps */
    ps_dpb_mgr->u1_num_gaps = 0;

    for(i = 0; i < MAX_FRAMES; i++)
    {
        ps_dpb_mgr->ai4_gaps_start_frm_num[i] = INVALID_FRAME_NUM;
        ps_dpb_mgr->ai4_gaps_end_frm_num[i] = 0;
        ps_dpb_mgr->ai1_gaps_per_seq[i] = 0;
    }
}

void imvcd_dpb_release_display_bufs(mvc_dpb_manager_t *ps_dpb_mgr)
{
    WORD32 i, j;

    display_buf_info_t *ps_display_buf_info = ps_dpb_mgr->as_display_buf_info;

    WORD32 i4_min_poc = 0x7fffffff;
    WORD32 i4_min_poc_buf_id = 0;
    WORD32 i4_min_index = 0;

    imvcd_dpb_delete_nonref_nondisplay_pics(ps_dpb_mgr);

    for(j = 0; j < ps_dpb_mgr->i1_poc_buf_id_entries; j++)
    {
        i4_min_poc = 0x7fffffff;

        for(i = 0; i < MAX_FRAMES; i++)
        {
            if(ps_display_buf_info[i].i4_poc_buf_id != -1)
            {
                /* Checking for <= is necessary to handle cases where there is one
                   valid buffer with poc set to 0x7FFFFFFF. */
                if(ps_display_buf_info[i].i4_poc <= i4_min_poc)
                {
                    i4_min_poc = ps_display_buf_info[i].i4_poc;
                    i4_min_poc_buf_id = ps_display_buf_info[i].i4_poc_buf_id;
                    i4_min_index = i;
                }
            }
        }

        if(DO_NOT_DISP != i4_min_poc_buf_id)
        {
            ps_dpb_mgr->i4_cur_display_seq++;

            ih264_disp_mgr_add(
                ps_dpb_mgr->ps_disp_buf_mgr, i4_min_poc_buf_id, ps_dpb_mgr->i4_cur_display_seq,
                ps_dpb_mgr->ps_mvc_au_buf_mgr->aps_buf_id_to_au_buf_map[i4_min_poc_buf_id]);

            ps_display_buf_info[i4_min_index].i4_poc_buf_id = -1;
            ps_display_buf_info[i4_min_index].i4_poc = 0x7fffffff;
            ps_display_buf_info[i4_min_index].i4_frame_num = 0;
        }
        else
        {
            ps_display_buf_info[i4_min_index].i4_poc_buf_id = -1;
            ps_display_buf_info[i4_min_index].i4_poc = 0x7fffffff;
            ps_display_buf_info[i4_min_index].i4_frame_num = 0;
        }
    }

    ps_dpb_mgr->i1_poc_buf_id_entries = 0;
}

void imvcd_assign_pic_num(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_max_frame_num,
                          WORD32 i4_cur_frame_num, bool b_are_gaps_in_frame_num_value_allowed)
{
    mvc_dpb_info_t *ps_next_dpb;

    WORD32 i;
    WORD32 i4_ref_frame_num;

    /* Start from ST head */
    ps_next_dpb = ps_dpb_mgr->ps_dpb_st_head;

    for(i = 0; i < ps_dpb_mgr->u1_num_st_ref_bufs; i++)
    {
        WORD32 i4_pic_num;

        i4_ref_frame_num = ps_next_dpb->ps_au_buf->i4_pic_num;

        if(i4_ref_frame_num > i4_cur_frame_num)
        {
            i4_pic_num = i4_ref_frame_num - i4_max_frame_num;
        }
        else
        {
            i4_pic_num = i4_ref_frame_num;
        }

        ps_next_dpb->ps_au_buf->i4_pic_num = i4_pic_num;

        ps_next_dpb = ps_next_dpb->ps_prev_short;
    }

    if(b_are_gaps_in_frame_num_value_allowed && ps_dpb_mgr->u1_num_gaps)
    {
        WORD32 i4_start_frm, i4_end_frm;

        /* Assign pic numbers for gaps */
        for(i = 0; i < MAX_FRAMES; i++)
        {
            i4_start_frm = ps_dpb_mgr->ai4_gaps_start_frm_num[i];

            if(i4_start_frm != INVALID_FRAME_NUM)
            {
                if(i4_start_frm > i4_cur_frame_num)
                {
                    /* gap's frame_num is before Current frame_num in
                     decode order */
                    i4_start_frm -= i4_max_frame_num;
                }

                ps_dpb_mgr->ai4_gaps_start_frm_num[i] = i4_start_frm;
                i4_end_frm = ps_dpb_mgr->ai4_gaps_end_frm_num[i];

                if(i4_end_frm > i4_cur_frame_num)
                {
                    /* gap's frame_num is before Current frame_num in
                     decode order */
                    i4_end_frm -= i4_max_frame_num;
                }

                ps_dpb_mgr->ai4_gaps_end_frm_num[i] = i4_end_frm;
            }
        }
    }
}

/* If there is common node in both lt_list and st_list, then delete it from */
/* st_list */
UWORD8 imvcd_dpb_st_lt_deduplicator(mvc_dpb_manager_t *ps_dpb_mgr)
{
    mvc_dpb_info_t *ps_dpb_lt_head = ps_dpb_mgr->ps_dpb_lt_head;
    mvc_dpb_info_t *ps_lt_curr_dpb = ps_dpb_mgr->ps_dpb_lt_head;
    mvc_dpb_info_t *ps_dpb_st_head = ps_dpb_mgr->ps_dpb_st_head;

    UWORD8 u1_no_of_nodes_deleted = 0;
    UWORD8 u1_num_lt_refs = ps_dpb_mgr->u1_num_lt_ref_bufs;
    UWORD8 u1_num_st_refs = ps_dpb_mgr->u1_num_st_ref_bufs;

    while(u1_num_lt_refs && ps_dpb_lt_head)
    {
        if(ps_dpb_st_head &&
           ((ps_dpb_lt_head->s_bot_field.u1_reference_info |
             ps_dpb_lt_head->s_top_field.u1_reference_info) == (IS_SHORT_TERM | IS_LONG_TERM)))
        {
            mvc_dpb_info_t *ps_st_next_dpb = ps_dpb_st_head;
            mvc_dpb_info_t *ps_st_curr_dpb = ps_dpb_st_head;

            while(u1_num_st_refs && ps_st_curr_dpb)
            {
                if(ps_st_curr_dpb == ps_lt_curr_dpb)
                {
                    if(u1_num_st_refs == ps_dpb_mgr->u1_num_st_ref_bufs)
                    {
                        ps_dpb_mgr->ps_dpb_st_head = ps_dpb_mgr->ps_dpb_st_head->ps_prev_short;
                        ps_st_curr_dpb = ps_dpb_mgr->ps_dpb_st_head;
                    }
                    else
                    {
                        ps_st_next_dpb->ps_prev_short = ps_st_curr_dpb->ps_prev_short;
                    }

                    ps_dpb_mgr->u1_num_st_ref_bufs--;
                    u1_no_of_nodes_deleted++;

                    break;
                }

                ps_st_next_dpb = ps_st_curr_dpb;
                ps_st_curr_dpb = ps_st_curr_dpb->ps_prev_short;
                u1_num_st_refs--;
            }
        }

        ps_lt_curr_dpb = ps_lt_curr_dpb->ps_prev_long;
        u1_num_lt_refs--;
    }

    return u1_no_of_nodes_deleted;
}

static int qsort_pic_num_compare(const void *pv_au1, const void *pv_au2)
{
    return ((mvc_au_buffer_t **) pv_au1)[0]->i4_pic_num -
           ((mvc_au_buffer_t **) pv_au2)[0]->i4_pic_num;
}

static int qsort_poc_compare(const void *pv_au1, const void *pv_au2)
{
    return ((mvc_au_buffer_t **) pv_au1)[0]->i4_poc - ((mvc_au_buffer_t **) pv_au2)[0]->i4_poc;
}

static int qsort_lt_idx_compare(const void *pv_au1, const void *pv_au2)
{
    return ((mvc_au_buffer_t **) pv_au1)[0]->u1_long_term_frm_idx -
           ((mvc_au_buffer_t **) pv_au2)[0]->u1_long_term_frm_idx;
}

WORD32 imvcd_init_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr, nalu_mvc_ext_t *ps_cur_nalu_mvc_ext,
                               mvc_au_buffer_t *ps_cur_au, UWORD16 u2_view_order_id)
{
    mvc_dpb_info_t *ps_ref_au_data;
    mvc_ivp_ref_data_t *ps_mvc_ivp_ref_data;

    WORD32 i, j;

    mvc_au_buffer_t *aps_st_au_bufs[MVC_MAX_REF_PICS] = {NULL};
    mvc_au_buffer_t *aps_lt_au_bufs[MVC_MAX_REF_PICS] = {NULL};
    mvc_au_buffer_t *aps_ref_pic_buf_lx[2] = {ps_dpb_mgr->as_init_dpb[0],
                                              ps_dpb_mgr->as_init_dpb[1]};
    sps_mvc_ext_t *ps_sps_mvc_ext = ps_dpb_mgr->s_dpb_ivp_ctxt.ps_sps_mvc_ext;

    UWORD16 u2_view_id = ps_cur_nalu_mvc_ext->u2_view_id;
    WORD32 i4_cur_poc = ps_cur_au->i4_poc;
    bool b_is_b_pic = !!(ps_cur_au->au4_pack_slc_typ[u2_view_order_id] & B_SLC_BIT);
    UWORD8 *pu1_num_short_term_refs = ps_dpb_mgr->au1_num_active_st_refs;
    UWORD8 *pu1_num_long_term_refs = ps_dpb_mgr->au1_num_active_lt_refs;
    UWORD8 au1_total_num_refs[2] = {0};

    memset(pu1_num_short_term_refs, 0, 2 * sizeof(pu1_num_short_term_refs[0]));

    memset(pu1_num_long_term_refs, 0, 2 * sizeof(pu1_num_long_term_refs[0]));

    ps_ref_au_data = ps_dpb_mgr->ps_dpb_st_head;

    for(i = 0; i < ps_dpb_mgr->u1_num_st_ref_bufs; i++)
    {
        aps_st_au_bufs[i] = ps_ref_au_data->ps_au_buf;
        ps_ref_au_data = ps_ref_au_data->ps_prev_short;
    }

    qsort(aps_st_au_bufs, ps_dpb_mgr->u1_num_st_ref_bufs, sizeof(aps_st_au_bufs[0]),
          b_is_b_pic ? qsort_poc_compare : qsort_pic_num_compare);

    ps_ref_au_data = ps_dpb_mgr->ps_dpb_lt_head;

    for(i = 0; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
    {
        aps_lt_au_bufs[i] = ps_ref_au_data->ps_au_buf;
        ps_ref_au_data = ps_ref_au_data->ps_prev_long;
    }

    qsort(aps_lt_au_bufs, ps_dpb_mgr->u1_num_lt_ref_bufs, sizeof(aps_lt_au_bufs[0]),
          qsort_lt_idx_compare);

    if(b_is_b_pic)
    {
        for(i = 0; i < ps_dpb_mgr->u1_num_st_ref_bufs; i++)
        {
            if(aps_st_au_bufs[i]->i4_poc >= i4_cur_poc)
            {
                break;
            }
        }

        for(j = i - 1; j >= 0; j--)
        {
            aps_ref_pic_buf_lx[0][0] = aps_st_au_bufs[j][0];
            aps_ref_pic_buf_lx[0]++;
            pu1_num_short_term_refs[0]++;
        }

        for(j = i; j < ps_dpb_mgr->u1_num_st_ref_bufs; j++)
        {
            aps_ref_pic_buf_lx[1][0] = aps_st_au_bufs[j][0];
            aps_ref_pic_buf_lx[1]++;
            pu1_num_short_term_refs[1]++;
        }
    }
    else
    {
        for(i = ps_dpb_mgr->u1_num_st_ref_bufs - 1; i >= 0; i--)
        {
            aps_ref_pic_buf_lx[0][0] = aps_st_au_bufs[i][0];
            aps_ref_pic_buf_lx[0]++;
            pu1_num_short_term_refs[0]++;
        }
    }

    for(i = 0; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
    {
        for(j = 0; j < 1 + ((WORD32) b_is_b_pic); j++)
        {
            aps_ref_pic_buf_lx[j][0] = aps_lt_au_bufs[i][0];
            aps_ref_pic_buf_lx[j]->u1_long_term_pic_num =
                aps_ref_pic_buf_lx[j]->u1_long_term_frm_idx;
            aps_ref_pic_buf_lx[j]++;
            pu1_num_long_term_refs[j]++;
        }
    }

    if(0 != u2_view_order_id)
    {
        ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs = 0;

        for(i = 0; i < 1 + ((WORD32) b_is_b_pic); i++)
        {
            WORD32 i4_num_refs;

            if(ps_cur_nalu_mvc_ext->u1_anchor_pic_flag)
            {
                ps_mvc_ivp_ref_data = &ps_sps_mvc_ext->as_anchor_ref_data[i][u2_view_order_id];
            }
            else
            {
                ps_mvc_ivp_ref_data = &ps_sps_mvc_ext->as_non_anchor_ref_data[i][u2_view_order_id];
            }

            i4_num_refs = ps_mvc_ivp_ref_data->u1_num_refs;

            for(j = 0; j < i4_num_refs; j++)
            {
                mvc_au_buffer_t *ps_au_buf;
                mvc_au_mv_pred_t *ps_au_mv_data;
                nalu_mvc_ext_t *ps_ref_nalu_mvc_ext;

                WORD32 i4_pic_buf_id;
                WORD32 i4_mv_buf_id;

                UWORD16 u2_ref_view_id = ps_mvc_ivp_ref_data->au2_ref_view_ids[j];

                ps_ref_nalu_mvc_ext = imvcd_get_nalu_mvc_ext(
                    ps_dpb_mgr->s_dpb_ivp_ctxt.ps_nalu_mvc_exts, u2_view_order_id, u2_ref_view_id);

                if(!ps_ref_nalu_mvc_ext->u1_inter_view_flag)
                {
                    continue;
                }

                ps_au_buf = ih264_buf_mgr_get_next_free(
                    ps_dpb_mgr->ps_mvc_au_buf_mgr->ps_buf_mgr_ctxt, &i4_pic_buf_id);

                if(NULL == ps_au_buf)
                {
                    return ERROR_UNAVAIL_PICBUF_T;
                }
                else
                {
                    ih264_buf_mgr_set_status(ps_dpb_mgr->ps_mvc_au_buf_mgr->ps_buf_mgr_ctxt,
                                             i4_pic_buf_id, BUF_MGR_REF);
                }

                ps_au_mv_data = ih264_buf_mgr_get_next_free(
                    ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr->ps_buf_mgr_ctxt, &i4_mv_buf_id);

                if(NULL == ps_au_mv_data)
                {
                    return ERROR_UNAVAIL_PICBUF_T;
                }
                else
                {
                    ih264_buf_mgr_set_status(ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr->ps_buf_mgr_ctxt,
                                             i4_mv_buf_id, BUF_MGR_REF);
                }

                ps_au_buf->i4_pic_buf_id = i4_pic_buf_id;
                ps_au_buf->i4_mv_buf_id = i4_mv_buf_id;
                ps_au_buf->s_ivp_data.b_is_ivp_ref = true;
                ps_au_buf->s_ivp_data.u2_ref_view_id = u2_ref_view_id;

                imvcd_ivp_buf_copier(ps_cur_au, ps_au_buf, ps_cur_au->ps_au_mv_data, ps_au_mv_data,
                                     u2_ref_view_id, u2_view_id);

                ps_dpb_mgr->ps_mvc_au_buf_mgr->au1_au_buf_id_to_mv_buf_id_map[i4_pic_buf_id] =
                    i4_mv_buf_id;
                ps_dpb_mgr->ps_mvc_au_buf_mgr->aps_buf_id_to_au_buf_map[i4_pic_buf_id] = ps_au_buf;
                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr->aps_buf_id_to_mv_pred_buf_map[i4_mv_buf_id] =
                    ps_au_mv_data;

                ps_dpb_mgr->s_dpb_ivp_ctxt
                    .au1_au_buf_ids[ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs] = i4_pic_buf_id;
                ps_dpb_mgr->s_dpb_ivp_ctxt
                    .au1_mv_buf_ids[ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs] = i4_mv_buf_id;

                aps_ref_pic_buf_lx[i][0] = ps_au_buf[0];

                aps_ref_pic_buf_lx[i]++;
                pu1_num_short_term_refs[i]++;
                ps_dpb_mgr->s_dpb_ivp_ctxt.u4_num_ivp_refs++;
            }
        }
    }

    for(i = 0; i < 1 + ((WORD32) b_is_b_pic); i++)
    {
        au1_total_num_refs[i] = pu1_num_short_term_refs[i] + pu1_num_long_term_refs[i];

        if(pu1_num_short_term_refs[i] > MVC_MAX_REF_PICS)
        {
            return ERROR_NUM_REF;
        }

        if(au1_total_num_refs[i] > MVC_MAX_REF_PICS)
        {
            return ERROR_NUM_REF;
        }

        if(0 == au1_total_num_refs[i])
        {
            return ERROR_NUM_REF;
        }
    }

    /* If list0 and list1 entries are same then swap the 0th and 1st entry */
    /* of list 1 */
    {
        mvc_au_buffer_t *aps_ref_pic_bufs_lx[2] = {ps_dpb_mgr->as_init_dpb[0],
                                                   ps_dpb_mgr->as_init_dpb[1]};

        if((au1_total_num_refs[0] == au1_total_num_refs[1]) && (au1_total_num_refs[0] > 1))
        {
            bool b_swap;

            b_swap = true;

            for(i = 0; i < au1_total_num_refs[0]; i++)
            {
                if(aps_ref_pic_bufs_lx[0][i]
                       .as_view_buffers[u2_view_id]
                       .as_component_bufs[Y]
                       .pv_data != aps_ref_pic_bufs_lx[1][i]
                                       .as_view_buffers[u2_view_id]
                                       .as_component_bufs[Y]
                                       .pv_data)
                {
                    b_swap = false;

                    break;
                }
            }

            if(b_swap)
            {
                SWAP(aps_ref_pic_bufs_lx[1][0], aps_ref_pic_bufs_lx[1][1], mvc_au_buffer_t);
            }
        }
    }

    return OK;
}

void imvcd_dpb_set_missing_refs_to_default(mvc_dpb_manager_t *ps_dpb_mgr,
                                           ref_pic_list_mod_data_t *ps_ref_pic_list_mod_data,
                                           mvc_au_buffer_t *ps_cur_au, UWORD8 u1_pred_lx)
{
    UWORD8 u1_num_refs = ps_dpb_mgr->au1_num_active_st_refs[u1_pred_lx] +
                         ps_dpb_mgr->au1_num_active_lt_refs[u1_pred_lx];

    while(u1_num_refs < ps_ref_pic_list_mod_data->au1_num_active_refs[u1_pred_lx])
    {
        ps_dpb_mgr->as_init_dpb[u1_pred_lx][u1_num_refs] = ps_cur_au[0];
        ps_dpb_mgr->au1_num_active_st_refs[u1_pred_lx]++;
        u1_num_refs++;
    }
}

void imvcd_dpb_normalise_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr, UWORD16 u2_buf_mod_bitfield,
                                      UWORD8 u1_num_bufs_modified, UWORD8 u1_pred_lx)
{
    WORD32 i;

    UWORD8 u1_num_ref_bufs = ps_dpb_mgr->au1_num_active_st_refs[u1_pred_lx] +
                             ps_dpb_mgr->au1_num_active_lt_refs[u1_pred_lx];

    for(i = 0; i < u1_num_ref_bufs; i++)
    {
        if(!(u2_buf_mod_bitfield & (1 << i)))
        {
            ps_dpb_mgr->aps_mod_dpb[u1_pred_lx][u1_num_bufs_modified++] =
                &ps_dpb_mgr->as_init_dpb[u1_pred_lx][i];
        }
    }
}

WORD32 imvcd_dpb_reorder_ref_pic_list(mvc_dpb_manager_t *ps_dpb_mgr,
                                      nalu_mvc_ext_t *ps_cur_nalu_mvc_ext,
                                      mvc_au_buffer_t *ps_cur_au,
                                      ref_pic_list_mod_data_t *ps_ref_pic_list_mod_data,
                                      UWORD16 u2_view_order_id)
{
    WORD32 i, j;

    sps_mvc_ext_t *ps_sps_mvc_ext = ps_dpb_mgr->s_dpb_ivp_ctxt.ps_sps_mvc_ext;

    UWORD8 u1_anchor_pic_flag = ps_cur_nalu_mvc_ext->u1_anchor_pic_flag;
    WORD32 i4_cur_pic_num = ps_cur_au->i4_pic_num;
    WORD32 i4_max_pic_num = ps_dpb_mgr->i4_max_pic_num;
    bool b_is_b_pic = !!(ps_cur_au->au4_pack_slc_typ[u2_view_order_id] & B_SLC_BIT);

    for(i = 0; i < 1 + ((WORD32) b_is_b_pic); i++)
    {
        mvc_ivp_ref_data_t *ps_mvc_ivp_ref_data;

        UWORD16 u2_max_view_idx;

        WORD32 i4_pred_pic_num = i4_cur_pic_num;
        WORD16 i2_pred_view_order_id = 0;  // -1; Need to check spec and JMVM implementation match
        UWORD8 *pu1_modification_of_pic_nums_idc =
            ps_ref_pic_list_mod_data->au1_modification_of_pic_nums_idc[i];
        WORD32 *pi4_abs_diff_pic_num_minus1 =
            ps_ref_pic_list_mod_data->ai4_abs_diff_pic_num_minus1[i];
        WORD32 *pi4_long_term_pic_num = ps_ref_pic_list_mod_data->ai4_long_term_pic_num[i];
        WORD32 *pi4_abs_diff_view_idx_minus1 =
            ps_ref_pic_list_mod_data->ai4_abs_diff_view_idx_minus1[i];
        UWORD8 u1_num_ref_bufs =
            ps_dpb_mgr->au1_num_active_st_refs[i] + ps_dpb_mgr->au1_num_active_lt_refs[i];
        UWORD16 u2_buf_mod_bitfield = 0;
        UWORD8 u1_num_bufs_modified = 0;

        if(!ps_ref_pic_list_mod_data->au1_ref_pic_list_modification_flag_lx[i] ||
           (3 == pu1_modification_of_pic_nums_idc[0]))
        {
            imvcd_dpb_set_missing_refs_to_default(ps_dpb_mgr, ps_ref_pic_list_mod_data, ps_cur_au,
                                                  i);

            imvcd_dpb_normalise_ref_pic_list(ps_dpb_mgr, u2_buf_mod_bitfield, u1_num_bufs_modified,
                                             i);

            continue;
        }

        ps_mvc_ivp_ref_data =
            (0 == u2_view_order_id)
                ? NULL
                : (u1_anchor_pic_flag
                       ? &ps_sps_mvc_ext->as_anchor_ref_data[i][u2_view_order_id]
                       : &ps_sps_mvc_ext->as_non_anchor_ref_data[i][u2_view_order_id]);
        u2_max_view_idx = (0 == u2_view_order_id) ? 0 : ps_mvc_ivp_ref_data->u1_num_refs;

        do
        {
            if((0 == pu1_modification_of_pic_nums_idc[0]) ||
               (1 == pu1_modification_of_pic_nums_idc[0]))
            {
                WORD32 i4_mod_pic_num = 1 + pi4_abs_diff_pic_num_minus1[0];
                UWORD8 u1_mod_buf_idx = u1_num_ref_bufs;

                if(pi4_abs_diff_pic_num_minus1[0] > i4_max_pic_num)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                if(0 == pu1_modification_of_pic_nums_idc[0])
                {
                    i4_mod_pic_num = i4_pred_pic_num - i4_mod_pic_num;

                    if(i4_mod_pic_num < 0)
                    {
                        i4_mod_pic_num += i4_max_pic_num;
                    }
                }
                else
                {
                    i4_mod_pic_num = i4_pred_pic_num + i4_mod_pic_num;

                    if(i4_mod_pic_num >= i4_max_pic_num)
                    {
                        i4_mod_pic_num -= i4_max_pic_num;
                    }
                }

                if(i4_mod_pic_num > i4_cur_pic_num)
                {
                    i4_mod_pic_num -= i4_max_pic_num;
                }

                for(j = 0; j < u1_num_ref_bufs; j++)
                {
                    if(ps_dpb_mgr->as_init_dpb[i][j].i4_pic_num == i4_mod_pic_num)
                    {
                        u1_mod_buf_idx = j;

                        break;
                    }
                }

                if(u1_mod_buf_idx == u1_num_ref_bufs)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                ps_dpb_mgr->aps_mod_dpb[i][u1_num_bufs_modified++] =
                    &ps_dpb_mgr->as_init_dpb[i][u1_mod_buf_idx];

                u2_buf_mod_bitfield |= (1 << u1_mod_buf_idx);
                i4_pred_pic_num = i4_mod_pic_num;
            }
            else if(2 == pu1_modification_of_pic_nums_idc[0])
            {
                WORD32 i4_mod_lt_pic_num = pi4_long_term_pic_num[0];
                UWORD8 u1_mod_buf_idx = u1_num_ref_bufs;

                if(pi4_long_term_pic_num[0] > (MAX_REF_BUFS + 1))
                {
                    return ERROR_DBP_MANAGER_T;
                }

                for(j = 0; j < u1_num_ref_bufs; j++)
                {
                    if(!ps_dpb_mgr->as_init_dpb[i][j].b_is_short_term_ref &&
                       (i4_mod_lt_pic_num == ps_dpb_mgr->as_init_dpb[i][j].u1_long_term_pic_num))
                    {
                        u1_mod_buf_idx = j;

                        break;
                    }
                }

                if(u1_mod_buf_idx == u1_num_ref_bufs)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                ps_dpb_mgr->aps_mod_dpb[i][u1_num_bufs_modified++] =
                    &ps_dpb_mgr->as_init_dpb[i][u1_mod_buf_idx];

                u2_buf_mod_bitfield |= (1 << u1_mod_buf_idx);
            }
            else if((4 == pu1_modification_of_pic_nums_idc[0]) ||
                    (5 == pu1_modification_of_pic_nums_idc[0]))
            {
                WORD32 i4_target_view_id;

                WORD32 i4_mod_view_order_id = pi4_abs_diff_view_idx_minus1[0] + 1;
                UWORD8 u1_mod_buf_idx = u1_num_ref_bufs;

                if(4 == pu1_modification_of_pic_nums_idc[0])
                {
                    i4_mod_view_order_id = i2_pred_view_order_id - i4_mod_view_order_id;

                    if(i4_mod_view_order_id < 0)
                    {
                        i4_mod_view_order_id += u2_max_view_idx;
                    }
                }
                else
                {
                    i4_mod_view_order_id = i2_pred_view_order_id + i4_mod_view_order_id;

                    if(i4_mod_view_order_id >= u2_max_view_idx)
                    {
                        i4_mod_view_order_id -= u2_max_view_idx;
                    }
                }

                if((0 == u2_view_order_id) ||
                   !((i4_mod_view_order_id >= 0) && (i4_mod_view_order_id <= u2_max_view_idx)) ||
                   (NULL == ps_mvc_ivp_ref_data))
                {
                    return ERROR_DBP_MANAGER_T;
                }

                i4_target_view_id = ps_mvc_ivp_ref_data->au2_ref_view_ids[i4_mod_view_order_id];

                for(j = 0; j < u1_num_ref_bufs; j++)
                {
                    if(ps_dpb_mgr->as_init_dpb[i][j].s_ivp_data.b_is_ivp_ref)
                    {
                        if((ps_dpb_mgr->as_init_dpb[i][j].i4_pic_num == i4_cur_pic_num) &&
                           (ps_dpb_mgr->as_init_dpb[i][j].s_ivp_data.u2_ref_view_id ==
                            i4_target_view_id))
                        {
                            u1_mod_buf_idx = j;

                            break;
                        }
                    }
                }

                if(u1_mod_buf_idx == u1_num_ref_bufs)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                u2_buf_mod_bitfield |= (1 << u1_mod_buf_idx);
                ps_dpb_mgr->aps_mod_dpb[i][u1_num_bufs_modified++] =
                    &ps_dpb_mgr->as_init_dpb[i][u1_mod_buf_idx];
                i2_pred_view_order_id = i4_mod_view_order_id;
            }
            else if(3 != pu1_modification_of_pic_nums_idc[0])
            {
                return ERROR_REFIDX_ORDER_T;
            }
            else
            {
                break;
            }

            pu1_modification_of_pic_nums_idc++;
            pi4_abs_diff_pic_num_minus1++;
            pi4_long_term_pic_num++;
            pi4_abs_diff_view_idx_minus1++;
        } while(true);

        imvcd_dpb_set_missing_refs_to_default(ps_dpb_mgr, ps_ref_pic_list_mod_data, ps_cur_au, i);

        imvcd_dpb_normalise_ref_pic_list(ps_dpb_mgr, u2_buf_mod_bitfield, u1_num_bufs_modified, i);
    }

    return OK;
}

WORD32 imvcd_dpb_insert_st_node(mvc_dpb_manager_t *ps_dpb_mgr, mvc_au_buffer_t *ps_au_buf)
{
    WORD32 i;

    mvc_dpb_info_t *ps_dpb_info = ps_dpb_mgr->as_dpb_info;
    UWORD8 u1_picture_type = ps_au_buf->u1_picturetype;

    for(i = 0; i < MVC_MAX_REF_PICS; i++)
    {
        if((ps_dpb_info[i].ps_au_buf == ps_au_buf) && ps_dpb_info[i].b_used_as_ref)
        {
            /* Can occur only for field bottom pictures */
            if(ps_dpb_info[i].ps_au_buf->u1_pic_type == FRM_PIC)
            {
                return ERROR_DBP_MANAGER_T;
            }
            else
            {
                ps_dpb_info[i].s_bot_field.u1_reference_info = IS_SHORT_TERM;

                return OK;
            }
        }

        if(!ps_dpb_info[i].b_used_as_ref &&
           (ps_dpb_info[i].s_top_field.u1_reference_info == UNUSED_FOR_REF) &&
           (ps_dpb_info[i].s_bot_field.u1_reference_info == UNUSED_FOR_REF))
        {
            break;
        }
    }

    if(i == MVC_MAX_REF_PICS)
    {
        return ERROR_DBP_MANAGER_T;
    }

    ps_dpb_info[i].ps_au_buf = ps_au_buf;
    ps_dpb_info[i].ps_prev_short = ps_dpb_mgr->ps_dpb_st_head;
    ps_dpb_info[i].b_used_as_ref = true;

    ps_dpb_mgr->ps_dpb_st_head = ps_dpb_info + i;

    ps_dpb_mgr->u1_num_st_ref_bufs++;

    ps_au_buf->b_is_short_term_ref = true;

    if((u1_picture_type & 0x03) == FRM_PIC)
    {
        ps_dpb_info[i].s_top_field.u1_reference_info = IS_SHORT_TERM;
        ps_dpb_info[i].s_bot_field.u1_reference_info = IS_SHORT_TERM;
    }
    else if((u1_picture_type & 0x03) == TOP_FLD)
    {
        ps_dpb_info[i].s_top_field.u1_reference_info = IS_SHORT_TERM;
    }
    else if((u1_picture_type & 0x03) == BOT_FLD)
    {
        ps_dpb_info[i].s_bot_field.u1_reference_info = IS_SHORT_TERM;
    }

    return OK;
}

static WORD32 imvcd_dpb_delete_gap_frm_mmco(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_frame_num,
                                            UWORD8 *pu1_del_node)
{
    WORD8 i, j;
    WORD32 *pi4_start, *pi4_end;
    WORD32 i4_start_frm_num, i4_end_frm_num, i4_max_pic_num;

    /* find the least frame num from gaps and current DPB node    */
    /* Delete the gaps                                            */
    *pu1_del_node = 1;
    pi4_start = ps_dpb_mgr->ai4_gaps_start_frm_num;
    pi4_end = ps_dpb_mgr->ai4_gaps_end_frm_num;
    i4_max_pic_num = ps_dpb_mgr->i4_max_pic_num;

    if(0 == ps_dpb_mgr->u1_num_gaps)
    {
        return OK;
    }

    if(i4_frame_num < 0)
    {
        i4_frame_num += i4_max_pic_num;
    }

    for(i = 0; i < MAX_FRAMES; i++)
    {
        i4_start_frm_num = pi4_start[i];

        if(i4_start_frm_num < 0)
        {
            i4_start_frm_num += i4_max_pic_num;
        }

        if(INVALID_FRAME_NUM != i4_start_frm_num)
        {
            i4_end_frm_num = pi4_end[i];

            if(i4_end_frm_num < 0)
            {
                i4_end_frm_num += i4_max_pic_num;
            }

            if((i4_frame_num >= i4_start_frm_num) && (i4_frame_num <= i4_end_frm_num))
            {
                break;
            }
            else
            {
                if(((i4_frame_num + i4_max_pic_num) >= i4_start_frm_num) &&
                   ((i4_frame_num + i4_max_pic_num) <= i4_end_frm_num))
                {
                    return ERROR_DBP_MANAGER_T;
                }
            }
        }
    }

    /* find frame_num index, in the poc_map which needs to be deleted */
    for(j = 0; j < MAX_FRAMES; j++)
    {
        if(i4_frame_num == ps_dpb_mgr->as_display_buf_info[j].i4_frame_num)
        {
            break;
        }
    }

    if(MAX_FRAMES != i)
    {
        if(j == MAX_FRAMES)
        {
            return ERROR_DBP_MANAGER_T;
        }

        ps_dpb_mgr->as_display_buf_info[j].i4_poc_buf_id = -1;
        ps_dpb_mgr->as_display_buf_info[j].i4_poc = 0x7fffffff;
        ps_dpb_mgr->as_display_buf_info[j].i4_frame_num = GAP_FRAME_NUM;

        ps_dpb_mgr->i1_gaps_deleted++;
        ps_dpb_mgr->ai1_gaps_per_seq[i]--;
        ps_dpb_mgr->u1_num_gaps--;
        *pu1_del_node = 0;

        if(0 == ps_dpb_mgr->ai1_gaps_per_seq[i])
        {
            ps_dpb_mgr->ai4_gaps_start_frm_num[i] = INVALID_FRAME_NUM;
            ps_dpb_mgr->ai4_gaps_end_frm_num[i] = 0;
        }
    }
    else
    {
        return ERROR_DBP_MANAGER_T;
    }

    return OK;
}

static WORD32 imvcd_dpb_insert_lt_node(mvc_dpb_manager_t *ps_dpb_mgr, mvc_dpb_info_t *ps_new_node,
                                       UWORD32 u4_lt_idx)
{
    ps_new_node->s_top_field.u1_reference_info = IS_LONG_TERM;
    ps_new_node->s_bot_field.u1_reference_info = IS_LONG_TERM;
    ps_new_node->s_top_field.u1_long_term_frame_idx = u4_lt_idx;
    ps_new_node->s_bot_field.u1_long_term_frame_idx = u4_lt_idx;
    ps_new_node->ps_au_buf->u1_long_term_frm_idx = u4_lt_idx;
    ps_new_node->b_used_as_ref = true;

    if(ps_dpb_mgr->u1_num_lt_ref_bufs > 0)
    {
        WORD32 i;

        mvc_dpb_info_t **pps_next_node = &ps_dpb_mgr->ps_dpb_lt_head;

        for(i = 0; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
        {
            if((*pps_next_node)->ps_au_buf->u1_long_term_frm_idx > u4_lt_idx)
            {
                ps_new_node->ps_prev_long = *pps_next_node;
                *pps_next_node = ps_new_node;

                break;
            }
            else if(NULL == (*pps_next_node)->ps_prev_long)
            {
                (*pps_next_node)->ps_prev_long = ps_new_node;
                ps_new_node->ps_prev_long = NULL;

                break;
            }

            pps_next_node = &(*pps_next_node)->ps_prev_long;
        }
    }
    else
    {
        ps_dpb_mgr->ps_dpb_lt_head = ps_new_node;
        ps_new_node->ps_prev_long = NULL;
    }

    ps_new_node->ps_au_buf->b_is_short_term_ref = false;

    ps_dpb_mgr->u1_num_lt_ref_bufs++;

    return OK;
}

static WORD32 imvcd_dpb_delete_lt_node(mvc_dpb_manager_t *ps_dpb_mgr, UWORD32 u4_lt_idx)
{
    mvc_dpb_info_t *ps_next_dpb;
    mvc_dpb_info_t *ps_unmark_node;

    WORD32 i;

    if(ps_dpb_mgr->u1_num_lt_ref_bufs > 0)
    {
        ps_next_dpb = ps_dpb_mgr->ps_dpb_lt_head;

        if(ps_next_dpb->ps_au_buf->u1_long_term_frm_idx == u4_lt_idx)
        {
            ps_unmark_node = ps_next_dpb;
        }
        else
        {
            for(i = 1; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
            {
                if(ps_next_dpb->ps_prev_long->ps_au_buf->u1_long_term_frm_idx == u4_lt_idx)
                {
                    break;
                }

                ps_next_dpb = ps_next_dpb->ps_prev_long;
            }

            if(i < ps_dpb_mgr->u1_num_lt_ref_bufs)
            {
                ps_unmark_node = ps_next_dpb->ps_prev_long;
            }
            else
            {
                return OK;
            }
        }

        ps_unmark_node->b_used_as_ref = false;

        if(ps_unmark_node == ps_dpb_mgr->ps_dpb_lt_head)
        {
            ps_dpb_mgr->ps_dpb_lt_head = ps_next_dpb->ps_prev_long;
        }

        ps_unmark_node->s_top_field.u1_reference_info = UNUSED_FOR_REF;
        ps_unmark_node->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

        imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr, ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                            ps_unmark_node->ps_au_buf->i4_pic_buf_id);

        ps_next_dpb->ps_prev_long = ps_unmark_node->ps_prev_long;
        ps_unmark_node->ps_prev_long = NULL;
        ps_dpb_mgr->u1_num_lt_ref_bufs--;
    }

    return OK;
}

WORD32 imvcd_dpb_delete_st_node_or_make_lt(mvc_dpb_manager_t *ps_dpb_mgr, WORD32 i4_pic_num,
                                           UWORD32 u4_lt_idx)
{
    WORD32 i4_error_code;

    mvc_dpb_info_t *ps_next_dpb = ps_dpb_mgr->ps_dpb_st_head;
    mvc_dpb_info_t *ps_unmark_node = NULL;

    UWORD8 u1_del_node = 0, u1_del_st = 0;
    WORD32 i = 0;

    if(ps_next_dpb->ps_au_buf->i4_pic_num == i4_pic_num)
    {
        ps_unmark_node = ps_next_dpb;
    }
    else
    {
        for(i = 1; i < ps_dpb_mgr->u1_num_st_ref_bufs; i++)
        {
            if(ps_next_dpb->ps_prev_short->ps_au_buf->i4_pic_num == i4_pic_num)
            {
                ps_unmark_node = ps_next_dpb->ps_prev_short;

                break;
            }

            ps_next_dpb = ps_next_dpb->ps_prev_short;
        }
    }

    if(i == ps_dpb_mgr->u1_num_st_ref_bufs)
    {
        if(ps_dpb_mgr->u1_num_gaps)
        {
            i4_error_code = imvcd_dpb_delete_gap_frm_mmco(ps_dpb_mgr, i4_pic_num, &u1_del_st);

            if(i4_error_code != OK)
            {
                return i4_error_code;
            }
        }
        else
        {
            return ERROR_DBP_MANAGER_T;
        }

        if(u1_del_st)
        {
            return ERROR_DBP_MANAGER_T;
        }
        else
        {
            return 0;
        }
    }

    ps_unmark_node->b_used_as_ref = false;
    ps_unmark_node->s_top_field.u1_reference_info = UNUSED_FOR_REF;
    ps_unmark_node->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

    if(ps_unmark_node == ps_dpb_mgr->ps_dpb_st_head)
    {
        ps_dpb_mgr->ps_dpb_st_head = ps_next_dpb->ps_prev_short;
    }
    else
    {
        ps_next_dpb->ps_prev_short = ps_unmark_node->ps_prev_short;
    }

    ps_dpb_mgr->u1_num_st_ref_bufs--;
    u1_del_node = 1;

    if(u4_lt_idx == (MAX_REF_BUFS + 1))
    {
        if(u1_del_node)
        {
            imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                ps_unmark_node->ps_au_buf->i4_pic_buf_id);

            ps_unmark_node->ps_prev_short = NULL;
        }
    }
    else
    {
        i4_error_code = imvcd_dpb_delete_lt_node(ps_dpb_mgr, u4_lt_idx);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        i4_error_code = imvcd_dpb_insert_lt_node(ps_dpb_mgr, ps_unmark_node, u4_lt_idx);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }
    }

    return OK;
}

WORD32 imvcd_dpb_do_mmco(dpb_commands_t *ps_dpb_cmds, mvc_dpb_manager_t *ps_dpb_mgr,
                         mvc_au_buffer_t *ps_cur_au, UWORD8 u1_max_num_ref_frames,
                         UWORD8 u1_curr_pic_in_err)
{
    mvc_dpb_info_t *ps_next_dpb;

    WORD32 i, j;
    UWORD8 u1_buf_mode, u1_marked_lt;
    UWORD8 u1_num_gaps;
    WORD32 i4_error_code;

    UWORD8 u1_del_node = 1;
    UWORD8 u1_insert_st_pic = 1;

    // 0 - sliding window; 1 - Adaptive
    u1_buf_mode = ps_dpb_cmds->u1_buf_mode;
    u1_marked_lt = 0;
    u1_num_gaps = ps_dpb_mgr->u1_num_gaps;

    if(!u1_buf_mode)
    {
        // Sliding window - implements 8.2.5.3
        if((ps_dpb_mgr->u1_num_st_ref_bufs + ps_dpb_mgr->u1_num_lt_ref_bufs + u1_num_gaps) ==
           u1_max_num_ref_frames)
        {
            UWORD8 u1_new_node_flag = 1;

            if((0 == ps_dpb_mgr->u1_num_st_ref_bufs) && (0 == u1_num_gaps))
            {
                return ERROR_DBP_MANAGER_T;
            }

            // Chase the links to reach the last but one picNum, if available
            ps_next_dpb = ps_dpb_mgr->ps_dpb_st_head;

            if(ps_dpb_mgr->u1_num_st_ref_bufs > 1)
            {
                if(ps_next_dpb->ps_au_buf->i4_pic_num == ps_cur_au->i4_pic_num)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                for(i = 1; i < (ps_dpb_mgr->u1_num_st_ref_bufs - 1); i++)
                {
                    if(ps_next_dpb == NULL)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    if(ps_next_dpb->ps_au_buf->i4_pic_num == ps_cur_au->i4_pic_num)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    ps_next_dpb = ps_next_dpb->ps_prev_short;
                }

                if(ps_next_dpb->ps_prev_short->ps_prev_short != NULL)
                {
                    return ERROR_DBP_MANAGER_T;
                }

                if(u1_new_node_flag)
                {
                    if(u1_num_gaps)
                    {
                        i4_error_code = imvcd_dpb_delete_gap_frm_sliding(
                            ps_dpb_mgr, ps_next_dpb->ps_prev_short->ps_au_buf->i4_pic_num,
                            &u1_del_node);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }

                    if(u1_del_node)
                    {
                        ps_dpb_mgr->u1_num_st_ref_bufs--;
                        ps_next_dpb->ps_prev_short->b_used_as_ref = false;
                        ps_next_dpb->ps_prev_short->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                        ps_next_dpb->ps_prev_short->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                        imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                            ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                            ps_next_dpb->ps_prev_short->ps_au_buf->i4_pic_buf_id);

                        ps_next_dpb->ps_prev_short->ps_au_buf = NULL;
                        ps_next_dpb->ps_prev_short = NULL;
                    }
                }
            }
            else
            {
                if(ps_dpb_mgr->u1_num_st_ref_bufs)
                {
                    i4_error_code = imvcd_dpb_delete_gap_frm_sliding(
                        ps_dpb_mgr, ps_next_dpb->ps_au_buf->i4_pic_num, &u1_del_node);

                    if(i4_error_code != OK)
                    {
                        return i4_error_code;
                    }

                    if((ps_next_dpb->ps_au_buf->i4_pic_num != ps_cur_au->i4_pic_num) && u1_del_node)
                    {
                        ps_dpb_mgr->u1_num_st_ref_bufs--;
                        ps_next_dpb->b_used_as_ref = false;
                        ps_next_dpb->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                        ps_next_dpb->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                        imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                            ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                            ps_next_dpb->ps_au_buf->i4_pic_buf_id);

                        ps_next_dpb->ps_au_buf = NULL;
                        ps_next_dpb->ps_prev_short = NULL;
                        ps_dpb_mgr->ps_dpb_st_head = NULL;
                        ps_next_dpb = NULL;
                    }
                    else if(ps_next_dpb->ps_au_buf->i4_pic_num == ps_cur_au->i4_pic_num)
                    {
                        if(u1_curr_pic_in_err)
                        {
                            u1_insert_st_pic = 0;
                        }
                        else if(ps_dpb_mgr->u1_num_st_ref_bufs > 0)
                        {
                            ps_dpb_mgr->u1_num_st_ref_bufs--;
                            ps_next_dpb->b_used_as_ref = false;
                            ps_next_dpb->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                            ps_next_dpb->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                            imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                                ps_next_dpb->ps_au_buf->i4_pic_buf_id);

                            ps_next_dpb->ps_au_buf = NULL;
                            ps_next_dpb = NULL;
                        }
                    }
                }
                else
                {
                    i4_error_code = imvcd_dpb_delete_gap_frm_sliding(ps_dpb_mgr, INVALID_FRAME_NUM,
                                                                     &u1_del_node);

                    if(i4_error_code != OK)
                    {
                        return i4_error_code;
                    }

                    if(u1_del_node)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }
                }
            }
        }
    }
    else
    {
        // Adaptive memory control - implements 8.2.5.4
        struct MMCParams *ps_mmc_params;

        UWORD32 u4_mmco;
        UWORD32 u4_diff_pic_num;
        UWORD32 u4_lt_idx;

        UWORD32 au4_num_mmco_cmds[NUM_MMCO_CMD_IDS] = {0};

        for(j = 0; j < ps_dpb_cmds->u1_num_of_commands; j++)
        {
            ps_mmc_params = &ps_dpb_cmds->as_mmc_params[j];
            u4_mmco = ps_mmc_params->u4_mmco;

            switch(u4_mmco)
            {
                case MARK_ST_PICNUM_AS_NONREF:
                {
                    WORD64 i8_pic_num;

                    u4_diff_pic_num = ps_mmc_params->u4_diff_pic_num;
                    i8_pic_num =
                        ((WORD64) ps_cur_au->i4_pic_num) - ((WORD64) (u4_diff_pic_num + 1));

                    if(IS_OUT_OF_RANGE_S32(i8_pic_num))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    if(ps_dpb_mgr->u1_num_st_ref_bufs > 0)
                    {
                        i4_error_code = imvcd_dpb_delete_st_node_or_make_lt(
                            ps_dpb_mgr, (WORD32) i8_pic_num, MAX_REF_BUFS + 1);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }
                    else
                    {
                        UWORD8 u1_dummy;

                        i4_error_code = imvcd_dpb_delete_gap_frm_mmco(
                            ps_dpb_mgr, (WORD32) i8_pic_num, &u1_dummy);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }

                    break;
                }
                case MARK_LT_INDEX_AS_NONREF:
                {
                    u4_lt_idx = ps_mmc_params->u4_lt_idx;

                    i4_error_code = imvcd_dpb_delete_lt_node(ps_dpb_mgr, u4_lt_idx);

                    if(i4_error_code != OK)
                    {
                        return i4_error_code;
                    }

                    break;
                }
                case MARK_ST_PICNUM_AS_LT_INDEX:
                {
                    WORD64 i8_pic_num;

                    u4_diff_pic_num = ps_mmc_params->u4_diff_pic_num;

                    i8_pic_num =
                        ((WORD64) ps_cur_au->i4_pic_num) - ((WORD64) (u4_diff_pic_num + 1));

                    if(IS_OUT_OF_RANGE_S32(i8_pic_num))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    u4_lt_idx = ps_mmc_params->u4_lt_idx;

                    if((ps_dpb_mgr->u1_max_lt_frame_idx == NO_LONG_TERM_INDICIES) ||
                       (u4_lt_idx > ps_dpb_mgr->u1_max_lt_frame_idx))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    if(ps_dpb_mgr->u1_num_st_ref_bufs > 0)
                    {
                        i4_error_code = imvcd_dpb_delete_st_node_or_make_lt(
                            ps_dpb_mgr, (WORD32) i8_pic_num, u4_lt_idx);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }

                    break;
                }
                case SET_MAX_LT_INDEX:
                {
                    if(au4_num_mmco_cmds[SET_MAX_LT_INDEX] > 0)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    u4_lt_idx =
                        ps_mmc_params->u4_max_lt_idx_plus1;  // Get Max_long_term_index_plus1

                    if((u4_lt_idx <= ps_dpb_mgr->u1_max_lt_frame_idx) &&
                       (ps_dpb_mgr->u1_num_lt_ref_bufs > 0))
                    {
                        mvc_dpb_info_t *ps_nxtDPB;

                        // Set all LT buffers with index >= u4_lt_idx to nonreference
                        ps_nxtDPB = ps_dpb_mgr->ps_dpb_lt_head;
                        ps_next_dpb = ps_nxtDPB->ps_prev_long;

                        if(ps_nxtDPB->ps_au_buf->u1_long_term_frm_idx >= u4_lt_idx)
                        {
                            i = 0;
                            ps_dpb_mgr->ps_dpb_lt_head = NULL;
                        }
                        else
                        {
                            for(i = 1; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
                            {
                                if(ps_next_dpb->ps_au_buf->u1_long_term_frm_idx >= u4_lt_idx)
                                {
                                    break;
                                }

                                ps_nxtDPB = ps_next_dpb;
                                ps_next_dpb = ps_next_dpb->ps_prev_long;
                            }

                            ps_nxtDPB->ps_prev_long = NULL;  // Terminate the link of the
                                                             // closest LTIndex that is <=Max
                        }

                        ps_dpb_mgr->u1_num_lt_ref_bufs = i;

                        if(i == 0)
                        {
                            ps_next_dpb = ps_nxtDPB;
                        }

                        for(; i < ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
                        {
                            ps_nxtDPB = ps_next_dpb;
                            ps_nxtDPB->b_used_as_ref = false;
                            ps_nxtDPB->s_top_field.u1_reference_info = UNUSED_FOR_REF;
                            ps_nxtDPB->s_bot_field.u1_reference_info = UNUSED_FOR_REF;

                            imvcd_free_ref_bufs(ps_dpb_mgr->ps_mvc_au_buf_mgr,
                                                ps_dpb_mgr->ps_mvc_au_mv_pred_buf_mgr,
                                                ps_nxtDPB->ps_au_buf->i4_pic_buf_id);

                            ps_nxtDPB->ps_au_buf = NULL;

                            ps_next_dpb = ps_nxtDPB->ps_prev_long;
                            ps_nxtDPB->ps_prev_long = NULL;
                        }
                    }

                    if(u4_lt_idx == 0)
                    {
                        ps_dpb_mgr->u1_max_lt_frame_idx = NO_LONG_TERM_INDICIES;
                    }
                    else
                    {
                        ps_dpb_mgr->u1_max_lt_frame_idx = u4_lt_idx - 1;
                    }

                    break;
                }
                case SET_LT_INDEX:
                {
                    if(au4_num_mmco_cmds[SET_LT_INDEX] > 0)
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    u4_lt_idx = ps_mmc_params->u4_lt_idx;  // Get long term index

                    if((ps_dpb_mgr->u1_max_lt_frame_idx == NO_LONG_TERM_INDICIES) ||
                       (u4_lt_idx > ps_dpb_mgr->u1_max_lt_frame_idx))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    i4_error_code = imvcd_dpb_insert_st_node(ps_dpb_mgr, ps_cur_au);

                    if(i4_error_code != OK)
                    {
                        return i4_error_code;
                    }

                    if(ps_dpb_mgr->u1_num_st_ref_bufs > 0)
                    {
                        i4_error_code = imvcd_dpb_delete_st_node_or_make_lt(
                            ps_dpb_mgr, ps_cur_au->i4_pic_num, u4_lt_idx);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }
                    else
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    u1_marked_lt = 1;

                    break;
                }
                case RESET_REF_PICTURES:
                {
                    if((au4_num_mmco_cmds[RESET_REF_PICTURES] > 0) ||
                       (au4_num_mmco_cmds[MARK_ST_PICNUM_AS_NONREF] > 0) ||
                       (au4_num_mmco_cmds[MARK_LT_INDEX_AS_NONREF] > 0) ||
                       (au4_num_mmco_cmds[MARK_ST_PICNUM_AS_LT_INDEX] > 0))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    if((j > 0) && (ps_dpb_cmds->as_mmc_params[j - 1].u4_mmco == SET_LT_INDEX))
                    {
                        return ERROR_DBP_MANAGER_T;
                    }

                    __attribute__((fallthrough));
                }
                case RESET_ALL_PICTURES:
                {
                    WORD32 i4_pic_num = ps_cur_au->i4_frame_num;

                    imvcd_reset_dpb(ps_dpb_mgr);

                    ps_cur_au->i4_frame_num = 0;

                    if(!u1_marked_lt && u1_insert_st_pic)
                    {
                        i4_error_code = imvcd_dpb_insert_st_node(ps_dpb_mgr, ps_cur_au);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }

                    ps_cur_au->i4_frame_num = i4_pic_num;

                    return OK;
                }
                default:
                {
                    return ERROR_DBP_MANAGER_T;
                }
            }

            au4_num_mmco_cmds[u4_mmco]++;
        }
    }

    if(!u1_marked_lt && u1_insert_st_pic)
    {
        i4_error_code = imvcd_dpb_insert_st_node(ps_dpb_mgr, ps_cur_au);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }
    }

    return OK;
}

WORD32 imvcd_dpb_update_default_index_list(mvc_dpb_manager_t *ps_dpb_mgr)
{
    WORD32 i;

    mvc_dpb_info_t *ps_next_dpb = ps_dpb_mgr->ps_dpb_st_head;

    for(i = 0; i < ps_dpb_mgr->u1_num_st_ref_bufs; i++)
    {
        ps_dpb_mgr->aps_def_dpb[i] = ps_next_dpb->ps_au_buf;
        ps_next_dpb = ps_next_dpb->ps_prev_short;
    }

    ps_next_dpb = ps_dpb_mgr->ps_dpb_lt_head;

    for(; i < ps_dpb_mgr->u1_num_st_ref_bufs + ps_dpb_mgr->u1_num_lt_ref_bufs; i++)
    {
        ps_dpb_mgr->aps_def_dpb[i] = ps_next_dpb->ps_au_buf;
        ps_next_dpb = ps_next_dpb->ps_prev_long;
    }

    return OK;
}
