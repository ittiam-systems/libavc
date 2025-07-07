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
/*  File Name         : imvcd_slice_functions.c                              */
/*                                                                           */
/*  Description       : Functions for MVC Slice parsing, etc.                */
/*                                                                           */
/*****************************************************************************/

#include "ih264_typedefs.h"
#include "ih264_error.h"
#include "ih264_buf_mgr.h"
#include "ih264d_bitstrm.h"
#include "ih264d_deblocking.h"
#include "ih264d_debug.h"
#include "ih264d_defs.h"
#include "ih264d_error_handler.h"
#include "ih264d_inter_pred.h"
#include "ih264d_mb_utils.h"
#include "ih264d_mvpred.h"
#include "ih264d_parse_slice.h"
#include "ih264d_parse_islice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_pslice.h"
#include "ih264d_quant_scaling.h"
#include "ih264d_tables.h"
#include "ih264d_thread_compute_bs.h"
#include "ih264d_thread_parse_decode.h"
#include "ih264d_structs.h"
#include "ih264d_utils.h"
#include "ih264d_api_utils.h"
#include "ithread.h"
#include "imvc_defs.h"
#include "imvcd_dpb_manager.h"
#include "imvcd_error_handler.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

static WORD32 imvcd_set_first_mb_in_slice(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u2_first_mb_in_slice = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(ps_cur_slice->u2_first_mb_in_slice >=
       (ps_view_ctxt->u2_frm_ht_in_mbs * ps_view_ctxt->u2_frm_wd_in_mbs))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    if(((ps_cur_slice->u2_first_mb_in_slice << ps_cur_slice->u1_mbaff_frame_flag) <=
        ps_view_ctxt->u4_cur_mb_addr) &&
       (ps_view_ctxt->u4_first_slice_in_pic == 0))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    COPYTHECONTEXT("SH: first_mb_in_slice", ps_cur_slice->u2_first_mb_in_slice);

    return OK;
}

static WORD32 imvcd_set_slice_type(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u1_slice_type = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(ps_cur_slice->u1_slice_type > 9)
    {
        return ERROR_INV_SLC_TYPE_T;
    }

    if(ps_cur_slice->u1_slice_type > 4)
    {
        ps_cur_slice->u1_slice_type -= 5;
    }

    COPYTHECONTEXT("SH: slice_type", ps_cur_slice->u1_slice_type);

    return OK;
}

static WORD32 imvcd_set_cur_pps(dec_struct_t *ps_view_ctxt, UWORD8 *pu1_pps_id)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pu1_pps_id[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    ps_view_ctxt->ps_cur_pps = &ps_view_ctxt->ps_pps[pu1_pps_id[0]];
    ps_view_ctxt->ps_cur_sps = ps_view_ctxt->ps_pps[pu1_pps_id[0]].ps_sps;

    if(!ps_view_ctxt->ps_cur_pps->u1_is_valid || !ps_view_ctxt->ps_cur_pps->ps_sps->u1_is_valid)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    COPYTHECONTEXT("SH: pps_id", pu1_pps_id[0]);

    return OK;
}

static WORD32 imvcd_set_frame_num(dec_struct_t *ps_view_ctxt, UWORD8 u1_bits_in_frm_num)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u2_frame_num = ih264d_get_bits_h264(ps_bitstrm, u1_bits_in_frm_num);

    COPYTHECONTEXT("SH: frame_num", ps_cur_slice->u2_frame_num);

    return OK;
}

static WORD32 imvcd_set_idr_pic_id(dec_struct_t *ps_view_ctxt, UWORD32 *pu4_idr_pic_id)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pu4_idr_pic_id[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(pu4_idr_pic_id[0] > 65535)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    COPYTHECONTEXT("SH: idr_pic_id", pu4_idr_pic_id[0]);

    return OK;
}

static WORD32 imvcd_set_poc_lsb(dec_struct_t *ps_view_ctxt, WORD32 *pi4_pic_order_cnt_lsb,
                                WORD32 i4_max_poc_lsb, UWORD8 u1_log2_max_poc_lsb)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pi4_pic_order_cnt_lsb[0] = ih264d_get_bits_h264(ps_bitstrm, u1_log2_max_poc_lsb);

    if((pi4_pic_order_cnt_lsb[0] < 0) || (pi4_pic_order_cnt_lsb[0] > i4_max_poc_lsb))
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    COPYTHECONTEXT("SH: pic_order_cnt_lsb", pi4_pic_order_cnt_lsb[0]);

    return OK;
}

static WORD32 imvcd_set_delta_poc(dec_struct_t *ps_view_ctxt, WORD32 *pi4_delta_poc)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pi4_delta_poc[0] = ih264d_sev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    COPYTHECONTEXT("SH: delta_pic_order_cnt", pi4_delta_poc[0]);

    return OK;
}

static WORD32 imvcd_set_redundant_pic_cnt(dec_struct_t *ps_view_ctxt, UWORD8 *pu1_redundant_pic_cnt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pu1_redundant_pic_cnt[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(pu1_redundant_pic_cnt[0] > MAX_REDUNDANT_PIC_CNT)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    COPYTHECONTEXT("SH: redundant_pic_cnt", pu1_redundant_pic_cnt[0]);

    return OK;
}

static WORD32 imvcd_set_direct_spatial_mv_pred_flag(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u1_direct_spatial_mv_pred_flag = ih264d_get_bit_h264(ps_bitstrm);

    COPYTHECONTEXT("SH: direct_spatial_mv_pred_flag", ps_cur_slice->u1_direct_spatial_mv_pred_flag);

    return OK;
}

static WORD32 imvcd_set_ref_idx_override_flag(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u1_num_ref_idx_active_override_flag = ih264d_get_bit_h264(ps_bitstrm);

    COPYTHECONTEXT("SH: num_ref_idx_override_flag",
                   ps_cur_slice->u1_num_ref_idx_active_override_flag);

    return OK;
}

static WORD32 imvcd_set_num_ref_idx_active(dec_struct_t *ps_view_ctxt, UWORD8 *pu1_num_ref_idx)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    UWORD32 u4_num_ref_idx_m1 = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(u4_num_ref_idx_m1 >= H264_MAX_REF_PICS)
    {
        return ERROR_NUM_REF;
    }

    pu1_num_ref_idx[0] = 1 + u4_num_ref_idx_m1;

    COPYTHECONTEXT("SH: num_ref_idx_lx_active_minus1", u4_num_ref_idx_m1);

    return OK;
}

static WORD32 imvcd_set_ref_pic_list_reordering_flag(dec_struct_t *ps_view_ctxt,
                                                     UWORD8 *pu1_ref_idx_reorder_flag)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pu1_ref_idx_reorder_flag[0] = ih264d_get_bit_h264(ps_bitstrm);

    COPYTHECONTEXT("SH: ref_pic_list_reordering_flag_lx", pu1_ref_idx_reorder_flag[0]);

    return OK;
}

static WORD32 imvcd_set_modification_of_pic_nums_idc(dec_struct_t *ps_view_ctxt,
                                                     UWORD8 *pu1_modification_of_pic_nums_idc)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pu1_modification_of_pic_nums_idc[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    COPYTHECONTEXT("SH: modification_of_pic_nums_idc", pu1_modification_of_pic_nums_idc[0]);

    return OK;
}

static WORD32 imvcd_set_abs_diff_pic_num_minus1(dec_struct_t *ps_view_ctxt,
                                                WORD32 *pi4_abs_diff_pic_num_minus1)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pi4_abs_diff_pic_num_minus1[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    COPYTHECONTEXT("SH: abs_diff_pic_num_minus1", pi4_abs_diff_pic_num_minus1[0]);

    return OK;
}

static WORD32 imvcd_set_abs_diff_view_idx_minus1(dec_struct_t *ps_view_ctxt,
                                                 WORD32 *pi4_abs_diff_view_idx_minus1)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pi4_abs_diff_view_idx_minus1[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    COPYTHECONTEXT("SH: abs_diff_view_idx_minus1", pi4_abs_diff_view_idx_minus1[0]);

    return OK;
}

static WORD32 imvcd_set_long_term_pic_num(dec_struct_t *ps_view_ctxt, WORD32 *pi4_long_term_pic_num)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    pi4_long_term_pic_num[0] = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    COPYTHECONTEXT("SH: long_term_pic_num", pi4_long_term_pic_num[0]);

    return OK;
}

static WORD32 imvcd_set_cabac_init_idc(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;

    ps_cur_slice->u1_cabac_init_idc = ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

    if(ps_cur_slice->u1_cabac_init_idc > MAX_CABAC_INIT_IDC)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    COPYTHECONTEXT("SH: cabac_init_idc", ps_cur_slice->u1_cabac_init_idc);

    return OK;
}

static WORD32 imvcd_set_slice_qp(dec_struct_t *ps_view_ctxt)
{
    WORD8 i1_slice_qp_delta;

    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    dec_pic_params_t *ps_cur_pps = ps_view_ctxt->ps_cur_pps;

    i1_slice_qp_delta = ih264d_sev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);
    ps_cur_slice->u1_slice_qp = i1_slice_qp_delta + ps_cur_pps->u1_pic_init_qp;

    if(ps_cur_slice->u1_slice_qp > MAX_H264_QP)
    {
        return ERROR_INV_RANGE_QP_T;
    }

    COPYTHECONTEXT("SH: slice_qp_delta", i1_slice_qp_delta);

    return OK;
}

static WORD32 imvcd_set_slice_deblk_params(dec_struct_t *ps_view_ctxt)
{
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    dec_pic_params_t *ps_cur_pps = ps_view_ctxt->ps_cur_pps;

    if(ps_cur_pps->u1_deblocking_filter_parameters_present_flag)
    {
        ps_cur_slice->u1_disable_dblk_filter_idc =
            ih264d_uev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer);

        if(ps_cur_slice->u1_disable_dblk_filter_idc > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }

        COPYTHECONTEXT("SH: disable_deblocking_filter_idc",
                       ps_cur_slice->u1_disable_dblk_filter_idc);

        if(ps_cur_slice->u1_disable_dblk_filter_idc != 1)
        {
            ps_cur_slice->i1_slice_alpha_c0_offset =
                ih264d_sev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer) << 1;

            if((MIN_DBLK_FIL_OFF > ps_cur_slice->i1_slice_alpha_c0_offset) ||
               (ps_cur_slice->i1_slice_alpha_c0_offset > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            COPYTHECONTEXT("SH: slice_alpha_c0_offset_div2",
                           ps_cur_slice->i1_slice_alpha_c0_offset >> 1);

            ps_cur_slice->i1_slice_beta_offset =
                ih264d_sev(&ps_bitstrm->u4_ofst, ps_bitstrm->pu4_buffer) << 1;

            if((MIN_DBLK_FIL_OFF > ps_cur_slice->i1_slice_beta_offset) ||
               (ps_cur_slice->i1_slice_beta_offset > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            COPYTHECONTEXT("SH: slice_beta_offset_div2", ps_cur_slice->i1_slice_beta_offset >> 1);
        }
        else
        {
            ps_cur_slice->i1_slice_alpha_c0_offset = 0;
            ps_cur_slice->i1_slice_beta_offset = 0;
        }
    }
    else
    {
        ps_cur_slice->u1_disable_dblk_filter_idc = 0;
        ps_cur_slice->i1_slice_alpha_c0_offset = 0;
        ps_cur_slice->i1_slice_beta_offset = 0;
    }

    return OK;
}

static WORD32 imvcd_set_ref_pic_list_mod_data(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i4_error_code;
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    ref_pic_list_mod_data_t *ps_ref_pic_list_mod_data =
        imvcd_get_cur_ref_pic_list_mod_data(ps_mvcd_ctxt);

    bool b_is_b_pic = ps_cur_slice->u1_slice_type == BSLICE;

    for(i = 0; i < 1 + ((WORD32) b_is_b_pic); i++)
    {
        ps_ref_pic_list_mod_data->au1_num_active_refs[i] =
            ps_cur_slice->u1_num_ref_idx_lx_active[i];

        i4_error_code = imvcd_set_ref_pic_list_reordering_flag(
            ps_view_ctxt, &ps_ref_pic_list_mod_data->au1_ref_pic_list_modification_flag_lx[i]);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        if(ps_ref_pic_list_mod_data->au1_ref_pic_list_modification_flag_lx[i])
        {
            UWORD8 *pu1_modification_of_pic_nums_idc =
                ps_ref_pic_list_mod_data->au1_modification_of_pic_nums_idc[i];
            WORD32 *pi4_abs_diff_pic_num_minus1 =
                ps_ref_pic_list_mod_data->ai4_abs_diff_pic_num_minus1[i];
            WORD32 *pi4_long_term_pic_num = ps_ref_pic_list_mod_data->ai4_long_term_pic_num[i];
            WORD32 *pi4_abs_diff_view_idx_minus1 =
                ps_ref_pic_list_mod_data->ai4_abs_diff_view_idx_minus1[i];
            UWORD32 u4_pic_num_mod_count = 0;

            do
            {
                i4_error_code = imvcd_set_modification_of_pic_nums_idc(
                    ps_view_ctxt, pu1_modification_of_pic_nums_idc);

                if(OK != i4_error_code)
                {
                    return i4_error_code;
                }

                if((0 == pu1_modification_of_pic_nums_idc[0]) ||
                   (1 == pu1_modification_of_pic_nums_idc[0]))
                {
                    i4_error_code = imvcd_set_abs_diff_pic_num_minus1(ps_view_ctxt,
                                                                      pi4_abs_diff_pic_num_minus1);

                    if(OK != i4_error_code)
                    {
                        return i4_error_code;
                    }
                }
                else if(2 == pu1_modification_of_pic_nums_idc[0])
                {
                    i4_error_code =
                        imvcd_set_long_term_pic_num(ps_view_ctxt, pi4_long_term_pic_num);

                    if(OK != i4_error_code)
                    {
                        return i4_error_code;
                    }
                }
                else if((4 == pu1_modification_of_pic_nums_idc[0]) ||
                        (5 == pu1_modification_of_pic_nums_idc[0]))
                {
                    i4_error_code = imvcd_set_abs_diff_view_idx_minus1(
                        ps_view_ctxt, pi4_abs_diff_view_idx_minus1);

                    if(OK != i4_error_code)
                    {
                        return i4_error_code;
                    }
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
                u4_pic_num_mod_count++;

                if(u4_pic_num_mod_count > ps_ref_pic_list_mod_data->au1_num_active_refs[i])
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
            } while(true);
        }
    }

    return OK;
}

static WORD32 imvcd_decode_gaps_in_frame_num(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    pocstruct_t s_tmp_poc;

    UWORD32 u4_start_frm_num;
    WORD32 i4_poc;
    WORD8 i1_gap_idx;
    WORD8 *pi1_gaps_per_seq;
    WORD32 i4_error_code;
    WORD64 i8_display_poc;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    dec_pic_params_t *ps_pps = ps_view_ctxt->ps_cur_pps;
    mvc_dpb_manager_t *ps_dpb_mgr = ps_mvcd_ctxt->ps_dpb_mgr;

    UWORD16 u2_frame_num = ps_cur_slice->u2_frame_num;
    UWORD32 u4_next_frm_num = ps_view_ctxt->u2_prev_ref_frame_num + 1;
    UWORD32 u4_max_frm_num = ps_view_ctxt->ps_cur_sps->u2_u4_max_pic_num_minus1 + 1;
    WORD32 *pi4_gaps_start_frm_num = ps_dpb_mgr->ai4_gaps_start_frm_num;
    bool b_is_idr_slice = imvcd_is_idr_au(ps_mvcd_ctxt);

    if(ps_cur_slice->u1_field_pic_flag)
    {
        if(ps_view_ctxt->u2_prev_ref_frame_num == u2_frame_num)
        {
            return OK;
        }
    }

    if(u4_next_frm_num >= u4_max_frm_num)
    {
        u4_next_frm_num -= u4_max_frm_num;
    }

    if(u4_next_frm_num == u2_frame_num)
    {
        return OK;
    }

    if(b_is_idr_slice && (u4_next_frm_num >= u2_frame_num))
    {
        return OK;
    }

    u4_start_frm_num = u4_next_frm_num;

    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[0] = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[1] = 0;

    for(i1_gap_idx = 0; i1_gap_idx < MAX_FRAMES; i1_gap_idx++)
    {
        if(INVALID_FRAME_NUM == pi4_gaps_start_frm_num[i1_gap_idx])
        {
            break;
        }
    }

    if(MAX_FRAMES == i1_gap_idx)
    {
        return ERROR_DBP_MANAGER_T;
    }

    i4_poc = 0;
    pi4_gaps_start_frm_num[i1_gap_idx] = u4_start_frm_num;
    ps_dpb_mgr->ai4_gaps_end_frm_num[i1_gap_idx] = u2_frame_num - 1;
    pi1_gaps_per_seq = ps_dpb_mgr->ai1_gaps_per_seq;
    pi1_gaps_per_seq[i1_gap_idx] = 0;

    while(u4_next_frm_num != u2_frame_num)
    {
        imvcd_dpb_delete_nonref_nondisplay_pics(ps_dpb_mgr);

        if(ps_pps->ps_sps->u1_pic_order_cnt_type)
        {
            /* allocate a picture buffer and insert it as ST node */
            i4_error_code =
                ih264d_decode_pic_order_cnt(0, u4_next_frm_num, &ps_view_ctxt->s_prev_pic_poc,
                                            &s_tmp_poc, ps_cur_slice, ps_pps, 1, 0, 0, &i4_poc);

            if(i4_error_code != OK)
            {
                return i4_error_code;
            }

            /* Display seq no calculations */
            if(i4_poc >= ps_view_ctxt->i4_max_poc)
            {
                ps_view_ctxt->i4_max_poc = i4_poc;
            }

            /* IDR Picture or POC wrap around */
            if(i4_poc == 0)
            {
                imvcd_modulate_max_disp_seq(ps_view_ctxt);
            }

            ps_cur_slice->u1_mmco_equalto5 = 0;
            ps_cur_slice->u2_frame_num = u4_next_frm_num;
        }

        if(ps_dpb_mgr->i1_poc_buf_id_entries >= ps_view_ctxt->u1_max_dec_frame_buffering)
        {
            i4_error_code = imvcd_dpb_assign_display_seq(ps_mvcd_ctxt->ps_dpb_mgr);

            if(i4_error_code != OK)
            {
                return i4_error_code;
            }
        }

        i8_display_poc = ((WORD64) ps_view_ctxt->i4_prev_max_display_seq) + ((WORD64) i4_poc);

        if(IS_OUT_OF_RANGE_S32(i8_display_poc))
        {
            ps_view_ctxt->i4_prev_max_display_seq = 0;
            i8_display_poc = i4_poc;
        }

        i4_error_code = imvcd_dpb_insert_pic_in_display_list(ps_dpb_mgr, (WORD32) i8_display_poc,
                                                             u4_next_frm_num, DO_NOT_DISP);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        pi1_gaps_per_seq[i1_gap_idx]++;

        i4_error_code =
            imvcd_dpb_do_mmco_for_gaps(ps_dpb_mgr, ps_view_ctxt->ps_cur_sps->u1_num_ref_frames);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        imvcd_dpb_delete_nonref_nondisplay_pics(ps_dpb_mgr);

        u4_next_frm_num++;

        if(u4_next_frm_num >= u4_max_frm_num)
        {
            u4_next_frm_num -= u4_max_frm_num;
        }
    }

    return OK;
}

static void imvcd_pocstruct_init(dec_struct_t *ps_view_ctxt)
{
    pocstruct_t *ps_prev_poc = &ps_view_ctxt->s_prev_pic_poc;
    pocstruct_t *ps_cur_poc = &ps_view_ctxt->s_cur_pic_poc;

    ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst;
    ps_prev_poc->u2_frame_num = ps_cur_poc->u2_frame_num;
    ps_prev_poc->u1_mmco_equalto5 = ps_cur_poc->u1_mmco_equalto5;

    if(ps_view_ctxt->ps_cur_slice->u1_nal_ref_idc)
    {
        ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb;
        ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb;
        ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom;
        ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0];
        ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1];
        ps_prev_poc->u1_bot_field = ps_cur_poc->u1_bot_field;
    }
}

static WORD32 imvcd_pic_init(mvc_dec_ctxt_t *ps_mvcd_ctxt, pocstruct_t *ps_cur_poc, WORD32 i4_poc,
                             bool b_is_idr_slice)
{
    WORD32 i4_error_code;
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    pocstruct_t *ps_prev_poc = &ps_view_ctxt->s_cur_pic_poc;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    dec_pic_params_t *ps_pps = ps_view_ctxt->ps_cur_pps;
    dec_seq_params_t *ps_sps = ps_pps->ps_sps;
    subset_sps_t *ps_subset_sps = imvcd_get_valid_subset_sps(ps_mvcd_ctxt);
    nalu_mvc_ext_t *ps_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);
    dec_err_status_t *ps_err = ps_view_ctxt->ps_dec_err_status;
    prev_seq_params_t *ps_prev_seq_params = &ps_view_ctxt->s_prev_seq_params;

    UWORD16 u2_num_views = ps_mvcd_ctxt->u2_num_views;
    UWORD16 u2_view_order_id = ps_mvcd_ctxt->u2_num_views_decoded;
    UWORD16 u2_view_id = ps_nalu_mvc_ext->u2_view_id;
    UWORD16 u2_frame_num = ps_cur_slice->u2_frame_num;

    ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb;
    ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb;
    ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom;
    ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0];
    ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1];
    ps_prev_poc->u1_bot_field = ps_view_ctxt->ps_cur_slice->u1_bottom_field_flag;
    ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst;
    ps_prev_poc->u2_frame_num = u2_frame_num;

    ps_view_ctxt->i1_prev_mb_qp_delta = 0;
    ps_view_ctxt->i1_next_ctxt_idx = 0;
    ps_view_ctxt->u4_use_intrapred_line_copy = 1;

    if(ps_view_ctxt->u4_num_cores == 1)
    {
        ps_view_ctxt->u4_nmb_deblk = 1;
    }
    else
    {
        ps_view_ctxt->u4_nmb_deblk = 0;
    }

    ps_view_ctxt->u4_app_disable_deblk_frm = 0;
    if(ps_view_ctxt->i4_degrade_type && ps_view_ctxt->i4_degrade_pics)
    {
        WORD32 i4_degrade_pic = 0;

        ps_view_ctxt->i4_degrade_pic_cnt++;

        /* If degrade is to be done in all frames, then do not check further */
        switch(ps_view_ctxt->i4_degrade_pics)
        {
            case 4:
            {
                i4_degrade_pic = 1;

                break;
            }
            case 3:
            {
                if(ps_cur_slice->u1_slice_type != I_SLICE)
                {
                    i4_degrade_pic = 1;
                }

                break;
            }
            case 2:
            {
                if((ps_cur_slice->u1_slice_type != I_SLICE) &&
                   (ps_view_ctxt->i4_degrade_pic_cnt != ps_view_ctxt->i4_nondegrade_interval))
                {
                    i4_degrade_pic = 1;
                }

                break;
            }
            case 1:
            {
                if(0 == ps_cur_slice->u1_nal_ref_idc)
                {
                    i4_degrade_pic = 1;
                }

                break;
            }
        }

        if(i4_degrade_pic)
        {
            if(ps_view_ctxt->i4_degrade_type & 0x2)
            {
                ps_view_ctxt->u4_app_disable_deblk_frm = 1;
            }

            if(0 == ps_cur_slice->u1_nal_ref_idc)
            {
                if(ps_view_ctxt->i4_degrade_type & 0x4)
                {
                    ps_view_ctxt->i4_mv_frac_mask = 0;
                }

                if(ps_view_ctxt->i4_degrade_type & 0x8)
                {
                    ps_view_ctxt->i4_mv_frac_mask = 0;
                }
            }
        }
        else
        {
            ps_view_ctxt->i4_degrade_pic_cnt = 0;
        }
    }

    if((ps_cur_slice->u1_slice_type == I_SLICE) || (ps_cur_slice->u1_slice_type == SI_SLICE))
    {
        ps_err->u1_cur_pic_type = PIC_TYPE_I;
    }
    else
    {
        ps_err->u1_cur_pic_type = PIC_TYPE_UNKNOWN;
    }

    if(ps_err->u1_pic_aud_i == PIC_TYPE_I)
    {
        ps_err->u1_cur_pic_type = PIC_TYPE_I;
        ps_err->u1_pic_aud_i = PIC_TYPE_UNKNOWN;
    }

    if(b_is_idr_slice)
    {
        if(ps_err->u1_err_flag)
        {
            imvcd_reset_dpb(ps_mvcd_ctxt->ps_dpb_mgr);
        }

        ps_err->u1_err_flag = ACCEPT_ALL_PICS;
    }

    if(ps_view_ctxt->u1_init_dec_flag && ps_view_ctxt->s_prev_seq_params.u1_eoseq_pending &&
       (u2_view_order_id == (u2_num_views - 1)))
    {
        imvcd_release_all_ref_and_io_bufs(ps_mvcd_ctxt, MAX_DISP_BUFS_NEW);

        ps_view_ctxt->u1_second_field = 0;
        ps_view_ctxt->i4_cur_display_seq = 0;
        ps_view_ctxt->s_prev_seq_params.u1_eoseq_pending = 0;

        imvcd_dpb_set_display_num(ps_mvcd_ctxt->ps_dpb_mgr, 0);
    }

    if(0 == u2_view_order_id)
    {
        imvcd_dpb_set_max_pic_num(ps_mvcd_ctxt->ps_dpb_mgr, ps_sps->u2_u4_max_pic_num_minus1 + 1);
        imvcd_dpb_set_num_views(ps_mvcd_ctxt->ps_dpb_mgr, u2_num_views);
    }

    ps_view_ctxt->i4_pic_type = NA_SLICE;
    ps_view_ctxt->i4_frametype = IV_NA_FRAME;
    ps_view_ctxt->i4_content_type = IV_CONTENTTYPE_NA;

    ps_sps->u4_max_mb_addr = ps_sps->u2_frm_wd_in_mbs * ps_sps->u2_frm_ht_in_mbs - 1;
    ps_view_ctxt->u2_frm_ht_in_mbs = ps_sps->u2_frm_ht_in_mbs;

    if(!ps_view_ctxt->u1_init_dec_flag)
    {
        ps_view_ctxt->u1_max_dec_frame_buffering = ih264d_get_dpb_size(ps_sps);

        ps_view_ctxt->i4_display_delay = ps_view_ctxt->u1_max_dec_frame_buffering;

        if(ps_sps->u1_vui_parameters_present_flag && ps_sps->s_vui.u1_bitstream_restriction_flag)
        {
            if(ps_sps->u1_frame_mbs_only_flag)
            {
                ps_view_ctxt->i4_display_delay = ps_sps->s_vui.u4_num_reorder_frames + 1;
            }
            else
            {
                ps_view_ctxt->i4_display_delay = ps_sps->s_vui.u4_num_reorder_frames * 2 + 2;
            }
        }

        if(IVD_DECODE_FRAME_OUT == ps_view_ctxt->e_frm_out_mode)
        {
            ps_view_ctxt->i4_display_delay = 0;
        }

        imvcd_dpb_set_display_delay(ps_mvcd_ctxt->ps_dpb_mgr, ps_view_ctxt->i4_display_delay);

        ps_view_ctxt->u1_pic_bufs = ps_view_ctxt->i4_display_delay + ps_sps->u1_num_ref_frames + 1;
        ps_view_ctxt->u1_pic_bufs += imvcd_get_max_num_ivp_refs(ps_mvcd_ctxt);
        ps_view_ctxt->u1_pic_bufs = CLIP3(2, MVC_MAX_REF_PICS, ps_view_ctxt->u1_pic_bufs);

        ps_view_ctxt->u1_max_dec_frame_buffering =
            MIN(ps_view_ctxt->u1_max_dec_frame_buffering, ps_view_ctxt->u1_pic_bufs);

        /*********************************************************************/
        /* Configuring decoder parameters based on level and then            */
        /* fresh pointer initialisation in decoder scratch and state buffers */
        /*********************************************************************/
        i4_error_code = ih264d_init_dec_mb_grp(ps_view_ctxt);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        i4_error_code = imvcd_allocate_dynamic_bufs(ps_mvcd_ctxt);

        if(i4_error_code != OK)
        {
            imvcd_free_dynamic_bufs(ps_mvcd_ctxt);

            return IVD_MEM_ALLOC_FAILED;
        }

        i4_error_code = imvcd_init_au_buffers(ps_mvcd_ctxt);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        i4_error_code = imvcd_init_au_mv_pred_bufs(ps_mvcd_ctxt);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        ps_view_ctxt->u1_init_dec_flag = 1;
        ps_prev_seq_params->u2_frm_wd_in_mbs = ps_sps->u2_frm_wd_in_mbs;
        ps_prev_seq_params->u1_level_idc = ps_sps->u1_level_idc;
        ps_prev_seq_params->u1_profile_idc = ps_sps->u1_profile_idc;
        ps_prev_seq_params->u2_frm_ht_in_mbs = ps_sps->u2_frm_ht_in_mbs;
        ps_prev_seq_params->u1_frame_mbs_only_flag = ps_sps->u1_frame_mbs_only_flag;
        ps_prev_seq_params->u1_direct_8x8_inference_flag = ps_sps->u1_direct_8x8_inference_flag;

        ps_view_ctxt->i4_cur_display_seq = 0;
        ps_view_ctxt->i4_prev_max_display_seq = 0;
        ps_view_ctxt->i4_max_poc = 0;

        imvcd_dpb_set_display_num(ps_mvcd_ctxt->ps_dpb_mgr, 0);

        {
            /* 0th entry of CtxtIncMbMap will be always be containing default values
             for CABAC context representing MB not available */
            ctxt_inc_mb_info_t *p_DefCtxt = ps_view_ctxt->p_ctxt_inc_mb_map - 1;
            UWORD8 *pu1_temp;

            p_DefCtxt->u1_mb_type = CAB_SKIP;

            p_DefCtxt->u1_cbp = 0x0f;
            p_DefCtxt->u1_intra_chroma_pred_mode = 0;

            p_DefCtxt->u1_yuv_dc_csbp = 0x7;

            p_DefCtxt->u1_transform8x8_ctxt = 0;

            pu1_temp = (UWORD8 *) p_DefCtxt->i1_ref_idx;
            for(i = 0; i < 4; i++, pu1_temp++)
            {
                (*pu1_temp) = 0;
            }

            pu1_temp = (UWORD8 *) p_DefCtxt->u1_mv;
            for(i = 0; i < 16; i++, pu1_temp++)
            {
                (*pu1_temp) = 0;
            }

            ps_view_ctxt->ps_def_ctxt_mb_info = p_DefCtxt;
        }
    }

    /* reset DBP commands read u4_flag */
    ps_view_ctxt->ps_dpb_cmds->u1_dpb_commands_read = 0;

    ps_view_ctxt->pv_parse_tu_coeff_data = ps_view_ctxt->pv_pic_tu_coeff_data;
    ps_view_ctxt->pv_proc_tu_coeff_data = ps_view_ctxt->pv_pic_tu_coeff_data;
    ps_view_ctxt->ps_nmb_info = ps_view_ctxt->ps_frm_mb_info;

    if(ps_view_ctxt->u1_separate_parse)
    {
        UWORD32 num_mbs;

        num_mbs = ps_view_ctxt->ps_cur_sps->u4_total_num_of_mbs;

        if(ps_view_ctxt->pu1_dec_mb_map)
        {
            memset((void *) ps_view_ctxt->pu1_dec_mb_map, 0, num_mbs);
        }

        if(ps_view_ctxt->pu1_recon_mb_map)
        {
            memset((void *) ps_view_ctxt->pu1_recon_mb_map, 0, num_mbs);
        }

        if(ps_view_ctxt->pu2_slice_num_map)
        {
            memset((void *) ps_view_ctxt->pu2_slice_num_map, 0, (num_mbs * sizeof(UWORD16)));
        }
    }

    ps_view_ctxt->ps_parse_cur_slice = &(ps_view_ctxt->ps_dec_slice_buf[0]);
    ps_view_ctxt->ps_decode_cur_slice = &(ps_view_ctxt->ps_dec_slice_buf[0]);
    ps_view_ctxt->ps_computebs_cur_slice = &(ps_view_ctxt->ps_dec_slice_buf[0]);
    ps_view_ctxt->u2_cur_slice_num = 0;

    ps_view_ctxt->s_high_profile.u1_scaling_present = 0;
    ps_view_ctxt->s_high_profile.u1_transform8x8_present = 0;

    if(0 == u2_view_order_id)
    {
        mvc_au_buffer_t *ps_cur_au;
        mvc_au_mv_pred_t *ps_au_mv_data;

        WORD32 i4_pic_buf_id, i4_mv_buf_id;

        ps_cur_au = (mvc_au_buffer_t *) ih264_buf_mgr_get_next_free(
            ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt, &i4_pic_buf_id);

        if(NULL == ps_cur_au)
        {
            return ERROR_UNAVAIL_PICBUF_T;
        }
        else
        {
            /* Buf will alwys be marked as REF here to ensure IVP works */
            /* If AU nalRefIdc=0, REF status will be removed during endOfAU processing
             */
            ih264_buf_mgr_set_status(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt, i4_pic_buf_id,
                                     BUF_MGR_IO | BUF_MGR_REF);
        }

        ps_au_mv_data = (mvc_au_mv_pred_t *) ih264_buf_mgr_get_next_free(
            ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt, &i4_mv_buf_id);

        if(ps_au_mv_data == NULL)
        {
            return ERROR_UNAVAIL_MVBUF_T;
        }
        else
        {
            /* Buf will alwys be marked as REF here to ensure IVP works */
            /* If AU nalRefIdc=0, REF status will be removed during endOfAU processing
             */
            ih264_buf_mgr_set_status(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                                     i4_mv_buf_id, BUF_MGR_REF);
        }

        ps_mvcd_ctxt->ps_cur_au = ps_cur_au;

        ps_cur_au->s_sei_pic = ps_view_ctxt->ps_sei[0];

        ps_cur_au->i4_mv_buf_id = i4_mv_buf_id;
        ps_cur_au->ps_au_mv_data = ps_au_mv_data;
        ps_cur_au->i4_poc = i4_poc;
        ps_cur_au->i4_avg_poc = i4_poc;
        ps_cur_au->i4_frame_num = u2_frame_num;
        ps_cur_au->i4_pic_num = u2_frame_num;
        ps_cur_au->u4_time_stamp = ps_view_ctxt->u4_ts;
        ps_cur_au->u1_picturetype = FRM_PIC;
        ps_cur_au->u2_disp_width = ps_view_ctxt->u2_disp_width;
        ps_cur_au->u2_disp_height = ps_view_ctxt->u2_disp_height;

        memset(ps_cur_au->au4_pack_slc_typ, 0, sizeof(ps_cur_au->au4_pack_slc_typ));

        ps_mvcd_ctxt->s_mvc_au_buf_mgr.au1_au_buf_id_to_mv_buf_id_map[i4_pic_buf_id] = i4_mv_buf_id;
        ps_mvcd_ctxt->s_mvc_au_buf_mgr.aps_buf_id_to_au_buf_map[i4_pic_buf_id] = ps_cur_au;
        ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.aps_buf_id_to_mv_pred_buf_map[i4_mv_buf_id] =
            ps_au_mv_data;

        ps_view_ctxt->au1_pic_buf_ref_flag[i4_pic_buf_id] = 0;

        ps_cur_au->s_ivp_data.b_is_ivp_ref = false;

        imvcd_dpb_init_au_bufs(ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_au);
    }

    if(u2_view_order_id > 0)
    {
        ps_mvcd_ctxt->ps_cur_au->as_disp_offsets[u2_view_id] =
            ps_mvcd_ctxt->aps_pps_id_to_subset_sps_map[ps_pps->u1_pic_parameter_set_id]
                ->s_disp_offsets;
    }
    else
    {
        /* Accounting for lihbavc's idiocy */
        ps_mvcd_ctxt->ps_cur_au->as_disp_offsets[u2_view_id].u2_left_offset =
            ps_view_ctxt->u2_crop_offset_y;
        ps_mvcd_ctxt->ps_cur_au->as_disp_offsets[u2_view_id].u2_right_offset = 0;
        ps_mvcd_ctxt->ps_cur_au->as_disp_offsets[u2_view_id].u2_top_offset = 0;
        ps_mvcd_ctxt->ps_cur_au->as_disp_offsets[u2_view_id].u2_bottom_offset = 0;
    }

    for(i = 0; i < 2; i++)
    {
        ps_view_ctxt->ps_ref_pic_buf_lx[i] = imvcd_dpb_get_view_ref_pic_list(
            ps_mvcd_ctxt->ps_dpb_mgr, u2_view_order_id, u2_view_id, i);

        imvcd_set_view_buf_id_to_buf_map(ps_view_ctxt);
    }

    if(ps_mvcd_ctxt->u2_num_views > 1)
    {
        imvcd_dpb_init_view_bufs(ps_mvcd_ctxt->ps_dpb_mgr, u2_view_order_id, u2_view_id);

        imvcd_dpb_init_ivp_ctxt(ps_mvcd_ctxt->ps_dpb_mgr, &ps_subset_sps->s_sps_mvc_ext,
                                ps_mvcd_ctxt->as_nalu_mvc_ext);
    }

    ps_view_ctxt->u4_pic_buf_got = 1;
    ps_cur_slice->u1_mbaff_frame_flag = 0;

    ps_view_ctxt->ps_cur_mb_row = ps_view_ctxt->ps_nbr_mb_row;
    // Increment by 2 ,so that left mb (mbaff decrements by 2)  will always be
    // valid
    ps_view_ctxt->ps_cur_mb_row += 2;
    ps_view_ctxt->ps_top_mb_row = ps_view_ctxt->ps_nbr_mb_row;
    ps_view_ctxt->ps_top_mb_row += ps_view_ctxt->u2_frm_wd_in_mbs + 2;
    // Increment by 2 ,so that left mb (mbaff decrements by 2)  will always be
    // valid
    ps_view_ctxt->ps_top_mb_row += 2;
    ps_view_ctxt->u4_mb_idx = 0;
    ps_view_ctxt->u4_total_mbs_coded = 0;
    ps_view_ctxt->i4_submb_ofst = -(SUB_BLK_SIZE);
    ps_view_ctxt->i2_prev_slice_mbx = -1;
    ps_view_ctxt->i2_prev_slice_mby = 0;

    ps_view_ctxt->u4_pred_info_idx = 0;
    ps_view_ctxt->u4_pred_info_pkd_idx = 0;
    ps_view_ctxt->ps_part = ps_view_ctxt->ps_parse_part_params;

    ps_view_ctxt->u4_dma_buf_idx = 0;

    ps_view_ctxt->ps_mv_cur = ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->aps_mvs[u2_view_id];
    ps_view_ctxt->ps_mv_top = ps_view_ctxt->ps_mv_top_p[0];
    ps_view_ctxt->u1_mv_top_p = 0;
    ps_view_ctxt->ps_mv_left = ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->aps_mvs[u2_view_id];
    ps_view_ctxt->ps_mv = ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->aps_mvs[u2_view_id];
    ps_view_ctxt->ps_mv_bank_cur = ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->aps_mvs[u2_view_id];
    ps_view_ctxt->pu1_col_zero_flag =
        ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->apu1_mode_descriptors[u2_view_id];
    ps_view_ctxt->u2_mv_2mb[0] = 0;
    ps_view_ctxt->u2_mv_2mb[1] = 0;

    ps_view_ctxt->u1_last_pic_not_decoded = 0;
    ps_view_ctxt->u2_cur_slice_num_dec_thread = 0;
    ps_view_ctxt->u2_cur_slice_num_bs = 0;

    ps_view_ctxt->u4_intra_pred_line_ofst = 0;
    ps_view_ctxt->pu1_cur_y_intra_pred_line = ps_view_ctxt->pu1_y_intra_pred_line;
    ps_view_ctxt->pu1_cur_u_intra_pred_line = ps_view_ctxt->pu1_u_intra_pred_line;
    ps_view_ctxt->pu1_cur_v_intra_pred_line = ps_view_ctxt->pu1_v_intra_pred_line;
    ps_view_ctxt->pu1_cur_y_intra_pred_line_base = ps_view_ctxt->pu1_y_intra_pred_line;
    ps_view_ctxt->pu1_cur_u_intra_pred_line_base = ps_view_ctxt->pu1_u_intra_pred_line;
    ps_view_ctxt->pu1_cur_v_intra_pred_line_base = ps_view_ctxt->pu1_v_intra_pred_line;
    ps_view_ctxt->pu1_prev_y_intra_pred_line =
        ps_view_ctxt->pu1_y_intra_pred_line + (ps_view_ctxt->u2_frm_wd_in_mbs * MB_SIZE);
    ps_view_ctxt->pu1_prev_u_intra_pred_line =
        ps_view_ctxt->pu1_u_intra_pred_line +
        ps_view_ctxt->u2_frm_wd_in_mbs * BLK8x8SIZE * YUV420SP_FACTOR;
    ps_view_ctxt->pu1_prev_v_intra_pred_line =
        ps_view_ctxt->pu1_v_intra_pred_line + ps_view_ctxt->u2_frm_wd_in_mbs * BLK8x8SIZE;

    ps_view_ctxt->ps_deblk_mbn = ps_view_ctxt->ps_deblk_pic;

    ps_view_ctxt->pf_compute_bs = ih264d_compute_bs_non_mbaff;
    ps_view_ctxt->u1_cur_mb_fld_dec_flag = ps_cur_slice->u1_field_pic_flag;

    if(0 == u2_view_order_id)
    {
        imvcd_assign_pic_num(ps_mvcd_ctxt->ps_dpb_mgr, ps_sps->u2_u4_max_pic_num_minus1 + 1,
                             ps_mvcd_ctxt->ps_cur_au->i4_frame_num,
                             ps_sps->u1_gaps_in_frame_num_value_allowed_flag);

        ps_view_ctxt->s_tran_addrecon.u2_mv_top_left_inc = (ps_view_ctxt->u4_recon_mb_grp << 2) - 1;
        ps_view_ctxt->s_tran_addrecon.u2_mv_left_inc = (ps_view_ctxt->u4_recon_mb_grp - 1) << 4;
    }

    if((ps_sps->u1_profile_idc == HIGH_PROFILE_IDC) ||
       (ps_sps->u1_profile_idc == MULTIVIEW_HIGH_PROFILE_IDC))
    {
        if((ps_sps->i4_seq_scaling_matrix_present_flag) ||
           (ps_pps->i4_pic_scaling_matrix_present_flag))
        {
            i4_error_code = ih264d_form_scaling_matrix_picture(ps_sps, ps_pps, ps_view_ctxt);
            ps_view_ctxt->s_high_profile.u1_scaling_present = 1;
        }
        else
        {
            i4_error_code = ih264d_form_default_scaling_matrix(ps_view_ctxt);
        }

        if(ps_pps->i4_transform_8x8_mode_flag)
        {
            ps_view_ctxt->s_high_profile.u1_transform8x8_present = 1;
        }
    }
    else
    {
        i4_error_code = ih264d_form_default_scaling_matrix(ps_view_ctxt);
    }

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    ps_view_ctxt->s_high_profile.u1_direct_8x8_inference_flag =
        ps_sps->u1_direct_8x8_inference_flag;
    ps_view_ctxt->s_high_profile.s_cavlc_ctxt = ps_view_ctxt->s_cavlc_ctxt;

    ps_view_ctxt->i1_recon_in_thread3_flag = 1;

    ps_view_ctxt->ps_cur_pic = &ps_view_ctxt->s_cur_pic;
    imvcd_convert_au_buf_to_view_buf(ps_mvcd_ctxt->ps_cur_au, &ps_view_ctxt->s_cur_pic,
                                     u2_view_order_id, u2_view_id);

    ih264d_init_deblk_tfr_ctxt(ps_view_ctxt, &ps_view_ctxt->s_pad_mgr,
                               &ps_view_ctxt->s_tran_addrecon, ps_view_ctxt->u2_frm_wd_in_mbs, 0);

    ps_view_ctxt->ps_frame_buf_ip_recon = &ps_view_ctxt->s_tran_addrecon;

    if(ps_view_ctxt->u1_separate_parse)
    {
        ps_view_ctxt->s_tran_addrecon_parse = ps_view_ctxt->s_tran_addrecon;

        if((ps_view_ctxt->u4_num_cores >= 3) && ps_view_ctxt->i1_recon_in_thread3_flag)
        {
            ps_view_ctxt->s_tran_iprecon = ps_view_ctxt->s_tran_addrecon;
            ps_view_ctxt->ps_frame_buf_ip_recon = &ps_view_ctxt->s_tran_iprecon;
        }
    }

    ps_view_ctxt->ps_cur_deblk_mb = ps_view_ctxt->ps_deblk_pic;
    ps_view_ctxt->u4_cur_deblk_mb_num = 0;

    ps_view_ctxt->u4_deblk_mb_x = 0;
    ps_view_ctxt->u4_deblk_mb_y = 0;
    ps_view_ctxt->pu4_wt_ofsts = ps_view_ctxt->pu4_wts_ofsts_mat;

    ps_view_ctxt->u4_first_slice_in_pic = 0;

    return OK;
}

static WORD32 imvcd_corrupted_slice_handler(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_mb_info_t *ps_cur_mb_info;
    parse_pmbarams_t *ps_parse_mb_data;
    deblk_mb_t *ps_cur_deblk_mb;
    parse_part_params_t *ps_part_info;

    UWORD32 u4_num_mbs_next;
    bool b_is_end_of_row;
    bool b_is_slice_end;
    bool b_tfr_n_mb;
    bool b_decode_nmb;
    UWORD8 u1_inter_mb_type;
    UWORD8 u1_deblk_mb_type;
    UWORD16 i2_cur_mb_addr;
    UWORD32 u4_mb_skip_run;
    WORD32 i, j;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_slice_params_t *ps_slice = ps_view_ctxt->ps_cur_slice;
    nalu_mvc_ext_t *ps_cur_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);

    UWORD32 u4_num_mbs = 0;
    UWORD32 u4_mb_idx = ps_view_ctxt->u4_mb_idx;
    UWORD32 u4_remaining_mbs =
        (ps_view_ctxt->ps_cur_sps->u4_max_mb_addr + 1) - ps_view_ctxt->u4_total_mbs_coded;

    if(ps_view_ctxt->ps_dec_err_status->u1_err_flag & REJECT_CUR_PIC)
    {
        imvcd_free_ref_and_io_bufs(&ps_mvcd_ctxt->s_mvc_au_buf_mgr,
                                   &ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr,
                                   ps_mvcd_ctxt->ps_cur_au->i4_pic_buf_id);

        return OK;
    }

    if((ISLICE == ps_slice->u1_slice_type) || (0 == ps_view_ctxt->u4_total_mbs_coded))
    {
        yuv_buf_props_t *ps_view_buf =
            &ps_mvcd_ctxt->ps_cur_au->as_view_buffers[ps_cur_nalu_mvc_ext->u2_view_id];

        for(i = 0; i < NUM_SP_COMPONENTS; i++)
        {
            buffer_container_t *ps_component_buf = &ps_view_buf->as_component_bufs[i];

            bool b_is_chroma = ((COMPONENT_TYPES_T) i) != Y;
            UWORD16 u2_height = ps_view_buf->u2_height >> b_is_chroma;
            UWORD16 u2_width = ps_view_buf->u2_width;

            for(j = 0; j < u2_height; j++)
            {
                UWORD8 *pu1_data =
                    ((UWORD8 *) ps_component_buf->pv_data) + j * ps_component_buf->i4_data_stride;

                memset(pu1_data, 128, u2_width * sizeof(pu1_data[0]));
            }
        }

        memset(ps_view_ctxt->apv_buf_id_pic_buf_map, 0,
               sizeof(ps_view_ctxt->apv_buf_id_pic_buf_map));

        ps_view_ctxt->apv_buf_id_pic_buf_map[ps_mvcd_ctxt->ps_cur_au->i4_pic_buf_id] =
            &ps_view_ctxt->s_cur_pic;
        ps_view_ctxt->ps_ref_pic_buf_lx[0] = &ps_view_ctxt->ps_cur_pic;
        (ps_view_ctxt->ppv_map_ref_idx_to_poc + FRM_LIST_L0)[0] =
            ps_view_ctxt->ps_cur_pic->pu1_buf1;
        (ps_view_ctxt->ppv_map_ref_idx_to_poc + FRM_LIST_L1)[0] = NULL;
    }

    ps_view_ctxt->ps_dpb_cmds->u1_long_term_reference_flag = 0;

    if(ps_view_ctxt->u4_total_mbs_coded > 0)
    {
        ps_view_ctxt->u4_total_mbs_coded -=
            ps_view_ctxt->u4_total_mbs_coded % ps_view_ctxt->ps_cur_sps->u2_frm_wd_in_mbs;
        u4_remaining_mbs =
            (ps_view_ctxt->ps_cur_sps->u4_max_mb_addr + 1) - ps_view_ctxt->u4_total_mbs_coded;

        while(ps_view_ctxt->u4_dec_thread_created &&
              (ps_view_ctxt->cur_dec_mb_num < ps_view_ctxt->u4_total_mbs_coded))
        {
            NOP(1 << 10);
        }

        while(ps_view_ctxt->u4_bs_deblk_thread_created &&
              (ps_view_ctxt->cur_recon_mb_num < ps_view_ctxt->u4_total_mbs_coded))
        {
            NOP(1 << 10);
        }

        while(ps_view_ctxt->u4_bs_deblk_thread_created &&
              (ps_view_ctxt->u4_cur_deblk_mb_num < ps_view_ctxt->u4_total_mbs_coded))
        {
            NOP(1 << 10);
        }

        ps_view_ctxt->ps_nmb_info = ps_view_ctxt->ps_frm_mb_info + ps_view_ctxt->u4_total_mbs_coded;
        ps_view_ctxt->ps_deblk_mbn = ps_view_ctxt->ps_cur_deblk_mb =
            ps_view_ctxt->ps_deblk_pic + ps_view_ctxt->u4_total_mbs_coded;
    }

    u4_num_mbs = ps_view_ctxt->u4_num_mbs_cur_nmb = 0;

    if(ps_view_ctxt->u1_separate_parse)
    {
        ps_cur_mb_info = ps_view_ctxt->ps_nmb_info;
    }
    else
    {
        ps_cur_mb_info = ps_view_ctxt->ps_nmb_info + ps_view_ctxt->u4_num_mbs_prev_nmb - 1;
    }

    ps_view_ctxt->u2_mby = ps_cur_mb_info->u2_mby;
    ps_view_ctxt->u2_mbx = ps_cur_mb_info->u2_mbx;

    ps_view_ctxt->u1_mb_ngbr_availablity = ps_cur_mb_info->u1_mb_ngbr_availablity;

    if(ps_view_ctxt->u4_total_mbs_coded >= (ps_view_ctxt->ps_cur_sps->u4_max_mb_addr + 1))
    {
        ps_view_ctxt->u1_pic_decode_done = 1;

        return OK;
    }

    /******************************************************/
    /* Initializations to new slice                       */
    /******************************************************/
    ps_view_ctxt->ps_parse_cur_slice->ppv_map_ref_idx_to_poc =
        (volatile void **) ps_view_ctxt->pv_map_ref_idx_to_poc_buf;
    ps_slice->i1_slice_alpha_c0_offset = 0;
    ps_slice->i1_slice_beta_offset = 0;
    ps_slice->u2_first_mb_in_slice = ps_view_ctxt->u4_total_mbs_coded;
    ps_view_ctxt->ps_parse_cur_slice->u4_first_mb_in_slice = ps_view_ctxt->u4_total_mbs_coded;
    ps_view_ctxt->ps_parse_cur_slice->u2_log2Y_crwd = ps_slice->u2_log2Y_crwd;

    if(ps_view_ctxt->u1_separate_parse)
    {
        ps_view_ctxt->ps_parse_cur_slice->pv_tu_coeff_data_start =
            ps_view_ctxt->pv_parse_tu_coeff_data;
    }
    else
    {
        ps_view_ctxt->pv_proc_tu_coeff_data = ps_view_ctxt->pv_parse_tu_coeff_data;
    }

    /******************************************************/
    /* Initializations specific to P slice                */
    /******************************************************/
    u1_inter_mb_type = P_MB;
    u1_deblk_mb_type = D_INTER_MB;

    ps_slice->u1_slice_type = P_SLICE;
    ps_view_ctxt->ps_parse_cur_slice->slice_type = P_SLICE;
    ps_view_ctxt->pf_mvpred_ref_tfr_nby2mb = ih264d_mv_pred_ref_tfr_nby2_pmb;
    ps_view_ctxt->ps_part = ps_view_ctxt->ps_parse_part_params;
    ps_view_ctxt->u2_mbx =
        MOD((WORD32)ps_view_ctxt->u4_total_mbs_coded - 1, ps_view_ctxt->u2_frm_wd_in_mbs);
    ps_view_ctxt->u2_mby =
        DIV((WORD32)ps_view_ctxt->u4_total_mbs_coded - 1, ps_view_ctxt->u2_frm_wd_in_mbs);

    /******************************************************/
    /* Parsing / decoding the slice                       */
    /******************************************************/
    ps_view_ctxt->u1_qp = ps_slice->u1_slice_qp;
    ih264d_update_qp(ps_view_ctxt, 0);
    u4_mb_idx = ps_view_ctxt->u4_mb_idx;
    ps_parse_mb_data = ps_view_ctxt->ps_parse_mb_data;
    u4_num_mbs = u4_mb_idx;

    b_is_slice_end = false;
    b_tfr_n_mb = false;
    b_decode_nmb = false;
    i2_cur_mb_addr = ps_view_ctxt->u4_total_mbs_coded;
    u4_mb_skip_run = u4_remaining_mbs;

    while(!b_is_slice_end)
    {
        if(i2_cur_mb_addr > ps_view_ctxt->ps_cur_sps->u4_max_mb_addr)
        {
            break;
        }

        ps_cur_mb_info = ps_view_ctxt->ps_nmb_info + u4_num_mbs;
        ps_view_ctxt->u4_num_mbs_cur_nmb = u4_num_mbs;

        ps_cur_mb_info->u1_Mux = 0;
        ps_cur_mb_info->u1_end_of_slice = 0;

        ps_view_ctxt->u4_num_pmbair = u4_num_mbs;
        ps_cur_deblk_mb = ps_view_ctxt->ps_deblk_mbn + u4_num_mbs;

        ps_parse_mb_data->u1_num_part = 1;
        ps_parse_mb_data->u4_isI_mb = 0;

        /**************************************************************/
        /* Get the required information for decoding of MB            */
        /**************************************************************/
        /* mb_x, mb_y, neighbor availablity, */
        ih264d_get_mb_info_cavlc_nonmbaff(ps_view_ctxt, i2_cur_mb_addr, ps_cur_mb_info,
                                          u4_mb_skip_run);

        if(ps_view_ctxt->u4_app_disable_deblk_frm == 0)
        {
            ih264d_set_deblocking_parameters(ps_cur_deblk_mb, ps_slice,
                                             ps_view_ctxt->u1_mb_ngbr_availablity,
                                             ps_view_ctxt->u1_cur_mb_fld_dec_flag);
        }

        ps_view_ctxt->i1_prev_mb_qp_delta = 0;
        ps_view_ctxt->u1_sub_mb_num = 0;
        ps_cur_mb_info->u1_mb_type = MB_SKIP;
        ps_cur_mb_info->u1_mb_mc_mode = PRED_16x16;
        ps_cur_mb_info->u1_cbp = 0;

        /* Storing Skip partition info */
        ps_part_info = ps_view_ctxt->ps_part;
        ps_part_info->u1_is_direct = PART_DIRECT_16x16;
        ps_part_info->u1_sub_mb_num = 0;
        ps_view_ctxt->ps_part++;

        /* Update Nnzs */
        ih264d_update_nnz_for_skipmb(ps_view_ctxt, ps_cur_mb_info, CAVLC);

        ps_cur_mb_info->ps_curmb->u1_mb_type = u1_inter_mb_type;
        ps_cur_deblk_mb->u1_mb_type |= u1_deblk_mb_type;

        u4_mb_skip_run--;

        ps_cur_deblk_mb->u1_mb_qp = ps_view_ctxt->u1_qp;
        ps_cur_deblk_mb->u1_deblocking_mode = MB_DISABLE_FILTERING;

        i2_cur_mb_addr++;
        u4_num_mbs++;
        ps_parse_mb_data++;

        /****************************************************************/
        /* Check for End Of Row and other flags that determine when to  */
        /* do DMA setup for N/2-Mb, Decode for N-Mb, and Transfer for   */
        /* N-Mb                                                         */
        /****************************************************************/
        u4_num_mbs_next = ps_view_ctxt->ps_cur_sps->u2_frm_wd_in_mbs - 1 - ps_view_ctxt->u2_mbx;
        b_is_end_of_row = (0 == u4_num_mbs_next);
        b_is_slice_end = !u4_mb_skip_run;
        ps_cur_mb_info->u1_end_of_slice = !u4_mb_skip_run;
        b_tfr_n_mb =
            (u4_num_mbs == ps_view_ctxt->u4_recon_mb_grp) || b_is_end_of_row || b_is_slice_end;
        b_decode_nmb = b_tfr_n_mb || b_is_slice_end;

        if(b_decode_nmb)
        {
            ps_view_ctxt->pf_mvpred_ref_tfr_nby2mb(ps_view_ctxt, u4_mb_idx, u4_num_mbs);

            ps_parse_mb_data = ps_view_ctxt->ps_parse_mb_data;
            ps_view_ctxt->ps_part = ps_view_ctxt->ps_parse_part_params;

            if(ps_view_ctxt->u1_separate_parse)
            {
                ih264d_parse_tfr_nmb(ps_view_ctxt, u4_mb_idx, u4_num_mbs, u4_num_mbs_next,
                                     b_tfr_n_mb, b_is_end_of_row);

                ps_view_ctxt->ps_nmb_info += u4_num_mbs;
            }
            else
            {
                ih264d_decode_recon_tfr_nmb(ps_view_ctxt, u4_mb_idx, u4_num_mbs, u4_num_mbs_next,
                                            b_tfr_n_mb, b_is_end_of_row);
            }

            ps_view_ctxt->u4_total_mbs_coded += u4_num_mbs;

            if(b_tfr_n_mb)
            {
                u4_num_mbs = 0;
            }

            u4_mb_idx = u4_num_mbs;
            ps_view_ctxt->u4_mb_idx = u4_num_mbs;
        }
    }

    ps_view_ctxt->u4_num_mbs_cur_nmb = 0;
    ps_view_ctxt->i2_prev_slice_mbx = ps_view_ctxt->u2_mbx;
    ps_view_ctxt->i2_prev_slice_mby = ps_view_ctxt->u2_mby;

    if(ps_view_ctxt->u4_total_mbs_coded >= (ps_view_ctxt->ps_cur_sps->u4_max_mb_addr + 1))
    {
        ps_view_ctxt->u1_pic_decode_done = 1;
    }

    return 0;
}

static WORD32 imvcd_parse_pslice(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    UWORD8 u1_num_ref_idx_l0, u1_num_ref_idx_l1;
    WORD32 i4_error_code;
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    nalu_mvc_ext_t *ps_cur_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);
    dec_pic_params_t *ps_pps = ps_view_ctxt->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_view_ctxt->ps_cur_slice;
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    ps_view_ctxt->s_default_mv_pred = imvcd_get_default_mv_pred();

    i4_error_code = imvcd_set_ref_idx_override_flag(ps_view_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    if(ps_slice->u1_num_ref_idx_active_override_flag)
    {
        i4_error_code = imvcd_set_num_ref_idx_active(ps_view_ctxt, &u1_num_ref_idx_l0);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }
    }
    else
    {
        u1_num_ref_idx_l0 = ps_view_ctxt->ps_cur_pps->u1_num_ref_idx_lx_active[0];
    }

    u1_num_ref_idx_l1 = 0;

    ps_slice->u1_num_ref_idx_lx_active[0] = u1_num_ref_idx_l0;
    ps_slice->u1_num_ref_idx_lx_active[1] = u1_num_ref_idx_l1;
    ps_view_ctxt->u1_num_ref_idx_lx_active_prev = ps_slice->u1_num_ref_idx_lx_active[0];

    i4_error_code =
        imvcd_init_ref_pic_list(ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_nalu_mvc_ext,
                                ps_mvcd_ctxt->ps_cur_au, ps_mvcd_ctxt->u2_num_views_decoded);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_ref_pic_list_mod_data(ps_mvcd_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_dpb_reorder_ref_pic_list(
        ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_nalu_mvc_ext, ps_mvcd_ctxt->ps_cur_au,
        imvcd_get_cur_ref_pic_list_mod_data(ps_mvcd_ctxt), ps_mvcd_ctxt->u2_num_views_decoded);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    ps_view_ctxt->ps_ref_pic_buf_lx[0] = imvcd_dpb_get_view_ref_pic_list(
        ps_mvcd_ctxt->ps_dpb_mgr, ps_mvcd_ctxt->u2_num_views_decoded,
        ps_cur_nalu_mvc_ext->u2_view_id, 0);

    for(i = 0; i < u1_num_ref_idx_l0; i++)
    {
        if(NULL == ps_view_ctxt->ps_ref_pic_buf_lx[0][i]->pu1_buf1)
        {
            return ERROR_FEATURE_UNAVAIL;
        }
    }

    imvcd_set_view_buf_id_to_buf_map(ps_view_ctxt);

    imvcd_init_ref_idx_to_ref_buf_map(ps_mvcd_ctxt);

    if(ps_pps->u1_wted_pred_flag)
    {
        i4_error_code = ih264d_parse_pred_weight_table(ps_slice, ps_bitstrm);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        ih264d_form_pred_weight_matrix(ps_view_ctxt);
    }
    else
    {
        ps_view_ctxt->ps_cur_slice->u2_log2Y_crwd = 0;
    }

    ps_view_ctxt->pu4_wt_ofsts = ps_view_ctxt->pu4_wts_ofsts_mat;
    ps_view_ctxt->ps_parse_cur_slice->u2_log2Y_crwd = ps_view_ctxt->ps_cur_slice->u2_log2Y_crwd;

    if(ps_slice->u1_nal_ref_idc != 0)
    {
        if(!ps_view_ctxt->ps_dpb_cmds->u1_dpb_commands_read)
        {
            WORD32 i4_bit_offset = ih264d_read_mmco_commands(ps_view_ctxt);

            if(i4_bit_offset < 0)
            {
                return ERROR_DBP_MANAGER_T;
            }

            ps_view_ctxt->u4_bitoffset = i4_bit_offset;
        }
        else
        {
            ps_bitstrm->u4_ofst += ps_view_ctxt->u4_bitoffset;
        }
    }

    if(ps_pps->u1_entropy_coding_mode == CABAC)
    {
        i4_error_code = imvcd_set_cabac_init_idc(ps_view_ctxt);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }
    }

    i4_error_code = imvcd_set_slice_qp(ps_view_ctxt);

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_slice_deblk_params(ps_view_ctxt);

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    ps_view_ctxt->u1_slice_header_done = 1;

    if(ps_pps->u1_entropy_coding_mode)
    {
        ps_view_ctxt->pf_parse_inter_slice = ih264d_parse_inter_slice_data_cabac;
        ps_view_ctxt->pf_parse_inter_mb = ih264d_parse_pmb_cabac;
        ps_view_ctxt->pf_get_mb_info = ih264d_get_mb_info_cabac_nonmbaff;

        ih264d_init_cabac_contexts(P_SLICE, ps_view_ctxt);
    }
    else
    {
        ps_view_ctxt->pf_parse_inter_slice = ih264d_parse_inter_slice_data_cavlc;
        ps_view_ctxt->pf_parse_inter_mb = ih264d_parse_pmb_cavlc;
        ps_view_ctxt->pf_get_mb_info = ih264d_get_mb_info_cavlc_nonmbaff;
    }

    ps_view_ctxt->pf_mvpred_ref_tfr_nby2mb = ih264d_mv_pred_ref_tfr_nby2_pmb;

    ps_view_ctxt->u1_B = 0;

    i4_error_code =
        ps_view_ctxt->pf_parse_inter_slice(ps_view_ctxt, ps_slice, ps_slice->u2_first_mb_in_slice);

    return i4_error_code;
}

static WORD32 imvcd_parse_bslice(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    UWORD8 u1_num_ref_idx_l0, u1_num_ref_idx_l1;
    WORD32 i4_error_code;
    WORD32 i, j;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    nalu_mvc_ext_t *ps_cur_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);
    dec_pic_params_t *ps_pps = ps_view_ctxt->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_view_ctxt->ps_cur_slice;
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    ps_view_ctxt->s_default_mv_pred = imvcd_get_default_mv_pred();

    i4_error_code = imvcd_set_ref_idx_override_flag(ps_view_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    if(ps_slice->u1_num_ref_idx_active_override_flag)
    {
        i4_error_code = imvcd_set_num_ref_idx_active(ps_view_ctxt, &u1_num_ref_idx_l0);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        i4_error_code = imvcd_set_num_ref_idx_active(ps_view_ctxt, &u1_num_ref_idx_l1);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }
    }
    else
    {
        u1_num_ref_idx_l0 = ps_view_ctxt->ps_cur_pps->u1_num_ref_idx_lx_active[0];
        u1_num_ref_idx_l1 = ps_view_ctxt->ps_cur_pps->u1_num_ref_idx_lx_active[1];
    }

    if((0 == u1_num_ref_idx_l0) || (0 == u1_num_ref_idx_l1))
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ps_slice->u1_num_ref_idx_lx_active[0] = u1_num_ref_idx_l0;
    ps_slice->u1_num_ref_idx_lx_active[1] = u1_num_ref_idx_l1;
    ps_view_ctxt->u1_num_ref_idx_lx_active_prev =
        ps_view_ctxt->ps_cur_slice->u1_num_ref_idx_lx_active[0];

    i4_error_code =
        imvcd_init_ref_pic_list(ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_nalu_mvc_ext,
                                ps_mvcd_ctxt->ps_cur_au, ps_mvcd_ctxt->u2_num_views_decoded);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_ref_pic_list_mod_data(ps_mvcd_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_dpb_reorder_ref_pic_list(
        ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_nalu_mvc_ext, ps_mvcd_ctxt->ps_cur_au,
        imvcd_get_cur_ref_pic_list_mod_data(ps_mvcd_ctxt), ps_mvcd_ctxt->u2_num_views_decoded);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    for(i = 0; i < 2; i++)
    {
        ps_view_ctxt->ps_ref_pic_buf_lx[i] = imvcd_dpb_get_view_ref_pic_list(
            ps_mvcd_ctxt->ps_dpb_mgr, ps_mvcd_ctxt->u2_num_views_decoded,
            ps_cur_nalu_mvc_ext->u2_view_id, i);

        for(j = 0; j < ps_slice->u1_num_ref_idx_lx_active[i]; j++)
        {
            if(NULL == ps_view_ctxt->ps_ref_pic_buf_lx[i][j]->pu1_buf1)
            {
                return ERROR_FEATURE_UNAVAIL;
            }
        }
    }

    imvcd_set_view_buf_id_to_buf_map(ps_view_ctxt);

    imvcd_init_ref_idx_to_ref_buf_map(ps_mvcd_ctxt);

    if(ps_pps->u1_wted_bipred_idc == 1)
    {
        i4_error_code = ih264d_parse_pred_weight_table(ps_slice, ps_bitstrm);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        ih264d_form_pred_weight_matrix(ps_view_ctxt);

        ps_view_ctxt->pu4_wt_ofsts = ps_view_ctxt->pu4_wts_ofsts_mat;
    }
    else if(ps_pps->u1_wted_bipred_idc == 2)
    {
        /* Implicit Weighted prediction */
        ps_slice->u2_log2Y_crwd = 0x0505;
        ps_view_ctxt->pu4_wt_ofsts = ps_view_ctxt->pu4_wts_ofsts_mat;

        ih264d_get_implicit_weights(ps_view_ctxt);
    }
    else
    {
        ps_view_ctxt->ps_cur_slice->u2_log2Y_crwd = 0;
    }

    ps_view_ctxt->ps_parse_cur_slice->u2_log2Y_crwd = ps_view_ctxt->ps_cur_slice->u2_log2Y_crwd;

    if(ps_slice->u1_nal_ref_idc != 0)
    {
        if(!ps_view_ctxt->ps_dpb_cmds->u1_dpb_commands_read)
        {
            WORD32 i4_bit_offset = ih264d_read_mmco_commands(ps_view_ctxt);

            if(i4_bit_offset < 0)
            {
                return ERROR_DBP_MANAGER_T;
            }

            ps_view_ctxt->u4_bitoffset = i4_bit_offset;
        }
        else
        {
            ps_bitstrm->u4_ofst += ps_view_ctxt->u4_bitoffset;
        }
    }

    if(ps_pps->u1_entropy_coding_mode == CABAC)
    {
        i4_error_code = imvcd_set_cabac_init_idc(ps_view_ctxt);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }
    }

    i4_error_code = imvcd_set_slice_qp(ps_view_ctxt);

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_slice_deblk_params(ps_view_ctxt);

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    ps_view_ctxt->u1_slice_header_done = 1;

    if(ps_pps->u1_entropy_coding_mode)
    {
        ps_view_ctxt->pf_parse_inter_slice = ih264d_parse_inter_slice_data_cabac;
        ps_view_ctxt->pf_parse_inter_mb = ih264d_parse_bmb_cabac;
        ps_view_ctxt->pf_get_mb_info = ih264d_get_mb_info_cabac_nonmbaff;

        ih264d_init_cabac_contexts(B_SLICE, ps_view_ctxt);
    }
    else
    {
        ps_view_ctxt->pf_parse_inter_slice = ih264d_parse_inter_slice_data_cavlc;
        ps_view_ctxt->pf_parse_inter_mb = ih264d_parse_bmb_cavlc;
        ps_view_ctxt->pf_get_mb_info = ih264d_get_mb_info_cavlc_nonmbaff;
    }

    i4_error_code = ih264d_cal_col_pic(ps_view_ctxt);

    if(i4_error_code != OK)
    {
        return i4_error_code;
    }

    ps_view_ctxt->u1_B = 1;

    ps_view_ctxt->pf_mvpred_ref_tfr_nby2mb = ih264d_mv_pred_ref_tfr_nby2_bmb;

    i4_error_code =
        ps_view_ctxt->pf_parse_inter_slice(ps_view_ctxt, ps_slice, ps_slice->u2_first_mb_in_slice);

    return i4_error_code;
}

static WORD32 imvcd_parse_islice(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i4_error_code;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    i4_error_code =
        ih264d_parse_islice(ps_view_ctxt, ps_view_ctxt->ps_cur_slice->u2_first_mb_in_slice);

    return i4_error_code;
}

static WORD32 imvcd_finish_slice_decode(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_err_status_t *ps_err = ps_view_ctxt->ps_dec_err_status;

    UWORD16 u2_view_order_id = ps_mvcd_ctxt->u2_num_views_decoded;
    UWORD16 u2_num_views = ps_mvcd_ctxt->u2_num_views;

    imvcd_dpb_reset_ivp_ctxt(ps_mvcd_ctxt->ps_dpb_mgr);

    /* End of Picture detection */
    if(ps_view_ctxt->u4_total_mbs_coded >= (ps_view_ctxt->ps_cur_sps->u4_max_mb_addr + 1))
    {
        ps_view_ctxt->u1_pic_decode_done = 1;
    }
    else
    {
        imvcd_corrupted_slice_handler(ps_mvcd_ctxt);

        return ERROR_CORRUPTED_SLICE;
    }

    if((ps_view_ctxt->u1_slice_header_done) && (u2_view_order_id == (u2_num_views - 1)))
    {
        ps_view_ctxt->u1_first_slice_in_stream = 0;
    }

    if((ps_mvcd_ctxt->au1_nal_ref_idc[u2_view_order_id] != 0) && (0 == u2_view_order_id))
    {
        if(!ps_view_ctxt->ps_dpb_cmds->u1_dpb_commands_read)
        {
            ps_view_ctxt->ps_dpb_cmds[0] = ps_view_ctxt->s_dpb_cmds_scratch;
        }
    }

    /* storing last Mb X and MbY of the slice */
    ps_view_ctxt->i2_prev_slice_mbx = ps_view_ctxt->u2_mbx;
    ps_view_ctxt->i2_prev_slice_mby = ps_view_ctxt->u2_mby;

    if((ps_err->u1_err_flag & REJECT_PB_PICS) && (ps_err->u1_cur_pic_type == PIC_TYPE_I))
    {
        ps_err->u1_err_flag = ACCEPT_ALL_PICS;
    }

    /* Accounting for idiocy in 'ih264d_parse_sps' */
    if(u2_view_order_id > 0)
    {
        for(i = 0; i < MAX_NUM_SEQ_PARAMS; i++)
        {
            if(ps_view_ctxt->ps_sps->u1_is_valid)
            {
                ps_view_ctxt->ps_cur_sps = ps_view_ctxt->ps_sps;

                break;
            }
        }
    }

    return OK;
}

WORD32 imvcd_parse_decode_slice(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_pic_params_t *ps_pps;
    dec_seq_params_t *ps_sps;
    dec_slice_params_t *ps_cur_slice;

    WORD32 i4_error_code;
    UWORD8 u1_pps_id;
    UWORD8 u1_pic_order_cnt_type;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    nalu_mvc_ext_t *ps_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;
    pocstruct_t s_tmp_poc = {0};
    dec_err_status_t *ps_err = ps_view_ctxt->ps_dec_err_status;

    WORD32 ai4_delta_poc[2] = {0};
    WORD32 i4_poc = 0;
    UWORD32 u4_idr_pic_id = 0;
    UWORD16 u2_view_id = ps_nalu_mvc_ext->u2_view_id;
    UWORD16 u2_view_order_id = ps_mvcd_ctxt->u2_num_views_decoded;
    bool b_is_idr_slice = imvcd_is_idr_au(ps_mvcd_ctxt);
    UWORD16 u2_num_views = ps_mvcd_ctxt->u2_num_views;
    UWORD8 u1_redundant_pic_cnt = 0;
    const UWORD8 u1_field_pic_flag = 0;
    const UWORD8 u1_bottom_field_flag = 0;

    ps_view_ctxt->ps_cur_slice = ps_cur_slice = &ps_mvcd_ctxt->as_slices[u2_view_id];
    ps_view_ctxt->ps_dpb_cmds->u1_dpb_commands_read_slc = 0;

    ps_cur_slice->u1_nal_unit_type = ps_mvcd_ctxt->ae_nalu_id[u2_view_order_id];
    ps_cur_slice->u1_nal_ref_idc = ps_mvcd_ctxt->au1_nal_ref_idc[u2_view_order_id];

    i4_error_code = imvcd_set_first_mb_in_slice(ps_view_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_slice_type(ps_view_ctxt);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    i4_error_code = imvcd_set_cur_pps(ps_view_ctxt, &u1_pps_id);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    ps_pps = ps_view_ctxt->ps_cur_pps;
    ps_sps = ps_view_ctxt->ps_cur_sps;

    i4_error_code = imvcd_set_frame_num(ps_view_ctxt, ps_sps->u1_bits_in_frm_num);

    if(OK != i4_error_code)
    {
        return i4_error_code;
    }

    if(!ps_view_ctxt->u1_first_slice_in_stream && ps_view_ctxt->u4_first_slice_in_pic)
    {
        ps_view_ctxt->u2_mbx = 0xffff;
        ps_view_ctxt->u2_mby = 0;
        ps_view_ctxt->u4_total_mbs_coded = 0;

        if(0 == u2_view_order_id)
        {
            if(b_is_idr_slice || ps_cur_slice->u1_mmco_equalto5)
            {
                ps_view_ctxt->u2_prev_ref_frame_num = 0;
            }

            if(ps_view_ctxt->ps_cur_sps->u1_gaps_in_frame_num_value_allowed_flag)
            {
                i4_error_code = imvcd_decode_gaps_in_frame_num(ps_mvcd_ctxt);

                if(OK != i4_error_code)
                {
                    return i4_error_code;
                }
            }

            if(!b_is_idr_slice && ps_cur_slice->u1_nal_ref_idc)
            {
                ps_view_ctxt->u2_prev_ref_frame_num = ps_cur_slice->u2_frame_num;
            }

            imvcd_pocstruct_init(ps_view_ctxt);
        }
    }

    if(b_is_idr_slice)
    {
        i4_error_code = imvcd_set_idr_pic_id(ps_view_ctxt, &u4_idr_pic_id);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        /* 'ih264d_read_mmco_commands' asssumes AVC semantics */
        ps_view_ctxt->u1_nal_unit_type = SLICE_IDR;
    }
    else
    {
        ps_view_ctxt->u1_nal_unit_type = SLICE_NON_IDR;
    }

    u1_pic_order_cnt_type = ps_sps->u1_pic_order_cnt_type;

    if(0 == u1_pic_order_cnt_type)
    {
        i4_error_code = imvcd_set_poc_lsb(ps_view_ctxt, &s_tmp_poc.i4_pic_order_cnt_lsb,
                                          ps_sps->i4_max_pic_order_cntLsb,
                                          ps_sps->u1_log2_max_pic_order_cnt_lsb_minus);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        if(ps_pps->u1_pic_order_present_flag)
        {
            i4_error_code =
                imvcd_set_delta_poc(ps_view_ctxt, &s_tmp_poc.i4_delta_pic_order_cnt_bottom);

            if(OK != i4_error_code)
            {
                return i4_error_code;
            }
        }
    }

    if((1 == u1_pic_order_cnt_type) && !ps_sps->u1_delta_pic_order_always_zero_flag)
    {
        i4_error_code = imvcd_set_delta_poc(ps_view_ctxt, &s_tmp_poc.i4_delta_pic_order_cnt[0]);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        if(ps_pps->u1_pic_order_present_flag)
        {
            i4_error_code = imvcd_set_delta_poc(ps_view_ctxt, &s_tmp_poc.i4_delta_pic_order_cnt[1]);

            if(OK != i4_error_code)
            {
                return i4_error_code;
            }
        }
    }

    if(ps_pps->u1_redundant_pic_cnt_present_flag)
    {
        i4_error_code = imvcd_set_redundant_pic_cnt(ps_view_ctxt, &u1_redundant_pic_cnt);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }
    }

    ps_view_ctxt->ps_dec_err_status->u1_err_flag &= MASK_REJECT_CUR_PIC;

    ps_view_ctxt->u1_slice_header_done = 0;

    if(ps_view_ctxt->u4_first_slice_in_pic)
    {
        i4_error_code = ih264d_decode_pic_order_cnt(
            b_is_idr_slice, ps_cur_slice->u2_frame_num, &ps_view_ctxt->s_prev_pic_poc, &s_tmp_poc,
            ps_cur_slice, ps_pps, ps_mvcd_ctxt->au1_nal_ref_idc[u2_view_order_id],
            u1_bottom_field_flag, u1_field_pic_flag, &i4_poc);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }

        /* Display seq no calculations */
        if(i4_poc >= ps_view_ctxt->i4_max_poc)
        {
            ps_view_ctxt->i4_max_poc = i4_poc;
        }

        /* IDR Picture or POC wrap around */
        if(i4_poc == 0)
        {
            imvcd_modulate_max_disp_seq(ps_view_ctxt);
        }
    }

    if((0 == i4_poc) && (ps_mvcd_ctxt->ae_nalu_id[u2_view_order_id] == SLICE_IDR) &&
       (ps_cur_slice->u1_slice_type != ISLICE))
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    /*--------------------------------------------------------------------*/
    /* Copy the values read from the bitstream to the slice header and then*/
    /* If the slice is first slice in picture, then do Start of Picture   */
    /* processing.                                                        */
    /*--------------------------------------------------------------------*/
    ps_cur_slice->i4_delta_pic_order_cnt[0] = ai4_delta_poc[0];
    ps_cur_slice->i4_delta_pic_order_cnt[1] = ai4_delta_poc[1];
    ps_cur_slice->u4_idr_pic_id = u4_idr_pic_id;
    ps_cur_slice->u1_field_pic_flag = u1_field_pic_flag;
    ps_cur_slice->u1_bottom_field_flag = u1_bottom_field_flag;
    ps_cur_slice->i4_pic_order_cnt_lsb = s_tmp_poc.i4_pic_order_cnt_lsb;
    ps_cur_slice->u1_redundant_pic_cnt = u1_redundant_pic_cnt;
    ps_cur_slice->u1_pic_order_cnt_type = u1_pic_order_cnt_type;
    ps_cur_slice->i4_poc = i4_poc;

    ps_cur_slice->u1_direct_8x8_inference_flag = ps_sps->u1_direct_8x8_inference_flag;

    if(IV_SUCCESS != imvcd_view_error_checks(ps_mvcd_ctxt))
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    if(ps_cur_slice->u1_slice_type == B_SLICE)
    {
        i4_error_code = imvcd_set_direct_spatial_mv_pred_flag(ps_view_ctxt);

        if(OK != i4_error_code)
        {
            return i4_error_code;
        }

        if(ps_cur_slice->u1_direct_spatial_mv_pred_flag)
        {
            ps_cur_slice->pf_decodeDirect = ih264d_decode_spatial_direct;
        }
        else
        {
            ps_cur_slice->pf_decodeDirect = ih264d_decode_temporal_direct;
        }

        ps_view_ctxt->pf_mvpred = ih264d_mvpred_nonmbaffB;
    }
    else
    {
        ps_view_ctxt->pf_mvpred = ih264d_mvpred_nonmbaff;
    }

    if(ps_view_ctxt->u4_first_slice_in_pic)
    {
        if(0 == ps_cur_slice->u2_first_mb_in_slice)
        {
            i4_error_code = imvcd_pic_init(ps_mvcd_ctxt, &s_tmp_poc, i4_poc, b_is_idr_slice);

            if(i4_error_code != OK)
            {
                return i4_error_code;
            }
        }
        else
        {
            return ERROR_INV_SLICE_HDR_T;
        }

        ps_view_ctxt->u4_output_present = 0;

        if(u2_view_order_id == (u2_num_views - 1))
        {
            if(IV_SUCCESS == imvcd_get_next_display_au_buf(ps_mvcd_ctxt))
            {
                ps_view_ctxt->u4_output_present = 1;
            }
        }

        if(!imvcd_dpb_is_diff_poc_valid(ps_mvcd_ctxt->ps_dpb_mgr, ps_cur_slice->i4_poc))
        {
            return ERROR_INV_SLICE_HDR_T;
        }

        if(ps_view_ctxt->u1_separate_parse == 1)
        {
            if(!ps_view_ctxt->u4_dec_thread_created)
            {
                ithread_create(ps_view_ctxt->pv_dec_thread_handle, NULL,
                               ih264d_decode_picture_thread, ps_view_ctxt);

                ps_view_ctxt->u4_dec_thread_created = 1;
            }

            if((3 == ps_view_ctxt->u4_num_cores) &&
               (!ps_view_ctxt->u4_app_disable_deblk_frm ||
                ps_view_ctxt->i1_recon_in_thread3_flag) &&
               !ps_view_ctxt->u4_bs_deblk_thread_created)
            {
                ps_view_ctxt->u4_start_recon_deblk = 0;

                ithread_create(ps_view_ctxt->pv_bs_deblk_thread_handle, NULL,
                               ih264d_recon_deblk_thread, ps_view_ctxt);

                ps_view_ctxt->u4_bs_deblk_thread_created = 1;
            }
        }
    }

    if((ps_cur_slice->u1_slice_type != B_SLICE) &&
       (ps_view_ctxt->ps_cur_pps->u1_wted_pred_flag == 0))
    {
        ps_view_ctxt->p_form_mb_part_info = ih264d_form_mb_part_info_bp;
        ps_view_ctxt->p_motion_compensate = ih264d_motion_compensate_bp;
    }
    else
    {
        ps_view_ctxt->p_form_mb_part_info = ih264d_form_mb_part_info_mp;
        ps_view_ctxt->p_motion_compensate = ih264d_motion_compensate_mp;
    }

    if(ps_err->u4_frm_sei_sync == ps_cur_slice->u2_frame_num)
    {
        ps_err->u1_err_flag = ACCEPT_ALL_PICS;
        ps_err->u4_frm_sei_sync = SYNC_FRM_DEFAULT;
    }

    ps_err->u4_cur_frm = ps_cur_slice->u2_frame_num;

    ps_view_ctxt->i4_submb_ofst = -SUB_BLK_SIZE;

    ps_view_ctxt->u4_cur_mb_addr = 0;
    ps_view_ctxt->ps_deblk_mbn = ps_view_ctxt->ps_deblk_pic;
    ps_view_ctxt->ps_mv_cur = ps_mvcd_ctxt->ps_cur_au->ps_au_mv_data->aps_mvs[u2_view_id];

    ps_view_ctxt->s_tran_addrecon.pu1_dest_y =
        ps_mvcd_ctxt->ps_cur_au->as_view_buffers[u2_view_id].as_component_bufs[Y].pv_data;
    ps_view_ctxt->s_tran_addrecon.pu1_dest_u =
        ps_mvcd_ctxt->ps_cur_au->as_view_buffers[u2_view_id].as_component_bufs[UV].pv_data;
    ps_view_ctxt->s_tran_addrecon.pu1_dest_v = NULL;

    ps_view_ctxt->s_tran_addrecon.pu1_mb_y =
        ps_mvcd_ctxt->ps_cur_au->as_view_buffers[u2_view_id].as_component_bufs[Y].pv_data;
    ps_view_ctxt->s_tran_addrecon.pu1_mb_u =
        ps_mvcd_ctxt->ps_cur_au->as_view_buffers[u2_view_id].as_component_bufs[UV].pv_data;
    ps_view_ctxt->s_tran_addrecon.pu1_mb_v = NULL;

    ps_view_ctxt->ps_part = ps_view_ctxt->ps_parse_part_params;

    ps_view_ctxt->u2_mbx = (MOD(ps_cur_slice->u2_first_mb_in_slice - 1, ps_sps->u2_frm_wd_in_mbs));
    ps_view_ctxt->u2_mby = (DIV(ps_cur_slice->u2_first_mb_in_slice - 1, ps_sps->u2_frm_wd_in_mbs));
    ps_view_ctxt->i2_prev_slice_mbx = ps_view_ctxt->u2_mbx;
    ps_view_ctxt->i2_prev_slice_mby = ps_view_ctxt->u2_mby;

    /* RBSP stop bit is used for CABAC decoding*/
    ps_bitstrm->u4_max_ofst += ps_view_ctxt->ps_cur_pps->u1_entropy_coding_mode;

    ps_view_ctxt->u1_B = (ps_cur_slice->u1_slice_type == B_SLICE);
    ps_view_ctxt->u4_next_mb_skip = 0;

    ps_view_ctxt->ps_parse_cur_slice->u4_first_mb_in_slice = ps_cur_slice->u2_first_mb_in_slice;
    ps_view_ctxt->ps_parse_cur_slice->slice_type = ps_cur_slice->u1_slice_type;

    ps_view_ctxt->u4_start_recon_deblk = 1;

    ps_view_ctxt->ps_parse_cur_slice->ppv_map_ref_idx_to_poc =
        ps_view_ctxt->pv_map_ref_idx_to_poc_buf;

    if(ps_view_ctxt->u1_separate_parse)
    {
        ps_view_ctxt->ps_parse_cur_slice->pv_tu_coeff_data_start =
            ps_view_ctxt->pv_parse_tu_coeff_data;
    }
    else
    {
        ps_view_ctxt->pv_proc_tu_coeff_data = ps_view_ctxt->pv_parse_tu_coeff_data;
    }

    if(0 == u2_view_order_id)
    {
        i4_error_code = imvcd_dpb_st_lt_deduplicator(ps_mvcd_ctxt->ps_dpb_mgr);

        if(i4_error_code < 0)
        {
            i4_error_code = ERROR_DBP_MANAGER_T;
        }
    }

    if(ps_cur_slice->u1_slice_type == I_SLICE)
    {
        ps_mvcd_ctxt->ps_cur_au->au4_pack_slc_typ[u2_view_order_id] |= I_SLC_BIT;

        i4_error_code = imvcd_parse_islice(ps_mvcd_ctxt);

        ps_view_ctxt->u1_pr_sl_type = ps_cur_slice->u1_slice_type;

        if(ps_view_ctxt->i4_pic_type != B_SLICE && ps_view_ctxt->i4_pic_type != P_SLICE)
        {
            ps_view_ctxt->i4_pic_type = I_SLICE;
        }
    }
    else if(ps_cur_slice->u1_slice_type == P_SLICE)
    {
        ps_mvcd_ctxt->ps_cur_au->au4_pack_slc_typ[u2_view_order_id] |= P_SLC_BIT;

        i4_error_code = imvcd_parse_pslice(ps_mvcd_ctxt);

        ps_view_ctxt->u1_pr_sl_type = ps_cur_slice->u1_slice_type;

        if(ps_view_ctxt->i4_pic_type != B_SLICE)
        {
            ps_view_ctxt->i4_pic_type = P_SLICE;
        }
    }
    else if(ps_cur_slice->u1_slice_type == B_SLICE)
    {
        ps_mvcd_ctxt->ps_cur_au->au4_pack_slc_typ[u2_view_order_id] |= B_SLC_BIT;

        i4_error_code = imvcd_parse_bslice(ps_mvcd_ctxt);

        ps_view_ctxt->u1_pr_sl_type = ps_cur_slice->u1_slice_type;

        ps_view_ctxt->i4_pic_type = B_SLICE;
    }
    else
    {
        i4_error_code = ERROR_INV_SLC_TYPE_T;
    }

    i4_error_code = imvcd_finish_slice_decode(ps_mvcd_ctxt);

    return i4_error_code;
}
