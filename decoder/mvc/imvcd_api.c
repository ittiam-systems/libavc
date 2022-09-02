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
/*  File Name         : imvcd_api.c                                          */
/*                                                                           */
/*  Description       : Has all MVC API functions                            */
/*                                                                           */
/*                                                                           */
/*  List of Functions :                                                      */
/*                                                                           */
/*****************************************************************************/
#include <string.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"
#include "imvcd.h"
#include "ih264_debug.h"
#include "ih264_disp_mgr.h"
#include "ih264_error.h"
#include "ih264_buf_mgr.h"
#include "ih264_platform_macros.h"
#include "ih264d_inter_pred.h"
#include "ih264d_structs.h"
#include "ih264d_deblocking.h"
#include "ih264d_error_handler.h"
#include "ih264d_function_selector.h"
#include "ih264d_nal.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_parse_headers.h"
#include "ih264d_tables.h"
#include "ih264d_thread_compute_bs.h"
#include "ih264d_utils.h"
#include "ih264d_api_utils.h"
#include "ithread.h"
#include "imvcd_api_utils.h"
#include "imvcd_dpb_manager.h"
#include "imvcd_error_handler.h"
#include "imvcd_nalu_parser.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

static void imvcd_free_static_bufs(iv_obj_t *ps_dec_hdl)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt;
    dec_struct_t *ps_view_ctxt;

    FT_ALIGNED_FREE *pf_aligned_free;

    void *pv_mem_ctxt;

    if(!ps_dec_hdl)
    {
        return;
    }

    ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;

    if(!ps_mvcd_ctxt)
    {
        return;
    }

    ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    pf_aligned_free = ps_view_ctxt->pf_aligned_free;
    pv_mem_ctxt = ps_view_ctxt->pv_mem_ctxt;

    imvcd_free_dynamic_bufs(ps_mvcd_ctxt);

    imvcd_bitsteam_buf_free(ps_view_ctxt);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_left_mvpred_addr);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu4_wts_ofsts_mat);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu4_mbaff_wt_mat);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_init_dpb_base);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_temp_mc_buffer);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pi2_pred1);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_ref_buff_base);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_left_mb_ctxt_info);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->p_cabac_ctxt_table_t);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ppv_map_ref_idx_to_poc_base);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_bits_buf_static);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pv_scratch_sps_pps);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_bitstrm);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_dpb_cmds);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_sei_parse);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_sei);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_dec_err_status);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_pred);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pv_bs_deblk_thread_handle);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pv_dec_thread_handle);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_pps);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_sps);

    ih264_buf_mgr_free(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_mem);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_mem);

    ih264_buf_mgr_free(ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_mem);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_mem);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_mvcd_ctxt->ps_dpb_mgr);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_mvcd_ctxt);

    if(ps_dec_hdl)
    {
        pf_aligned_free(pv_mem_ctxt, ps_dec_hdl);
    }
}

static IV_API_CALL_STATUS_T imvcd_view_ctxt_init(imvcd_create_ip_t *ps_ip,
                                                 dec_struct_t *ps_view_ctxt)
{
    pocstruct_t *ps_prev_poc, *ps_cur_poc;

    WORD32 i4_mem_size;
    void *pv_buf;

    FT_ALIGNED_ALLOC *pf_aligned_alloc = ps_ip->s_ivd_ip.pf_aligned_alloc;

    void *pv_mem_ctxt = ps_ip->s_ivd_ip.pv_mem_ctxt;
    const WORD32 i4_default_alignment = 128;

    ps_view_ctxt->u4_share_disp_buf = 0;
    ps_view_ctxt->u1_chroma_format = ps_ip->s_ivd_ip.e_output_format;

    ps_view_ctxt->pf_aligned_alloc = pf_aligned_alloc;
    ps_view_ctxt->pf_aligned_free = ps_ip->s_ivd_ip.pf_aligned_free;
    ps_view_ctxt->pv_mem_ctxt = ps_ip->s_ivd_ip.pv_mem_ctxt;

    i4_mem_size = ((sizeof(dec_seq_params_t)) * MAX_NUM_SEQ_PARAMS);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_sps = pv_buf;

    i4_mem_size = (sizeof(dec_pic_params_t)) * MAX_NUM_PIC_PARAMS;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_pps = pv_buf;

    i4_mem_size = ithread_get_handle_size();
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pv_dec_thread_handle = pv_buf;

    i4_mem_size = ithread_get_handle_size();
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pv_bs_deblk_thread_handle = pv_buf;

    i4_mem_size = sizeof(pred_info_t) * 2 * 32;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_pred = pv_buf;

    ps_view_ctxt->pv_disp_buf_mgr = NULL;

    ps_view_ctxt->pv_pic_buf_mgr = NULL;

    ps_view_ctxt->ps_pic_buf_base = NULL;

    i4_mem_size = sizeof(dec_err_status_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_dec_err_status = (dec_err_status_t *) pv_buf;

    i4_mem_size = sizeof(sei);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_sei = (sei *) pv_buf;

    i4_mem_size = sizeof(sei);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_sei_parse = (sei *) pv_buf;

    i4_mem_size = sizeof(dpb_commands_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_dpb_cmds = (dpb_commands_t *) pv_buf;

    i4_mem_size = sizeof(dec_bit_stream_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_bitstrm = (dec_bit_stream_t *) pv_buf;

    i4_mem_size = MAX(sizeof(dec_seq_params_t), sizeof(dec_pic_params_t));
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pv_scratch_sps_pps = pv_buf;

    ps_view_ctxt->u4_static_bits_buf_size = MIN_BITSTREAMS_BUF_SIZE;
    pv_buf =
        pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, ps_view_ctxt->u4_static_bits_buf_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, ps_view_ctxt->u4_static_bits_buf_size);
    ps_view_ctxt->pu1_bits_buf_static = pv_buf;

    i4_mem_size = (TOTAL_LIST_ENTRIES + PAD_MAP_IDX_POC) * sizeof(void *);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->ppv_map_ref_idx_to_poc_base = pv_buf;
    ps_view_ctxt->ppv_map_ref_idx_to_poc =
        ps_view_ctxt->ppv_map_ref_idx_to_poc_base + OFFSET_MAP_IDX_POC;
    memset(ps_view_ctxt->ppv_map_ref_idx_to_poc_base, 0, i4_mem_size);

    i4_mem_size = (sizeof(bin_ctxt_model_t) * NUM_CABAC_CTXTS);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->p_cabac_ctxt_table_t = pv_buf;

    i4_mem_size = sizeof(ctxt_inc_mb_info_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_left_mb_ctxt_info = pv_buf;

    i4_mem_size = MAX_REF_BUF_SIZE * 2;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu1_ref_buff_base = pv_buf;
    ps_view_ctxt->pu1_ref_buff = ps_view_ctxt->pu1_ref_buff_base + MAX_REF_BUF_SIZE;

    i4_mem_size = sizeof(WORD16) * PRED_BUFFER_WIDTH * PRED_BUFFER_HEIGHT * 2;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pi2_pred1 = pv_buf;

    i4_mem_size = sizeof(UWORD8) * (MB_LUM_SIZE);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu1_temp_mc_buffer = pv_buf;

    i4_mem_size = (sizeof(UWORD32) * 2 * 3 * ((MAX_FRAMES << 1) * (MAX_FRAMES << 1)) * 2);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu4_mbaff_wt_mat = pv_buf;

    i4_mem_size = sizeof(UWORD32) * 2 * 3 * ((MAX_FRAMES << 1) * (MAX_FRAMES << 1));
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu4_wts_ofsts_mat = pv_buf;

    i4_mem_size = (sizeof(neighbouradd_t) << 2);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_left_mvpred_addr = pv_buf;

    ps_view_ctxt->pv_mv_buf_mgr = NULL;

    ps_view_ctxt->ps_col_mv_base = NULL;

    ps_view_ctxt->init_done = 0;
    ps_view_ctxt->u4_num_cores = 1;
    ps_view_ctxt->u2_pic_ht = ps_view_ctxt->u2_pic_wd = 0;
    ps_view_ctxt->u1_separate_parse = DEFAULT_SEPARATE_PARSE;
    ps_view_ctxt->u4_app_disable_deblk_frm = 0;
    ps_view_ctxt->i4_degrade_type = 0;
    ps_view_ctxt->i4_degrade_pics = 0;

    memset(ps_view_ctxt->ps_pps, 0, ((sizeof(dec_pic_params_t)) * MAX_NUM_PIC_PARAMS));
    memset(ps_view_ctxt->ps_sps, 0, ((sizeof(dec_seq_params_t)) * MAX_NUM_SEQ_PARAMS));

    ps_view_ctxt->p_DeblockPicture[0] = ih264d_deblock_picture_non_mbaff;
    ps_view_ctxt->p_DeblockPicture[1] = ih264d_deblock_picture_mbaff;
    ps_view_ctxt->s_cab_dec_env.pv_codec_handle = ps_view_ctxt;
    ps_view_ctxt->u4_num_fld_in_frm = 0;
    ps_view_ctxt->ps_sei->u1_is_valid = 0;
    ps_view_ctxt->ps_cur_pps = NULL;
    ps_view_ctxt->ps_cur_sps = NULL;
    ps_view_ctxt->ps_cur_slice = NULL;
    ps_view_ctxt->u1_init_dec_flag = 0;
    ps_view_ctxt->u1_first_slice_in_stream = 1;
    ps_view_ctxt->u1_last_pic_not_decoded = 0;
    ps_view_ctxt->u4_app_disp_width = 0;
    ps_view_ctxt->i4_header_decoded = 0;
    ps_view_ctxt->u4_total_frames_decoded = 0;
    ps_view_ctxt->i4_error_code = 0;
    ps_view_ctxt->i4_content_type = IV_CONTENTTYPE_NA;
    ps_view_ctxt->ps_dec_err_status->u1_err_flag = ACCEPT_ALL_PICS;
    ps_view_ctxt->ps_dec_err_status->u1_cur_pic_type = PIC_TYPE_UNKNOWN;
    ps_view_ctxt->ps_dec_err_status->u4_frm_sei_sync = SYNC_FRM_DEFAULT;
    ps_view_ctxt->ps_dec_err_status->u4_cur_frm = INIT_FRAME;
    ps_view_ctxt->ps_dec_err_status->u1_pic_aud_i = PIC_TYPE_UNKNOWN;
    ps_view_ctxt->u1_pr_sl_type = 0xFF;
    ps_view_ctxt->u2_mbx = 0xffff;
    ps_view_ctxt->u2_mby = 0;
    ps_view_ctxt->u2_total_mbs_coded = 0;

    ps_prev_poc = &ps_view_ctxt->s_prev_pic_poc;
    ps_cur_poc = &ps_view_ctxt->s_cur_pic_poc;
    ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb = 0;
    ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb = 0;
    ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom = 0;
    ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0] = 0;
    ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1] = 0;
    ps_prev_poc->u1_mmco_equalto5 = ps_cur_poc->u1_mmco_equalto5 = 0;
    ps_prev_poc->i4_top_field_order_count = ps_cur_poc->i4_top_field_order_count = 0;
    ps_prev_poc->i4_bottom_field_order_count = ps_cur_poc->i4_bottom_field_order_count = 0;
    ps_prev_poc->u1_bot_field = ps_cur_poc->u1_bot_field = 0;
    ps_prev_poc->u1_mmco_equalto5 = ps_cur_poc->u1_mmco_equalto5 = 0;
    ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst = 0;

    ps_view_ctxt->i4_max_poc = 0;
    ps_view_ctxt->i4_prev_max_display_seq = 0;
    ps_view_ctxt->u1_recon_mb_grp = 4;
    ps_view_ctxt->i4_reorder_depth = -1;
    ps_view_ctxt->u1_second_field = 0;
    ps_view_ctxt->s_prev_seq_params.u1_eoseq_pending = 0;
    ps_view_ctxt->u2_crop_offset_y = 0;
    ps_view_ctxt->u2_crop_offset_uv = 0;
    ps_view_ctxt->i4_vui_frame_rate = -1;
    ps_view_ctxt->i4_pic_type = NA_SLICE;
    ps_view_ctxt->i4_frametype = IV_NA_FRAME;
    ps_view_ctxt->i4_content_type = IV_CONTENTTYPE_NA;
    ps_view_ctxt->u1_res_changed = 0;
    ps_view_ctxt->u1_frame_decoded_flag = 0;
    ps_view_ctxt->u4_skip_frm_mask = SKIP_NONE;

    ps_view_ctxt->pf_cavlc_4x4res_block[0] = ih264d_cavlc_4x4res_block_totalcoeff_1;
    ps_view_ctxt->pf_cavlc_4x4res_block[1] = ih264d_cavlc_4x4res_block_totalcoeff_2to10;
    ps_view_ctxt->pf_cavlc_4x4res_block[2] = ih264d_cavlc_4x4res_block_totalcoeff_11to16;
    ps_view_ctxt->pf_cavlc_parse4x4coeff[0] = ih264d_cavlc_parse4x4coeff_n0to7;
    ps_view_ctxt->pf_cavlc_parse4x4coeff[1] = ih264d_cavlc_parse4x4coeff_n8;
    ps_view_ctxt->pf_cavlc_parse_8x8block[0] = ih264d_cavlc_parse_8x8block_none_available;
    ps_view_ctxt->pf_cavlc_parse_8x8block[1] = ih264d_cavlc_parse_8x8block_left_available;
    ps_view_ctxt->pf_cavlc_parse_8x8block[2] = ih264d_cavlc_parse_8x8block_top_available;
    ps_view_ctxt->pf_cavlc_parse_8x8block[3] = ih264d_cavlc_parse_8x8block_both_available;

    ps_view_ctxt->pf_fill_bs1[0][0] = ih264d_fill_bs1_16x16mb_pslice;
    ps_view_ctxt->pf_fill_bs1[0][1] = ih264d_fill_bs1_non16x16mb_pslice;
    ps_view_ctxt->pf_fill_bs1[1][0] = ih264d_fill_bs1_16x16mb_bslice;
    ps_view_ctxt->pf_fill_bs1[1][1] = ih264d_fill_bs1_non16x16mb_bslice;
    ps_view_ctxt->pf_fill_bs_xtra_left_edge[0] = ih264d_fill_bs_xtra_left_edge_cur_frm;
    ps_view_ctxt->pf_fill_bs_xtra_left_edge[1] = ih264d_fill_bs_xtra_left_edge_cur_fld;

    ps_view_ctxt->u2_prv_frame_num = 0;
    ps_view_ctxt->u1_top_bottom_decoded = 0;
    ps_view_ctxt->u1_dangling_field = 0;
    ps_view_ctxt->s_cab_dec_env.cabac_table = gau4_ih264d_cabac_table;
    ps_view_ctxt->pu1_left_mv_ctxt_inc = ps_view_ctxt->u1_left_mv_ctxt_inc_arr[0];
    ps_view_ctxt->pi1_left_ref_idx_ctxt_inc = &ps_view_ctxt->i1_left_ref_idx_ctx_inc_arr[0][0];
    ps_view_ctxt->pu1_left_yuv_dc_csbp = &ps_view_ctxt->u1_yuv_dc_csbp_topmb;
    ps_view_ctxt->u1_flushfrm = 0;
    ps_view_ctxt->s_cab_dec_env.pv_codec_handle = ps_view_ctxt;
    ps_view_ctxt->ps_bitstrm->pv_codec_handle = ps_view_ctxt;

    memset(ps_view_ctxt->disp_bufs, 0, (MAX_DISP_BUFS_NEW) * sizeof(disp_buf_t));
    memset(ps_view_ctxt->u4_disp_buf_mapping, 0, (MAX_DISP_BUFS_NEW) * sizeof(UWORD32));
    memset(ps_view_ctxt->u4_disp_buf_to_be_freed, 0, (MAX_DISP_BUFS_NEW) * sizeof(UWORD32));

    ih264d_init_arch(ps_view_ctxt);
    ih264d_init_function_ptr(ps_view_ctxt);
    ps_view_ctxt->e_frm_out_mode = IVD_DISPLAY_FRAME_OUT;
    ps_view_ctxt->init_done = 1;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctxt_init(imvcd_create_ip_t *ps_ip, mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i4_mem_size;
    void *pv_buf;

    FT_ALIGNED_ALLOC *pf_aligned_alloc = ps_ip->s_ivd_ip.pf_aligned_alloc;

    void *pv_mem_ctxt = ps_ip->s_ivd_ip.pv_mem_ctxt;
    const WORD32 i4_default_alignment = 128;

    memset(ps_mvcd_ctxt, 0, sizeof(ps_mvcd_ctxt[0]));

    i4_mem_size = sizeof(mvc_dpb_manager_t);
    ps_mvcd_ctxt->ps_dpb_mgr = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == ps_mvcd_ctxt->ps_dpb_mgr), IV_FAIL);
    memset(ps_mvcd_ctxt->ps_dpb_mgr, 0, i4_mem_size);

    imvcd_init_dpb_mgr(ps_mvcd_ctxt->ps_dpb_mgr, &ps_mvcd_ctxt->s_mvc_au_buf_mgr,
                       &ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr, &ps_mvcd_ctxt->s_mvc_disp_buf_mgr);

    ih264_disp_mgr_init(&ps_mvcd_ctxt->s_mvc_disp_buf_mgr);

    i4_mem_size = sizeof(buf_mgr_t) + ithread_get_mutex_lock_size();
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_mem = pv_buf;
    ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt = pv_buf;

    ih264_buf_mgr_init(ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_mem);

    i4_mem_size = sizeof(buf_mgr_t) + ithread_get_mutex_lock_size();
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_mem = pv_buf;
    ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt = pv_buf;

    ih264_buf_mgr_init(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_mem);

    if(IV_SUCCESS != imvcd_view_ctxt_init(ps_ip, &ps_mvcd_ctxt->s_view_dec_ctxt))
    {
        return IV_FAIL;
    }

    ps_mvcd_ctxt->u2_num_views = 0;
    ps_mvcd_ctxt->u2_num_views_decoded = 0;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_allocate_static_bufs(imvcd_create_ip_t *ps_ip,
                                                       imvcd_create_op_t *ps_op)
{
    iv_obj_t *ps_dec_hdl;
    mvc_dec_ctxt_t *ps_mvcd_ctxt;

    WORD32 i4_mem_size;

    FT_ALIGNED_ALLOC *pf_aligned_alloc = ps_ip->s_ivd_ip.pf_aligned_alloc;

    void *pv_mem_ctxt = ps_ip->s_ivd_ip.pv_mem_ctxt;
    const WORD32 i4_default_alignment = 128;

    i4_mem_size = sizeof(ps_dec_hdl[0]);
    ps_op->s_ivd_op.pv_handle = ps_dec_hdl =
        (iv_obj_t *) pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);

    if(NULL == ps_dec_hdl)
    {
        return IV_FAIL;
    }

    i4_mem_size = sizeof(ps_mvcd_ctxt[0]);
    ps_dec_hdl->pv_codec_handle = ps_mvcd_ctxt =
        (mvc_dec_ctxt_t *) pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);

    if(NULL == ps_mvcd_ctxt)
    {
        return IV_FAIL;
    }

    if(IV_SUCCESS != imvcd_ctxt_init(ps_ip, ps_mvcd_ctxt))
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

/* Description - 'Create' API for MVC Decoder */
static IV_API_CALL_STATUS_T imvcd_create(imvcd_create_ip_t *ps_ip, imvcd_create_op_t *ps_op)
{
    if(IV_SUCCESS != imvcd_check_create_structs(ps_ip, ps_op))
    {
        return IV_FAIL;
    }

    if(IV_SUCCESS != imvcd_allocate_static_bufs(ps_ip, ps_op))
    {
        imvcd_free_static_bufs((iv_obj_t *) ps_op->s_ivd_op.pv_handle);

        return IV_FAIL;
    }

    return IV_SUCCESS;
}

/* Description - 'Delete' API for MVC Decoder */
static IV_API_CALL_STATUS_T imvcd_delete(iv_obj_t *ps_dec_hdl)
{
    if(IV_SUCCESS != imvcd_check_dec_handle(ps_dec_hdl))
    {
        return IV_FAIL;
    }

    imvcd_free_static_bufs(ps_dec_hdl);

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_flush_mode_decode(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                                    imvcd_video_decode_op_t *ps_op)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    if(!ps_view_ctxt->u1_init_dec_flag)
    {
        ps_view_ctxt->u1_flushfrm = 0;
        ps_mvcd_ctxt->b_flush_enabled = false;
        ps_op->s_ivd_op.u4_output_present = 0;

        return IV_FAIL;
    }

    if(IV_SUCCESS != imvcd_get_next_display_au_buf(ps_mvcd_ctxt))
    {
        ps_view_ctxt->u1_flushfrm = 0;
        ps_mvcd_ctxt->b_flush_enabled = false;
        ps_op->s_ivd_op.u4_output_present = false;

        return IV_SUCCESS;
    }

    ih264d_export_sei_params(&ps_op->s_ivd_op.s_sei_decode_op, ps_view_ctxt);

    ps_op->s_ivd_op.u4_pic_wd = ps_view_ctxt->u2_disp_width;
    ps_op->s_ivd_op.u4_pic_ht = ps_view_ctxt->u2_disp_height;
    ps_op->s_ivd_op.u4_ts = ps_view_ctxt->s_disp_op.u4_ts;
    ps_op->s_ivd_op.u4_output_present = 1;
    ps_op->s_ivd_op.e_output_format = IV_YUV_420P;

    imvcd_convert_to_app_disp_buf(ps_mvcd_ctxt, ps_op->ps_view_disp_bufs);

    return IV_SUCCESS;
}

static void imvcd_fill_output_struct_from_context(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                                  imvcd_video_decode_op_t *ps_op)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    if((ps_op->s_ivd_op.u4_error_code & 0xff) != ERROR_DYNAMIC_RESOLUTION_NOT_SUPPORTED)
    {
        ps_op->s_ivd_op.u4_pic_wd = ps_view_ctxt->u2_disp_width;
        ps_op->s_ivd_op.u4_pic_ht = ps_view_ctxt->u2_disp_height;
    }

    ps_op->s_ivd_op.u4_output_present = ps_view_ctxt->u4_output_present;

    ps_op->s_ivd_op.e_output_format = IV_YUV_420P;

    imvcd_convert_to_app_disp_buf(ps_mvcd_ctxt, ps_op->ps_view_disp_bufs);

    ih264d_export_sei_params(&ps_op->s_ivd_op.s_sei_decode_op, ps_view_ctxt);
}

static void imvcd_video_decode_clean_return(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                            imvcd_video_decode_ip_t *ps_ip,
                                            imvcd_video_decode_op_t *ps_op)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ih264d_signal_decode_thread(ps_view_ctxt);
    ih264d_signal_bs_deblk_thread(ps_view_ctxt);

    imvcd_fill_output_struct_from_context(ps_mvcd_ctxt, ps_op);

    ps_op->s_ivd_op.u4_frame_decoded_flag = 0;
    ps_op->s_ivd_op.u4_num_bytes_consumed = ps_ip->s_ivd_ip.u4_num_Bytes;
}

static FORCEINLINE void imvcd_update_num_pps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_mvcd_ctxt->u1_num_pps = 0;

    for(i = 0; i < MAX_NUM_PIC_PARAMS; i++)
    {
        if(ps_view_ctxt->ps_pps[i].u1_is_valid)
        {
            UWORD8 u1_sps_id = ps_view_ctxt->ps_pps[i].ps_sps->u1_seq_parameter_set_id;

            if(ps_mvcd_ctxt->as_subset_sps[u1_sps_id].s_sps_data.u1_is_valid)
            {
                ps_mvcd_ctxt->aps_pps_id_to_subset_sps_map[i] =
                    &ps_mvcd_ctxt->as_subset_sps[u1_sps_id];
            }

            ps_mvcd_ctxt->u1_num_pps++;
        }
    }
}

static FORCEINLINE void imvcd_update_num_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_mvcd_ctxt->u1_num_sps = 0;

    for(i = 0; i < MAX_NUM_SEQ_PARAMS; i++)
    {
        if(ps_view_ctxt->ps_sps[i].u1_is_valid)
        {
            ps_mvcd_ctxt->u1_num_sps++;
        }
    }
}

static FORCEINLINE void imvcd_update_num_subset_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    ps_mvcd_ctxt->u1_num_subset_sps = 0;

    for(i = 0; i < MAX_NUM_SEQ_PARAMS; i++)
    {
        if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_is_valid)
        {
            ps_mvcd_ctxt->u1_num_subset_sps++;
        }
    }
}

static IV_API_CALL_STATUS_T imvcd_view_decode(iv_obj_t *ps_dec_hdl, imvcd_video_decode_ip_t *ps_ip,
                                              imvcd_video_decode_op_t *ps_op)
{
    UWORD8 *pu1_input_buffer;
    UWORD8 *pu1_bitstream_buf;
    UWORD32 u4_bitstream_buf_size;
    WORD32 i4_nalu_length;
    UWORD32 u4_length_of_start_code;
    WORD32 i4_error_code;

    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UWORD32 u4_num_bytes_consumed = 0;
    UWORD32 u4_num_bytes_remaining = ps_ip->s_ivd_ip.u4_num_Bytes;
    bool b_first_start_code_found = false;
    bool b_frame_data_left = true;
    bool b_header_data_left = true;
    UWORD32 u4_next_is_aud = 0;

    ASSERT(u4_num_bytes_remaining > 0);

    imvcd_view_init(ps_mvcd_ctxt);

    do
    {
        pu1_input_buffer = ((UWORD8 *) ps_ip->s_ivd_ip.pv_stream_buffer) + u4_num_bytes_consumed;

        if(!ps_view_ctxt->pu1_bits_buf_dynamic &&
           is_header_decoded(ps_view_ctxt->i4_header_decoded, PPS))
        {
            if(IV_SUCCESS !=
               imvcd_bitstream_buf_alloc(
                   ps_view_ctxt, is_header_decoded(ps_view_ctxt->i4_header_decoded, SUBSET_SPS)
                                     ? ps_mvcd_ctxt->u2_num_views
                                     : 1))
            {
                return IV_FAIL;
            }
        }

        if(ps_view_ctxt->pu1_bits_buf_dynamic)
        {
            pu1_bitstream_buf = ps_view_ctxt->pu1_bits_buf_dynamic;
            u4_bitstream_buf_size = ps_view_ctxt->u4_dynamic_bits_buf_size;
        }
        else
        {
            pu1_bitstream_buf = ps_view_ctxt->pu1_bits_buf_static;
            u4_bitstream_buf_size = ps_view_ctxt->u4_static_bits_buf_size;
        }

        i4_nalu_length = ih264d_find_start_code(pu1_input_buffer, 0, u4_num_bytes_remaining,
                                                &u4_length_of_start_code, &u4_next_is_aud);

        if(i4_nalu_length == -1)
        {
            i4_nalu_length = 0;
        }

        if((0 != u4_next_is_aud) && (1 != u4_next_is_aud))
        {
            return IV_FAIL;
        }

        /* Ignore bytes beyond the allocated size of intermediate buffer */
        /* Since 8 bytes are read ahead, ensure 8 bytes are free at the
        end of the buffer, which will be memset to 0 after emulation prevention */
        i4_nalu_length = MIN((UWORD32) i4_nalu_length, u4_bitstream_buf_size - 8);

        if(i4_nalu_length)
        {
            memcpy(pu1_bitstream_buf, pu1_input_buffer + u4_length_of_start_code, i4_nalu_length);

            /* Decoder may read extra 8 bytes near end of the frame */
            if(((UWORD32) (i4_nalu_length + 8)) < u4_bitstream_buf_size)
            {
                memset(pu1_bitstream_buf + i4_nalu_length, 0, 8 * sizeof(pu1_bitstream_buf[0]));
            }

            b_first_start_code_found = true;
        }
        else
        {
            if(!b_first_start_code_found)
            {
                ps_view_ctxt->i4_error_code = ERROR_START_CODE_NOT_FOUND;
                ps_op->s_ivd_op.u4_error_code |= 1 << IVD_INSUFFICIENTDATA;

                if(ps_view_ctxt->u4_pic_buf_got == 0)
                {
                    imvcd_fill_output_struct_from_context(ps_mvcd_ctxt, ps_op);

                    ps_op->s_ivd_op.u4_error_code = ps_view_ctxt->i4_error_code;

                    imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

                    return IV_FAIL;
                }
                else
                {
                    ps_view_ctxt->u1_pic_decode_done = 1;

                    continue;
                }
            }
            else
            {
                /* a start code has already been found earlier in the same process
                 * call*/
                b_frame_data_left = false;
                b_header_data_left = false;

                if(!ps_view_ctxt->i4_decode_header && !ps_view_ctxt->u4_pic_buf_got)
                {
                    ps_op->s_ivd_op.u4_error_code = ih264d_map_error(ERROR_UNKNOWN_NAL);

                    imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

                    return IV_FAIL;
                }

                continue;
            }
        }

        ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded] =
            NAL_UNIT_TYPE(pu1_bitstream_buf[0]);
        ps_mvcd_ctxt->au1_nal_ref_idc[ps_mvcd_ctxt->u2_num_views_decoded] =
            NAL_REF_IDC(pu1_bitstream_buf[0]);

        if(ps_view_ctxt->u4_dec_thread_created &&
           !is_slice_nalu_type(ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded]))
        {
            ps_op->s_ivd_op.u4_error_code = ERROR_FEATURE_UNAVAIL;

            imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

            return IV_FAIL;
        }

        if(!is_mvc_nalu(ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded]))
        {
            ivd_video_decode_op_t s_avc_op;

            i4_error_code =
                ih264d_parse_nal_unit(ps_dec_hdl, &s_avc_op, pu1_bitstream_buf, i4_nalu_length);
        }
        else
        {
            i4_error_code = imvcd_nalu_parser(ps_mvcd_ctxt, pu1_bitstream_buf, i4_nalu_length);
        }

        if(OK != i4_error_code)
        {
            ps_op->s_ivd_op.u4_error_code = i4_error_code;

            imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

            return IV_FAIL;
        }
        else if(PPS == ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded])
        {
            imvcd_update_num_pps(ps_mvcd_ctxt);
        }
        else if(SPS == ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded])
        {
            imvcd_update_num_sps(ps_mvcd_ctxt);
        }
        else if(SUBSET_SPS == ps_mvcd_ctxt->ae_nalu_id[ps_mvcd_ctxt->u2_num_views_decoded])
        {
            imvcd_update_num_subset_sps(ps_mvcd_ctxt);
        }

        b_header_data_left = ps_view_ctxt->i4_decode_header &&
                             (!is_header_decoded(ps_view_ctxt->i4_header_decoded, SPS) ||
                              !is_header_decoded(ps_view_ctxt->i4_header_decoded, PPS)) &&
                             (u4_num_bytes_consumed < ps_ip->s_ivd_ip.u4_num_Bytes);
        b_frame_data_left = (!ps_view_ctxt->i4_decode_header &&
                             (!ps_view_ctxt->u1_pic_decode_done || u4_next_is_aud)) &&
                            (u4_num_bytes_consumed < ps_ip->s_ivd_ip.u4_num_Bytes);

        u4_num_bytes_consumed += i4_nalu_length + u4_length_of_start_code;
        u4_num_bytes_remaining -= i4_nalu_length + u4_length_of_start_code;

    } while(b_header_data_left || b_frame_data_left);

    if((i4_error_code == IVD_RES_CHANGED) || (i4_error_code == IVD_MEM_ALLOC_FAILED) ||
       (i4_error_code == ERROR_UNAVAIL_PICBUF_T) || (i4_error_code == ERROR_UNAVAIL_MVBUF_T) ||
       (i4_error_code == ERROR_INV_SPS_PPS_T))
    {
        ih264d_signal_decode_thread(ps_view_ctxt);

        if(ps_view_ctxt->u4_num_cores == 3)
        {
            ih264d_signal_bs_deblk_thread(ps_view_ctxt);
        }

        /* dont consume bitstream for change in resolution case */
        if(i4_error_code == IVD_RES_CHANGED)
        {
            ps_op->s_ivd_op.u4_num_bytes_consumed -= u4_num_bytes_consumed;
        }

        imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

        return IV_FAIL;
    }

    if(ps_view_ctxt->u1_separate_parse)
    {
        if(ps_view_ctxt->u4_num_cores == 2)
        {
            if((ps_view_ctxt->u4_nmb_deblk == 0) && (ps_view_ctxt->u4_start_recon_deblk == 1))
            {
                tfr_ctxt_t s_tfr_ctxt;

                UWORD32 u4_num_mbs, u4_max_addr;

                tfr_ctxt_t *ps_tfr_cxt = &s_tfr_ctxt;
                pad_mgr_t *ps_pad_mgr = &ps_view_ctxt->s_pad_mgr;
                nalu_mvc_ext_t *ps_cur_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);

                /*BS is done for all mbs while parsing*/
                u4_max_addr = (ps_view_ctxt->u2_frm_wd_in_mbs * ps_view_ctxt->u2_frm_ht_in_mbs) - 1;
                ps_view_ctxt->u4_cur_bs_mb_num = u4_max_addr + 1;

                ps_view_ctxt->ps_cur_pic = &ps_view_ctxt->s_cur_pic;
                imvcd_convert_au_buf_to_view_buf(ps_mvcd_ctxt->ps_cur_au, &ps_view_ctxt->s_cur_pic,
                                                 ps_mvcd_ctxt->u2_num_views_decoded,
                                                 ps_cur_nalu_mvc_ext->u2_view_id);

                ih264d_init_deblk_tfr_ctxt(ps_view_ctxt, ps_pad_mgr, ps_tfr_cxt,
                                           ps_view_ctxt->u2_frm_wd_in_mbs, 0);

                u4_num_mbs = u4_max_addr - ps_view_ctxt->u4_cur_deblk_mb_num + 1;

                if(u4_num_mbs != 0)
                {
                    ih264d_check_mb_map_deblk(ps_view_ctxt, u4_num_mbs, ps_tfr_cxt, 1);
                }

                ps_view_ctxt->u4_start_recon_deblk = 0;
            }
        }

        ih264d_signal_decode_thread(ps_view_ctxt);

        if(ps_view_ctxt->u4_num_cores == 3)
        {
            ih264d_signal_bs_deblk_thread(ps_view_ctxt);
        }
    }

    DATA_SYNC();

    // Report if header (sps and pps) has not been decoded yet
    if(ps_view_ctxt->i4_decode_header &&
       (!is_header_decoded(ps_view_ctxt->i4_header_decoded, SPS) &&
        !is_header_decoded(ps_view_ctxt->i4_header_decoded, PPS)))
    {
        ps_op->s_ivd_op.u4_error_code |= (1 << IVD_INSUFFICIENTDATA);

        imvcd_video_decode_clean_return(ps_mvcd_ctxt, ps_ip, ps_op);

        return IV_FAIL;
    }

    if(ps_view_ctxt->u4_pic_buf_got)
    {
        ps_view_ctxt->u1_top_bottom_decoded = TOP_FIELD_ONLY | BOT_FIELD_ONLY;

        if(((ps_view_ctxt->ps_dec_err_status->u1_err_flag & REJECT_CUR_PIC) == 0) &&
           ps_view_ctxt->u1_pic_decode_done)
        {
            nalu_mvc_ext_t *ps_cur_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);

            if(!ps_mvcd_ctxt->au1_nal_ref_idc[ps_mvcd_ctxt->u2_num_views_decoded] &&
               ps_cur_nalu_mvc_ext->u1_inter_view_flag)
            {
                ps_view_ctxt->ps_cur_slice->u1_nal_ref_idc = 1;
            }

            /* Padding only. Deblk has happened already. */
            ih264d_deblock_picture_progressive(ps_view_ctxt);

            if(!ps_mvcd_ctxt->au1_nal_ref_idc[ps_mvcd_ctxt->u2_num_views_decoded] &&
               ps_cur_nalu_mvc_ext->u1_inter_view_flag)
            {
                ps_view_ctxt->ps_cur_slice->u1_nal_ref_idc = 0;
            }
        }

        /*Update the i4_frametype at the end of picture*/
        if(imvcd_is_idr_au(ps_mvcd_ctxt))
        {
            ps_view_ctxt->i4_frametype = IV_IDR_FRAME;
        }
        else if(ps_view_ctxt->i4_pic_type == B_SLICE)
        {
            ps_view_ctxt->i4_frametype = IV_B_FRAME;
        }
        else if(ps_view_ctxt->i4_pic_type == P_SLICE)
        {
            ps_view_ctxt->i4_frametype = IV_P_FRAME;
        }
        else if(ps_view_ctxt->i4_pic_type == I_SLICE)
        {
            ps_view_ctxt->i4_frametype = IV_I_FRAME;
        }

        ps_view_ctxt->i4_content_type = ps_view_ctxt->ps_cur_slice->u1_field_pic_flag;
    }

    /* close deblock thread if it is not closed yet*/
    if(ps_view_ctxt->u4_num_cores == 3)
    {
        ih264d_signal_bs_deblk_thread(ps_view_ctxt);
    }

    if(ps_view_ctxt->u4_dec_thread_created)
    {
        ih264d_signal_decode_thread(ps_view_ctxt);
    }

    if(ps_view_ctxt->u4_bs_deblk_thread_created)
    {
        ih264d_signal_bs_deblk_thread(ps_view_ctxt);
    }

    ps_op->s_ivd_op.u4_num_bytes_consumed = u4_num_bytes_consumed;

    DATA_SYNC();

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_finish_au_decode(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                                   imvcd_video_decode_op_t *ps_op)
{
    WORD32 i4_error_code;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_slice_params_t *ps_cur_slice = ps_view_ctxt->ps_cur_slice;
    mvc_au_buffer_t *ps_cur_au = ps_mvcd_ctxt->ps_cur_au;
    mvc_dpb_manager_t *ps_dpb_mgr = ps_mvcd_ctxt->ps_dpb_mgr;

    bool b_is_idr = imvcd_is_idr_au(ps_mvcd_ctxt);
    bool b_is_ref_au = !!ps_mvcd_ctxt->au1_nal_ref_idc[ps_mvcd_ctxt->u2_num_views - 1];
    WORD64 i8_display_poc =
        ((WORD64) ps_view_ctxt->i4_prev_max_display_seq) + ((WORD64) ps_cur_au->i4_poc);

    imvcd_dpb_delete_nonref_nondisplay_pics(ps_dpb_mgr);

    if(ps_cur_slice->u1_mmco_equalto5 || b_is_idr)
    {
        ps_cur_au->i4_poc = 0;
        ps_cur_au->i4_avg_poc = 0;

        if(ps_view_ctxt->u2_total_mbs_coded == (ps_view_ctxt->ps_cur_sps->u2_max_mb_addr + 1))
        {
            imvcd_reset_dpb(ps_dpb_mgr);
        }

        imvcd_dpb_release_display_bufs(ps_dpb_mgr);
    }

    if(IVD_DECODE_FRAME_OUT != ps_view_ctxt->e_frm_out_mode)
    {
        i4_error_code = imvcd_dpb_assign_display_seq(ps_dpb_mgr);

        if(OK != i4_error_code)
        {
            return IV_FAIL;
        }
    }

    if(b_is_ref_au)
    {
        ih264_buf_mgr_set_status(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt,
                                 ps_cur_au->i4_pic_buf_id, BUF_MGR_REF);

        ih264_buf_mgr_set_status(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                                 ps_cur_au->i4_mv_buf_id, BUF_MGR_REF);

        ps_view_ctxt->au1_pic_buf_ref_flag[ps_cur_au->i4_pic_buf_id] = 1;
    }
    else
    {
        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt,
                              ps_cur_au->i4_pic_buf_id, BUF_MGR_REF);

        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                              ps_cur_au->i4_mv_buf_id, BUF_MGR_REF | BUF_MGR_IO);

        ps_view_ctxt->au1_pic_buf_ref_flag[ps_cur_au->i4_pic_buf_id] = 0;
    }

    if((!ps_view_ctxt->u1_last_pic_not_decoded &&
        (0 == (ps_view_ctxt->ps_cur_pic->u4_pack_slc_typ & ps_view_ctxt->u4_skip_frm_mask))) ||
       b_is_idr)
    {
        ih264_buf_mgr_set_status(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt,
                                 ps_cur_au->i4_pic_buf_id, BUF_MGR_IO);
    }

    if(IS_OUT_OF_RANGE_S32(i8_display_poc))
    {
        ps_view_ctxt->i4_prev_max_display_seq = 0;
    }

    i4_error_code = imvcd_dpb_insert_pic_in_display_list(
        ps_dpb_mgr, i8_display_poc, ps_cur_au->i4_frame_num, ps_cur_au->i4_pic_buf_id);

    if(i4_error_code != OK)
    {
        return IV_FAIL;
    }

    if(IVD_DECODE_FRAME_OUT == ps_view_ctxt->e_frm_out_mode)
    {
        i4_error_code = imvcd_dpb_assign_display_seq(ps_dpb_mgr);

        if(i4_error_code != OK)
        {
            return IV_FAIL;
        }
    }

    ps_view_ctxt->u4_total_frames_decoded++;

    /* In case the decoder is configured to run in low delay mode,
     * then get display buffer and then format convert.
     * Note in this mode, format conversion does not run paralelly in a thread
     * and adds to the codec cycles
     */
    if((IVD_DECODE_FRAME_OUT == ps_view_ctxt->e_frm_out_mode) && ps_view_ctxt->u1_init_dec_flag)
    {
        i4_error_code = imvcd_get_next_display_au_buf(ps_mvcd_ctxt);

        if(i4_error_code != OK)
        {
            return IV_FAIL;
        }

        ps_op->s_ivd_op.u4_output_present = 1;
    }

    ps_cur_au->u1_pic_type |= TOP_REF | BOT_REF;

    if(ps_view_ctxt->u4_pic_buf_got)
    {
        if(ps_view_ctxt->u1_last_pic_not_decoded)
        {
            return IV_FAIL;
        }
        else if(b_is_ref_au)
        {
            if(b_is_idr)
            {
                ps_dpb_mgr->u1_mmco_error_in_seq = 0;

                if(!ps_view_ctxt->ps_dpb_cmds->u1_long_term_reference_flag)
                {
                    imvcd_reset_dpb(ps_dpb_mgr);

                    i4_error_code = imvcd_dpb_insert_st_node(ps_dpb_mgr, ps_cur_au);

                    if(i4_error_code != OK)
                    {
                        return IV_FAIL;
                    }

                    ps_dpb_mgr->u1_max_lt_frame_idx = NO_LONG_TERM_INDICIES;
                }
                else
                {
                    i4_error_code = imvcd_dpb_insert_st_node(ps_dpb_mgr, ps_cur_au);

                    if(i4_error_code != OK)
                    {
                        return IV_FAIL;
                    }

                    imvcd_dpb_delete_st_node_or_make_lt(ps_dpb_mgr, ps_cur_au->i4_pic_num, 0);

                    ps_dpb_mgr->u1_max_lt_frame_idx = 0;
                }
            }
            else if(!ps_dpb_mgr->u1_mmco_error_in_seq)
            {
                i4_error_code = imvcd_dpb_do_mmco(ps_view_ctxt->ps_dpb_cmds, ps_dpb_mgr, ps_cur_au,
                                                  ps_view_ctxt->ps_cur_sps->u1_num_ref_frames,
                                                  ps_view_ctxt->e_dec_status);

                ps_dpb_mgr->u1_mmco_error_in_seq = i4_error_code != OK;
            }

            i4_error_code = imvcd_dpb_update_default_index_list(ps_dpb_mgr);

            if(i4_error_code != OK)
            {
                return IV_FAIL;
            }
        }
    }

    ps_op->s_ivd_op.u4_frame_decoded_flag = 1;

    return IV_SUCCESS;
}

/* Description - 'AU Decode' API for MVC Decoder */
static IV_API_CALL_STATUS_T imvcd_decode(iv_obj_t *ps_dec_hdl, imvcd_video_decode_ip_t *ps_ip,
                                         imvcd_video_decode_op_t *ps_op)
{
    IV_API_CALL_STATUS_T e_retval;

    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    imvcd_video_decode_ip_t s_view_ip = ps_ip[0];
    imvcd_video_decode_op_t s_view_op = ps_op[0];

    UWORD16 u2_num_views_decoded = 0;
    UWORD16 u2_num_views = (ps_mvcd_ctxt->b_flush_enabled || ps_mvcd_ctxt->b_header_only_decode)
                               ? 1
                               : ps_mvcd_ctxt->u2_num_views;

    ps_mvcd_ctxt->u2_num_views_decoded = 0;

    if(IV_SUCCESS != imvcd_check_dec_handle(ps_dec_hdl))
    {
        return IV_FAIL;
    }

    if(IV_SUCCESS != imvcd_check_decode_structs(ps_dec_hdl, ps_ip, ps_op))
    {
        return IV_FAIL;
    }

    if(!ps_mvcd_ctxt->b_header_only_decode)
    {
        if(IV_SUCCESS != imvcd_au_error_checks(ps_mvcd_ctxt, ps_ip))
        {
            return IV_FAIL;
        }
    }

    /*Data memory barries instruction,so that bitstream write by the application
     * is complete*/
    DATA_SYNC();

    imvcd_au_init(ps_dec_hdl, ps_ip, ps_op);

    if(ps_mvcd_ctxt->b_flush_enabled)
    {
        return imvcd_flush_mode_decode(ps_mvcd_ctxt, ps_op);
    }

    while(u2_num_views_decoded < u2_num_views)
    {
        e_retval = imvcd_view_decode(ps_dec_hdl, &s_view_ip, &s_view_op);

        if(IV_SUCCESS != e_retval)
        {
            ps_op->s_ivd_op.u4_error_code = s_view_op.s_ivd_op.u4_error_code;

            return IV_FAIL;
        }

        s_view_ip.s_ivd_ip.pv_stream_buffer = ((UWORD8 *) s_view_ip.s_ivd_ip.pv_stream_buffer) +
                                              s_view_op.s_ivd_op.u4_num_bytes_consumed;
        s_view_ip.s_ivd_ip.u4_num_Bytes -= s_view_op.s_ivd_op.u4_num_bytes_consumed;
        ps_op->s_ivd_op.u4_num_bytes_consumed += s_view_op.s_ivd_op.u4_num_bytes_consumed;

        u2_num_views_decoded++;
        ps_mvcd_ctxt->u2_num_views_decoded++;
    }

    if(!ps_mvcd_ctxt->b_header_only_decode)
    {
        e_retval = imvcd_finish_au_decode(ps_mvcd_ctxt, ps_op);

        if(IV_SUCCESS != e_retval)
        {
            return IV_FAIL;
        }
    }

    ps_op->s_ivd_op.u4_pic_wd = ps_view_ctxt->u2_disp_width;
    ps_op->s_ivd_op.u4_pic_ht = ps_view_ctxt->u2_disp_height;
    ps_op->s_ivd_op.u4_output_present = ps_view_ctxt->u4_output_present;
    ps_op->s_ivd_op.u4_ts = ps_view_ctxt->s_disp_op.u4_ts;
    ps_op->s_ivd_op.i4_reorder_depth = ps_view_ctxt->i4_reorder_depth;
    ps_op->s_ivd_op.e_output_format = IV_YUV_420P;

    if(ps_op->s_ivd_op.u4_output_present)
    {
        imvcd_convert_to_app_disp_buf(ps_mvcd_ctxt, ps_op->ps_view_disp_bufs);
    }

    return e_retval;
}

static IV_API_CALL_STATUS_T imvcd_ctl_set_dec_mode(iv_obj_t *ps_dec_hdl,
                                                   imvcd_set_config_ip_t *ps_ip,
                                                   imvcd_set_config_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_view_ctxt->u4_skip_frm_mask = SKIP_NONE;

    ps_op->s_ivd_op.u4_error_code = 0;

    ps_view_ctxt->u4_app_disp_width = 0;

    if(ps_ip->s_ivd_ip.e_frm_skip_mode != IVD_SKIP_NONE)
    {
        ps_op->s_ivd_op.u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);

        return IV_FAIL;
    }

    if(ps_ip->s_ivd_ip.e_vid_dec_mode == IVD_DECODE_FRAME)
    {
        ps_view_ctxt->i4_decode_header = 0;
        ps_mvcd_ctxt->b_header_only_decode = false;
    }
    else if(ps_ip->s_ivd_ip.e_vid_dec_mode == IVD_DECODE_HEADER)
    {
        ps_view_ctxt->i4_decode_header = 1;
        ps_mvcd_ctxt->b_header_only_decode = true;
    }
    else
    {
        ps_op->s_ivd_op.u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);

        return IV_FAIL;
    }

    if((ps_ip->s_ivd_ip.e_frm_out_mode != IVD_DECODE_FRAME_OUT) &&
       (ps_ip->s_ivd_ip.e_frm_out_mode != IVD_DISPLAY_FRAME_OUT))
    {
        ps_op->s_ivd_op.u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);

        return IV_FAIL;
    }

    ps_mvcd_ctxt->b_flush_enabled = false;
    ps_view_ctxt->e_frm_out_mode = ps_ip->s_ivd_ip.e_frm_out_mode;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_set_num_cores(iv_obj_t *ps_dec_hdl,
                                                    imvcd_set_num_cores_ip_t *ps_ip,
                                                    imvcd_set_num_cores_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_view_ctxt->u4_num_cores = ps_ip->u4_num_cores;

    ps_op->u4_error_code = 0;

    if(ps_view_ctxt->u4_num_cores == 1)
    {
        ps_view_ctxt->u1_separate_parse = 0;
    }
    else
    {
        ps_view_ctxt->u1_separate_parse = 1;
    }

    /*using only upto three threads currently*/
    if(ps_view_ctxt->u4_num_cores > 3)
    {
        ps_view_ctxt->u4_num_cores = 3;
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_set_arch(iv_obj_t *ps_dec_hdl, imvcd_set_arch_ip_t *ps_ip,
                                               imvcd_set_arch_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_view_ctxt->e_processor_arch = ps_ip->e_arch;
    ps_view_ctxt->e_processor_soc = ps_ip->e_soc;

    ps_op->u4_error_code = 0;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_set_degrade_mode(iv_obj_t *ps_dec_hdl,
                                                       imvcd_set_degrade_mode_ip_t *ps_ip,
                                                       imvcd_set_degrade_mode_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_view_ctxt->i4_degrade_type = ps_ip->i4_degrade_type;
    ps_view_ctxt->i4_nondegrade_interval = ps_ip->i4_nondegrade_interval;
    ps_view_ctxt->i4_degrade_pics = ps_ip->i4_degrade_pics;
    ps_view_ctxt->i4_degrade_pic_cnt = 0;

    ps_op->u4_error_code = 0;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_flush_dec(iv_obj_t *ps_dec_hdl, imvcd_flush_dec_ip_t *ps_ip,
                                                imvcd_flush_dec_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UNUSED(ps_ip);

    ps_op->s_ivd_op.u4_error_code = 0;

    ps_mvcd_ctxt->b_flush_enabled = true;
    ps_view_ctxt->u1_flushfrm = 1;

    if(ps_view_ctxt->u1_init_dec_flag)
    {
        imvcd_release_all_ref_bufs(ps_mvcd_ctxt, ps_view_ctxt->u1_pic_bufs);
        imvcd_dpb_release_display_bufs(ps_mvcd_ctxt->ps_dpb_mgr);
    }

    /* Ignore dangling fields during flush */
    ps_view_ctxt->u1_top_bottom_decoded = 0;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_get_buf_info(iv_obj_t *ps_dec_hdl,
                                                   imvcd_get_buf_info_ip_t *ps_ip,
                                                   imvcd_get_buf_info_op_t *ps_op)
{
    UWORD32 au4_min_out_buf_size[IVD_VIDDEC_MAX_IO_BUFFERS];
    UWORD32 u4_pic_wd, u4_pic_ht;
    UWORD32 i;

    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UNUSED(ps_ip);

    ps_op->s_ivd_op.u4_error_code = 0;

    ps_op->s_ivd_op.u4_num_disp_bufs = 0;
    ps_op->s_ivd_op.u4_min_num_in_bufs = MIN_IN_BUFS;

    u4_pic_wd = 0;
    u4_pic_ht = 0;

    if(is_header_decoded(ps_view_ctxt->i4_header_decoded, SPS))
    {
        u4_pic_wd = ps_view_ctxt->u2_disp_width;
        u4_pic_ht = ps_view_ctxt->u2_disp_height;
    }

    ps_op->s_mvc_buf_info.u2_num_views = ps_mvcd_ctxt->u2_num_views;

    for(i = 0; i < ps_op->s_ivd_op.u4_min_num_in_bufs; i++)
    {
        ps_op->s_ivd_op.u4_min_in_buf_size[i] =
            MAX(256000, u4_pic_wd * u4_pic_ht * ps_mvcd_ctxt->u2_num_views * 3 / 2);
    }

    ps_op->s_ivd_op.u4_min_num_out_bufs = ih264d_get_outbuf_size(
        u4_pic_wd, u4_pic_ht, ps_view_ctxt->u1_chroma_format, &au4_min_out_buf_size[0]);
    ps_op->s_ivd_op.u4_min_num_out_bufs *= ps_mvcd_ctxt->u2_num_views;

    for(i = 0; i < ps_op->s_ivd_op.u4_min_num_out_bufs; i++)
    {
        ps_op->s_ivd_op.u4_min_out_buf_size[i] = au4_min_out_buf_size[i % NUM_COMPONENTS];
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_ctl_get_vui(iv_obj_t *ps_dec_hdl, imvcd_get_vui_ip_t *ps_ip,
                                              imvcd_get_vui_op_t *ps_op)
{
    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UNUSED(ps_ip);

    ps_op->u4_error_code = 0;
    ps_op->b_is_vui_available = false;

    if((ps_mvcd_ctxt->u1_num_sps > 0) && ps_view_ctxt->ps_cur_sps)
    {
        ps_op->b_is_vui_available = ps_view_ctxt->ps_cur_sps->u1_vui_parameters_present_flag;

        if(ps_op->b_is_vui_available)
        {
            ps_op->u1_aspect_ratio_idc = ps_view_ctxt->ps_cur_sps->s_vui.u1_aspect_ratio_idc;
            ps_op->u2_sar_width = ps_view_ctxt->ps_cur_sps->s_vui.u2_sar_width;
            ps_op->u2_sar_height = ps_view_ctxt->ps_cur_sps->s_vui.u2_sar_height;
            ps_op->u1_overscan_appropriate_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_overscan_appropriate_flag;
            ps_op->u1_video_format = ps_view_ctxt->ps_cur_sps->s_vui.u1_video_format;
            ps_op->u1_video_full_range_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_video_full_range_flag;
            ps_op->u1_colour_primaries = ps_view_ctxt->ps_cur_sps->s_vui.u1_colour_primaries;
            ps_op->u1_tfr_chars = ps_view_ctxt->ps_cur_sps->s_vui.u1_tfr_chars;
            ps_op->u1_matrix_coeffs = ps_view_ctxt->ps_cur_sps->s_vui.u1_matrix_coeffs;
            ps_op->u1_cr_top_field = ps_view_ctxt->ps_cur_sps->s_vui.u1_cr_top_field;
            ps_op->u1_cr_bottom_field = ps_view_ctxt->ps_cur_sps->s_vui.u1_cr_bottom_field;
            ps_op->u4_num_units_in_tick = ps_view_ctxt->ps_cur_sps->s_vui.u4_num_units_in_tick;
            ps_op->u4_time_scale = ps_view_ctxt->ps_cur_sps->s_vui.u4_time_scale;
            ps_op->u1_fixed_frame_rate_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_fixed_frame_rate_flag;
            ps_op->u1_nal_hrd_params_present =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_nal_hrd_params_present;
            ps_op->u1_vcl_hrd_params_present =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_vcl_hrd_params_present;
            ps_op->u1_low_delay_hrd_flag = ps_view_ctxt->ps_cur_sps->s_vui.u1_low_delay_hrd_flag;
            ps_op->u1_pic_struct_present_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_pic_struct_present_flag;
            ps_op->u1_bitstream_restriction_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_bitstream_restriction_flag;
            ps_op->u1_mv_over_pic_boundaries_flag =
                ps_view_ctxt->ps_cur_sps->s_vui.u1_mv_over_pic_boundaries_flag;
            ps_op->u4_max_bytes_per_pic_denom =
                ps_view_ctxt->ps_cur_sps->s_vui.u4_max_bytes_per_pic_denom;
            ps_op->u4_max_bits_per_mb_denom =
                ps_view_ctxt->ps_cur_sps->s_vui.u4_max_bits_per_mb_denom;
            ps_op->u4_log2_max_mv_length_horz =
                ps_view_ctxt->ps_cur_sps->s_vui.u4_log2_max_mv_length_horz;
            ps_op->u4_log2_max_mv_length_vert =
                ps_view_ctxt->ps_cur_sps->s_vui.u4_log2_max_mv_length_vert;
            ps_op->u4_num_reorder_frames = ps_view_ctxt->ps_cur_sps->s_vui.u4_num_reorder_frames;
            ps_op->u4_max_dec_frame_buffering =
                ps_view_ctxt->ps_cur_sps->s_vui.u4_max_dec_frame_buffering;
        }
    }

    return IV_SUCCESS;
}

/* Description - 'Control Cmd' API for MVC Decoder */
static IV_API_CALL_STATUS_T imvcd_ctl_cmd_handler(iv_obj_t *ps_dec_hdl, void *pv_ip, void *pv_op)
{
    ivd_ctl_set_config_ip_t *ps_ip = (ivd_ctl_set_config_ip_t *) pv_ip;

    WORD32 i4_sub_cmd = ps_ip->e_sub_cmd;

    if(IV_SUCCESS != imvcd_check_dec_handle(ps_dec_hdl))
    {
        return IV_FAIL;
    }

    if(IV_SUCCESS != imvcd_check_ctl_structs(pv_ip, pv_op))
    {
        return IV_FAIL;
    }

    switch(i4_sub_cmd)
    {
        case IVD_CMD_CTL_SETPARAMS:
        {
            return imvcd_ctl_set_dec_mode(ps_dec_hdl, pv_ip, pv_op);
        }
        case IMVCD_CTL_SET_NUM_CORES:
        {
            return imvcd_ctl_set_num_cores(ps_dec_hdl, pv_ip, pv_op);
        }
        case IMVCD_CTL_SET_PROCESSOR:
        {
            return imvcd_ctl_set_arch(ps_dec_hdl, pv_ip, pv_op);
        }
        case IMVCD_CTL_DEGRADE:
        {
            return imvcd_ctl_set_degrade_mode(ps_dec_hdl, pv_ip, pv_op);
        }
        case IVD_CMD_CTL_FLUSH:
        {
            return imvcd_ctl_flush_dec(ps_dec_hdl, pv_ip, pv_op);
        }
        case IVD_CMD_CTL_GETBUFINFO:
        {
            return imvcd_ctl_get_buf_info(ps_dec_hdl, pv_ip, pv_op);
        }
        case IMVCD_CTL_GET_VUI_PARAMS:
        {
            return imvcd_ctl_get_vui(ps_dec_hdl, pv_ip, pv_op);
        }
        default:
        {
            return IV_FAIL;
        }
    }
}

IV_API_CALL_STATUS_T imvcd_api_function(iv_obj_t *ps_dec_hdl, void *pv_ip, void *pv_op)
{
    IVD_API_COMMAND_TYPE_T e_cmd = ((WORD32 *) pv_ip)[1];

    switch(e_cmd)
    {
        case IVD_CMD_CREATE:
        {
            return imvcd_create(pv_ip, pv_op);
        }
        case IVD_CMD_DELETE:
        {
            return imvcd_delete(ps_dec_hdl);
        }
        case IVD_CMD_VIDEO_CTL:
        {
            return imvcd_ctl_cmd_handler(ps_dec_hdl, pv_ip, pv_op);
        }
        case IVD_CMD_VIDEO_DECODE:
        {
            return imvcd_decode(ps_dec_hdl, pv_ip, pv_op);
        }
        default:
        {
            return IV_FAIL;
        }
    }
}
