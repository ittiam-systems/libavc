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
 *  isvcd_parse_slice.c
 *
 * @brief
 *  Contains routines that decodes a slice NAL unit
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include <assert.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ithread.h"
#include "isvcd_structs.h"
#include "ih264d_debug.h"
#include "ih264d_bitstrm.h"
#include "ih264d_parse_mb_header.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_pslice.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_utils.h"
#include "isvcd_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_defs.h"
#include "ih264d_error_handler.h"
#include "ih264d_tables.h"
#include "ih264d_defs.h"
#include "ih264d_mem_request.h"
#include "ih264d_parse_islice.h"
#include "ih264d_parse_slice.h"
#include "ih264d_mvpred.h"
#include "ih264d_mb_utils.h"
#include "ih264d_defs.h"
#include "ih264d_quant_scaling.h"
#include "ih264d_inter_pred.h"
#include "ih264d_sei.h"
#include "ih264_error.h"
#include "ih264_disp_mgr.h"
#include "ih264_buf_mgr.h"
#include "ih264d_thread_parse_decode.h"
#include "ih264d_thread_compute_bs.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_parse_islice.h"
#include "isvcd_parse_slice.h"
#include "isvcd_process_epslice.h"
#include "isvcd_process_ebslice.h"
#include "isvcd_thread_compute_bs.h"
#include "isvcd_thread_parse_decode.h"
#include "isvcd_deblocking.h"

#define RET_LAST_SKIP 0x80000000

WORD32 check_app_out_buf_size(dec_struct_t *ps_dec);
/*!
**************************************************************************
* \if Function name :  isvcd_verify_level \endif
*
* \brief
*    Initialize the Parameter required for all the slices for a picture
*
* \return           : Nothing
*
**************************************************************************
*/
WORD32 isvcd_verify_level(UWORD8 u1_level_idc)
{
    switch(u1_level_idc)
    {
        case H264_LEVEL_1_0:
        case H264_LEVEL_1_1:
        case H264_LEVEL_1_2:
        case H264_LEVEL_1_3:
        case H264_LEVEL_2_0:
        case H264_LEVEL_2_1:
        case H264_LEVEL_2_2:
        case H264_LEVEL_3_0:
        case H264_LEVEL_3_1:
        case H264_LEVEL_3_2:
        case H264_LEVEL_4_0:
        case H264_LEVEL_4_1:
        case H264_LEVEL_4_2:
        case H264_LEVEL_5_0:
        case H264_LEVEL_5_1:
            return OK;
        default:
            return NOT_OK;
    }
}

/*!
**************************************************************************
* \if Function name :  isvcd_start_of_pic \endif
*
* \brief
*    Initialize the Parameter required for all the slices for a picture
*
* \return           : Nothing
*
**************************************************************************
*/

WORD32 isvcd_start_of_pic(svc_dec_lyr_struct_t *ps_svc_lyr_dec, WORD32 i4_poc,
                          pocstruct_t *ps_temp_poc, UWORD16 u2_frame_num, dec_pic_params_t *ps_pps)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    pocstruct_t *ps_prev_poc = &ps_dec->s_cur_pic_poc;
    pocstruct_t *ps_cur_poc = ps_temp_poc;

    dec_slice_params_t *ps_cur_slice = ps_dec->ps_cur_slice;
    dec_seq_params_t *ps_seq = ps_dec->ps_cur_sps;
    UWORD8 u1_bottom_field_flag = ps_cur_slice->u1_bottom_field_flag;
    UWORD8 u1_field_pic_flag = ps_cur_slice->u1_field_pic_flag;
    /* high profile related declarations */
    WORD32 ret;

    H264_MUTEX_LOCK(&ps_dec->process_disp_mutex);

    if(u1_field_pic_flag == 1)
    {
        ps_dec->i4_error_code = ERROR_SVC_FIELD_PIC_UNSUPPORTED;
        return ERROR_SVC_FIELD_PIC_UNSUPPORTED;
    }

    /* check output buffer size given by the application */
    if(check_app_out_buf_size(ps_dec) != IV_SUCCESS) return IVD_DISP_FRM_ZERO_OP_BUF_SIZE;

    ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb;
    ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb;
    ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom;
    ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0];
    ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1];
    ps_prev_poc->u1_bot_field = ps_dec->ps_cur_slice->u1_bottom_field_flag;
    ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst;
    ps_prev_poc->u2_frame_num = u2_frame_num;
    ps_dec->i1_prev_mb_qp_delta = 0;
    ps_dec->i1_next_ctxt_idx = 0;

    ps_dec->u4_nmb_deblk = 0;
    if(ps_dec->u4_num_cores == 1) ps_dec->u4_nmb_deblk = 1;

    if(ps_seq->u1_mb_aff_flag == 1)
    {
        ps_dec->u4_nmb_deblk = 0;
        if(ps_dec->u4_num_cores > 2) ps_dec->u4_num_cores = 2;
    }

    ps_dec->u4_use_intrapred_line_copy = 0;

    if(ps_seq->u1_mb_aff_flag == 0)
    {
        ps_dec->u4_use_intrapred_line_copy = 1;
    }

    ps_dec->u4_app_disable_deblk_frm = 0;
    /* If degrade is enabled, set the degrade flags appropriately */
    if(ps_dec->i4_degrade_type && ps_dec->i4_degrade_pics)
    {
        WORD32 degrade_pic;
        ps_dec->i4_degrade_pic_cnt++;
        degrade_pic = 0;

        /* If degrade is to be done in all frames, then do not check further */
        switch(ps_dec->i4_degrade_pics)
        {
            case 4:
            {
                degrade_pic = 1;
                break;
            }
            case 3:
            {
                if(ps_cur_slice->u1_slice_type != I_SLICE) degrade_pic = 1;

                break;
            }
            case 2:
            {
                /* If pic count hits non-degrade interval or it is an islice, then do not
                 * degrade */
                if((ps_cur_slice->u1_slice_type != I_SLICE) &&
                   (ps_dec->i4_degrade_pic_cnt != ps_dec->i4_nondegrade_interval))
                    degrade_pic = 1;

                break;
            }
            case 1:
            {
                /* Check if the current picture is non-ref */
                if(0 == ps_cur_slice->u1_nal_ref_idc)
                {
                    degrade_pic = 1;
                }
                break;
            }
        }
        if(degrade_pic)
        {
            if(ps_dec->i4_degrade_type & 0x2) ps_dec->u4_app_disable_deblk_frm = 1;

            /* MC degrading is done only for non-ref pictures */
            if(0 == ps_cur_slice->u1_nal_ref_idc)
            {
                if(ps_dec->i4_degrade_type & 0x4) ps_dec->i4_mv_frac_mask = 0;

                if(ps_dec->i4_degrade_type & 0x8) ps_dec->i4_mv_frac_mask = 0;
            }
        }
        else
            ps_dec->i4_degrade_pic_cnt = 0;
    }

    {
        dec_err_status_t *ps_err = ps_dec->ps_dec_err_status;
        if((ps_cur_slice->u1_slice_type == I_SLICE) || (ps_cur_slice->u1_slice_type == SI_SLICE))
            ps_err->u1_cur_pic_type = PIC_TYPE_I;
        else
            ps_err->u1_cur_pic_type = PIC_TYPE_UNKNOWN;

        if(ps_err->u1_pic_aud_i == PIC_TYPE_I)
        {
            ps_err->u1_cur_pic_type = PIC_TYPE_I;
            ps_err->u1_pic_aud_i = PIC_TYPE_UNKNOWN;
        }

        if(ps_cur_slice->u1_nal_unit_type == IDR_SLICE_NAL)
        {
            if(ps_err->u1_err_flag) ih264d_reset_ref_bufs(ps_dec->ps_dpb_mgr);
            ps_err->u1_err_flag = ACCEPT_ALL_PICS;
        }
    }

    if(ps_dec->u1_init_dec_flag && ps_dec->s_prev_seq_params.u1_eoseq_pending)
    {
        /* Reset the decoder picture buffers */
        WORD32 j;
        for(j = 0; j < MAX_DISP_BUFS_NEW; j++)
        {
            ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, j, BUF_MGR_REF);
            ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_mv_buf_mgr,
                                  ps_dec->as_buf_id_info_map[j].mv_buf_id, BUF_MGR_REF);
            ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, j, BUF_MGR_IO);
        }

        /* reset the decoder structure parameters related to buffer handling */
        ps_dec->u1_second_field = 0;
        ps_dec->i4_cur_display_seq = 0;

        /********************************************************************/
        /* indicate in the decoder output i4_status that some frames are being */
        /* dropped, so that it resets timestamp and wait for a new sequence */
        /********************************************************************/
        ps_dec->s_prev_seq_params.u1_eoseq_pending = 0;
    }
    ret = isvcd_init_pic(ps_svc_lyr_dec, u2_frame_num, i4_poc, ps_pps);
    if(ret != OK) return ret;

    ps_dec->pv_parse_tu_coeff_data = ps_dec->pv_pic_tu_coeff_data;
    ps_dec->pv_proc_tu_coeff_data = ps_dec->pv_pic_tu_coeff_data;
    ps_dec->ps_nmb_info = ps_dec->ps_frm_mb_info;
    ps_svc_lyr_dec->ps_svc_nmb_info = ps_svc_lyr_dec->ps_svc_frm_mb_info;
    if(ps_dec->u1_separate_parse)
    {
        UWORD32 num_mbs;
        num_mbs = ps_dec->ps_cur_sps->u2_total_num_of_mbs
                  << (1 - ps_dec->ps_cur_sps->u1_frame_mbs_only_flag);

        if(ps_dec->pu1_dec_mb_map)
        {
            memset((void *) ps_dec->pu1_dec_mb_map, 0, num_mbs);
        }

        if(ps_dec->pu1_recon_mb_map)
        {
            memset((void *) ps_dec->pu1_recon_mb_map, 0, num_mbs);
        }

        if(ps_dec->pu2_slice_num_map)
        {
            memset((void *) ps_dec->pu2_slice_num_map, 0, (num_mbs * sizeof(UWORD16)));
        }
    }

    ps_dec->ps_parse_cur_slice = &(ps_dec->ps_dec_slice_buf[0]);
    ps_dec->ps_decode_cur_slice = &(ps_dec->ps_dec_slice_buf[0]);
    ps_dec->ps_computebs_cur_slice = &(ps_dec->ps_dec_slice_buf[0]);
    ps_dec->u2_cur_slice_num = 0;

    /* Initialize all the HP toolsets to zero */
    ps_dec->s_high_profile.u1_scaling_present = 0;
    ps_dec->s_high_profile.u1_transform8x8_present = 0;

    /* Get Next Free Picture */
    if(1 == ps_dec->u4_share_disp_buf)
    {
        UWORD32 i;
        /* Free any buffer that is in the queue to be freed */
        for(i = 0; i < MAX_DISP_BUFS_NEW; i++)
        {
            if(0 == ps_dec->u4_disp_buf_to_be_freed[i]) continue;
            ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, i, BUF_MGR_IO);
            ps_dec->u4_disp_buf_to_be_freed[i] = 0;
            ps_dec->u4_disp_buf_mapping[i] = 0;
        }
    }
    if(!(u1_field_pic_flag && 0 != ps_dec->u1_top_bottom_decoded))
    {
        pic_buffer_t *ps_cur_pic;
        WORD32 cur_pic_buf_id, cur_mv_buf_id;
        col_mv_buf_t *ps_col_mv;
        while(1)
        {
            ps_cur_pic = (pic_buffer_t *) ih264_buf_mgr_get_next_free(
                (buf_mgr_t *) ps_dec->pv_pic_buf_mgr, &cur_pic_buf_id);

            /* In case of IDR slices, if there is no free picture buffer, then release
             * all buffers from display and reference
             */
            if((ps_cur_pic == NULL) && (ps_cur_slice->u1_nal_unit_type == IDR_SLICE_NAL))
            {
                WORD32 j;

                for(j = 0; j < MAX_DISP_BUFS_NEW; j++)
                {
                    ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, j, BUF_MGR_REF);
                    ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_mv_buf_mgr,
                                          ps_dec->as_buf_id_info_map[j].mv_buf_id, BUF_MGR_REF);

                    ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, j, BUF_MGR_IO);
                }
                ps_cur_pic = (pic_buffer_t *) ih264_buf_mgr_get_next_free(
                    (buf_mgr_t *) ps_dec->pv_pic_buf_mgr, &cur_pic_buf_id);
            }
            if(ps_cur_pic == NULL)
            {
                ps_dec->i4_error_code = ERROR_UNAVAIL_PICBUF_T;
                ps_dec->ps_dec_err_status->u1_err_flag |= REJECT_CUR_PIC;
                return ERROR_UNAVAIL_PICBUF_T;
            }
            if(0 == ps_dec->u4_disp_buf_mapping[cur_pic_buf_id])
            {
                break;
            }
        }
        ps_col_mv = (col_mv_buf_t *) ih264_buf_mgr_get_next_free(
            (buf_mgr_t *) ps_dec->pv_mv_buf_mgr, &cur_mv_buf_id);
        if(ps_col_mv == NULL)
        {
            ps_dec->i4_error_code = ERROR_UNAVAIL_MVBUF_T;
            ps_dec->ps_dec_err_status->u1_err_flag |= REJECT_CUR_PIC;
            return ERROR_UNAVAIL_MVBUF_T;
        }

        ps_dec->ps_cur_pic = ps_cur_pic;
        ps_dec->u1_pic_buf_id = cur_pic_buf_id;
        ps_cur_pic->u4_ts = ps_dec->u4_ts;
        memcpy(&ps_cur_pic->s_sei_pic, ps_dec->ps_sei, sizeof(sei));

        ps_cur_pic->u1_mv_buf_id = cur_mv_buf_id;
        ps_dec->as_buf_id_info_map[cur_pic_buf_id].mv_buf_id = cur_mv_buf_id;

        if(ps_dec->u1_enable_mb_info)
        {
            UWORD32 mb_info_map_size = ps_dec->u4_total_mbs << 2;
            ps_dec->as_buf_id_info_map[cur_pic_buf_id].pu1_qp_map =
                ps_dec->pu1_qp_map_base + cur_pic_buf_id * mb_info_map_size;
            ps_dec->as_buf_id_info_map[cur_pic_buf_id].pu1_mb_type_map =
                ps_dec->pu1_mb_type_map_base + cur_pic_buf_id * mb_info_map_size;
            memset(ps_dec->as_buf_id_info_map[cur_pic_buf_id].pu1_qp_map, 0, mb_info_map_size);
            memset(ps_dec->as_buf_id_info_map[cur_pic_buf_id].pu1_mb_type_map, 0, mb_info_map_size);
        }
        ps_cur_pic->pu1_col_zero_flag = (UWORD8 *) ps_col_mv->pv_col_zero_flag;
        ps_cur_pic->ps_mv = (mv_pred_t *) ps_col_mv->pv_mv;
        ps_dec->au1_pic_buf_ref_flag[cur_pic_buf_id] = 0;

        {
            /*make first entry of list0 and list1 point to cur pic,
             *so that if first slice is in error, ref pic struct will have valid
             *entries*/
            ps_dec->ps_ref_pic_buf_lx[0] = ps_dec->ps_dpb_mgr->ps_init_dpb[0];
            ps_dec->ps_ref_pic_buf_lx[1] = ps_dec->ps_dpb_mgr->ps_init_dpb[1];
            *(ps_dec->ps_dpb_mgr->ps_init_dpb[0][0]) = *ps_cur_pic;
            /* Initialize for field reference as well */
            *(ps_dec->ps_dpb_mgr->ps_init_dpb[0][MAX_REF_BUFS]) = *ps_cur_pic;

            *(ps_dec->ps_dpb_mgr->ps_mod_dpb[0][0]) = *ps_cur_pic;
            /* Initialize for field reference as well */
            *(ps_dec->ps_dpb_mgr->ps_mod_dpb[0][MAX_REF_BUFS]) = *ps_cur_pic;
            *(ps_dec->ps_dpb_mgr->ps_init_dpb[1][0]) = *ps_cur_pic;
            /* Initialize for field reference as well */
            *(ps_dec->ps_dpb_mgr->ps_init_dpb[1][MAX_REF_BUFS]) = *ps_cur_pic;
            *(ps_dec->ps_dpb_mgr->ps_mod_dpb[1][0]) = *ps_cur_pic;
            /* Initialize for field reference as well */
            *(ps_dec->ps_dpb_mgr->ps_mod_dpb[1][MAX_REF_BUFS]) = *ps_cur_pic;
        }

        ps_dec->ps_cur_pic->u1_picturetype = u1_field_pic_flag;
        ps_dec->ps_cur_pic->u4_pack_slc_typ = SKIP_NONE;
        H264_DEC_DEBUG_PRINT("got a buffer\n");
    }
    else
    {
        H264_DEC_DEBUG_PRINT("did not get a buffer\n");
    }

    ps_dec->u4_pic_buf_got = 1;

    ps_dec->ps_cur_pic->i4_poc = i4_poc;
    ps_dec->ps_cur_pic->i4_frame_num = u2_frame_num;
    ps_dec->ps_cur_pic->i4_pic_num = u2_frame_num;
    ps_dec->ps_cur_pic->i4_top_field_order_cnt = ps_pps->i4_top_field_order_cnt;
    ps_dec->ps_cur_pic->i4_bottom_field_order_cnt = ps_pps->i4_bottom_field_order_cnt;
    ps_dec->ps_cur_pic->i4_avg_poc = ps_pps->i4_avg_poc;
    ps_dec->ps_cur_pic->u4_time_stamp = ps_dec->u4_pts;

    ps_dec->s_cur_pic = *(ps_dec->ps_cur_pic);
    if(u1_field_pic_flag && u1_bottom_field_flag)
    {
        WORD32 i4_temp_poc;
        WORD32 i4_top_field_order_poc, i4_bot_field_order_poc;
        /* Point to odd lines, since it's bottom field */
        ps_dec->s_cur_pic.pu1_buf1 += ps_dec->s_cur_pic.u2_frm_wd_y;
        ps_dec->s_cur_pic.pu1_buf2 += ps_dec->s_cur_pic.u2_frm_wd_uv;
        ps_dec->s_cur_pic.pu1_buf3 += ps_dec->s_cur_pic.u2_frm_wd_uv;
        ps_dec->s_cur_pic.ps_mv += ((ps_dec->u2_pic_ht * ps_dec->u2_pic_wd) >> 5);
        ps_dec->s_cur_pic.pu1_col_zero_flag += ((ps_dec->u2_pic_ht * ps_dec->u2_pic_wd) >> 5);
        ps_dec->ps_cur_pic->u1_picturetype |= BOT_FLD;
        i4_top_field_order_poc = ps_dec->ps_cur_pic->i4_top_field_order_cnt;
        i4_bot_field_order_poc = ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;
        i4_temp_poc = MIN(i4_top_field_order_poc, i4_bot_field_order_poc);
        ps_dec->ps_cur_pic->i4_avg_poc = i4_temp_poc;
    }

    ps_cur_slice->u1_mbaff_frame_flag = ps_seq->u1_mb_aff_flag && (!u1_field_pic_flag);
    ps_dec->ps_cur_pic->u1_picturetype |= (ps_cur_slice->u1_mbaff_frame_flag << 2);

    ps_dec->ps_cur_mb_row = ps_dec->ps_nbr_mb_row;
    ps_dec->ps_cur_mb_row += 2;
    ps_dec->ps_top_mb_row = ps_dec->ps_nbr_mb_row;
    ps_dec->ps_top_mb_row +=
        ((ps_dec->u2_frm_wd_in_mbs + 2) << (1 - ps_dec->ps_cur_sps->u1_frame_mbs_only_flag));
    // Increment by 2 ,so that left mb (mbaff decrements by 2)  will always be valid
    ps_dec->ps_top_mb_row += 2;
    ps_dec->ps_mv_cur = ps_dec->s_cur_pic.ps_mv;
    ps_dec->ps_mv_top = ps_dec->ps_mv_top_p[0];
    ps_dec->u1_mv_top_p = 0;
    ps_dec->u1_mb_idx = 0;
    ps_dec->ps_mv_left = ps_dec->s_cur_pic.ps_mv;
    ps_dec->u2_total_mbs_coded = 0;
    ps_dec->i4_submb_ofst = -(SUB_BLK_SIZE);
    ps_dec->u4_pred_info_idx = 0;
    ps_dec->u4_pred_info_pkd_idx = 0;
    ps_dec->u4_dma_buf_idx = 0;
    ps_dec->ps_mv = ps_dec->s_cur_pic.ps_mv;
    ps_dec->ps_mv_bank_cur = ps_dec->s_cur_pic.ps_mv;
    ps_dec->pu1_col_zero_flag = ps_dec->s_cur_pic.pu1_col_zero_flag;
    ps_dec->ps_part = ps_dec->ps_parse_part_params;
    ps_dec->i2_prev_slice_mbx = -1;
    ps_dec->i2_prev_slice_mby = 0;
    ps_dec->u2_mv_2mb[0] = 0;
    ps_dec->u2_mv_2mb[1] = 0;
    ps_dec->u1_last_pic_not_decoded = 0;

    ps_dec->u2_cur_slice_num_dec_thread = 0;
    ps_dec->u2_cur_slice_num_bs = 0;
    ps_dec->u4_intra_pred_line_ofst = 0;
    ps_dec->pu1_cur_y_intra_pred_line = ps_dec->pu1_y_intra_pred_line;
    ps_dec->pu1_cur_u_intra_pred_line = ps_dec->pu1_u_intra_pred_line;
    ps_dec->pu1_cur_v_intra_pred_line = ps_dec->pu1_v_intra_pred_line;

    ps_dec->pu1_cur_y_intra_pred_line_base = ps_dec->pu1_y_intra_pred_line;
    ps_dec->pu1_cur_u_intra_pred_line_base = ps_dec->pu1_u_intra_pred_line;
    ps_dec->pu1_cur_v_intra_pred_line_base = ps_dec->pu1_v_intra_pred_line;

    ps_dec->pu1_prev_y_intra_pred_line =
        ps_dec->pu1_y_intra_pred_line + (ps_dec->u2_frm_wd_in_mbs * MB_SIZE);

    ps_dec->pu1_prev_u_intra_pred_line =
        ps_dec->pu1_u_intra_pred_line + ps_dec->u2_frm_wd_in_mbs * BLK8x8SIZE * YUV420SP_FACTOR;
    ps_dec->pu1_prev_v_intra_pred_line =
        ps_dec->pu1_v_intra_pred_line + ps_dec->u2_frm_wd_in_mbs * BLK8x8SIZE;

    ps_dec->ps_deblk_mbn = ps_dec->ps_deblk_pic;
    /* Initialize The Function Pointer Depending Upon the Entropy and MbAff Flag
     */
    {
        if(ps_cur_slice->u1_mbaff_frame_flag)
        {
            ps_dec->pf_compute_bs = ih264d_compute_bs_mbaff;
            ps_dec->pf_mvpred = ih264d_mvpred_mbaff;
            ps_svc_lyr_dec->pf_svc_compute_bs = isvcd_compute_bs_non_mbaff;
        }
        else
        {
            ps_dec->pf_compute_bs = ih264d_compute_bs_non_mbaff;
            ps_svc_lyr_dec->pf_svc_compute_bs = isvcd_compute_bs_non_mbaff;
            ps_dec->u1_cur_mb_fld_dec_flag = ps_cur_slice->u1_field_pic_flag;

            if((ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) &&
               (0 == ps_svc_lyr_dec->u1_base_res_flag))
            {
                ps_svc_lyr_dec->pf_svc_compute_bs = isvcd_compute_bs_non_mbaff_target_lyr;
            }

            if((ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) &&
               (1 == ps_svc_lyr_dec->u1_base_res_flag))
            {
                ps_svc_lyr_dec->pf_svc_compute_bs =
                    isvcd_compute_bs_non_mbaff_target_lyr_no_inter_layer;
            }

            if((ps_svc_lyr_dec->u1_layer_identifier == MEDIAL_ENHANCEMENT_LAYER) &&
               (0 == ps_svc_lyr_dec->u1_base_res_flag))
            {
                ps_svc_lyr_dec->pf_svc_compute_bs = isvcd_compute_bs_non_mbaff_medial_lyr;
            }
        }
    }
    /* Set up the Parameter for DMA transfer */
    {
        UWORD8 u1_field_pic_flag = ps_dec->ps_cur_slice->u1_field_pic_flag;
        UWORD8 u1_mbaff = ps_cur_slice->u1_mbaff_frame_flag;
        UWORD8 uc_lastmbs = (((ps_dec->u2_pic_wd) >> 4) % (ps_dec->u1_recon_mb_grp >> u1_mbaff));
        UWORD16 ui16_lastmbs_widthY =
            (uc_lastmbs ? (uc_lastmbs << 4) : ((ps_dec->u1_recon_mb_grp >> u1_mbaff) << 4));
        UWORD16 ui16_lastmbs_widthUV =
            uc_lastmbs ? (uc_lastmbs << 3) : ((ps_dec->u1_recon_mb_grp >> u1_mbaff) << 3);

        ps_dec->s_tran_addrecon.pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1;
        ps_dec->s_tran_addrecon.pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2;
        ps_dec->s_tran_addrecon.pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3;

        ps_dec->s_tran_addrecon.u2_frm_wd_y = ps_dec->u2_frm_wd_y << u1_field_pic_flag;
        ps_dec->s_tran_addrecon.u2_frm_wd_uv = ps_dec->u2_frm_wd_uv << u1_field_pic_flag;

        if(u1_field_pic_flag)
        {
            ui16_lastmbs_widthY += ps_dec->u2_frm_wd_y;
            ui16_lastmbs_widthUV += ps_dec->u2_frm_wd_uv;
        }

        /* Normal Increment of Pointer */
        ps_dec->s_tran_addrecon.u4_inc_y[0] = ((ps_dec->u1_recon_mb_grp << 4) >> u1_mbaff);
        ps_dec->s_tran_addrecon.u4_inc_uv[0] = ((ps_dec->u1_recon_mb_grp << 4) >> u1_mbaff);

        /* End of Row Increment */
        ps_dec->s_tran_addrecon.u4_inc_y[1] =
            (ui16_lastmbs_widthY + (PAD_LEN_Y_H << 1) +
             ps_dec->s_tran_addrecon.u2_frm_wd_y * ((15 << u1_mbaff) + u1_mbaff));
        ps_dec->s_tran_addrecon.u4_inc_uv[1] =
            (ui16_lastmbs_widthUV + (PAD_LEN_UV_H << 2) +
             ps_dec->s_tran_addrecon.u2_frm_wd_uv * ((15 << u1_mbaff) + u1_mbaff));

        /* Assign picture numbers to each frame/field  */
        /* only once per picture.                      */
        ih264d_assign_pic_num(ps_dec);
        ps_dec->s_tran_addrecon.u2_mv_top_left_inc =
            (ps_dec->u1_recon_mb_grp << 2) - 1 - (u1_mbaff << 2);
        ps_dec->s_tran_addrecon.u2_mv_left_inc = ((ps_dec->u1_recon_mb_grp >> u1_mbaff) - 1)
                                                 << (4 + u1_mbaff);
    }
    /**********************************************************************/
    /* High profile related initialization at pictrue level               */
    /**********************************************************************/
    if((ps_seq->u1_profile_idc == HIGH_PROFILE_IDC) ||
       (ps_seq->u1_profile_idc == SCALABLE_HIGH_PROFILE_IDC) ||
       (ps_seq->u1_profile_idc == SCALABLE_BASELINE_PROFILE_IDC))
    {
        if((ps_seq->i4_seq_scaling_matrix_present_flag) ||
           (ps_pps->i4_pic_scaling_matrix_present_flag))
        {
            ret = ih264d_form_scaling_matrix_picture(ps_seq, ps_pps, ps_dec);
            ps_dec->s_high_profile.u1_scaling_present = 1;
        }
        else
        {
            ret = ih264d_form_default_scaling_matrix(ps_dec);
        }

        if(ps_pps->i4_transform_8x8_mode_flag)
        {
            ps_dec->s_high_profile.u1_transform8x8_present = 1;
        }
    }
    else
    {
        ret = ih264d_form_default_scaling_matrix(ps_dec);
    }

    if(ret != OK) return ret;

    /* required while reading the transform_size_8x8 u4_flag */
    ps_dec->s_high_profile.u1_direct_8x8_inference_flag = ps_seq->u1_direct_8x8_inference_flag;
    ps_dec->s_high_profile.s_cavlc_ctxt = ps_dec->s_cavlc_ctxt;

    ps_dec->i1_recon_in_thread3_flag = 1;
    ps_dec->ps_frame_buf_ip_recon = &ps_dec->s_tran_addrecon;
    if(ps_dec->u1_separate_parse)
    {
        memcpy(&ps_dec->s_tran_addrecon_parse, &ps_dec->s_tran_addrecon, sizeof(tfr_ctxt_t));
    }

    ih264d_init_deblk_tfr_ctxt(ps_dec, &(ps_dec->s_pad_mgr), &(ps_dec->s_tran_addrecon),
                               ps_dec->u2_frm_wd_in_mbs, 0);

    ps_dec->ps_cur_deblk_mb = ps_dec->ps_deblk_pic;
    ps_dec->u4_cur_deblk_mb_num = 0;
    ps_dec->u4_deblk_mb_x = 0;
    ps_dec->u4_deblk_mb_y = 0;
    ps_dec->pu4_wt_ofsts = ps_dec->pu4_wts_ofsts_mat;

    ps_dec->u4_first_slice_in_pic = 0;
    H264_MUTEX_UNLOCK(&ps_dec->process_disp_mutex);
    return OK;
}
/*!
**************************************************************************
* \if Function name : isvcd_parse_decode_slice_ext_nal \endif
*
* \brief
*    Parses a slice extension NAL
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_decode_slice_ext_nal(UWORD8 u1_is_idr_slice, UWORD8 u1_nal_ref_idc,
                                        svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    dec_pic_params_t *ps_pps;
    dec_seq_params_t *ps_seq;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_params_t *ps_cur_slice = NULL;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    pocstruct_t s_tmp_poc = {0};
    WORD32 i_delta_poc[2] = {0};
    WORD32 i4_poc = 0;
    UWORD16 u2_first_mb_in_slice, u2_frame_num;
    UWORD8 u1_field_pic_flag, u1_redundant_pic_cnt = 0, u1_slice_type;
    UWORD32 u4_idr_pic_id = 0;
    UWORD8 u1_bottom_field_flag, u1_pic_order_cnt_type;

    UWORD8 u1_nal_unit_type;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    WORD8 i1_is_end_of_poc;
    WORD32 ret;
    WORD32 prev_slice_err, num_mb_skipped;
    UWORD8 u1_mbaff;
    pocstruct_t *ps_cur_poc;
    UWORD32 u4_temp;
    WORD32 i_temp;
    svc_dec_ctxt_t *psvcd_dec_ctxt;
    dec_struct_t *ps_dec_cur_lyr_minus_1;
    svc_dec_lyr_struct_t *ps_svc_cur_lyr_dec_minus_1;

    ps_cur_slice = ps_dec->ps_cur_slice;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    /* read FirstMbInSlice  and slice type*/
    ps_dec->ps_dpb_cmds->u1_dpb_commands_read_slc = 0;
    u2_first_mb_in_slice = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u2_first_mb_in_slice > (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    /*we currently don not support ASO*/
    if(((u2_first_mb_in_slice << ps_cur_slice->u1_mbaff_frame_flag) <= ps_dec->u2_cur_mb_addr) &&
       (ps_dec->u4_first_slice_in_pic == 0))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    if(ps_dec->u4_first_slice_in_pic == 1)
    {
        if(u2_first_mb_in_slice != 0)
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }

    COPYTHECONTEXT("Slice Header SVC ext: first_mb_in_slice", u2_first_mb_in_slice);

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(u4_temp > 9) return ERROR_INV_SLC_TYPE_T;

    u1_slice_type = u4_temp;
    COPYTHECONTEXT("Slice Header SVC ext: slice_type", (u1_slice_type));
    /* Find Out the Slice Type is 5 to 9 or not then Set the Flag   */
    /* u1_sl_typ_5_9 = 1 .Which tells that all the slices in the Pic*/
    /* will be of same type of current                            */
    if(u1_slice_type > 4)
    {
        u1_slice_type -= 5;
    }

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_PIC_SET_ID) return ERROR_INV_SLICE_HDR_T;
    /* discard slice if pic param is invalid */
    COPYTHECONTEXT("Slice Header SVC ext: pic_parameter_set_id", u4_temp);
    ps_pps = &ps_dec->ps_pps[u4_temp];
    if(FALSE == ps_pps->u1_is_valid)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    /* slices in a layer should have same PPS id*/
    if(UINT32_MAX == ps_svc_lyr_dec->u4_pps_id_for_layer)
    {
        ps_svc_lyr_dec->u4_pps_id_for_layer = u4_temp;
    }
    else if(u4_temp != ps_svc_lyr_dec->u4_pps_id_for_layer)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    ps_seq = ps_pps->ps_sps;
    ps_seq += MAX_NUM_SEQ_PARAMS;
    ps_subset_seq =
        &ps_svc_lyr_dec->ps_subset_sps[MAX_NUM_SEQ_PARAMS + ps_seq->u1_seq_parameter_set_id];

    ps_dec->ps_cur_sps = ps_seq;
    ps_svc_lyr_dec->ps_cur_subset_sps = ps_subset_seq;

    if(!ps_seq) return ERROR_INV_SLICE_HDR_T;
    if(FALSE == ps_seq->u1_is_valid) return ERROR_INV_SLICE_HDR_T;
    if(ps_seq->u1_mb_aff_flag) return ERROR_INV_SLICE_HDR_T;
    if(ps_seq->u1_level_idc > H264_LEVEL_4_2) return ERROR_INV_SLICE_HDR_T;
    if(!ps_seq->u1_frame_mbs_only_flag) return ERROR_INV_SLICE_HDR_T;
    if(OK != isvcd_verify_level(ps_seq->u1_level_idc)) return ERROR_INV_SLICE_HDR_T;

    if(ps_dec->u1_init_dec_flag == 1)
    {
        if(ps_dec->u2_frm_wd_in_mbs != ps_seq->u2_frm_wd_in_mbs) return ERROR_INV_SLICE_HDR_T;
        if(ps_dec->u2_frm_ht_in_mbs != ps_seq->u2_frm_ht_in_mbs) return ERROR_INV_SLICE_HDR_T;
    }

    ps_dec->i4_reorder_depth = ps_subset_seq->i4_reorder_depth;

    ps_dec->u2_disp_height = ps_subset_seq->u2_disp_height;
    ps_dec->u2_disp_width = ps_subset_seq->u2_disp_width;

    if(ps_svc_lyr_dec->u1_layer_id > 0)
    {
        psvcd_dec_ctxt = ps_svc_lyr_dec->ps_svcd_ctxt;
        ps_svc_cur_lyr_dec_minus_1 =
            &psvcd_dec_ctxt->ps_svc_dec_lyr[ps_svc_lyr_dec->u1_layer_id - 1];

        ps_dec_cur_lyr_minus_1 = &ps_svc_cur_lyr_dec_minus_1->s_dec;

        if((ps_dec_cur_lyr_minus_1->u2_pic_wd > ps_subset_seq->u2_pic_wd) ||
           (ps_dec_cur_lyr_minus_1->u2_pic_ht > ps_subset_seq->u2_pic_ht))
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }

    ps_dec->u2_pic_wd = ps_subset_seq->u2_pic_wd;
    ps_dec->u2_pic_ht = ps_subset_seq->u2_pic_ht;
    ps_dec->u4_total_mbs = ps_seq->u2_total_num_of_mbs << (1 - ps_seq->u1_frame_mbs_only_flag);

    /* Determining the Width and Height of Frame from that of Picture */
    ps_dec->u2_frm_wd_y = ps_subset_seq->u2_frm_wd_y;
    ps_dec->u2_frm_ht_y = ps_subset_seq->u2_frm_ht_y;

    ps_dec->u2_frm_wd_uv = ps_subset_seq->u2_frm_wd_uv;
    ps_dec->u2_frm_ht_uv = ps_subset_seq->u2_frm_ht_uv;

    ps_dec->s_pad_mgr.u1_pad_len_y_v = ps_subset_seq->u1_pad_len_y_v;
    ps_dec->s_pad_mgr.u1_pad_len_cr_v = ps_subset_seq->u1_pad_len_cr_v;

    ps_dec->u2_frm_wd_in_mbs = ps_seq->u2_frm_wd_in_mbs;
    ps_dec->u2_frm_ht_in_mbs = ps_seq->u2_frm_ht_in_mbs;

    ps_dec->u2_crop_offset_y = ps_subset_seq->u2_crop_offset_y;
    ps_dec->u2_crop_offset_uv = ps_subset_seq->u2_crop_offset_uv;

    /* Get the frame num */
    u2_frame_num = ih264d_get_bits_h264(ps_bitstrm, ps_seq->u1_bits_in_frm_num);

    COPYTHECONTEXT("Slice Header SVC ext: frame_num", u2_frame_num);
    if(!ps_dec->u1_first_slice_in_stream && ps_dec->u4_first_slice_in_pic)
    {
        pocstruct_t *ps_prev_poc = &ps_dec->s_prev_pic_poc;
        pocstruct_t *ps_cur_poc = &ps_dec->s_cur_pic_poc;

        ps_dec->u2_mbx = 0xffff;
        ps_dec->u2_mby = 0;

        if((0 == u1_is_idr_slice) && ps_cur_slice->u1_nal_ref_idc)
            ps_dec->u2_prev_ref_frame_num = ps_cur_slice->u2_frame_num;

        if(u1_is_idr_slice || ps_cur_slice->u1_mmco_equalto5) ps_dec->u2_prev_ref_frame_num = 0;

        if(ps_dec->ps_cur_sps->u1_gaps_in_frame_num_value_allowed_flag)
        {
            isvcd_decode_gaps_in_frame_num(ps_dec, u2_frame_num);
        }

        ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst;
        ps_prev_poc->u2_frame_num = ps_cur_poc->u2_frame_num;
        ps_prev_poc->u1_mmco_equalto5 = ps_cur_slice->u1_mmco_equalto5;
        if(ps_cur_slice->u1_nal_ref_idc)
        {
            ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb;
            ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb;
            ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom;
            ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0];
            ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1];
            ps_prev_poc->u1_bot_field = ps_cur_poc->u1_bot_field;
        }

        ps_dec->u2_total_mbs_coded = 0;
    }
    /* Get the field related flags  */
    if(!ps_seq->u1_frame_mbs_only_flag)
    {
        u1_field_pic_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: field_pic_flag", u1_field_pic_flag);
        u1_bottom_field_flag = 0;

        if(u1_field_pic_flag)
        {
            ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan_fld;
            u1_bottom_field_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("Slice Header SVC ext: bottom_field_flag", u1_bottom_field_flag);
        }
        else
        {
            ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan;
        }
    }
    else
    {
        u1_field_pic_flag = 0;
        u1_bottom_field_flag = 0;
        ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan;
    }

    u1_nal_unit_type = SLICE_NAL;
    if(u1_is_idr_slice)
    {
        u1_nal_unit_type = IDR_SLICE_NAL;
        u4_idr_pic_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_idr_pic_id > 65535) return ERROR_INV_SLICE_HDR_T;
        COPYTHECONTEXT("Slice Header SVC ext:  ", u4_idr_pic_id);
    }

    /* read delta pic order count information*/
    i_delta_poc[0] = i_delta_poc[1] = 0;
    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    u1_pic_order_cnt_type = ps_seq->u1_pic_order_cnt_type;
    if(u1_pic_order_cnt_type == 0)
    {
        i_temp = ih264d_get_bits_h264(ps_bitstrm, ps_seq->u1_log2_max_pic_order_cnt_lsb_minus);
        if(i_temp < 0 || i_temp >= ps_seq->i4_max_pic_order_cntLsb) return ERROR_INV_SLICE_HDR_T;
        s_tmp_poc.i4_pic_order_cnt_lsb = i_temp;
        COPYTHECONTEXT("Slice Header SVC ext: pic_order_cnt_lsb", s_tmp_poc.i4_pic_order_cnt_lsb);

        if((ps_pps->u1_pic_order_present_flag == 1) && (!u1_field_pic_flag))
        {
            s_tmp_poc.i4_delta_pic_order_cnt_bottom = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: delta_pic_order_cnt_bottom",
                           s_tmp_poc.i4_delta_pic_order_cnt_bottom);
        }
    }

    s_tmp_poc.i4_delta_pic_order_cnt[0] = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[1] = 0;
    if(u1_pic_order_cnt_type == 1 && (!ps_seq->u1_delta_pic_order_always_zero_flag))
    {
        s_tmp_poc.i4_delta_pic_order_cnt[0] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("Slice Header SVC ext: delta_pic_order_cnt[0]",
                       s_tmp_poc.i4_delta_pic_order_cnt[0]);

        if(ps_pps->u1_pic_order_present_flag && !u1_field_pic_flag)
        {
            s_tmp_poc.i4_delta_pic_order_cnt[1] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: delta_pic_order_cnt[1]",
                           s_tmp_poc.i4_delta_pic_order_cnt[1]);
        }
    }

    if(ps_pps->u1_redundant_pic_cnt_present_flag)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > MAX_REDUNDANT_PIC_CNT) return ERROR_INV_SLICE_HDR_T;
        u1_redundant_pic_cnt = u4_temp;
        COPYTHECONTEXT("Slice Header SVC ext: redundant_pic_cnt", u1_redundant_pic_cnt);
    }

    /*--------------------------------------------------------------------*/
    /* Check if the slice is part of new picture                          */
    /*--------------------------------------------------------------------*/
    /* First slice of a picture is always considered as part of new picture */
    i1_is_end_of_poc = 1;
    ps_dec->ps_dec_err_status->u1_err_flag &= MASK_REJECT_CUR_PIC;

    if(ps_dec->u4_first_slice_in_pic == 0)
    {
        i1_is_end_of_poc =
            ih264d_is_end_of_pic(u2_frame_num, u1_nal_ref_idc, &s_tmp_poc, &ps_dec->s_cur_pic_poc,
                                 ps_cur_slice, u1_pic_order_cnt_type, u1_nal_unit_type,
                                 u4_idr_pic_id, u1_field_pic_flag, u1_bottom_field_flag);
        if(i1_is_end_of_poc)
        {
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_INCOMPLETE_FRAME;
        }
    }

    /*--------------------------------------------------------------------*/
    /* Check for error in slice and parse the missing/corrupted MB's      */
    /* as skip-MB's in an inserted P-slice                                */
    /*--------------------------------------------------------------------*/
    u1_mbaff = ps_seq->u1_mb_aff_flag && (!u1_field_pic_flag);
    prev_slice_err = 0;

    if(i1_is_end_of_poc || ps_dec->u1_first_slice_in_stream)
    {
        /* If the current slice is not a field or frame number of the current
         * slice doesn't match with previous slice, and decoder is expecting
         * to decode a field i.e. ps_dec->u1_top_bottom_decoded is not 0 and
         * is not (TOP_FIELD_ONLY | BOT_FIELD_ONLY), treat it as a dangling
         * field */
        if((u1_field_pic_flag == 0 || u2_frame_num != ps_dec->u2_prv_frame_num) &&
           ps_dec->u1_top_bottom_decoded != 0 &&
           ps_dec->u1_top_bottom_decoded != (TOP_FIELD_ONLY | BOT_FIELD_ONLY))
        {
            ps_dec->u1_dangling_field = 1;
            if(ps_dec->u4_first_slice_in_pic)
            {
                // first slice - dangling field
                prev_slice_err = 1;
            }
            else
            {
                // last slice - dangling field
                prev_slice_err = 2;
            }

            if(ps_dec->u1_top_bottom_decoded == TOP_FIELD_ONLY)
                ps_cur_slice->u1_bottom_field_flag = 1;
            else
                ps_cur_slice->u1_bottom_field_flag = 0;

            num_mb_skipped =
                (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &ps_dec->s_cur_pic_poc;

            u1_is_idr_slice = ps_cur_slice->u1_nal_unit_type == IDR_SLICE_NAL;
        }
        else if(ps_dec->u4_first_slice_in_pic)
        {
            if(u2_first_mb_in_slice > 0)
            {
                // first slice - missing/header corruption
                prev_slice_err = 1;
                num_mb_skipped = u2_first_mb_in_slice << u1_mbaff;
                ps_cur_poc = &s_tmp_poc;

                // initializing slice parameters
                ps_cur_slice->u4_idr_pic_id = u4_idr_pic_id;
                ps_cur_slice->u1_field_pic_flag = u1_field_pic_flag;
                ps_cur_slice->u1_bottom_field_flag = u1_bottom_field_flag;
                ps_cur_slice->i4_pic_order_cnt_lsb = s_tmp_poc.i4_pic_order_cnt_lsb;
                ps_cur_slice->u1_nal_unit_type = u1_nal_unit_type;
                ps_cur_slice->u1_redundant_pic_cnt = u1_redundant_pic_cnt;
                ps_cur_slice->u1_nal_ref_idc = u1_nal_ref_idc;
                ps_cur_slice->u1_pic_order_cnt_type = u1_pic_order_cnt_type;
                ps_cur_slice->u1_mbaff_frame_flag = ps_seq->u1_mb_aff_flag && (!u1_field_pic_flag);
            }
        }
        else
        {
            /* since i1_is_end_of_poc is set ,means new frame num is encountered. so
             * conceal the current frame completely */
            prev_slice_err = 2;
            num_mb_skipped =
                (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &s_tmp_poc;
        }
    }
    else
    {
        if((u2_first_mb_in_slice << u1_mbaff) > ps_dec->u2_total_mbs_coded)
        {
            // previous slice - missing/corruption
            prev_slice_err = 2;
            num_mb_skipped = (u2_first_mb_in_slice << u1_mbaff) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &s_tmp_poc;
        }
        else if((u2_first_mb_in_slice << u1_mbaff) < ps_dec->u2_total_mbs_coded)
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }
    if(prev_slice_err)
    {
        ret = isvcd_mark_err_slice_skip((svc_dec_lyr_struct_t *) ps_dec, num_mb_skipped,
                                        u1_is_idr_slice, u2_frame_num, ps_cur_poc, prev_slice_err);

        if(ps_dec->u1_dangling_field == 1)
        {
            ps_dec->u1_second_field = 1 - ps_dec->u1_second_field;
            ps_dec->u1_first_slice_in_stream = 0;
            ps_dec->u1_top_bottom_decoded = TOP_FIELD_ONLY | BOT_FIELD_ONLY;
            return ERROR_DANGLING_FIELD_IN_PIC;
        }

        if(prev_slice_err == 2)
        {
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_INCOMPLETE_FRAME;
        }

        if(ps_dec->u2_total_mbs_coded >= ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs)
        {
            /* return if all MBs in frame are parsed*/
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_IN_LAST_SLICE_OF_PIC;
        }

        if(ps_dec->ps_dec_err_status->u1_err_flag & REJECT_CUR_PIC)
        {
            ih264d_err_pic_dispbuf_mgr(ps_dec);
            return ERROR_NEW_FRAME_EXPECTED;
        }

        if(ret != OK) return ret;

        i1_is_end_of_poc = 0;
    }

    if(u1_field_pic_flag)
    {
        ps_dec->u2_prv_frame_num = u2_frame_num;
    }

    if(ps_cur_slice->u1_mmco_equalto5)
    {
        WORD32 i4_temp_poc;
        WORD32 i4_top_field_order_poc, i4_bot_field_order_poc;
        WORD64 i8_result;
        if(!ps_cur_slice->u1_field_pic_flag)  // or a complementary field pair
        {
            i4_top_field_order_poc = ps_dec->ps_cur_pic->i4_top_field_order_cnt;
            i4_bot_field_order_poc = ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;
            i4_temp_poc = MIN(i4_top_field_order_poc, i4_bot_field_order_poc);
        }
        else if(!ps_cur_slice->u1_bottom_field_flag)
            i4_temp_poc = ps_dec->ps_cur_pic->i4_top_field_order_cnt;
        else
            i4_temp_poc = ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;

        i8_result = (WORD64) i4_temp_poc - ps_dec->ps_cur_pic->i4_top_field_order_cnt;
        if(IS_OUT_OF_RANGE_S32(i8_result))
        {
            return ERROR_INV_POC;
        }
        ps_dec->ps_cur_pic->i4_top_field_order_cnt = (WORD32) i8_result;
        i8_result = (WORD64) i4_temp_poc - ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;
        if(IS_OUT_OF_RANGE_S32(i8_result))
        {
            return ERROR_INV_POC;
        }
        ps_dec->ps_cur_pic->i4_bottom_field_order_cnt = (WORD32) i8_result;
        ps_dec->ps_cur_pic->i4_poc = i4_temp_poc;
        ps_dec->ps_cur_pic->i4_avg_poc = i4_temp_poc;
    }
    if(ps_dec->u4_first_slice_in_pic)
    {
        ret = isvcd_decode_pic_order_cnt(u1_is_idr_slice, u2_frame_num, &ps_dec->s_prev_pic_poc,
                                         &s_tmp_poc, ps_cur_slice, ps_pps, u1_nal_ref_idc,
                                         u1_bottom_field_flag, u1_field_pic_flag, &i4_poc, ps_dec);
        if(ret != OK) return ret;
        /* Display seq no calculations */
        if(i4_poc >= ps_dec->i4_max_poc) ps_dec->i4_max_poc = i4_poc;
        /* IDR Picture or POC wrap around */
        if(i4_poc == 0)
        {
            WORD64 i8_temp;
            i8_temp = (WORD64) ps_dec->i4_prev_max_display_seq + ps_dec->i4_max_poc +
                      ps_dec->u1_max_dec_frame_buffering + 1;
            /*If i4_prev_max_display_seq overflows integer range, reset it */
            ps_dec->i4_prev_max_display_seq = IS_OUT_OF_RANGE_S32(i8_temp) ? 0 : i8_temp;
            ps_dec->i4_max_poc = 0;
        }
    }

    /* Increment only if the current slice has atleast 1 more MB */
    if(ps_dec->u4_first_slice_in_pic == 0 &&
       (ps_dec->ps_parse_cur_slice->u4_first_mb_in_slice <
        (UWORD32) (ps_dec->u2_total_mbs_coded >> ps_dec->ps_cur_slice->u1_mbaff_frame_flag)))
    {
        ps_dec->ps_parse_cur_slice++;
        ps_dec->u2_cur_slice_num++;
        // in the case of single core increment ps_decode_cur_slice
        if(ps_dec->u1_separate_parse == 0)
        {
            ps_dec->ps_decode_cur_slice++;
        }
    }

    ps_dec->u1_slice_header_done = 0;

    /*--------------------------------------------------------------------*/
    /* Copy the values read from the bitstream to the slice header and then*/
    /* If the slice is first slice in picture, then do Start of Picture   */
    /* processing.                                                        */
    /*--------------------------------------------------------------------*/
    ps_cur_slice->i4_delta_pic_order_cnt[0] = i_delta_poc[0];
    ps_cur_slice->i4_delta_pic_order_cnt[1] = i_delta_poc[1];
    ps_cur_slice->u4_idr_pic_id = u4_idr_pic_id;
    ps_cur_slice->u2_first_mb_in_slice = u2_first_mb_in_slice;
    ps_cur_slice->u1_field_pic_flag = u1_field_pic_flag;
    ps_cur_slice->u1_bottom_field_flag = u1_bottom_field_flag;
    ps_cur_slice->u1_slice_type = u1_slice_type;
    ps_cur_slice->i4_pic_order_cnt_lsb = s_tmp_poc.i4_pic_order_cnt_lsb;

    ps_cur_slice->u1_nal_unit_type = u1_nal_unit_type;
    ps_cur_slice->u1_redundant_pic_cnt = u1_redundant_pic_cnt;
    ps_cur_slice->u1_nal_ref_idc = u1_nal_ref_idc;
    ps_cur_slice->u1_pic_order_cnt_type = u1_pic_order_cnt_type;

    if(ps_seq->u1_frame_mbs_only_flag)
        ps_cur_slice->u1_direct_8x8_inference_flag = ps_seq->u1_direct_8x8_inference_flag;
    else
        ps_cur_slice->u1_direct_8x8_inference_flag = 1;

    if(0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id)
    {
        if(B_SLICE == u1_slice_type)
        {
            ps_cur_slice->u1_direct_spatial_mv_pred_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("Slice Header SVC ext: direct_spatial_mv_pred_flag",
                           ps_cur_slice->u1_direct_spatial_mv_pred_flag);

            if(ps_cur_slice->u1_direct_spatial_mv_pred_flag)
                ps_cur_slice->pf_decodeDirect = isvcd_decode_spatial_direct;
            else
                ps_cur_slice->pf_decodeDirect = ih264d_decode_temporal_direct;
            if(!((ps_seq->u1_mb_aff_flag) && (!u1_field_pic_flag)))
                ps_dec->pf_mvpred = ih264d_mvpred_nonmbaffB;
        }
        else
        {
            if(!((ps_seq->u1_mb_aff_flag) && (!u1_field_pic_flag))) /*check if this is valid here */
                ps_dec->pf_mvpred = ih264d_mvpred_nonmbaff;
        }
    }

    if(ps_dec->u4_first_slice_in_pic)
    {
        if(u2_first_mb_in_slice == 0)
        {
            ret = isvcd_start_of_pic(ps_svc_lyr_dec, i4_poc, &s_tmp_poc, u2_frame_num, ps_pps);
            if(ret != OK) return ret;
            /*inter layer buffer intialization */
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
                ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start;
            ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_cur_mb =
                ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_base;
        }

        ps_dec->u4_output_present = 0;

        {
            ih264d_get_next_display_field(ps_dec, ps_dec->ps_out_buffer, &(ps_dec->s_disp_op));
            /* If error code is non-zero then there is no buffer available for
            display, hence avoid format conversion */

            if(0 != ps_dec->s_disp_op.u4_error_code)
            {
                ps_dec->u4_output_present = 0;
                ps_dec->u4_fmt_conv_cur_row = ps_dec->s_disp_frame_info.u4_y_ht;
            }
            else
                ps_dec->u4_output_present = 1;
        }
        ret = isvcd_parse_interlayer_resamp_func_init(ps_svc_lyr_dec, u2_first_mb_in_slice);
        if(ret != OK)
        {
            return ERROR_CORRUPTED_SLICE;
        }
        if((ps_dec->u1_separate_parse == 1) &&
           (ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) && (ps_svc_lyr_dec->u1_res_init_done == 1))
        {
            if(ps_dec->u4_dec_thread_created == 0)
            {
                ithread_create(ps_dec->pv_dec_thread_handle, NULL,
                               (void *) isvcd_decode_picture_thread, (void *) ps_dec);

                ps_dec->u4_dec_thread_created = 1;
            }
#ifdef KEEP_THREADS_ACTIVE
            ret = ithread_mutex_lock(ps_dec->apv_proc_start_mutex[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);

            ps_dec->ai4_process_start[0] = PROC_START;
            ret = ithread_cond_signal(ps_dec->apv_proc_start_condition[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);

            ret = ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);
#endif
#ifdef KEEP_THREADS_ACTIVE
            if(ps_dec->u4_bs_deblk_thread_created)
            {
                ret = ithread_mutex_lock(ps_dec->apv_proc_start_mutex[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ps_dec->ai4_process_start[1] = PROC_START;
                ret = ithread_cond_signal(ps_dec->apv_proc_start_condition[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ret = ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);
            }
#endif
        }
    }

    /* INITIALIZATION of fn ptrs for MC and formMbPartInfo functions */
    {
        UWORD8 uc_nofield_nombaff;

        uc_nofield_nombaff =
            ((ps_dec->ps_cur_slice->u1_field_pic_flag == 0) &&
             (ps_dec->ps_cur_slice->u1_mbaff_frame_flag == 0) && (u1_slice_type != B_SLICE) &&
             (ps_dec->ps_cur_pps->u1_wted_pred_flag == 0));

        /* Initialise MC and formMbPartInfo fn ptrs one time based on profile_idc */

        if(uc_nofield_nombaff)
        {
            ps_dec->p_form_mb_part_info = ih264d_form_mb_part_info_bp;
            ps_dec->p_motion_compensate = ih264d_motion_compensate_bp;
        }
        else
        {
            ps_dec->p_form_mb_part_info = ih264d_form_mb_part_info_mp;
            ps_dec->p_motion_compensate = ih264d_motion_compensate_mp;
        }
    }

    /*
     * Decide whether to decode the current picture or not
     */
    {
        dec_err_status_t *ps_err = ps_dec->ps_dec_err_status;
        if(ps_err->u4_frm_sei_sync == u2_frame_num)
        {
            ps_err->u1_err_flag = ACCEPT_ALL_PICS;
            ps_err->u4_frm_sei_sync = SYNC_FRM_DEFAULT;
        }
        ps_err->u4_cur_frm = u2_frame_num;
    }

    /* Decision for decoding if the picture is to be skipped */
    {
        WORD32 i4_skip_b_pic, i4_skip_p_pic;

        i4_skip_b_pic = (ps_dec->u4_skip_frm_mask & B_SLC_BIT) && (B_SLICE == u1_slice_type) &&
                        (0 == u1_nal_ref_idc);

        i4_skip_p_pic = (ps_dec->u4_skip_frm_mask & P_SLC_BIT) && (P_SLICE == u1_slice_type) &&
                        (0 == u1_nal_ref_idc);

        /**************************************************************/
        /* Skip the B picture if skip mask is set for B picture and   */
        /* Current B picture is a non reference B picture or there is */
        /* no user for reference B picture                            */
        /**************************************************************/
        if(i4_skip_b_pic)
        {
            ps_dec->ps_cur_pic->u4_pack_slc_typ |= B_SLC_BIT;
            /* Don't decode the picture in SKIP-B mode if that picture is B */
            /* and also it is not to be used as a reference picture         */
            ps_dec->u1_last_pic_not_decoded = 1;

            return OK;
        }
        /**************************************************************/
        /* Skip the P picture if skip mask is set for P picture and   */
        /* Current P picture is a non reference P picture or there is */
        /* no user for reference P picture                            */
        /**************************************************************/
        if(i4_skip_p_pic)
        {
            ps_dec->ps_cur_pic->u4_pack_slc_typ |= P_SLC_BIT;
            /* Don't decode the picture in SKIP-P mode if that picture is P */
            /* and also it is not to be used as a reference picture         */
            ps_dec->u1_last_pic_not_decoded = 1;

            return OK;
        }
    }

    {
        UWORD16 u2_mb_x, u2_mb_y;

        ps_dec->i4_submb_ofst =
            ((u2_first_mb_in_slice << ps_cur_slice->u1_mbaff_frame_flag) * SUB_BLK_SIZE) -
            SUB_BLK_SIZE;
        if(u2_first_mb_in_slice)
        {
            UWORD8 u1_mb_aff;
            UWORD8 u1_field_pic;
            UWORD16 u2_frm_wd_in_mbs;
            u2_frm_wd_in_mbs = ps_seq->u2_frm_wd_in_mbs;
            u1_mb_aff = ps_cur_slice->u1_mbaff_frame_flag;
            u1_field_pic = ps_cur_slice->u1_field_pic_flag;

            {
                UWORD32 x_offset;
                UWORD32 y_offset;
                UWORD32 u4_frame_stride;
                tfr_ctxt_t *ps_trns_addr;  // = &ps_dec->s_tran_addrecon_parse;

                if(ps_dec->u1_separate_parse)
                {
                    ps_trns_addr = &ps_dec->s_tran_addrecon_parse;
                }
                else
                {
                    ps_trns_addr = &ps_dec->s_tran_addrecon;
                }
                u2_mb_x = MOD(u2_first_mb_in_slice, u2_frm_wd_in_mbs);
                u2_mb_y = DIV(u2_first_mb_in_slice, u2_frm_wd_in_mbs);

                u2_mb_y <<= u1_mb_aff;

                if((u2_mb_x > u2_frm_wd_in_mbs - 1) || (u2_mb_y > ps_dec->u2_frm_ht_in_mbs - 1))
                {
                    return ERROR_CORRUPTED_SLICE;
                }

                u4_frame_stride = ps_dec->u2_frm_wd_y << u1_field_pic;
                x_offset = u2_mb_x << 4;
                y_offset = (u2_mb_y * u4_frame_stride) << 4;

                ps_trns_addr->pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1 + x_offset + y_offset;

                u4_frame_stride = ps_dec->u2_frm_wd_uv << u1_field_pic;
                x_offset >>= 1;
                y_offset = (u2_mb_y * u4_frame_stride) << 3;

                x_offset *= YUV420SP_FACTOR;

                ps_trns_addr->pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2 + x_offset + y_offset;
                ps_trns_addr->pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3 + x_offset + y_offset;

                ps_trns_addr->pu1_mb_y = ps_trns_addr->pu1_dest_y;
                ps_trns_addr->pu1_mb_u = ps_trns_addr->pu1_dest_u;
                ps_trns_addr->pu1_mb_v = ps_trns_addr->pu1_dest_v;

                // assign the deblock structure pointers to start of slice
                if(ps_dec->u1_separate_parse == 1)
                {
                    ps_dec->ps_deblk_mbn =
                        ps_dec->ps_deblk_pic + (u2_first_mb_in_slice << u1_mb_aff);
                }
                else
                {
                    ps_dec->ps_deblk_mbn =
                        ps_dec->ps_deblk_pic + (u2_first_mb_in_slice << u1_mb_aff);
                }

                ps_dec->u2_cur_mb_addr = (u2_first_mb_in_slice << u1_mb_aff);

                ps_dec->ps_mv_cur =
                    ps_dec->s_cur_pic.ps_mv + ((u2_first_mb_in_slice << u1_mb_aff) << 4);
            }
        }
        else
        {
            tfr_ctxt_t *ps_trns_addr;

            if(ps_dec->u1_separate_parse)
            {
                ps_trns_addr = &ps_dec->s_tran_addrecon_parse;
            }
            else
            {
                ps_trns_addr = &ps_dec->s_tran_addrecon;
            }

            u2_mb_x = 0xffff;
            u2_mb_y = 0;
            // assign the deblock structure pointers to start of slice
            ps_dec->u2_cur_mb_addr = 0;
            ps_dec->ps_deblk_mbn = ps_dec->ps_deblk_pic;
            ps_dec->ps_mv_cur = ps_dec->s_cur_pic.ps_mv;
            ps_trns_addr->pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1;
            ps_trns_addr->pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2;
            ps_trns_addr->pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3;

            ps_trns_addr->pu1_mb_y = ps_trns_addr->pu1_dest_y;
            ps_trns_addr->pu1_mb_u = ps_trns_addr->pu1_dest_u;
            ps_trns_addr->pu1_mb_v = ps_trns_addr->pu1_dest_v;
        }

        ps_dec->ps_part = ps_dec->ps_parse_part_params;

        ps_dec->u2_mbx = (MOD(u2_first_mb_in_slice - 1, ps_seq->u2_frm_wd_in_mbs));
        ps_dec->u2_mby = (DIV(u2_first_mb_in_slice - 1, ps_seq->u2_frm_wd_in_mbs));
        ps_dec->u2_mby <<= ps_cur_slice->u1_mbaff_frame_flag;
        ps_dec->i2_prev_slice_mbx = ps_dec->u2_mbx;
        ps_dec->i2_prev_slice_mby = ps_dec->u2_mby;
    }

    /* RBSP stop bit is used for CABAC decoding*/
    ps_bitstrm->u4_max_ofst += ps_dec->ps_cur_pps->u1_entropy_coding_mode;

    ps_dec->u1_B = (u1_slice_type == B_SLICE);
    ps_dec->u4_next_mb_skip = 0;

    ps_dec->ps_parse_cur_slice->u4_first_mb_in_slice = ps_dec->ps_cur_slice->u2_first_mb_in_slice;
    ps_dec->ps_parse_cur_slice->slice_type = ps_dec->ps_cur_slice->u1_slice_type;

    ps_dec->u4_start_recon_deblk = 1;
    {
        WORD32 num_entries;
        WORD32 size;
        UWORD8 *pu1_buf;

        num_entries = MAX_FRAMES;
        if((1 >= ps_dec->ps_cur_sps->u1_num_ref_frames) && (0 == ps_dec->i4_display_delay))
        {
            num_entries = 1;
        }
        num_entries = ((2 * num_entries) + 1);
        num_entries *= 2;

        size = num_entries * sizeof(void *);
        size += PAD_MAP_IDX_POC * sizeof(void *);

        pu1_buf = (UWORD8 *) ps_dec->pv_map_ref_idx_to_poc_buf;
        pu1_buf += size * ps_dec->u2_cur_slice_num;
        ps_dec->ps_parse_cur_slice->ppv_map_ref_idx_to_poc = (void *) pu1_buf;
    }

    if(ps_dec->u1_separate_parse)
    {
        ps_dec->ps_parse_cur_slice->pv_tu_coeff_data_start = ps_dec->pv_parse_tu_coeff_data;
    }
    else
    {
        ps_dec->pv_proc_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;
    }

    ret = ih264d_fix_error_in_dpb(ps_dec);
    if(ret < 0) return ERROR_DBP_MANAGER_T;

    /*Default initializing default values for some parameters*/
    ps_svc_slice_params->u1_slice_skip_flag = 0;
    ps_svc_slice_params->u1_adaptive_base_mode_flag = 0;
    ps_svc_slice_params->u1_default_base_mode_flag = 0;
    ps_svc_slice_params->u1_adaptive_motion_prediction_flag = 0;
    ps_svc_slice_params->u1_default_motion_prediction_flag = 0;
    ps_svc_slice_params->u1_adaptive_residual_prediction_flag = 0;
    ps_svc_slice_params->u1_default_residual_prediction_flag = 0;

    if(u1_slice_type == I_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= I_SLC_BIT;

        ret = isvcd_parse_eislice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        if(ps_dec->i4_pic_type != B_SLICE && ps_dec->i4_pic_type != P_SLICE)
            ps_dec->i4_pic_type = I_SLICE;
    }
    else if(u1_slice_type == P_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= P_SLC_BIT;
        ret = isvcd_parse_epslice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        if(ps_dec->i4_pic_type != B_SLICE) ps_dec->i4_pic_type = P_SLICE;
    }
    else if(u1_slice_type == B_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= B_SLC_BIT;
        ret = isvcd_parse_ebslice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        ps_dec->i4_pic_type = B_SLICE;
    }
    else
        return ERROR_INV_SLC_TYPE_T;

    if(ps_dec->u1_slice_header_done)
    {
        /* set to zero to indicate a valid slice has been decoded */
        ps_dec->u1_first_slice_in_stream = 0;
    }

    if(ret != OK) return ret;

    if(u1_nal_ref_idc != 0)
    {
        if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
        {
            memcpy((void *) ps_dec->ps_dpb_cmds, (void *) (&(ps_dec->s_dpb_cmds_scratch)),
                   sizeof(dpb_commands_t));
        }
    }

    /* storing last Mb X and MbY of the slice */
    ps_dec->i2_prev_slice_mbx = ps_dec->u2_mbx;
    ps_dec->i2_prev_slice_mby = ps_dec->u2_mby;

    /* End of Picture detection */
    if(ps_dec->u2_total_mbs_coded >= (ps_seq->u2_max_mb_addr + 1))
    {
        ps_dec->u1_pic_decode_done = 1;
    }

    {
        dec_err_status_t *ps_err = ps_dec->ps_dec_err_status;
        if((ps_err->u1_err_flag & REJECT_PB_PICS) && (ps_err->u1_cur_pic_type == PIC_TYPE_I))
        {
            ps_err->u1_err_flag = ACCEPT_ALL_PICS;
        }
    }

    PRINT_BIN_BIT_RATIO(ps_dec)

    return ret;
}

/*!
**************************************************************************
* \if Function name : isvcd_set_default_slice_header_ext \endif
*
* \brief
*    sets the default values for the svc slice header attr
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_set_default_slice_header_ext(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD32 i_status = OK;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_seq_params_t *ps_seq;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    dec_subset_seq_params_t *ps_sps_svc_ext = NULL;
    ps_seq = ps_pps->ps_sps;
    ps_seq += MAX_NUM_SEQ_PARAMS;
    ps_subset_seq =
        &ps_svc_lyr_dec->ps_subset_sps[MAX_NUM_SEQ_PARAMS + ps_seq->u1_seq_parameter_set_id];
    ps_sps_svc_ext = &ps_subset_seq->s_sps_svc_ext;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    if(0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id)
    {
        ps_svc_slice_params->u1_ref_layer_chroma_phase_y_plus1 =
            ps_sps_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1;

        ps_svc_slice_params->u1_ref_layer_chroma_phase_x_plus1_flag =
            ps_sps_svc_ext->u1_seq_ref_layer_chroma_phase_x_plus1_flag;
    }

    ps_svc_slice_params->u4_ref_layer_dq_id = UINT32_MAX;
    ps_svc_slice_params->u4_disable_inter_layer_deblk_filter_idc = 0;
    ps_svc_slice_params->u1_scan_idx_start = 0;
    ps_svc_slice_params->u1_scan_idx_end = 15;
    ps_svc_slice_params->i4_inter_layer_slice_alpha_c0_offset_div2 = 0;
    ps_svc_slice_params->i4_inter_layer_slice_beta_offset_div2 = 0;
    ps_svc_slice_params->u1_constrained_intra_resampling_flag = 0;

    return i_status;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_slice_header \endif
*
* \brief
*    parses the svc slice header attr
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_slice_header(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    dec_seq_params_t *ps_seq;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    dec_subset_seq_params_t *ps_sps_svc_ext = NULL;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    UWORD32 *pu4_bitstrm_buf = ps_dec->ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_dec->ps_bitstrm->u4_ofst;
    ps_svcd_ctxt = ps_svc_lyr_dec->ps_svcd_ctxt;
    ps_seq = ps_pps->ps_sps;
    ps_seq += MAX_NUM_SEQ_PARAMS;
    ps_subset_seq =
        &ps_svc_lyr_dec->ps_subset_sps[MAX_NUM_SEQ_PARAMS + ps_seq->u1_seq_parameter_set_id];
    ps_sps_svc_ext = &ps_subset_seq->s_sps_svc_ext;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    if(!ps_svc_lyr_dec->ps_nal_svc_ext->u1_no_inter_layer_pred_flag &&
       (0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id))
    {
        ps_svc_slice_params->u4_ref_layer_dq_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("Slice Header SVC ext: u4_ref_layer_dq_id",
                       ps_svc_slice_params->u4_ref_layer_dq_id);
        if(ps_svc_slice_params->u4_ref_layer_dq_id > MAX_REF_DEP_ID)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        /* Reference layer id update is taken care during resolution init */
        /*
        ps_svc_lyr_dec->u1_ref_layer_id = ps_svc_slice_params->u4_ref_layer_dq_id >> 4;
        if(ps_svc_lyr_dec->u1_ref_layer_id >= ps_svc_lyr_dec->u1_layer_id)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        */
        ps_svc_lyr_dec->ps_dec_svc_ref_layer =
            &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svc_lyr_dec->u1_ref_layer_id];

        if(ps_sps_svc_ext->u1_inter_layer_deblocking_filter_control_present_flag)
        {
            ps_svc_slice_params->u4_disable_inter_layer_deblk_filter_idc =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: u4_disable_inter_layer_deblk_filter_idc",
                           ps_svc_slice_params->u4_disable_inter_layer_deblk_filter_idc);

            if(ps_svc_slice_params->u4_disable_inter_layer_deblk_filter_idc > 6)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            if(1 != ps_svc_slice_params->u4_disable_inter_layer_deblk_filter_idc)
            {
                ps_svc_slice_params->i4_inter_layer_slice_alpha_c0_offset_div2 =
                    ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                COPYTHECONTEXT("Slice Header SVC ext: i4_inter_layer_slice_alpha_c0_offset_div2",
                               ps_svc_slice_params->i4_inter_layer_slice_alpha_c0_offset_div2);

                if(ps_svc_slice_params->i4_inter_layer_slice_alpha_c0_offset_div2 > 6 ||
                   ps_svc_slice_params->i4_inter_layer_slice_alpha_c0_offset_div2 < -6)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }

                ps_svc_slice_params->i4_inter_layer_slice_beta_offset_div2 =
                    ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                COPYTHECONTEXT("Slice Header SVC ext: i4_inter_layer_slice_beta_offset_div2",
                               ps_svc_slice_params->i4_inter_layer_slice_beta_offset_div2);

                if(ps_svc_slice_params->i4_inter_layer_slice_beta_offset_div2 > 6 ||
                   ps_svc_slice_params->i4_inter_layer_slice_beta_offset_div2 < -6)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
            }
        }

        ps_svc_slice_params->u1_constrained_intra_resampling_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: u1_constrained_intra_resampling_flag",
                       ps_svc_slice_params->u1_constrained_intra_resampling_flag);

        ps_svc_lyr_dec->s_res_prms.i1_constrained_intra_rsmpl_flag =
            ps_svc_lyr_dec->s_svc_slice_params.u1_constrained_intra_resampling_flag;
        isvcd_intra_resamp_res_init_update_flags(ps_svc_lyr_dec);

        if(2 == ps_sps_svc_ext->u1_extended_spatial_scalability_idc)
        {
            /* ChromaArrayType = i4_chroma_format_idc  if  separate_colour_plane_flag
             * = 0 for all chroma format except 4:4:4 */
            if(ps_dec->ps_cur_sps->i4_chroma_format_idc >= 0)
            {
                ps_svc_slice_params->u1_ref_layer_chroma_phase_x_plus1_flag =
                    ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("Slice Header SVC ext: u1_ref_layer_chroma_phase_x_plus1_flag",
                               ps_svc_slice_params->u1_ref_layer_chroma_phase_x_plus1_flag);

                ps_svc_slice_params->u1_ref_layer_chroma_phase_y_plus1 =
                    ih264d_get_bits_h264(ps_bitstrm, 2);
                COPYTHECONTEXT("Slice Header SVC ext: u1_ref_layer_chroma_phase_y_plus1",
                               ps_svc_slice_params->u1_ref_layer_chroma_phase_y_plus1);

                if(ps_svc_slice_params->u1_ref_layer_chroma_phase_y_plus1 > 2)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
            }
            else
            {
                if(0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id)
                {
                    ps_svc_slice_params->u1_ref_layer_chroma_phase_y_plus1 =
                        ps_sps_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1;
                }
            }

            ps_svc_slice_params->i4_scaled_ref_layer_left_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: i4_scaled_ref_layer_left_offset",
                           ps_svc_slice_params->i4_scaled_ref_layer_left_offset);

            if(ps_svc_slice_params->i4_scaled_ref_layer_left_offset != 0)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            if(ps_svc_slice_params->i4_scaled_ref_layer_left_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_svc_slice_params->i4_scaled_ref_layer_left_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            ps_svc_slice_params->i4_scaled_ref_layer_top_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: i4_scaled_ref_layer_top_offset",
                           ps_svc_slice_params->i4_scaled_ref_layer_top_offset);

            if(ps_svc_slice_params->i4_scaled_ref_layer_top_offset != 0)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            if(ps_svc_slice_params->i4_scaled_ref_layer_top_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_svc_slice_params->i4_scaled_ref_layer_top_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            ps_svc_slice_params->i4_scaled_ref_layer_right_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: i4_scaled_ref_layer_right_offset",
                           ps_svc_slice_params->i4_scaled_ref_layer_right_offset);

            if(ps_svc_slice_params->i4_scaled_ref_layer_right_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_svc_slice_params->i4_scaled_ref_layer_right_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_INV_SLICE_HDR_T;
            }

            ps_svc_slice_params->i4_scaled_ref_layer_bottom_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: i4_scaled_ref_layer_bottom_offset",
                           ps_svc_slice_params->i4_scaled_ref_layer_bottom_offset);

            if(ps_svc_slice_params->i4_scaled_ref_layer_bottom_offset >=
                   MAX_SCLD_REF_LAYER_OFFSET ||
               ps_svc_slice_params->i4_scaled_ref_layer_bottom_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_INV_SLICE_HDR_T;
            }
        }
        else
        {
            ps_svc_slice_params->i4_scaled_ref_layer_left_offset =
                ps_sps_svc_ext->i4_seq_scaled_ref_layer_left_offset;
            ps_svc_slice_params->i4_scaled_ref_layer_top_offset =
                ps_sps_svc_ext->i4_seq_scaled_ref_layer_top_offset;
            ps_svc_slice_params->i4_scaled_ref_layer_right_offset =
                ps_sps_svc_ext->i4_seq_scaled_ref_layer_right_offset;
            ps_svc_slice_params->i4_scaled_ref_layer_bottom_offset =
                ps_sps_svc_ext->i4_seq_scaled_ref_layer_bottom_offset;
        }
    }

    if(!ps_svc_lyr_dec->ps_nal_svc_ext->u1_no_inter_layer_pred_flag)
    {
        ps_svc_slice_params->u1_slice_skip_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: u1_slice_skip_flag",
                       ps_svc_slice_params->u1_slice_skip_flag);

        if(ps_svc_slice_params->u1_slice_skip_flag)
        {
            ps_svc_slice_params->u4_num_mbs_in_slice_minus1 =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: u4_num_mbs_in_slice_minus1",
                           ps_svc_slice_params->u4_num_mbs_in_slice_minus1);
        }
        else
        {
            ps_svc_slice_params->u1_adaptive_base_mode_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("Slice Header SVC ext: u1_adaptive_base_mode_flag",
                           ps_svc_slice_params->u1_adaptive_base_mode_flag);

            if(!ps_svc_slice_params->u1_adaptive_base_mode_flag)
            {
                ps_svc_slice_params->u1_default_base_mode_flag = ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("Slice Header SVC ext: u1_default_base_mode_flag",
                               ps_svc_slice_params->u1_default_base_mode_flag);
            }
            if(!ps_svc_slice_params->u1_default_base_mode_flag)
            {
                ps_svc_slice_params->u1_adaptive_motion_prediction_flag =
                    ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("Slice Header SVC ext: u1_adaptive_motion_prediction_flag",
                               ps_svc_slice_params->u1_adaptive_motion_prediction_flag);

                if(!ps_svc_slice_params->u1_adaptive_motion_prediction_flag)
                {
                    ps_svc_slice_params->u1_default_motion_prediction_flag =
                        ih264d_get_bit_h264(ps_bitstrm);
                    COPYTHECONTEXT("Slice Header SVC ext: u1_default_motion_prediction_flag",
                                   ps_svc_slice_params->u1_default_motion_prediction_flag);
                }
            }
            ps_svc_slice_params->u1_adaptive_residual_prediction_flag =
                ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("Slice Header SVC ext: u1_adaptive_residual_prediction_flag",
                           ps_svc_slice_params->u1_adaptive_residual_prediction_flag);

            if(!ps_svc_slice_params->u1_adaptive_residual_prediction_flag)
            {
                ps_svc_slice_params->u1_default_residual_prediction_flag =
                    ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("Slice Header SVC ext: u1_default_residual_prediction_flag",
                               ps_svc_slice_params->u1_default_residual_prediction_flag);
            }
        }

        if(ps_sps_svc_ext->u1_adaptive_tcoeff_level_prediction_flag)
        {
            ps_svc_slice_params->u1_tcoeff_level_prediction_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("Slice Header SVC ext: u1_tcoeff_level_prediction_flag",
                           ps_svc_slice_params->u1_tcoeff_level_prediction_flag);

            if(ps_svc_slice_params->u1_tcoeff_level_prediction_flag != 0)
            {
                return ERROR_INV_SPS_PPS_T;
            }
        }
    }

    if(!ps_sps_svc_ext->u1_slice_header_restriction_flag &&
       !ps_svc_slice_params->u1_slice_skip_flag)
    {
        ps_svc_slice_params->u1_scan_idx_start = ih264d_get_bits_h264(ps_bitstrm, 4);
        COPYTHECONTEXT("Slice Header SVC ext: u1_scan_idx_start",
                       ps_svc_slice_params->u1_scan_idx_start);
        ps_svc_slice_params->u1_scan_idx_end = ih264d_get_bits_h264(ps_bitstrm, 4);
        COPYTHECONTEXT("Slice Header SVC ext: u1_scan_idx_end",
                       ps_svc_slice_params->u1_scan_idx_end);

        if(0 != ps_svc_slice_params->u1_scan_idx_start &&
           15 != ps_svc_slice_params->u1_scan_idx_end)
            return ERROR_SVC_INV_SCAN_IDX;
    }
    return OK;
}

/*!
**************************************************************************
* \if Function name : DecodeSlice \endif
*
* \brief
*    Parses a slice
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/

WORD32 isvcd_parse_decode_slice(UWORD8 u1_is_idr_slice, UWORD8 u1_nal_ref_idc,
                                svc_dec_lyr_struct_t *ps_svc_lyr_dec /* SVC Decoder parameters */
)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    dec_pic_params_t *ps_pps;
    dec_seq_params_t *ps_seq;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_params_t *ps_cur_slice = ps_dec->ps_cur_slice;
    pocstruct_t s_tmp_poc = {0};
    WORD32 i_delta_poc[2] = {0};
    WORD32 i4_poc = 0;
    UWORD16 u2_first_mb_in_slice, u2_frame_num;
    UWORD8 u1_field_pic_flag, u1_redundant_pic_cnt = 0, u1_slice_type;
    UWORD32 u4_idr_pic_id = 0;
    UWORD8 u1_bottom_field_flag, u1_pic_order_cnt_type;
    UWORD8 u1_nal_unit_type;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    WORD8 i1_is_end_of_poc;

    WORD32 ret;
    WORD32 prev_slice_err, num_mb_skipped;
    UWORD8 u1_mbaff;
    pocstruct_t *ps_cur_poc;

    UWORD32 u4_temp;
    WORD32 i_temp;
    svc_dec_ctxt_t *psvcd_dec_ctxt;
    dec_struct_t *ps_dec_cur_lyr_minus_1;
    svc_dec_lyr_struct_t *ps_svc_cur_lyr_dec_minus_1;

    /* read FirstMbInSlice  and slice type*/
    ps_dec->ps_dpb_cmds->u1_dpb_commands_read_slc = 0;
    u2_first_mb_in_slice = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u2_first_mb_in_slice > (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    /*we currently don not support ASO*/
    if(((u2_first_mb_in_slice << ps_cur_slice->u1_mbaff_frame_flag) <= ps_dec->u2_cur_mb_addr) &&
       (ps_dec->u4_first_slice_in_pic == 0))
    {
        return ERROR_CORRUPTED_SLICE;
    }

    if(ps_dec->u4_first_slice_in_pic == 1)
    {
        if(u2_first_mb_in_slice != 0)
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }

    COPYTHECONTEXT("SH: first_mb_in_slice", u2_first_mb_in_slice);

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp > 9) return ERROR_INV_SLC_TYPE_T;

    u1_slice_type = u4_temp;
    COPYTHECONTEXT("SH: slice_type", (u1_slice_type));
    /* Find Out the Slice Type is 5 to 9 or not then Set the Flag   */
    /* u1_sl_typ_5_9 = 1 .Which tells that all the slices in the Pic*/
    /* will be of same type of current                            */
    if(u1_slice_type > 4)
    {
        u1_slice_type -= 5;
    }

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_PIC_SET_ID) return ERROR_INV_SLICE_HDR_T;
    /* discard slice if pic param is invalid */
    COPYTHECONTEXT("SH: pic_parameter_set_id", u4_temp);
    ps_pps = &ps_dec->ps_pps[u4_temp];
    if(FALSE == ps_pps->u1_is_valid)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    /* slices in a layer should have same PPS id*/
    if(UINT32_MAX == ps_svc_lyr_dec->u4_pps_id_for_layer)
    {
        ps_svc_lyr_dec->u4_pps_id_for_layer = u4_temp;
    }
    else if(u4_temp != ps_svc_lyr_dec->u4_pps_id_for_layer)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    ps_seq = ps_pps->ps_sps;
    ps_dec->ps_cur_sps = ps_seq;
    ps_subset_seq = &ps_svc_lyr_dec->ps_subset_sps[ps_seq->u1_seq_parameter_set_id];
    ps_svc_lyr_dec->ps_cur_subset_sps = ps_subset_seq;
    if(!ps_seq) return ERROR_INV_SLICE_HDR_T;
    if(FALSE == ps_seq->u1_is_valid) return ERROR_INV_SLICE_HDR_T;
    if(ps_seq->u1_mb_aff_flag) return ERROR_INV_SLICE_HDR_T;
    if(ps_seq->u1_level_idc > H264_LEVEL_4_2) return ERROR_INV_SLICE_HDR_T;
    if(!ps_seq->u1_frame_mbs_only_flag) return ERROR_INV_SLICE_HDR_T;
    if(OK != isvcd_verify_level(ps_seq->u1_level_idc)) return ERROR_INV_SLICE_HDR_T;
    if(ps_dec->u1_init_dec_flag == 1)
    {
        if(ps_dec->u2_frm_wd_in_mbs != ps_seq->u2_frm_wd_in_mbs) return ERROR_INV_SLICE_HDR_T;
        if(ps_dec->u2_frm_ht_in_mbs != ps_seq->u2_frm_ht_in_mbs) return ERROR_INV_SLICE_HDR_T;
    }

    if(ps_seq->u1_profile_idc == BASE_PROFILE_IDC)
    {
        if(ps_pps->u1_entropy_coding_mode != 0)
        {
            return ERROR_INV_SPS_PPS_T;
        }
    }

    ps_dec->i4_reorder_depth = ps_subset_seq->i4_reorder_depth;
    ps_dec->u2_disp_height = ps_subset_seq->u2_disp_height;
    ps_dec->u2_disp_width = ps_subset_seq->u2_disp_width;

    if(ps_svc_lyr_dec->u1_layer_id > 0)
    {
        psvcd_dec_ctxt = ps_svc_lyr_dec->ps_svcd_ctxt;
        ps_svc_cur_lyr_dec_minus_1 =
            &psvcd_dec_ctxt->ps_svc_dec_lyr[ps_svc_lyr_dec->u1_layer_id - 1];

        ps_dec_cur_lyr_minus_1 = &ps_svc_cur_lyr_dec_minus_1->s_dec;

        if((ps_dec_cur_lyr_minus_1->u2_pic_wd > ps_subset_seq->u2_pic_wd) ||
           (ps_dec_cur_lyr_minus_1->u2_pic_ht > ps_subset_seq->u2_pic_ht))
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }

    ps_dec->u2_pic_wd = ps_subset_seq->u2_pic_wd;
    ps_dec->u2_pic_ht = ps_subset_seq->u2_pic_ht;
    ps_dec->u4_total_mbs = ps_seq->u2_total_num_of_mbs << (1 - ps_seq->u1_frame_mbs_only_flag);

    /* Determining the Width and Height of Frame from that of Picture */
    ps_dec->u2_frm_wd_y = ps_subset_seq->u2_frm_wd_y;
    ps_dec->u2_frm_ht_y = ps_subset_seq->u2_frm_ht_y;

    ps_dec->u2_frm_wd_uv = ps_subset_seq->u2_frm_wd_uv;
    ps_dec->u2_frm_ht_uv = ps_subset_seq->u2_frm_ht_uv;

    ps_dec->s_pad_mgr.u1_pad_len_y_v = ps_subset_seq->u1_pad_len_y_v;
    ps_dec->s_pad_mgr.u1_pad_len_cr_v = ps_subset_seq->u1_pad_len_cr_v;
    ps_dec->u2_frm_wd_in_mbs = ps_seq->u2_frm_wd_in_mbs;
    ps_dec->u2_frm_ht_in_mbs = ps_seq->u2_frm_ht_in_mbs;

    ps_dec->u2_crop_offset_y = ps_subset_seq->u2_crop_offset_y;
    ps_dec->u2_crop_offset_uv = ps_subset_seq->u2_crop_offset_uv;

    /* Get the frame num */
    u2_frame_num = ih264d_get_bits_h264(ps_bitstrm, ps_seq->u1_bits_in_frm_num);
    COPYTHECONTEXT("SH: frame_num", u2_frame_num);

    if(!ps_dec->u1_first_slice_in_stream && ps_dec->u4_first_slice_in_pic)
    {
        pocstruct_t *ps_prev_poc = &ps_dec->s_prev_pic_poc;
        pocstruct_t *ps_cur_poc = &ps_dec->s_cur_pic_poc;

        ps_dec->u2_mbx = 0xffff;
        ps_dec->u2_mby = 0;

        if((0 == u1_is_idr_slice) && ps_cur_slice->u1_nal_ref_idc)
            ps_dec->u2_prev_ref_frame_num = ps_cur_slice->u2_frame_num;

        if(u1_is_idr_slice || ps_cur_slice->u1_mmco_equalto5) ps_dec->u2_prev_ref_frame_num = 0;

        if(ps_dec->ps_cur_sps->u1_gaps_in_frame_num_value_allowed_flag)
        {
            isvcd_decode_gaps_in_frame_num(ps_dec, u2_frame_num);
        }

        ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst;
        ps_prev_poc->u2_frame_num = ps_cur_poc->u2_frame_num;
        ps_prev_poc->u1_mmco_equalto5 = ps_cur_slice->u1_mmco_equalto5;
        if(ps_cur_slice->u1_nal_ref_idc)
        {
            ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb;
            ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb;
            ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom;
            ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0];
            ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1];
            ps_prev_poc->u1_bot_field = ps_cur_poc->u1_bot_field;
        }

        ps_dec->u2_total_mbs_coded = 0;
    }
    /* Get the field related flags  */
    if(!ps_seq->u1_frame_mbs_only_flag)
    {
        u1_field_pic_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SH: field_pic_flag", u1_field_pic_flag);
        u1_bottom_field_flag = 0;

        if(u1_field_pic_flag)
        {
            ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan_fld;
            u1_bottom_field_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("SH: bottom_field_flag", u1_bottom_field_flag);
        }
        else
        {
            ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan;
        }
    }
    else
    {
        u1_field_pic_flag = 0;
        u1_bottom_field_flag = 0;

        ps_dec->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan;
    }

    u1_nal_unit_type = SLICE_NAL;
    if(u1_is_idr_slice)
    {
        u1_nal_unit_type = IDR_SLICE_NAL;
        u4_idr_pic_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_idr_pic_id > 65535) return ERROR_INV_SLICE_HDR_T;
        COPYTHECONTEXT("SH:  ", u4_idr_pic_id);
    }

    /* read delta pic order count information*/
    i_delta_poc[0] = i_delta_poc[1] = 0;
    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    u1_pic_order_cnt_type = ps_seq->u1_pic_order_cnt_type;
    if(u1_pic_order_cnt_type == 0)
    {
        i_temp = ih264d_get_bits_h264(ps_bitstrm, ps_seq->u1_log2_max_pic_order_cnt_lsb_minus);
        if(i_temp < 0 || i_temp >= ps_seq->i4_max_pic_order_cntLsb) return ERROR_INV_SLICE_HDR_T;
        s_tmp_poc.i4_pic_order_cnt_lsb = i_temp;
        COPYTHECONTEXT("SH: pic_order_cnt_lsb", s_tmp_poc.i4_pic_order_cnt_lsb);

        if((ps_pps->u1_pic_order_present_flag == 1) && (!u1_field_pic_flag))
        {
            s_tmp_poc.i4_delta_pic_order_cnt_bottom = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SH: delta_pic_order_cnt_bottom",
                           s_tmp_poc.i4_delta_pic_order_cnt_bottom);
        }
    }

    s_tmp_poc.i4_delta_pic_order_cnt[0] = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[1] = 0;
    if(u1_pic_order_cnt_type == 1 && (!ps_seq->u1_delta_pic_order_always_zero_flag))
    {
        s_tmp_poc.i4_delta_pic_order_cnt[0] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SH: delta_pic_order_cnt[0]", s_tmp_poc.i4_delta_pic_order_cnt[0]);

        if(ps_pps->u1_pic_order_present_flag && !u1_field_pic_flag)
        {
            s_tmp_poc.i4_delta_pic_order_cnt[1] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SH: delta_pic_order_cnt[1]", s_tmp_poc.i4_delta_pic_order_cnt[1]);
        }
    }

    if(ps_pps->u1_redundant_pic_cnt_present_flag)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > MAX_REDUNDANT_PIC_CNT) return ERROR_INV_SLICE_HDR_T;
        u1_redundant_pic_cnt = u4_temp;
        COPYTHECONTEXT("SH: redundant_pic_cnt", u1_redundant_pic_cnt);
    }

    /*--------------------------------------------------------------------*/
    /* Check if the slice is part of new picture                          */
    /*--------------------------------------------------------------------*/
    /* First slice of a picture is always considered as part of new picture */
    i1_is_end_of_poc = 1;
    ps_dec->ps_dec_err_status->u1_err_flag &= MASK_REJECT_CUR_PIC;

    if(ps_dec->u4_first_slice_in_pic == 0)
    {
        i1_is_end_of_poc =
            ih264d_is_end_of_pic(u2_frame_num, u1_nal_ref_idc, &s_tmp_poc, &ps_dec->s_cur_pic_poc,
                                 ps_cur_slice, u1_pic_order_cnt_type, u1_nal_unit_type,
                                 u4_idr_pic_id, u1_field_pic_flag, u1_bottom_field_flag);
        if(i1_is_end_of_poc)
        {
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_INCOMPLETE_FRAME;
        }
    }

    /*--------------------------------------------------------------------*/
    /* Check for error in slice and parse the missing/corrupted MB's      */
    /* as skip-MB's in an inserted P-slice                                */
    /*--------------------------------------------------------------------*/
    u1_mbaff = ps_seq->u1_mb_aff_flag && (!u1_field_pic_flag);
    prev_slice_err = 0;

    if(i1_is_end_of_poc || ps_dec->u1_first_slice_in_stream)
    {
        /* If the current slice is not a field or frame number of the current
         * slice doesn't match with previous slice, and decoder is expecting
         * to decode a field i.e. ps_dec->u1_top_bottom_decoded is not 0 and
         * is not (TOP_FIELD_ONLY | BOT_FIELD_ONLY), treat it as a dangling
         * field */
        if((u1_field_pic_flag == 0 || u2_frame_num != ps_dec->u2_prv_frame_num) &&
           ps_dec->u1_top_bottom_decoded != 0 &&
           ps_dec->u1_top_bottom_decoded != (TOP_FIELD_ONLY | BOT_FIELD_ONLY))
        {
            ps_dec->u1_dangling_field = 1;
            if(ps_dec->u4_first_slice_in_pic)
            {
                // first slice - dangling field
                prev_slice_err = 1;
            }
            else
            {
                // last slice - dangling field
                prev_slice_err = 2;
            }

            if(ps_dec->u1_top_bottom_decoded == TOP_FIELD_ONLY)
                ps_cur_slice->u1_bottom_field_flag = 1;
            else
                ps_cur_slice->u1_bottom_field_flag = 0;

            num_mb_skipped =
                (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &ps_dec->s_cur_pic_poc;

            u1_is_idr_slice = ps_cur_slice->u1_nal_unit_type == IDR_SLICE_NAL;
        }
        else if(ps_dec->u4_first_slice_in_pic)
        {
            if(u2_first_mb_in_slice > 0)
            {
                /* first slice - missing/header corruption */
                prev_slice_err = 1;
                num_mb_skipped = u2_first_mb_in_slice << u1_mbaff;
                ps_cur_poc = &s_tmp_poc;

                /* initializing slice parameters */
                ps_cur_slice->u4_idr_pic_id = u4_idr_pic_id;
                ps_cur_slice->u1_field_pic_flag = u1_field_pic_flag;
                ps_cur_slice->u1_bottom_field_flag = u1_bottom_field_flag;
                ps_cur_slice->i4_pic_order_cnt_lsb = s_tmp_poc.i4_pic_order_cnt_lsb;
                ps_cur_slice->u1_nal_unit_type = u1_nal_unit_type;
                ps_cur_slice->u1_redundant_pic_cnt = u1_redundant_pic_cnt;
                ps_cur_slice->u1_nal_ref_idc = u1_nal_ref_idc;
                ps_cur_slice->u1_pic_order_cnt_type = u1_pic_order_cnt_type;
                ps_cur_slice->u1_mbaff_frame_flag = ps_seq->u1_mb_aff_flag && (!u1_field_pic_flag);
            }
        }
        else
        {
            /* since i1_is_end_of_poc is set ,means new frame num is encountered. so
             * conceal the current frame completely */
            prev_slice_err = 2;
            num_mb_skipped =
                (ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &s_tmp_poc;
        }
    }
    else
    {
        if((u2_first_mb_in_slice << u1_mbaff) > ps_dec->u2_total_mbs_coded)
        {
            // previous slice - missing/corruption
            prev_slice_err = 2;
            num_mb_skipped = (u2_first_mb_in_slice << u1_mbaff) - ps_dec->u2_total_mbs_coded;
            ps_cur_poc = &s_tmp_poc;
        }
        else if((u2_first_mb_in_slice << u1_mbaff) < ps_dec->u2_total_mbs_coded)
        {
            return ERROR_CORRUPTED_SLICE;
        }
    }
    if(prev_slice_err)
    {
        ret = isvcd_mark_err_slice_skip((svc_dec_lyr_struct_t *) ps_dec, num_mb_skipped,
                                        u1_is_idr_slice, u2_frame_num, ps_cur_poc, prev_slice_err);

        if(ps_dec->u1_dangling_field == 1)
        {
            ps_dec->u1_second_field = 1 - ps_dec->u1_second_field;
            ps_dec->u1_first_slice_in_stream = 0;
            ps_dec->u1_top_bottom_decoded = TOP_FIELD_ONLY | BOT_FIELD_ONLY;
            return ERROR_DANGLING_FIELD_IN_PIC;
        }

        if(prev_slice_err == 2)
        {
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_INCOMPLETE_FRAME;
        }

        if(ps_dec->u2_total_mbs_coded >= ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs)
        {
            /* return if all MBs in frame are parsed*/
            ps_dec->u1_first_slice_in_stream = 0;
            return ERROR_IN_LAST_SLICE_OF_PIC;
        }

        if(ps_dec->ps_dec_err_status->u1_err_flag & REJECT_CUR_PIC)
        {
            ih264d_err_pic_dispbuf_mgr(ps_dec);
            return ERROR_NEW_FRAME_EXPECTED;
        }

        if(ret != OK) return ret;

        i1_is_end_of_poc = 0;
    }

    if(u1_field_pic_flag)
    {
        ps_dec->u2_prv_frame_num = u2_frame_num;
    }

    if(ps_cur_slice->u1_mmco_equalto5 && NULL != ps_dec->ps_cur_pic)
    {
        WORD32 i4_temp_poc;
        WORD32 i4_top_field_order_poc, i4_bot_field_order_poc;
        WORD64 i8_result;
        if(!ps_cur_slice->u1_field_pic_flag)
        {
            i4_top_field_order_poc = ps_dec->ps_cur_pic->i4_top_field_order_cnt;
            i4_bot_field_order_poc = ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;
            i4_temp_poc = MIN(i4_top_field_order_poc, i4_bot_field_order_poc);
        }
        else if(!ps_cur_slice->u1_bottom_field_flag)
            i4_temp_poc = ps_dec->ps_cur_pic->i4_top_field_order_cnt;
        else
            i4_temp_poc = ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;

        i8_result = (WORD64) i4_temp_poc - ps_dec->ps_cur_pic->i4_top_field_order_cnt;
        if(IS_OUT_OF_RANGE_S32(i8_result))
        {
            return ERROR_INV_POC;
        }
        ps_dec->ps_cur_pic->i4_top_field_order_cnt = (WORD32) i8_result;
        i8_result = (WORD64) i4_temp_poc - ps_dec->ps_cur_pic->i4_bottom_field_order_cnt;
        if(IS_OUT_OF_RANGE_S32(i8_result))
        {
            return ERROR_INV_POC;
        }
        ps_dec->ps_cur_pic->i4_bottom_field_order_cnt = (WORD32) i8_result;
        ps_dec->ps_cur_pic->i4_poc = i4_temp_poc;
        ps_dec->ps_cur_pic->i4_avg_poc = i4_temp_poc;
    }
    if(ps_dec->u4_first_slice_in_pic)
    {
        ret = isvcd_decode_pic_order_cnt(u1_is_idr_slice, u2_frame_num, &ps_dec->s_prev_pic_poc,
                                         &s_tmp_poc, ps_cur_slice, ps_pps, u1_nal_ref_idc,
                                         u1_bottom_field_flag, u1_field_pic_flag, &i4_poc, ps_dec);
        if(ret != OK) return ret;
        /* Display seq no calculations */
        if(i4_poc >= ps_dec->i4_max_poc) ps_dec->i4_max_poc = i4_poc;
        /* IDR Picture or POC wrap around */
        if(i4_poc == 0)
        {
            WORD64 i8_temp;
            i8_temp = (WORD64) ps_dec->i4_prev_max_display_seq + ps_dec->i4_max_poc +
                      ps_dec->u1_max_dec_frame_buffering + 1;
            /*If i4_prev_max_display_seq overflows integer range, reset it */
            ps_dec->i4_prev_max_display_seq = IS_OUT_OF_RANGE_S32(i8_temp) ? 0 : i8_temp;
            ps_dec->i4_max_poc = 0;
        }
    }

    /* Increment only if the current slice has atleast 1 more MB */
    if(ps_dec->u4_first_slice_in_pic == 0 &&
       (ps_dec->ps_parse_cur_slice->u4_first_mb_in_slice <
        (UWORD32) (ps_dec->u2_total_mbs_coded >> ps_dec->ps_cur_slice->u1_mbaff_frame_flag)))
    {
        ps_dec->ps_parse_cur_slice++;
        ps_dec->u2_cur_slice_num++;
        // in the case of single core increment ps_decode_cur_slice
        if(ps_dec->u1_separate_parse == 0)
        {
            ps_dec->ps_decode_cur_slice++;
        }
    }

    ps_dec->u1_slice_header_done = 0;

    /*--------------------------------------------------------------------*/
    /* Copy the values read from the bitstream to the slice header and then*/
    /* If the slice is first slice in picture, then do Start of Picture   */
    /* processing.                                                        */
    /*--------------------------------------------------------------------*/
    ps_cur_slice->i4_delta_pic_order_cnt[0] = i_delta_poc[0];
    ps_cur_slice->i4_delta_pic_order_cnt[1] = i_delta_poc[1];
    ps_cur_slice->u4_idr_pic_id = u4_idr_pic_id;
    ps_cur_slice->u2_first_mb_in_slice = u2_first_mb_in_slice;
    ps_cur_slice->u1_field_pic_flag = u1_field_pic_flag;
    ps_cur_slice->u1_bottom_field_flag = u1_bottom_field_flag;
    ps_cur_slice->u1_slice_type = u1_slice_type;
    ps_cur_slice->i4_pic_order_cnt_lsb = s_tmp_poc.i4_pic_order_cnt_lsb;

    ps_cur_slice->u1_nal_unit_type = u1_nal_unit_type;
    ps_cur_slice->u1_redundant_pic_cnt = u1_redundant_pic_cnt;
    ps_cur_slice->u1_nal_ref_idc = u1_nal_ref_idc;
    ps_cur_slice->u1_pic_order_cnt_type = u1_pic_order_cnt_type;

    if(ps_seq->u1_frame_mbs_only_flag)
        ps_cur_slice->u1_direct_8x8_inference_flag = ps_seq->u1_direct_8x8_inference_flag;
    else
        ps_cur_slice->u1_direct_8x8_inference_flag = 1;

    if(u1_slice_type == B_SLICE)
    {
        ps_cur_slice->u1_direct_spatial_mv_pred_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SH: direct_spatial_mv_pred_flag",
                       ps_cur_slice->u1_direct_spatial_mv_pred_flag);

        if(ps_cur_slice->u1_direct_spatial_mv_pred_flag)
            ps_cur_slice->pf_decodeDirect = ih264d_decode_spatial_direct;
        else
            ps_cur_slice->pf_decodeDirect = ih264d_decode_temporal_direct;
        if(!((ps_seq->u1_mb_aff_flag) && (!u1_field_pic_flag)))
            ps_dec->pf_mvpred = ih264d_mvpred_nonmbaffB;
    }
    else
    {
        if(!((ps_seq->u1_mb_aff_flag) && (!u1_field_pic_flag)))
            ps_dec->pf_mvpred = ih264d_mvpred_nonmbaff;
    }

    if(ps_dec->u4_first_slice_in_pic)
    {
        if(u2_first_mb_in_slice == 0)
        {
            ret = isvcd_start_of_pic(ps_svc_lyr_dec, i4_poc, &s_tmp_poc, u2_frame_num, ps_pps);
            if(ret != OK) return ret;
            /*inter layer buffer intialization */
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
                ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start;
            ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_cur_mb =
                ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_base;
        }

        ps_dec->u4_output_present = 0;

        {
            ih264d_get_next_display_field(ps_dec, ps_dec->ps_out_buffer, &(ps_dec->s_disp_op));
            /* If error code is non-zero then there is no buffer available for
            display, hence avoid format conversion */

            if(0 != ps_dec->s_disp_op.u4_error_code)
            {
                ps_dec->u4_output_present = 0;
                ps_dec->u4_fmt_conv_cur_row = ps_dec->s_disp_frame_info.u4_y_ht;
            }
            else
                ps_dec->u4_output_present = 1;
        }
        ret = isvcd_parse_interlayer_resamp_func_init(ps_svc_lyr_dec, u2_first_mb_in_slice);
        if(ret != OK)
        {
            return ERROR_CORRUPTED_SLICE;
        }
        if((ps_dec->u1_separate_parse == 1) && (ps_svc_lyr_dec->u1_res_init_done == 1))
        {
            if(ps_dec->u4_dec_thread_created == 0)
            {
                if(ps_svc_lyr_dec->u1_layer_identifier != TARGET_LAYER)
                {
                    ithread_create(ps_dec->pv_dec_thread_handle, NULL,
                                   (void *) isvcd_decode_picture_thread, (void *) ps_dec);

                    ps_dec->u4_dec_thread_created = 1;
                }
                else
                {
                    ithread_create(ps_dec->pv_dec_thread_handle, NULL,
                                   (void *) ih264d_decode_picture_thread, (void *) ps_dec);

                    ps_dec->u4_dec_thread_created = 1;
                }
            }
#ifdef KEEP_THREADS_ACTIVE
            ret = ithread_mutex_lock(ps_dec->apv_proc_start_mutex[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);

            ps_dec->ai4_process_start[0] = PROC_START;
            ret = ithread_cond_signal(ps_dec->apv_proc_start_condition[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);

            ret = ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[0]);
            RETURN_IF((ret != IV_SUCCESS), ret);
#endif
#ifdef KEEP_THREADS_ACTIVE
            if(ps_dec->u4_bs_deblk_thread_created)
            {
                ret = ithread_mutex_lock(ps_dec->apv_proc_start_mutex[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ps_dec->ai4_process_start[1] = PROC_START;
                ret = ithread_cond_signal(ps_dec->apv_proc_start_condition[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ret = ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[1]);
                RETURN_IF((ret != IV_SUCCESS), ret);
            }
#endif
        }
    }

    /* INITIALIZATION of fn ptrs for MC and formMbPartInfo functions */
    {
        UWORD8 uc_nofield_nombaff;

        uc_nofield_nombaff =
            ((ps_dec->ps_cur_slice->u1_field_pic_flag == 0) &&
             (ps_dec->ps_cur_slice->u1_mbaff_frame_flag == 0) && (u1_slice_type != B_SLICE) &&
             (ps_dec->ps_cur_pps->u1_wted_pred_flag == 0));

        /* Initialise MC and formMbPartInfo fn ptrs one time based on profile_idc */

        if(uc_nofield_nombaff)
        {
            ps_dec->p_form_mb_part_info = ih264d_form_mb_part_info_bp;
            ps_dec->p_motion_compensate = ih264d_motion_compensate_bp;
        }
        else
        {
            ps_dec->p_form_mb_part_info = ih264d_form_mb_part_info_mp;
            ps_dec->p_motion_compensate = ih264d_motion_compensate_mp;
        }
    }

    /*
     * Decide whether to decode the current picture or not
     */
    {
        dec_err_status_t *ps_err = ps_dec->ps_dec_err_status;
        if(ps_err->u4_frm_sei_sync == u2_frame_num)
        {
            ps_err->u1_err_flag = ACCEPT_ALL_PICS;
            ps_err->u4_frm_sei_sync = SYNC_FRM_DEFAULT;
        }
        ps_err->u4_cur_frm = u2_frame_num;
    }

    /* Decision for decoding if the picture is to be skipped */
    {
        WORD32 i4_skip_b_pic, i4_skip_p_pic;

        i4_skip_b_pic = (ps_dec->u4_skip_frm_mask & B_SLC_BIT) && (B_SLICE == u1_slice_type) &&
                        (0 == u1_nal_ref_idc);

        i4_skip_p_pic = (ps_dec->u4_skip_frm_mask & P_SLC_BIT) && (P_SLICE == u1_slice_type) &&
                        (0 == u1_nal_ref_idc);

        /**************************************************************/
        /* Skip the B picture if skip mask is set for B picture and   */
        /* Current B picture is a non reference B picture or there is */
        /* no user for reference B picture                            */
        /**************************************************************/
        if(i4_skip_b_pic)
        {
            ps_dec->ps_cur_pic->u4_pack_slc_typ |= B_SLC_BIT;
            /* Don't decode the picture in SKIP-B mode if that picture is B */
            /* and also it is not to be used as a reference picture         */
            ps_dec->u1_last_pic_not_decoded = 1;

            return OK;
        }
        /**************************************************************/
        /* Skip the P picture if skip mask is set for P picture and   */
        /* Current P picture is a non reference P picture or there is */
        /* no user for reference P picture                            */
        /**************************************************************/
        if(i4_skip_p_pic)
        {
            ps_dec->ps_cur_pic->u4_pack_slc_typ |= P_SLC_BIT;
            /* Don't decode the picture in SKIP-P mode if that picture is P */
            /* and also it is not to be used as a reference picture         */
            ps_dec->u1_last_pic_not_decoded = 1;

            return OK;
        }
    }

    {
        UWORD16 u2_mb_x, u2_mb_y;

        ps_dec->i4_submb_ofst =
            ((u2_first_mb_in_slice << ps_cur_slice->u1_mbaff_frame_flag) * SUB_BLK_SIZE) -
            SUB_BLK_SIZE;
        if(u2_first_mb_in_slice)
        {
            UWORD8 u1_mb_aff;
            UWORD8 u1_field_pic;
            UWORD16 u2_frm_wd_in_mbs;
            u2_frm_wd_in_mbs = ps_seq->u2_frm_wd_in_mbs;
            u1_mb_aff = ps_cur_slice->u1_mbaff_frame_flag;
            u1_field_pic = ps_cur_slice->u1_field_pic_flag;

            {
                UWORD32 x_offset;
                UWORD32 y_offset;
                UWORD32 u4_frame_stride;
                tfr_ctxt_t *ps_trns_addr;

                if(ps_dec->u1_separate_parse)
                {
                    ps_trns_addr = &ps_dec->s_tran_addrecon_parse;
                }
                else
                {
                    ps_trns_addr = &ps_dec->s_tran_addrecon;
                }
                u2_mb_x = MOD(u2_first_mb_in_slice, u2_frm_wd_in_mbs);
                u2_mb_y = DIV(u2_first_mb_in_slice, u2_frm_wd_in_mbs);

                u2_mb_y <<= u1_mb_aff;

                if((u2_mb_x > u2_frm_wd_in_mbs - 1) || (u2_mb_y > ps_dec->u2_frm_ht_in_mbs - 1))
                {
                    return ERROR_CORRUPTED_SLICE;
                }

                u4_frame_stride = ps_dec->u2_frm_wd_y << u1_field_pic;
                x_offset = u2_mb_x << 4;
                y_offset = (u2_mb_y * u4_frame_stride) << 4;

                ps_trns_addr->pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1 + x_offset + y_offset;

                u4_frame_stride = ps_dec->u2_frm_wd_uv << u1_field_pic;
                x_offset >>= 1;
                y_offset = (u2_mb_y * u4_frame_stride) << 3;

                x_offset *= YUV420SP_FACTOR;

                ps_trns_addr->pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2 + x_offset + y_offset;
                ps_trns_addr->pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3 + x_offset + y_offset;

                ps_trns_addr->pu1_mb_y = ps_trns_addr->pu1_dest_y;
                ps_trns_addr->pu1_mb_u = ps_trns_addr->pu1_dest_u;
                ps_trns_addr->pu1_mb_v = ps_trns_addr->pu1_dest_v;

                /* assign the deblock structure pointers to start of slice */
                if(ps_dec->u1_separate_parse == 1)
                {
                    ps_dec->ps_deblk_mbn =
                        ps_dec->ps_deblk_pic + (u2_first_mb_in_slice << u1_mb_aff);
                }
                else
                {
                    ps_dec->ps_deblk_mbn =
                        ps_dec->ps_deblk_pic + (u2_first_mb_in_slice << u1_mb_aff);
                }

                ps_dec->u2_cur_mb_addr = (u2_first_mb_in_slice << u1_mb_aff);

                ps_dec->ps_mv_cur =
                    ps_dec->s_cur_pic.ps_mv + ((u2_first_mb_in_slice << u1_mb_aff) << 4);
            }
        }
        else
        {
            tfr_ctxt_t *ps_trns_addr;

            if(ps_dec->u1_separate_parse)
            {
                ps_trns_addr = &ps_dec->s_tran_addrecon_parse;
            }
            else
            {
                ps_trns_addr = &ps_dec->s_tran_addrecon;
            }

            u2_mb_x = 0xffff;
            u2_mb_y = 0;
            // assign the deblock structure pointers to start of slice
            ps_dec->u2_cur_mb_addr = 0;
            ps_dec->ps_deblk_mbn = ps_dec->ps_deblk_pic;
            ps_dec->ps_mv_cur = ps_dec->s_cur_pic.ps_mv;
            ps_trns_addr->pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1;
            ps_trns_addr->pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2;
            ps_trns_addr->pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3;

            ps_trns_addr->pu1_mb_y = ps_trns_addr->pu1_dest_y;
            ps_trns_addr->pu1_mb_u = ps_trns_addr->pu1_dest_u;
            ps_trns_addr->pu1_mb_v = ps_trns_addr->pu1_dest_v;
        }

        ps_dec->ps_part = ps_dec->ps_parse_part_params;

        ps_dec->u2_mbx = (MOD(u2_first_mb_in_slice - 1, ps_seq->u2_frm_wd_in_mbs));
        ps_dec->u2_mby = (DIV(u2_first_mb_in_slice - 1, ps_seq->u2_frm_wd_in_mbs));
        ps_dec->u2_mby <<= ps_cur_slice->u1_mbaff_frame_flag;
        ps_dec->i2_prev_slice_mbx = (WORD16) ps_dec->u2_mbx;
        ps_dec->i2_prev_slice_mby = (WORD16) ps_dec->u2_mby;
    }

    /* RBSP stop bit is used for CABAC decoding*/
    ps_bitstrm->u4_max_ofst += ps_dec->ps_cur_pps->u1_entropy_coding_mode;

    ps_dec->u1_B = (u1_slice_type == B_SLICE);
    ps_dec->u4_next_mb_skip = 0;

    ps_dec->ps_parse_cur_slice->u4_first_mb_in_slice = ps_dec->ps_cur_slice->u2_first_mb_in_slice;
    ps_dec->ps_parse_cur_slice->slice_type = ps_dec->ps_cur_slice->u1_slice_type;

    ps_dec->u4_start_recon_deblk = 1;
    {
        WORD32 num_entries;
        WORD32 size;
        UWORD8 *pu1_buf;

        num_entries = MAX_FRAMES;
        if((1 >= ps_dec->ps_cur_sps->u1_num_ref_frames) && (0 == ps_dec->i4_display_delay))
        {
            num_entries = 1;
        }
        num_entries = ((2 * num_entries) + 1);
        num_entries *= 2;

        size = num_entries * sizeof(void *);
        size += PAD_MAP_IDX_POC * sizeof(void *);

        pu1_buf = (UWORD8 *) ps_dec->pv_map_ref_idx_to_poc_buf;
        pu1_buf += size * ps_dec->u2_cur_slice_num;
        ps_dec->ps_parse_cur_slice->ppv_map_ref_idx_to_poc = (void *) pu1_buf;
    }

    if(ps_dec->u1_separate_parse)
    {
        ps_dec->ps_parse_cur_slice->pv_tu_coeff_data_start = ps_dec->pv_parse_tu_coeff_data;
    }
    else
    {
        ps_dec->pv_proc_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;
    }

    ret = ih264d_fix_error_in_dpb(ps_dec);
    if(ret < 0) return ERROR_DBP_MANAGER_T;

    if(u1_slice_type == I_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= I_SLC_BIT;

        ret = isvcd_parse_islice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        if(ps_dec->i4_pic_type != B_SLICE && ps_dec->i4_pic_type != P_SLICE)
            ps_dec->i4_pic_type = I_SLICE;
    }
    else if(u1_slice_type == P_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= P_SLC_BIT;
        ret = isvcd_parse_pslice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        if(ps_dec->i4_pic_type != B_SLICE) ps_dec->i4_pic_type = P_SLICE;
    }
    else if(u1_slice_type == B_SLICE)
    {
        ps_dec->ps_cur_pic->u4_pack_slc_typ |= B_SLC_BIT;
        ret = isvcd_parse_bslice(ps_svc_lyr_dec, u2_first_mb_in_slice);
        ps_dec->u1_pr_sl_type = u1_slice_type;
        ps_dec->i4_pic_type = B_SLICE;
    }
    else
        return ERROR_INV_SLC_TYPE_T;

    if(ps_dec->u1_slice_header_done)
    {
        /* set to zero to indicate a valid slice has been decoded */
        ps_dec->u1_first_slice_in_stream = 0;
    }

    if(ret != OK) return ret;

    if(u1_nal_ref_idc != 0)
    {
        if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
        {
            memcpy((void *) ps_dec->ps_dpb_cmds, (void *) (&(ps_dec->s_dpb_cmds_scratch)),
                   sizeof(dpb_commands_t));
        }
    }

    /* storing last Mb X and MbY of the slice */
    ps_dec->i2_prev_slice_mbx = ps_dec->u2_mbx;
    ps_dec->i2_prev_slice_mby = ps_dec->u2_mby;

    /* End of Picture detection */

    if(ps_dec->u2_total_mbs_coded >= (ps_seq->u2_max_mb_addr + 1))
    {
        ps_dec->u1_pic_decode_done = 1;
    }

    {
        dec_err_status_t *ps_err = ps_dec->ps_dec_err_status;
        if((ps_err->u1_err_flag & REJECT_PB_PICS) && (ps_err->u1_cur_pic_type == PIC_TYPE_I))
        {
            ps_err->u1_err_flag = ACCEPT_ALL_PICS;
        }
    }

    PRINT_BIN_BIT_RATIO(ps_dec)

    return ret;
}