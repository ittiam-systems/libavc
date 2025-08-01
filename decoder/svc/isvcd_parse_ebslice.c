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
 *  isvcd_parse_ebslice.c
 *
 * @brief
 *  Contains routines that decode a I slice type
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_parse_bmb_cabac()
 *  - isvcd_mv_pred_ref_tfr_nby2_ebmb()
 *  - isvcd_parse_bmb_cavlc()
 *  - isvcd_parse_bmb_non_direct_cabac()
 *  - isvcd_parse_bmb_non_direct_cavlc()
 *  - isvcd_parse_ebslice()
 *  - isvcd_parse_bslice()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <assert.h>
#include <string.h>
#include "ih264_defs.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_parse_slice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264d_mvpred.h"
#include "ih264d_parse_islice.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "isvcd_process_epslice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_deblocking.h"
#include "ih264d_cabac.h"
#include "ih264d_parse_mb_header.h"
#include "ih264d_error_handler.h"
#include "ih264d_mvpred.h"
#include "ih264d_cabac.h"
#include "ih264d_utils.h"
#include "ih264d_parse_headers.h"
#include "isvcd_parse_headers.h"
#include "isvcd_process_ebslice.h"
#include "isvcd_parse_slice.h"
#include "ih264_debug.h"
#include "isvcd_parse_cavlc.h"
#include "isvcd_mb_utils.h"

void ih264d_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec);
void isvcd_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec);
void ih264d_get_implicit_weights(dec_struct_t *ps_dec);
WORD32 isvcd_interlyr_motion_mode_pred(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                       dec_mb_info_t *ps_cur_mb_info,
                                       dec_svc_mb_info_t *ps_svc_cur_mb_info,
                                       parse_pmbarams_t *ps_mb_part_info,
                                       parse_part_params_t *ps_part);

/*!
 **************************************************************************
 * \if Function name : isvcd_parse_bmb_cabac \endif
 *
 * \brief
 *    This function parses CABAC syntax of a B MB.
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_parse_bmb_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num,
                             UWORD8 u1_num_mbsNby2)
{
    UWORD8 u1_cbp = 0;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    deblk_mb_t *ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_mb_num;
    const UWORD8 *puc_mb_mc_mode = (const UWORD8 *) gau1_ih264d_mb_mc_mode;
    UWORD8 u1_mb_type = ps_cur_mb_info->u1_mb_type;
    ctxt_inc_mb_info_t *p_curr_ctxt = ps_dec->ps_curr_ctxt_mb_info;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    decoding_envirnoment_t *ps_cab_env = &ps_dec->s_cab_dec_env;
    WORD32 ret;
    UWORD8 u1_Bdirect_tranform_read = 1;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 1;
    ps_cur_mb_info->u1_mb_mc_mode = puc_mb_mc_mode[5 + u1_mb_type];
    ps_cur_mb_info->u1_yuv_dc_block_flag = 0;
    ps_cur_deblk_mb->u1_mb_type |= D_B_SLICE;

    if(u1_mb_type != B_DIRECT)
    {
        ret = isvcd_parse_bmb_non_direct_cabac(ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info,
                                               u1_mb_num, u1_num_mbsNby2);
        if(ret != OK) return ret;
    }
    else
    {
        /************ STORING PARTITION INFO ***********/
        parse_part_params_t *ps_part_info;
        ps_part_info = ps_dec->ps_part;
        ps_part_info->u1_is_direct = PART_DIRECT_16x16;
        ps_part_info->u1_sub_mb_num = 0;
        ps_dec->ps_part++;
        p_curr_ctxt->u1_mb_type = CAB_BD16x16;

        MEMSET_16BYTES(&ps_dec->pu1_left_mv_ctxt_inc[0][0], 0);
        memset(ps_dec->pi1_left_ref_idx_ctxt_inc, 0, 4);
        MEMSET_16BYTES(p_curr_ctxt->u1_mv, 0);
        memset(p_curr_ctxt->i1_ref_idx, 0, 4);

        /* check whether transform8x8 u4_flag to be read or not */
        u1_Bdirect_tranform_read = ps_dec->s_high_profile.u1_direct_8x8_inference_flag;
    }

    if(ps_svc_slice_params->u1_adaptive_residual_prediction_flag &&
       ps_svc_cur_mb_info->u1_crop_window_flag)
    {
        ps_svc_cur_mb_info->u1_residual_prediction_flag = ih264d_decode_bin(
            1, ps_svc_lyr_dec->ps_residual_prediction_flag, ps_bitstrm, ps_cab_env);
        COPYTHECONTEXT("SVC ext: u1_residual_prediction_flag",
                       ps_svc_cur_mb_info->u1_residual_prediction_flag);
    }
    else
    {
        /*residual flag inference code */
        if(1 == ps_svc_cur_mb_info->u1_crop_window_flag &&
           1 == ps_svc_cur_mb_info->u1_base_mode_flag)
        {
            ps_svc_cur_mb_info->u1_residual_prediction_flag =
                ps_svc_slice_params->u1_default_residual_prediction_flag;
        }
        else
        {
            ps_svc_cur_mb_info->u1_residual_prediction_flag = 0;
        }
    }

    if(ps_svc_slice_params->u1_scan_idx_end >= ps_svc_slice_params->u1_scan_idx_start)
    {
        /* Read the Coded block pattern */
        u1_cbp = (WORD8) ih264d_parse_ctx_cbp_cabac(ps_dec);
        p_curr_ctxt->u1_cbp = u1_cbp;
        ps_cur_mb_info->u1_cbp = u1_cbp;

        if(u1_cbp > 47) return ERROR_CBP;
        COPYTHECONTEXT("coded_block_pattern", u1_cbp);

        ps_cur_mb_info->u1_tran_form8x8 = 0;
        ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;

        if((ps_dec->s_high_profile.u1_transform8x8_present) && (u1_cbp & (0xf)) &&
           (ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag) && (u1_Bdirect_tranform_read))
        {
            ps_cur_mb_info->u1_tran_form8x8 =
                ih264d_parse_transform8x8flag_cabac(ps_dec, ps_cur_mb_info);
            COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);

            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
            p_curr_ctxt->u1_transform8x8_ctxt = ps_cur_mb_info->u1_tran_form8x8;
        }
        else
        {
            p_curr_ctxt->u1_transform8x8_ctxt = 0;
        }
    }

    p_curr_ctxt->u1_intra_chroma_pred_mode = 0;
    p_curr_ctxt->u1_yuv_dc_csbp &= 0xFE;
    ps_dec->pu1_left_yuv_dc_csbp[0] &= 0x6;

    /* Read mb_qp_delta */
    if(u1_cbp)
    {
        WORD8 c_temp;
        ret = ih264d_parse_mb_qp_delta_cabac(ps_dec, &c_temp);
        if(ret != OK) return ret;
        COPYTHECONTEXT("mb_qp_delta", c_temp);
        if(c_temp)
        {
            ret = ih264d_update_qp(ps_dec, c_temp);
            if(ret != OK) return ret;
        }
    }
    else
        ps_dec->i1_prev_mb_qp_delta = 0;

    /*RESIDUAL FOR Start to end idx*/
    ih264d_parse_residual4x4_cabac(ps_dec, ps_cur_mb_info, 0);
    if(EXCEED_OFFSET(ps_dec->ps_bitstrm)) return ERROR_EOB_TERMINATE_T;
    return OK;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_mv_pred_ref_tfr_nby2_ebmb \endif
 *
 * \brief
 *    This function computes the mv pred for b frame
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_mv_pred_ref_tfr_nby2_ebmb(dec_struct_t *ps_dec, UWORD32 u4_mb_idx, UWORD32 u4_num_mbs)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) ps_dec;
    parse_pmbarams_t *ps_mb_part_info;
    parse_part_params_t *ps_part;
    mv_pred_t *ps_mv_nmb, *ps_mv_nmb_start, *ps_mv_ntop, *ps_mv_ntop_start;
    pic_buffer_t *ps_ref_frame;
    UWORD8 u1_direct_mode_width;
    UWORD8 i, j;
    dec_mb_info_t *ps_cur_mb_info;
    dec_svc_mb_info_t *ps_svc_cur_mb_info;
    const UWORD8 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    UWORD8 u1_field;
    WORD32 ret = 0;
    WORD16 i2_mv_x, i2_mv_y;

    ps_dec->i4_submb_ofst -= (u4_num_mbs - u4_mb_idx) << 4;
    ps_mb_part_info = ps_dec->ps_parse_mb_data;
    ps_part = ps_dec->ps_parse_part_params;

    /* N/2 Mb MvPred and Transfer Setup Loop */
    for(i = u4_mb_idx; i < u4_num_mbs; i++, ps_mb_part_info++)
    {
        UWORD8 u1_colz = 0;
        ps_dec->i4_submb_ofst += SUB_BLK_SIZE;
        /* Restore the slice scratch MbX and MbY context */
        ps_cur_mb_info = ps_dec->ps_nmb_info + i;
        ps_svc_cur_mb_info = ps_svc_lyr_dec->ps_svc_nmb_info + i;
        u1_field = ps_cur_mb_info->u1_mb_field_decodingflag;
        ps_mv_nmb_start = ps_dec->ps_mv_cur + (i << 4);
        ps_dec->u2_mbx = ps_cur_mb_info->u2_mbx;
        ps_dec->u2_mby = ps_cur_mb_info->u2_mby;
        ps_dec->u1_currB_type = 0;
        ps_dec->u2_mv_2mb[i & 0x1] = 0;

        /* Look for MV Prediction and Reference Transfer in Non-I Mbs */
        if(!ps_mb_part_info->u4_isI_mb)
        {
            UWORD8 u1_blk_no;
            WORD16 i1_ref_idx, i1_ref_idx1;
            UWORD8 u1_pred_mode;
            UWORD8 u1_sub_mb_x, u1_sub_mb_y, u1_sub_mb_num;
            UWORD8 u1_lx, u1_lx_start, u1_lxend, u1_tmp_lx;
            UWORD8 u1_num_part, u1_num_ref, u1_wd, u1_ht;
            UWORD32 *pu4_wt_offst;
            UWORD8 u1_scale_ref, u4_bot_mb;
            deblk_mb_t *ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + i;
            WORD8(*pi1_ref_idx)[MAX_REFIDX_INFO_PER_MB] = ps_mb_part_info->i1_ref_idx;
            WORD8 *pi1_ref_idx0 = pi1_ref_idx[0], *pi1_ref_idx1 = pi1_ref_idx[1];
            UWORD32 **ppu4_wt_ofst = ps_mb_part_info->pu4_wt_offst;
            WORD32 i4_mb_mode_svc;
            UWORD8 u1_motion_pred_flag_l0 = ps_svc_cur_mb_info->au1_motion_pred_flag[0];
            UWORD8 u1_motion_pred_flag_l1 = ps_svc_cur_mb_info->au1_motion_pred_flag[1];

            /* MB Level initialisations */
            ps_dec->u4_num_pmbair = i >> u1_mbaff;
            ps_dec->u4_mb_idx_mv = i;

            i4_mb_mode_svc = isvcd_interlyr_motion_mode_pred(
                ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info, ps_mb_part_info, ps_part);

            if((-1 == i4_mb_mode_svc) || (SVC_INTER_MB == i4_mb_mode_svc))
            {
                ps_mv_ntop_start =
                    ps_mv_nmb_start - (ps_dec->u2_frm_wd_in_mbs << (4 + u1_mbaff)) + 12;

                u1_num_part = ps_mb_part_info->u1_num_part;
                ps_cur_deblk_mb->u1_mb_type |= (u1_num_part > 1) << 1;
                u1_direct_mode_width = (1 == ps_mb_part_info->u1_num_part) ? 16 : 8;

                ps_cur_mb_info->u4_pred_info_pkd_idx = ps_dec->u4_pred_info_pkd_idx;
                ps_cur_mb_info->u1_num_pred_parts = 0;

                /****************************************************/
                /* weighted u4_ofst pointer calculations, this loop  */
                /* runs maximum 4 times, even in direct cases       */
                /****************************************************/
                u1_scale_ref = u1_mbaff & ps_cur_mb_info->u1_mb_field_decodingflag;
                u4_bot_mb = 1 - ps_cur_mb_info->u1_topmb;
                if(ps_dec->ps_cur_pps->u1_wted_bipred_idc)
                {
                    u1_num_ref = MIN(u1_num_part, 4);
                    if(PART_DIRECT_16x16 != ps_part->u1_is_direct)
                    {
                        for(u1_blk_no = 0; u1_blk_no < u1_num_ref; u1_blk_no++)
                        {
                            i1_ref_idx = MAX(pi1_ref_idx0[u1_blk_no], 0);
                            if(u1_scale_ref) i1_ref_idx >>= 1;
                            i1_ref_idx *= ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1];
                            if(u1_scale_ref)
                                i1_ref_idx += (MAX(pi1_ref_idx1[u1_blk_no], 0) >> 1);
                            else
                                i1_ref_idx += MAX(pi1_ref_idx1[u1_blk_no], 0);
                            pu4_wt_offst = (UWORD32 *) &ps_dec->pu4_wt_ofsts[2 * X3(i1_ref_idx)];

                            if(pi1_ref_idx0[u1_blk_no] < 0) pu4_wt_offst += 1;

                            ppu4_wt_ofst[u1_blk_no] = pu4_wt_offst;
                            if(u1_scale_ref && (ps_dec->ps_cur_pps->u1_wted_bipred_idc == 2))
                            {
                                i1_ref_idx = MAX(pi1_ref_idx0[u1_blk_no], 0);
                                i1_ref_idx *=
                                    (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1] << 1);
                                i1_ref_idx += MAX(pi1_ref_idx1[u1_blk_no], 0);
                                if(u4_bot_mb)
                                {
                                    i1_ref_idx +=
                                        (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0] << 1) *
                                        (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1] << 1);
                                }
                                pu4_wt_offst =
                                    (UWORD32 *) &ps_dec->pu4_mbaff_wt_mat[2 * X3(i1_ref_idx)];
                                ppu4_wt_ofst[u1_blk_no] = pu4_wt_offst;
                            }
                        }
                    }
                }

                /**************************************************/
                /* Loop on Partitions                             */
                /* direct mode is reflected as a single partition */
                /**************************************************/
                for(j = 0; j < u1_num_part; j++, ps_part++)
                {
                    u1_sub_mb_num = ps_part->u1_sub_mb_num;
                    ps_dec->u1_sub_mb_num = u1_sub_mb_num;

                    if(PART_NOT_DIRECT != ps_part->u1_is_direct)
                    {
                        /**************************************************/
                        /* Direct Mode, Call DecodeSpatial/TemporalDirect */
                        /* only (those will in turn call FormMbPartInfo)  */
                        /**************************************************/
                        ret = isvcd_decode_spatial_direct(ps_dec, u1_direct_mode_width,
                                                          ps_cur_mb_info, i);
                        if(ret != OK) return ret;
                        ps_cur_deblk_mb->u1_mb_type |= (ps_dec->u1_currB_type << 1);
                    }
                    else
                    {
                        mv_pred_t s_mvPred = {0};
                        /**************************************************/
                        /* Non Direct Mode, Call Motion Vector Predictor  */
                        /* and FormMbpartInfo                             */
                        /**************************************************/
                        u1_sub_mb_x = u1_sub_mb_num & 0x03;
                        u1_sub_mb_y = u1_sub_mb_num >> 2;
                        u1_blk_no = (u1_num_part < 4)
                                        ? j
                                        : (((u1_sub_mb_y >> 1) << 1) + (u1_sub_mb_x >> 1));

                        ps_mv_ntop = ps_mv_ntop_start + u1_sub_mb_x;
                        ps_mv_nmb = ps_mv_nmb_start + u1_sub_mb_num;

                        /* Populate the colpic info and reference frames */
                        s_mvPred.i1_ref_frame[0] = pi1_ref_idx0[u1_blk_no];
                        s_mvPred.i1_ref_frame[1] = pi1_ref_idx1[u1_blk_no];
                        u1_pred_mode = ps_part->u1_pred_mode;
                        u1_wd = ps_part->u1_partwidth;
                        u1_ht = ps_part->u1_partheight;

                        if(1 != ps_svc_cur_mb_info->u1_base_mode_flag)
                        {
                            u1_lx_start = 0;
                            u1_lxend = 2;
                            if(PRED_L0 == u1_pred_mode)
                            {
                                s_mvPred.i2_mv[2] = 0;
                                s_mvPred.i2_mv[3] = 0;
                                if(0 == (u1_motion_pred_flag_l0 & (1 << u1_blk_no)))
                                {
                                    u1_lxend = 1;
                                }
                                else
                                {
                                    u1_lxend = 0;
                                }
                            }
                            else if(PRED_L1 == u1_pred_mode)
                            {
                                s_mvPred.i2_mv[0] = 0;
                                s_mvPred.i2_mv[1] = 0;
                                if(0 == (u1_motion_pred_flag_l1 & (1 << u1_blk_no)))
                                {
                                    u1_lx_start = 1;
                                }
                                else
                                {
                                    u1_lx_start = 2;
                                }
                            }
                            else  // Bi Pred
                            {
                                if(0 == (u1_motion_pred_flag_l0 & (1 << u1_blk_no)))
                                {
                                    u1_lxend = 1;
                                }
                                if(0 == (u1_motion_pred_flag_l1 & (1 << u1_blk_no)))
                                {
                                    u1_lx_start = 1;
                                }
                                if((0 != (u1_motion_pred_flag_l0 & (1 << u1_blk_no))) &&
                                   (0 != (u1_motion_pred_flag_l1 & (1 << u1_blk_no))))
                                {
                                    u1_lx_start = 0;
                                    u1_lxend = 0;
                                }
                                if((0 == (u1_motion_pred_flag_l0 & (1 << u1_blk_no))) &&
                                   (0 == (u1_motion_pred_flag_l1 & (1 << u1_blk_no))))
                                {
                                    u1_lx_start = 0;
                                    u1_lxend = 2;
                                }
                            }
                            ps_dec->pf_mvpred(ps_dec, ps_cur_mb_info, ps_mv_nmb, ps_mv_ntop,
                                              &s_mvPred, u1_sub_mb_num, u1_wd, u1_lx_start,
                                              u1_lxend, ps_cur_mb_info->u1_mb_mc_mode);
                        }

                        /* for generic case based on pred mode derived / signalled */
                        u1_lx_start = 0;
                        u1_lxend = 2;
                        if(PRED_L0 == u1_pred_mode)
                        {
                            s_mvPred.i2_mv[2] = 0;
                            s_mvPred.i2_mv[3] = 0;
                            u1_lxend = 1;
                        }
                        if(PRED_L1 == u1_pred_mode)
                        {
                            s_mvPred.i2_mv[0] = 0;
                            s_mvPred.i2_mv[1] = 0;
                            u1_lx_start = 1;
                        }

                        /**********************************************************/
                        /* Loop on number of predictors, 1 Each for Forw Backw    */
                        /* Loop 2 times for BiDirect mode                         */
                        /**********************************************************/
                        for(u1_lx = u1_lx_start; u1_lx < u1_lxend; u1_lx++)
                        {
                            UWORD8 u1_motion_pred_flag =
                                u1_lx ? u1_motion_pred_flag_l1 : u1_motion_pred_flag_l0;

                            if((0 != (u1_motion_pred_flag & (1 << u1_blk_no))) ||
                               (ps_svc_cur_mb_info->u1_base_mode_flag))
                            {
                                isvcd_retrive_infer_mode_mv(ps_svc_lyr_dec, &s_mvPred, u1_lx,
                                                            u1_sub_mb_num);
                            }
                            /********************************************************/
                            /* Predict Mv                                           */
                            /* Add Mv Residuals and store back                      */
                            /********************************************************/
                            u1_tmp_lx = (u1_lx << 1);
                            i1_ref_idx = s_mvPred.i1_ref_frame[u1_lx];
                            /********************************************************************/
                            /* If reference index is inferred from the base layer and it is     */
                            /* exceeding the number of active reference in the current layer.   */
                            /* Then reference index is clipped to the max in the current layer  */
                            /********************************************************************/
                            if(ps_svc_cur_mb_info->u1_base_mode_flag == 1)
                            {
                                if(i1_ref_idx > (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[u1_lx] - 1))
                                {
                                    i1_ref_idx = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[u1_lx] - 1;
                                }
                            }
                            if(0 == ps_svc_cur_mb_info->u1_base_mode_flag)
                            {
                                i2_mv_x = ps_mv_nmb->i2_mv[u1_tmp_lx];
                                i2_mv_y = ps_mv_nmb->i2_mv[u1_tmp_lx + 1];

                                i2_mv_x += s_mvPred.i2_mv[u1_tmp_lx];
                                i2_mv_y += s_mvPred.i2_mv[u1_tmp_lx + 1];

                                s_mvPred.i2_mv[u1_tmp_lx] = i2_mv_x;
                                s_mvPred.i2_mv[u1_tmp_lx + 1] = i2_mv_y;
                            }
                            else
                            {
                                i2_mv_x = s_mvPred.i2_mv[u1_tmp_lx];
                                i2_mv_y = s_mvPred.i2_mv[u1_tmp_lx + 1];
                            }

                            /********************************************************/
                            /* Transfer setup call                                  */
                            /* convert RefIdx if it is MbAff                        */
                            /* Pass Weight Offset and refFrame                      */
                            /********************************************************/
                            i1_ref_idx1 = i1_ref_idx >> u1_scale_ref;

                            if(-1 == i1_ref_idx1) return NOT_OK;
                            if(u1_scale_ref && ((i1_ref_idx & 0x01) != u4_bot_mb))
                                i1_ref_idx1 += MAX_REF_BUFS;
                            ps_ref_frame = ps_dec->ps_ref_pic_buf_lx[u1_lx][i1_ref_idx1];

                            /* Storing Colocated-Zero u4_flag */
                            if(u1_lx == u1_lx_start)
                            {
                                /* Fill colocated info in MvPred structure */
                                s_mvPred.u1_col_ref_pic_idx = ps_ref_frame->u1_mv_buf_id;
                                s_mvPred.u1_pic_type = ps_ref_frame->u1_pic_type;

                                /* Calculating colocated zero information */
                                u1_colz =
                                    (u1_field << 1) | ((i1_ref_idx == 0) && (ABS(i2_mv_x) <= 1) &&
                                                       (ABS(i2_mv_y) <= 1));
                                u1_colz |= ps_mb_part_info->u1_col_info[u1_blk_no];
                            }

                            pu4_wt_offst = ppu4_wt_ofst[u1_blk_no];
                            {
                                pred_info_pkd_t *ps_pred_pkd;
                                WORD16 i2_mv[2];

                                i2_mv[0] = i2_mv_x;
                                i2_mv[1] = i2_mv_y;

                                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                                ih264d_fill_pred_info(i2_mv, u1_wd, u1_ht, u1_sub_mb_num,
                                                      u1_pred_mode, ps_pred_pkd,
                                                      ps_ref_frame->u1_pic_buf_id, i1_ref_idx,
                                                      pu4_wt_offst, ps_ref_frame->u1_pic_type);
                                ps_dec->u4_pred_info_pkd_idx++;
                                ps_cur_mb_info->u1_num_pred_parts++;
                            }
                        }
                        if(ps_mv_nmb)
                        {
                            ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb, u1_sub_mb_num, u1_colz,
                                               u1_ht, u1_wd);
                        }
                        else
                        {
                            return NOT_OK;
                        }
                    }
                }
                /* to take care of 16 parttitions increment for base mode flag case*/
                if(1 == ps_svc_cur_mb_info->u1_base_mode_flag)
                {
                    ps_part += (MAX_NUM_MB_PART - u1_num_part);
                }
            }
            else
            {
                /* Set zero values in case of Intra Mbs */
                mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
                /* to take care of 16 parttitions increment for base mode flag case*/
                if(1 != ps_svc_cur_mb_info->u1_base_mode_flag)
                {
                    return NOT_OK;
                }

                ps_cur_deblk_mb->u1_mb_type |= D_INTRA_IBL;
                if((ps_svc_lyr_dec->u1_layer_identifier != TARGET_LAYER) &&
                   (DBLK_ENABLED == ps_dec->ps_cur_slice->u1_disable_dblk_filter_idc))
                {
                    ps_cur_deblk_mb->u1_deblocking_mode = MB_ENABLE_FILTERING;
                }

                ps_part += (MAX_NUM_MB_PART);
                /* Storing colocated zero information */
                if(ps_mv_nmb_start)
                {
                    ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0,
                                       (UWORD8) (u1_field << 1), 4, 4);
                }
                else
                {
                    return NOT_OK;
                }
            }
        }
        else
        {
            /* Set zero values in case of Intra Mbs */
            mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
            /* Storing colocated zero information */
            if(ps_mv_nmb_start)
            {
                ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0, (UWORD8) (u1_field << 1),
                                   4, 4);
            }
            else
            {
                return NOT_OK;
            }
        }
    }

    return OK;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_bmb_cavlc \endif
*
* \brief
*    This function parses CAVLC syntax of a B MB.
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_bmb_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num,
                             UWORD8 u1_num_mbsNby2)
{
    UWORD32 u4_cbp = 0;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    deblk_mb_t *ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_mb_num;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    const UWORD8 *puc_mb_mc_mode = (const UWORD8 *) gau1_ih264d_mb_mc_mode;
    UWORD8 u1_mb_type = ps_cur_mb_info->u1_mb_type;
    WORD32 ret;
    UWORD8 u1_Bdirect_tranform_read = 1;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 1;
    ps_cur_mb_info->u1_tran_form8x8 = 0;
    ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
    ps_cur_mb_info->u1_yuv_dc_block_flag = 0;
    ps_cur_mb_info->u1_mb_mc_mode = puc_mb_mc_mode[5 + u1_mb_type];
    ps_cur_deblk_mb->u1_mb_type |= D_B_SLICE;
    if(u1_mb_type != B_DIRECT)
    {
        ret = isvcd_parse_bmb_non_direct_cavlc(ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info,
                                               u1_mb_num, u1_num_mbsNby2);
        if(ret != OK) return ret;
    }
    else
    {
        /************ STORING PARTITION INFO ***********/
        parse_part_params_t *ps_part_info;
        ps_part_info = ps_dec->ps_part;
        ps_part_info->u1_is_direct = PART_DIRECT_16x16;
        ps_part_info->u1_sub_mb_num = 0;
        ps_dec->ps_part++;
        /* check whether transform8x8 u4_flag to be read or not */
        u1_Bdirect_tranform_read = ps_dec->s_high_profile.u1_direct_8x8_inference_flag;
    }

    if(ps_svc_slice_params->u1_adaptive_residual_prediction_flag &&
       ps_svc_cur_mb_info->u1_crop_window_flag)
    {
        ps_svc_cur_mb_info->u1_residual_prediction_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SVC ext: u1_residual_prediction_flag",
                       ps_svc_cur_mb_info->u1_residual_prediction_flag);
    }
    else
    {
        /*residual flag inference code */
        if(1 == ps_svc_cur_mb_info->u1_crop_window_flag &&
           1 == ps_svc_cur_mb_info->u1_base_mode_flag)
        {
            ps_svc_cur_mb_info->u1_residual_prediction_flag =
                ps_svc_slice_params->u1_default_residual_prediction_flag;
        }
        else
        {
            ps_svc_cur_mb_info->u1_residual_prediction_flag = 0;
        }
    }

    if(ps_svc_slice_params->u1_scan_idx_end >= ps_svc_slice_params->u1_scan_idx_start)
    {
        /* Read the Coded block pattern */
        const UWORD8 *puc_CbpInter = gau1_ih264d_cbp_inter;
        // Inlined ih264d_uev
        UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
        UWORD32 u4_word, u4_ldz;

        /***************************************************************/
        /* Find leading zeros in next 32 bits                          */
        /***************************************************************/
        NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
        u4_ldz = CLZ(u4_word);
        /* Flush the ps_bitstrm */
        u4_bitstream_offset += (u4_ldz + 1);
        /* Read the suffix from the ps_bitstrm */
        u4_word = 0;
        if(u4_ldz) GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
        *pu4_bitstrm_ofst = u4_bitstream_offset;
        u4_cbp = ((1 << u4_ldz) + u4_word - 1);
        // Inlined ih264d_uev
        if(u4_cbp > 47) return ERROR_CBP;
        u4_cbp = puc_CbpInter[u4_cbp];

        if((ps_dec->s_high_profile.u1_transform8x8_present) && (u4_cbp & (0xf)) &&
           (ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag) && (u1_Bdirect_tranform_read))
        {
            ps_cur_mb_info->u1_tran_form8x8 = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
        }

        COPYTHECONTEXT("coded_block_pattern", u4_cbp);
        ps_cur_mb_info->u1_cbp = u4_cbp;
    }

    /* Read mb_qp_delta */
    if(u4_cbp)
    {
        WORD32 i_temp;
        // inlining ih264d_sev
        UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
        UWORD32 u4_word, u4_ldz, u4_abs_val;

        /***************************************************************/
        /* Find leading zeros in next 32 bits                          */
        /***************************************************************/
        NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
        u4_ldz = CLZ(u4_word);

        /* Flush the ps_bitstrm */
        u4_bitstream_offset += (u4_ldz + 1);

        /* Read the suffix from the ps_bitstrm */
        u4_word = 0;
        if(u4_ldz) GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);

        *pu4_bitstrm_ofst = u4_bitstream_offset;
        u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

        if(u4_word & 0x1)
            i_temp = (-(WORD32) u4_abs_val);
        else
            i_temp = (u4_abs_val);

        if(i_temp < -26 || i_temp > 25) return ERROR_INV_RANGE_QP_T;
        // inlinined ih264d_sev
        COPYTHECONTEXT("mb_qp_delta", i_temp);
        if(i_temp)
        {
            ret = ih264d_update_qp(ps_dec, (WORD8) i_temp);
            if(ret != OK) return ret;
        }

        /*SVC residual from start to end idx*/
        ret = ih264d_parse_residual4x4_cavlc(ps_dec, ps_cur_mb_info, 0);
        if(ret != OK) return ret;
        if(EXCEED_OFFSET(ps_bitstrm)) return ERROR_EOB_TERMINATE_T;
    }
    else
    {
        ps_dec->i1_prev_mb_qp_delta = 0;
        ih264d_update_nnz_for_skipmb(ps_dec, ps_cur_mb_info, CAVLC);
    }

    return OK;
}

/*!
**************************************************************************
* \if Function name : ParseMb_SubMb_PredBCab_svc\endif
*
* \brief
*    Implements sub_mb_pred() of 7.3.5.2. & mb_pred() of 7.3.5.1
*
* \return
*    None.
*
**************************************************************************
*/

WORD32 isvcd_parse_bmb_non_direct_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                        dec_mb_info_t *ps_cur_mb_info,
                                        dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num,
                                        UWORD8 u1_num_mbsNby2)
{
    /* Loads from ps_dec */
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    decoding_envirnoment_t *ps_cab_env = &ps_dec->s_cab_dec_env;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    ctxt_inc_mb_info_t *p_curr_ctxt = ps_dec->ps_curr_ctxt_mb_info;
    parse_pmbarams_t *ps_parse_mb_data = ps_dec->ps_parse_mb_data + u1_num_mbsNby2;

    /* table pointer loads */
    const UWORD8 *pu1_sub_mb_pred_modes = (UWORD8 *) (gau1_ih264d_submb_pred_modes) + 4;
    const UWORD8(*pu1_mb_pred_modes)[32] = (const UWORD8(*)[32]) gau1_ih264d_mb_pred_modes;
    const UWORD8 *pu1_num_mb_part = (const UWORD8 *) gau1_ih264d_num_mb_part;
    const UWORD8 *pu1_sub_mb_mc_mode = (UWORD8 *) (gau1_ih264d_submb_mc_mode) + 4;

    const UWORD8 u1_mb_type = ps_cur_mb_info->u1_mb_type;
    UWORD8 *pu1_col_info = ps_parse_mb_data->u1_col_info;
    WORD8 *pi1_ref_idx_l0 = &ps_parse_mb_data->i1_ref_idx[0][0];
    WORD8 *pi1_ref_idx_l1 = &ps_parse_mb_data->i1_ref_idx[1][0];
    UWORD8 u1_dec_ref_l0, u1_dec_ref_l1;

    UWORD8 u1_num_mb_part, u1_mb_mc_mode, u1_sub_mb, u1_mbpred_mode = 5 + u1_mb_type;
    UWORD32 u4_mb_mc_mode = 0, u4_mb_pred_mode = 0;
    WORD32 ret;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    p_curr_ctxt->u1_mb_type = CAB_NON_BD16x16;
    u1_sub_mb = !(u1_mb_type ^ B_8x8);

    {
        UWORD8 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
        UWORD8 *pu1_num_ref_idx_lx_active = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active;
        UWORD8 uc_field = ps_cur_mb_info->u1_mb_field_decodingflag;
        UWORD8 u1_mbaff_field = (u1_mbaff & uc_field);
        u1_dec_ref_l0 = (pu1_num_ref_idx_lx_active[0] << u1_mbaff_field) - 1;
        u1_dec_ref_l1 = (pu1_num_ref_idx_lx_active[1] << u1_mbaff_field) - 1;
    }

    if(u1_sub_mb)
    {
        const UWORD8 u1_colz = ((PRED_8x8) << 6);
        UWORD8 uc_i;
        u1_mb_mc_mode = 0;
        u1_num_mb_part = 4;
        /* Reading the subMB type */
        for(uc_i = 0; uc_i < 4; uc_i++)
        {
            UWORD8 u1_sub_mb_mode, u1_subMbPredModes;
            u1_sub_mb_mode =
                ih264d_parse_submb_type_cabac(1, ps_cab_env, ps_bitstrm, ps_dec->p_sub_mb_type_t);

            if(u1_sub_mb_mode > 12) return ERROR_SUB_MB_TYPE;

            u1_subMbPredModes = pu1_sub_mb_pred_modes[u1_sub_mb_mode];
            u4_mb_mc_mode = (u4_mb_mc_mode << 8) | pu1_sub_mb_mc_mode[u1_sub_mb_mode];
            u4_mb_pred_mode = (u4_mb_pred_mode << 8) | u1_subMbPredModes;
            *pi1_ref_idx_l0++ = (u1_subMbPredModes & PRED_L0) ? u1_dec_ref_l0 : -1;
            *pi1_ref_idx_l1++ = (u1_subMbPredModes & PRED_L1) ? u1_dec_ref_l1 : -1;
            COPYTHECONTEXT("sub_mb_type", u1_sub_mb_mode);
            /* Storing collocated Mb and SubMb mode information */
            *pu1_col_info++ = (u1_colz | (pu1_sub_mb_mc_mode[u1_sub_mb_mode] << 4));
            if(u1_sub_mb_mode != B_DIRECT_8x8)
            {
                if(u1_sub_mb_mode > B_BI_8x8)
                {
                    ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 0;
                }
            }
            else if(!ps_dec->s_high_profile.u1_direct_8x8_inference_flag)
            {
                ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 0;
            }
        }
        pi1_ref_idx_l0 -= 4;
        pi1_ref_idx_l1 -= 4;
    }
    else
    {
        UWORD8 u1_mb_pred_mode_part0 = pu1_mb_pred_modes[0][u1_mbpred_mode];
        UWORD8 u1_mb_pred_mode_part1 = pu1_mb_pred_modes[1][u1_mbpred_mode];
        u1_mb_mc_mode = ps_cur_mb_info->u1_mb_mc_mode;
        u1_num_mb_part = pu1_num_mb_part[u1_mb_mc_mode];
        /* Storing collocated Mb and SubMb mode information */
        *pu1_col_info++ = (u1_mb_mc_mode << 6);
        if(u1_mb_mc_mode) *pu1_col_info++ = (u1_mb_mc_mode << 6);
        u4_mb_mc_mode = u1_mb_mc_mode | (u1_mb_mc_mode << 8);
        u4_mb_mc_mode <<= 16;
        u4_mb_pred_mode = ((u1_mb_pred_mode_part0 << 8) | u1_mb_pred_mode_part1) << 16;

        *pi1_ref_idx_l0++ = (u1_mb_pred_mode_part0 & PRED_L0) ? u1_dec_ref_l0 : -1;
        *pi1_ref_idx_l0-- = (u1_mb_pred_mode_part1 & PRED_L0) ? u1_dec_ref_l0 : -1;
        *pi1_ref_idx_l1++ = (u1_mb_pred_mode_part0 & PRED_L1) ? u1_dec_ref_l1 : -1;
        *pi1_ref_idx_l1-- = (u1_mb_pred_mode_part1 & PRED_L1) ? u1_dec_ref_l1 : -1;
    }

    /*Adding SVC extension code to get Motion Prediction Flags*/
    {
        UWORD8 uc_i, u1_mvp_l1, u1_mvp_l0;
        UWORD8 *pu1_motion_pred_flag_l0;
        UWORD8 *pu1_motion_pred_flag_l1;
        WORD8 *pi1_ref_idx;
        WORD8 *pi1_lft_cxt = ps_dec->pi1_left_ref_idx_ctxt_inc;
        WORD8 *pi1_top_cxt = p_curr_ctxt->i1_ref_idx;

        pu1_motion_pred_flag_l0 = &ps_svc_cur_mb_info->au1_motion_pred_flag[0];
        *pu1_motion_pred_flag_l0 = 0;
        pu1_motion_pred_flag_l1 = &ps_svc_cur_mb_info->au1_motion_pred_flag[1];
        *pu1_motion_pred_flag_l1 = 0;

        if(ps_svc_cur_mb_info->u1_crop_window_flag &&
           ps_svc_slice_params->u1_adaptive_motion_prediction_flag)
        {
            pi1_ref_idx = pi1_ref_idx_l0;
            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++, pi1_ref_idx++)
            {
                if(*pi1_ref_idx > 0)
                {
                    u1_mvp_l0 = ih264d_decode_bin(0, ps_svc_lyr_dec->ps_motion_prediction_flag_l0,
                                                  ps_bitstrm, ps_cab_env);
                    COPYTHECONTEXT("SVC ext: u1_motion_prediction_flag_l0", u1_mvp_l0);

                    *pu1_motion_pred_flag_l0 |= (u1_mvp_l0 << uc_i);
                    if((u1_mvp_l0 & 0x01))
                    {
                        *pi1_ref_idx = -1;
                    }
                }
            }
            pi1_ref_idx = pi1_ref_idx_l1;
            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++, pi1_ref_idx++)
            {
                if(*pi1_ref_idx > 0)
                {
                    u1_mvp_l1 = ih264d_decode_bin(0, ps_svc_lyr_dec->ps_motion_prediction_flag_l1,
                                                  ps_bitstrm, ps_cab_env);
                    COPYTHECONTEXT("SVC ext: u1_motion_prediction_flag_l1", u1_mvp_l1);

                    *pu1_motion_pred_flag_l1 |= (u1_mvp_l1 << uc_i);

                    if((u1_mvp_l1 & 0x01))
                    {
                        *pi1_ref_idx = -1;
                    }
                }
            }
        }
        else if(ps_svc_cur_mb_info->u1_crop_window_flag)
        {
            *pu1_motion_pred_flag_l0 = ps_svc_slice_params->u1_default_motion_prediction_flag
                                           ? ((1 << u1_num_mb_part) - 1)
                                           : 0;
            *pu1_motion_pred_flag_l1 = ps_svc_slice_params->u1_default_motion_prediction_flag
                                           ? ((1 << u1_num_mb_part) - 1)
                                           : 0;
            if(ps_svc_slice_params->u1_default_motion_prediction_flag)
            {
                pi1_ref_idx_l0[0] = -1;
                pi1_ref_idx_l0[1] = -1;
                pi1_ref_idx_l0[2] = -1;
                pi1_ref_idx_l0[3] = -1;

                pi1_ref_idx_l1[0] = -1;
                pi1_ref_idx_l1[1] = -1;
                pi1_ref_idx_l1[2] = -1;
                pi1_ref_idx_l1[3] = -1;
            }
        }

        ret = ih264d_parse_ref_idx_cabac(u1_num_mb_part, 0, u1_dec_ref_l0, u1_mb_mc_mode,
                                         pi1_ref_idx_l0, pi1_lft_cxt, pi1_top_cxt, ps_cab_env,
                                         ps_bitstrm, ps_dec->p_ref_idx_t);

        if(ret != OK)
        {
            return ret;
        }
        ret = ih264d_parse_ref_idx_cabac(u1_num_mb_part, 2, u1_dec_ref_l1, u1_mb_mc_mode,
                                         pi1_ref_idx_l1, pi1_lft_cxt, pi1_top_cxt, ps_cab_env,
                                         ps_bitstrm, ps_dec->p_ref_idx_t);

        if(ret != OK)
        {
            return ret;
        }
    }
    /* Read MotionVectors */
    {
        const UWORD8 *pu1_top_left_sub_mb_indx;
        UWORD8 uc_j, uc_lx;
        UWORD8 u1_mb_part_wd, u1_mb_part_ht;

        const UWORD8 *pu1_sub_mb_indx_mod =
            (const UWORD8 *) gau1_ih264d_submb_indx_mod + (u1_sub_mb * 6);
        const UWORD8 *pu1_sub_mb_partw = (const UWORD8 *) gau1_ih264d_submb_partw;
        const UWORD8 *pu1_sub_mb_parth = (const UWORD8 *) gau1_ih264d_submb_parth;
        const UWORD8 *pu1_num_sub_mb_part = (const UWORD8 *) gau1_ih264d_num_submb_part;
        const UWORD8 *pu1_mb_partw = (const UWORD8 *) gau1_ih264d_mb_partw;
        const UWORD8 *pu1_mb_parth = (const UWORD8 *) gau1_ih264d_mb_parth;

        UWORD8 u1_p_idx = 0;
        UWORD8 u1_num_submb_part;
        parse_part_params_t *ps_part;
        mv_pred_t *ps_mv_start = ps_dec->ps_mv_cur + (u1_mb_num << 4);
        ps_part = ps_dec->ps_part;

        /* Default initialization for non subMb case */
        u1_mb_part_wd = pu1_mb_partw[u1_mb_mc_mode];
        u1_mb_part_ht = pu1_mb_parth[u1_mb_mc_mode];
        u1_num_submb_part = 1;

        /* Decoding the MV for the subMB */
        for(uc_lx = 0; uc_lx < 2; uc_lx++)
        {
            UWORD8 u1_sub_mb_num = 0;
            UWORD32 u4_mb_pred_mode_tmp = u4_mb_pred_mode;
            UWORD32 u4_mb_mc_mode_tmp = u4_mb_mc_mode;
            UWORD8 u1_mb_mc_mode_1, u1_pred_mode, uc_i;
            UWORD16 u2_sub_mb_num = 0x028A;
            UWORD8 u1_b2 = uc_lx << 1;
            u1_pred_mode = (uc_lx) ? PRED_L1 : PRED_L0;
            /* Default for Cabac */
            pu1_top_left_sub_mb_indx = pu1_sub_mb_indx_mod + (u1_mb_mc_mode << 1);
            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++)
            {
                WORD8 i1_pred = (UWORD8) (u4_mb_pred_mode_tmp >> 24);
                u1_mb_mc_mode_1 = (UWORD8) (u4_mb_mc_mode_tmp >> 24);
                u4_mb_pred_mode_tmp <<= 8;
                u4_mb_mc_mode_tmp <<= 8;

                /* subMb prediction mode */
                if(u1_sub_mb)
                {
                    u1_mb_part_wd = pu1_sub_mb_partw[u1_mb_mc_mode_1];
                    u1_mb_part_ht = pu1_sub_mb_parth[u1_mb_mc_mode_1];
                    u1_sub_mb_num = u2_sub_mb_num >> 12;
                    pu1_top_left_sub_mb_indx = pu1_sub_mb_indx_mod + (u1_mb_mc_mode_1 << 1);
                    u1_num_submb_part = pu1_num_sub_mb_part[u1_mb_mc_mode_1];
                    u2_sub_mb_num = u2_sub_mb_num << 4;
                }

                for(uc_j = 0; uc_j < u1_num_submb_part; uc_j++, pu1_top_left_sub_mb_indx++)
                {
                    mv_pred_t *ps_mv;
                    u1_sub_mb_num = u1_sub_mb_num + *pu1_top_left_sub_mb_indx;
                    ps_mv = ps_mv_start + u1_sub_mb_num;

                    /* Storing Info for partitions, writing only once */
                    if(uc_lx)
                    {
                        ps_part->u1_is_direct = (!i1_pred);
                        ps_part->u1_pred_mode = i1_pred;
                        ps_part->u1_sub_mb_num = u1_sub_mb_num;
                        ps_part->u1_partheight = u1_mb_part_ht;
                        ps_part->u1_partwidth = u1_mb_part_wd;

                        /* Increment partition Index */
                        u1_p_idx++;
                        ps_part++;
                    }

                    ih264d_get_mvd_cabac(u1_sub_mb_num, u1_b2, u1_mb_part_wd, u1_mb_part_ht,
                                         (UWORD8) (i1_pred & u1_pred_mode), ps_dec, ps_mv);
                }
            }
        }
        /* write back to the scratch partition info */
        ps_dec->ps_part = ps_part;
        ps_parse_mb_data->u1_num_part = u1_sub_mb ? u1_p_idx : u1_num_mb_part;
    }

    return OK;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_bmb_non_direct_cavlc\endif
*
* \brief
*    Implements sub_mb_pred() of 7.3.5.2. & mb_pred() of 7.3.5.1
*
* \return
*    None.
*
**************************************************************************
*/
WORD32 isvcd_parse_bmb_non_direct_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                        dec_mb_info_t *ps_cur_mb_info,
                                        dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num,
                                        UWORD8 u1_num_mbsNby2)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD8 *pu1_sub_mb_pred_modes = (UWORD8 *) (gau1_ih264d_submb_pred_modes) + 4;
    const UWORD8(*pu1_mb_pred_modes)[32] = (const UWORD8(*)[32]) gau1_ih264d_mb_pred_modes;
    const UWORD8 *pu1_num_mb_part = (const UWORD8 *) gau1_ih264d_num_mb_part;
    const UWORD8 *pu1_sub_mb_mc_mode = (const UWORD8 *) (gau1_ih264d_submb_mc_mode) + 4;
    parse_pmbarams_t *ps_parse_mb_data = ps_dec->ps_parse_mb_data + u1_num_mbsNby2;
    UWORD8 *pu1_col_info = ps_parse_mb_data->u1_col_info;
    WORD8(*pi1_ref_idx)[MAX_REFIDX_INFO_PER_MB] = ps_parse_mb_data->i1_ref_idx;
    UWORD8 u1_mb_type = ps_cur_mb_info->u1_mb_type;
    UWORD8 u1_mb_mc_mode, u1_num_mb_part, u1_sub_mb = !(u1_mb_type ^ B_8x8);
    UWORD32 u4_mb_mc_mode = 0, u4_mb_pred_mode = 0;
    WORD32 ret = OK;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    if(u1_sub_mb)
    {
        UWORD8 uc_i;
        u1_mb_mc_mode = 0;
        u1_num_mb_part = 4;
        /* Reading the subMB type */
        for(uc_i = 0; uc_i < 4; uc_i++)
        {
            UWORD32 ui_sub_mb_mode;
            // Inlined ih264d_uev
            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
            UWORD32 u4_word, u4_ldz;

            /***************************************************************/
            /* Find leading zeros in next 32 bits                          */
            /***************************************************************/
            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
            u4_ldz = CLZ(u4_word);
            /* Flush the ps_bitstrm */
            u4_bitstream_offset += (u4_ldz + 1);
            /* Read the suffix from the ps_bitstrm */
            u4_word = 0;
            if(u4_ldz) GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
            *pu4_bitstrm_ofst = u4_bitstream_offset;
            ui_sub_mb_mode = ((1 << u4_ldz) + u4_word - 1);
            // Inlined ih264d_uev
            if(ui_sub_mb_mode > 12)
                return ERROR_SUB_MB_TYPE;
            else
            {
                UWORD8 u1_subMbPredMode = pu1_sub_mb_pred_modes[ui_sub_mb_mode];
                u4_mb_mc_mode = (u4_mb_mc_mode << 8) | pu1_sub_mb_mc_mode[ui_sub_mb_mode];
                u4_mb_pred_mode = (u4_mb_pred_mode << 8) | u1_subMbPredMode;
                pi1_ref_idx[0][uc_i] = ((u1_subMbPredMode & PRED_L0) - 1) >> 1;
                pi1_ref_idx[1][uc_i] = ((u1_subMbPredMode & PRED_L1) - 1) >> 1;
                COPYTHECONTEXT("sub_mb_type", u1_subMbPredMode);
            }
            /* Storing collocated Mb and SubMb mode information */
            *pu1_col_info++ = ((PRED_8x8) << 6) | ((pu1_sub_mb_mc_mode[ui_sub_mb_mode] << 4));
            if(ui_sub_mb_mode != B_DIRECT_8x8)
            {
                if(ui_sub_mb_mode > B_BI_8x8)
                {
                    ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 0;
                }
            }
            else if(!ps_dec->s_high_profile.u1_direct_8x8_inference_flag)
            {
                ps_dec->s_high_profile.u1_no_submb_part_size_lt8x8_flag = 0;
            }
        }
    }
    else
    {
        UWORD8 u1_mb_pred_mode_idx = 5 + u1_mb_type;
        UWORD8 u1_mb_pred_mode_part0 = pu1_mb_pred_modes[0][u1_mb_pred_mode_idx];
        UWORD8 u1_mb_pred_mode_part1 = pu1_mb_pred_modes[1][u1_mb_pred_mode_idx];
        u1_mb_mc_mode = ps_cur_mb_info->u1_mb_mc_mode;
        u1_num_mb_part = pu1_num_mb_part[u1_mb_mc_mode];

        pi1_ref_idx[0][0] = ((u1_mb_pred_mode_part0 & PRED_L0) - 1) >> 1;
        pi1_ref_idx[1][0] = ((u1_mb_pred_mode_part0 & PRED_L1) - 1) >> 1;
        pi1_ref_idx[0][1] = ((u1_mb_pred_mode_part1 & PRED_L0) - 1) >> 1;
        pi1_ref_idx[1][1] = ((u1_mb_pred_mode_part1 & PRED_L1) - 1) >> 1;

        u4_mb_pred_mode = (u1_mb_pred_mode_part0 << 8) | u1_mb_pred_mode_part1;
        u4_mb_mc_mode = u1_mb_mc_mode | (u1_mb_mc_mode << 8);
        u4_mb_mc_mode <<= 16;
        u4_mb_pred_mode <<= 16;

        /* Storing collocated Mb and SubMb mode information */
        *pu1_col_info++ = (u1_mb_mc_mode << 6);
        if(u1_mb_mc_mode) *pu1_col_info++ = (u1_mb_mc_mode << 6);
    }

    {
        UWORD8 uc_i;
        UWORD8 *pu1_motion_pred_flag_l0;
        UWORD8 *pu1_motion_pred_flag_l1;
        UWORD8 u1_mvp_l1;
        UWORD8 u1_mvp_l0;

        pu1_motion_pred_flag_l0 = &ps_svc_cur_mb_info->au1_motion_pred_flag[0];
        *pu1_motion_pred_flag_l0 = 0;
        pu1_motion_pred_flag_l1 = &ps_svc_cur_mb_info->au1_motion_pred_flag[1];
        *pu1_motion_pred_flag_l1 = 0;

        if(ps_svc_cur_mb_info->u1_crop_window_flag &&
           ps_svc_slice_params->u1_adaptive_motion_prediction_flag)
        {
            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++)
            {
                if(pi1_ref_idx[0][uc_i] > -1)
                {
                    u1_mvp_l0 = ih264d_get_bit_h264(ps_bitstrm);
                    COPYTHECONTEXT("SVC ext: ps_motion_prediction_flag_l0", u1_mvp_l0);
                    *pu1_motion_pred_flag_l0 |= (u1_mvp_l0 << uc_i);
                }
            }

            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++)
            {
                if(pi1_ref_idx[1][uc_i] > -1)
                {
                    u1_mvp_l1 = ih264d_get_bit_h264(ps_bitstrm);
                    COPYTHECONTEXT("SVC ext: ps_motion_prediction_flag_l1", u1_mvp_l1);
                    *pu1_motion_pred_flag_l1 |= (u1_mvp_l1 << uc_i);
                }
            }
        }
        else if(ps_svc_cur_mb_info->u1_crop_window_flag)
        {
            *pu1_motion_pred_flag_l0 = ps_svc_slice_params->u1_default_motion_prediction_flag
                                           ? ((1 << u1_num_mb_part) - 1)
                                           : 0;
            *pu1_motion_pred_flag_l1 = ps_svc_slice_params->u1_default_motion_prediction_flag
                                           ? ((1 << u1_num_mb_part) - 1)
                                           : 0;
        }

        {
            UWORD8 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
            UWORD8 uc_field = ps_cur_mb_info->u1_mb_field_decodingflag;
            UWORD8 *pu1_num_ref_idx_lx_active = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active;
            const UWORD8 u1_mbaff_field = (u1_mbaff & uc_field);
            UWORD8 u4_num_ref_idx_lx_active;

            u4_num_ref_idx_lx_active = (pu1_num_ref_idx_lx_active[0] << u1_mbaff_field) - 1;

            if(u4_num_ref_idx_lx_active)
            {
                if(1 == u4_num_ref_idx_lx_active)
                    isvcd_parse_bmb_ref_index_cavlc_range1(u1_num_mb_part, ps_bitstrm,
                                                           pi1_ref_idx[0], u4_num_ref_idx_lx_active,
                                                           pu1_motion_pred_flag_l0);
                else
                {
                    isvcd_parse_bmb_ref_index_cavlc(u1_num_mb_part, ps_bitstrm, pi1_ref_idx[0],
                                                    u4_num_ref_idx_lx_active,
                                                    pu1_motion_pred_flag_l0);
                    if(ret != OK) return ret;
                }
            }

            u4_num_ref_idx_lx_active = (pu1_num_ref_idx_lx_active[1] << u1_mbaff_field) - 1;

            if(u4_num_ref_idx_lx_active)
            {
                if(1 == u4_num_ref_idx_lx_active)
                    isvcd_parse_bmb_ref_index_cavlc_range1(u1_num_mb_part, ps_bitstrm,
                                                           pi1_ref_idx[1], u4_num_ref_idx_lx_active,
                                                           pu1_motion_pred_flag_l1);
                else
                {
                    ret = isvcd_parse_bmb_ref_index_cavlc(u1_num_mb_part, ps_bitstrm,
                                                          pi1_ref_idx[1], u4_num_ref_idx_lx_active,
                                                          pu1_motion_pred_flag_l1);
                    if(ret != OK) return ret;
                }
            }
        }
    }
    /* Read MotionVectors */
    {
        const UWORD8 *pu1_top_left_sub_mb_indx;
        const UWORD8 *pu1_sub_mb_indx_mod =
            (const UWORD8 *) (gau1_ih264d_submb_indx_mod) + (u1_sub_mb * 6);
        const UWORD8 *pu1_sub_mb_partw = (const UWORD8 *) gau1_ih264d_submb_partw;
        const UWORD8 *pu1_sub_mb_parth = (const UWORD8 *) gau1_ih264d_submb_parth;
        const UWORD8 *pu1_num_sub_mb_part = (const UWORD8 *) gau1_ih264d_num_submb_part;
        const UWORD8 *pu1_mb_partw = (const UWORD8 *) gau1_ih264d_mb_partw;
        const UWORD8 *pu1_mb_parth = (const UWORD8 *) gau1_ih264d_mb_parth;
        UWORD8 u1_p_idx = 0, u1_num_submb_part, uc_lx;
        parse_part_params_t *ps_part;
        mv_pred_t *ps_mv_start = ps_dec->ps_mv_cur + (u1_mb_num << 4);
        UWORD8 u1_mb_part_wd, u1_mb_part_ht;

        ps_part = ps_dec->ps_part;
        /* Default Initialization for Non subMb Case Mode */
        u1_mb_part_wd = pu1_mb_partw[u1_mb_mc_mode];
        u1_mb_part_ht = pu1_mb_parth[u1_mb_mc_mode];
        u1_num_submb_part = 1;

        /* Decoding the MV for the subMB */
        for(uc_lx = 0; uc_lx < 2; uc_lx++)
        {
            UWORD8 u1_sub_mb_num = 0, u1_pred_mode, uc_i;
            UWORD32 u4_mb_mc_mode_tmp = u4_mb_mc_mode;
            UWORD32 u4_mb_pred_mode_tmp = u4_mb_pred_mode;
            UWORD16 u2_sub_mb_num = 0x028A;  // for sub mb case
            UWORD8 u1_b2 = uc_lx << 1;
            u1_pred_mode = (uc_lx) ? PRED_L1 : PRED_L0;
            pu1_top_left_sub_mb_indx = pu1_sub_mb_indx_mod + (u1_mb_mc_mode << 1);

            for(uc_i = 0; uc_i < u1_num_mb_part; uc_i++)
            {
                UWORD8 u1_mb_mc_mode, uc_j;
                UWORD8 i1_pred = u4_mb_pred_mode_tmp >> 24;
                u1_mb_mc_mode = u4_mb_mc_mode_tmp >> 24;
                u4_mb_pred_mode_tmp <<= 8;
                u4_mb_mc_mode_tmp <<= 8;
                /* subMb prediction mode */
                if(u1_sub_mb)
                {
                    u1_mb_part_wd = pu1_sub_mb_partw[u1_mb_mc_mode];
                    u1_mb_part_ht = pu1_sub_mb_parth[u1_mb_mc_mode];
                    u1_sub_mb_num = u2_sub_mb_num >> 12;
                    u1_num_submb_part = pu1_num_sub_mb_part[u1_mb_mc_mode];
                    pu1_top_left_sub_mb_indx = pu1_sub_mb_indx_mod + (u1_mb_mc_mode << 1);
                    u2_sub_mb_num <<= 4;
                }
                for(uc_j = 0; uc_j < u1_num_submb_part; uc_j++, pu1_top_left_sub_mb_indx++)
                {
                    mv_pred_t *ps_mv;
                    u1_sub_mb_num = u1_sub_mb_num + *pu1_top_left_sub_mb_indx;
                    ps_mv = ps_mv_start + u1_sub_mb_num;

                    /* Storing Info for partitions, writing only once */
                    if(uc_lx)
                    {
                        ps_part->u1_is_direct = (!i1_pred);
                        ps_part->u1_pred_mode = i1_pred;
                        ps_part->u1_sub_mb_num = u1_sub_mb_num;
                        ps_part->u1_partheight = u1_mb_part_ht;
                        ps_part->u1_partwidth = u1_mb_part_wd;
                        /* Increment partition Index */
                        u1_p_idx++;
                        ps_part++;
                    }

                    if(i1_pred & u1_pred_mode)
                    {
                        WORD16 i2_mvx, i2_mvy;
                        // inlining ih264d_sev
                        {
                            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
                            UWORD32 u4_word, u4_ldz, u4_abs_val;

                            /***************************************************************/
                            /* Find leading zeros in next 32 bits                          */
                            /***************************************************************/
                            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
                            u4_ldz = CLZ(u4_word);

                            /* Flush the ps_bitstrm */
                            u4_bitstream_offset += (u4_ldz + 1);

                            /* Read the suffix from the ps_bitstrm */
                            u4_word = 0;
                            if(u4_ldz)
                                GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);

                            *pu4_bitstrm_ofst = u4_bitstream_offset;
                            u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

                            if(u4_word & 0x1)
                                i2_mvx = (-(WORD32) u4_abs_val);
                            else
                                i2_mvx = (u4_abs_val);
                        }
                        // inlinined ih264d_sev
                        {
                            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
                            UWORD32 u4_word, u4_ldz, u4_abs_val;

                            /***************************************************************/
                            /* Find leading zeros in next 32 bits                          */
                            /***************************************************************/
                            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
                            u4_ldz = CLZ(u4_word);

                            /* Flush the ps_bitstrm */
                            u4_bitstream_offset += (u4_ldz + 1);

                            /* Read the suffix from the ps_bitstrm */
                            u4_word = 0;
                            if(u4_ldz)
                                GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);

                            *pu4_bitstrm_ofst = u4_bitstream_offset;
                            u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

                            if(u4_word & 0x1)
                                i2_mvy = (-(WORD32) u4_abs_val);
                            else
                                i2_mvy = (u4_abs_val);
                        }
                        // inlinined ih264d_sev
                        /* Storing Mv residuals */
                        ps_mv->i2_mv[u1_b2] = i2_mvx;
                        ps_mv->i2_mv[u1_b2 + 1] = i2_mvy;
                    }
                }
            }
        }
        /* write back to the scratch partition info */
        ps_dec->ps_part = ps_part;
        ps_parse_mb_data->u1_num_part = u1_sub_mb ? u1_p_idx : u1_num_mb_part;
    }
    return OK;
}

/*!
**************************************************************************
* \if Function name : ih264d_decode_ebslice \endif
*
* \brief
*    Decodes an EB Slice
*
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_ebslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD32 i_status = OK;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_dec->ps_cur_slice;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    dec_subset_seq_params_t *ps_sps_svc_ext = NULL;
    dec_nal_unit_svc_ext_params_t *ps_nal_svc_ext = NULL;
    UWORD8 u1_ref_idx_re_flag_lx;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD64 u8_ref_idx_l0, u8_ref_idx_l1;
    UWORD32 u4_temp;
    WORD32 i_temp;
    WORD32 ret;

    ps_nal_svc_ext = ps_svc_lyr_dec->ps_nal_svc_ext;
    ps_subset_seq = ps_svc_lyr_dec->ps_cur_subset_sps;
    ps_sps_svc_ext = &ps_subset_seq->s_sps_svc_ext;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    /*--------------------------------------------------------------------*/
    /* Read remaining contents of the slice header                        */
    /*--------------------------------------------------------------------*/
    {
        WORD8 *pi1_buf;
        WORD16 *pi2_mv = ps_dec->s_default_mv_pred.i2_mv;
        WORD32 *pi4_mv = (WORD32 *) pi2_mv;
        WORD16 *pi16_refFrame;
        pi1_buf = ps_dec->s_default_mv_pred.i1_ref_frame;
        pi16_refFrame = (WORD16 *) pi1_buf;
        *pi4_mv = 0;
        *(pi4_mv + 1) = 0;
        *pi16_refFrame = OUT_OF_RANGE_REF;
        ps_dec->s_default_mv_pred.u1_col_ref_pic_idx = (UWORD8) -1;
        ps_dec->s_default_mv_pred.u1_pic_type = (UWORD8) -1;
    }

    if(0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id)
    {
        ps_slice->u1_num_ref_idx_active_override_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: num_ref_idx_override_flag",
                       ps_slice->u1_num_ref_idx_active_override_flag);

        u8_ref_idx_l0 = ps_dec->ps_cur_pps->u1_num_ref_idx_lx_active[0];
        u8_ref_idx_l1 = ps_dec->ps_cur_pps->u1_num_ref_idx_lx_active[1];
        if(ps_slice->u1_num_ref_idx_active_override_flag)
        {
            u8_ref_idx_l0 = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: num_ref_idx_l0_active_minus1", u8_ref_idx_l0 - 1);

            u8_ref_idx_l1 = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("Slice Header SVC ext: num_ref_idx_l1_active_minus1", u8_ref_idx_l1 - 1);
        }

        {
            UWORD8 u1_max_ref_idx = H264_MAX_REF_PICS;
            if(ps_slice->u1_field_pic_flag)
            {
                u1_max_ref_idx = H264_MAX_REF_PICS << 1;
            }
            if((u8_ref_idx_l0 >= u1_max_ref_idx) || (u8_ref_idx_l1 >= u1_max_ref_idx))
            {
                return ERROR_NUM_REF;
            }
            ps_slice->u1_num_ref_idx_lx_active[0] = (UWORD8) u8_ref_idx_l0;
            ps_slice->u1_num_ref_idx_lx_active[1] = (UWORD8) u8_ref_idx_l1;
        }

        ih264d_init_ref_idx_lx_b(ps_dec);
        /* Store the value for future slices in the same picture */
        ps_dec->u1_num_ref_idx_lx_active_prev = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0];

        u1_ref_idx_re_flag_lx = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: ref_pic_list_reordering_flag_l0",
                       u1_ref_idx_re_flag_lx);

        /* Modified temporarily */
        if(u1_ref_idx_re_flag_lx)
        {
            WORD8 ret;
            ps_dec->ps_ref_pic_buf_lx[0] = ps_dec->ps_dpb_mgr->ps_mod_dpb[0];
            ret = ih264d_ref_idx_reordering(ps_dec, 0);
            if(ret == -1) return ERROR_REFIDX_ORDER_T;
        }
        else
            ps_dec->ps_ref_pic_buf_lx[0] = ps_dec->ps_dpb_mgr->ps_init_dpb[0];

        u1_ref_idx_re_flag_lx = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("Slice Header SVC ext: ref_pic_list_reordering_flag_l1",
                       u1_ref_idx_re_flag_lx);

        /* Modified temporarily */
        if(u1_ref_idx_re_flag_lx)
        {
            WORD8 ret;
            ps_dec->ps_ref_pic_buf_lx[1] = ps_dec->ps_dpb_mgr->ps_mod_dpb[1];
            ret = ih264d_ref_idx_reordering(ps_dec, 1);
            if(ret == -1) return ERROR_REFIDX_ORDER_T;
        }
        else
            ps_dec->ps_ref_pic_buf_lx[1] = ps_dec->ps_dpb_mgr->ps_init_dpb[1];

        /* Create refIdx to POC mapping */
        {
            void **ppv_map_ref_idx_to_poc_lx;
            WORD8 idx;
            struct pic_buffer_t *ps_pic;

            ppv_map_ref_idx_to_poc_lx = ps_dec->ppv_map_ref_idx_to_poc + FRM_LIST_L0;
            ppv_map_ref_idx_to_poc_lx[0] = 0;
            ppv_map_ref_idx_to_poc_lx++;
            for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0]; idx++)
            {
                ps_pic = ps_dec->ps_ref_pic_buf_lx[0][idx];
                ppv_map_ref_idx_to_poc_lx[idx] = (ps_pic->pu1_buf1);
            }

            ppv_map_ref_idx_to_poc_lx = ps_dec->ppv_map_ref_idx_to_poc + FRM_LIST_L1;

            ppv_map_ref_idx_to_poc_lx[0] = 0;
            ppv_map_ref_idx_to_poc_lx++;
            for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1]; idx++)
            {
                ps_pic = ps_dec->ps_ref_pic_buf_lx[1][idx];
                ppv_map_ref_idx_to_poc_lx[idx] = (ps_pic->pu1_buf1);
            }

            if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
            {
                void **ppv_map_ref_idx_to_poc_lx_t, **ppv_map_ref_idx_to_poc_lx_b;

                ppv_map_ref_idx_to_poc_lx_t = ps_dec->ppv_map_ref_idx_to_poc + TOP_LIST_FLD_L0;
                ppv_map_ref_idx_to_poc_lx_b = ps_dec->ppv_map_ref_idx_to_poc + BOT_LIST_FLD_L0;

                ppv_map_ref_idx_to_poc_lx_t[0] = 0;
                ppv_map_ref_idx_to_poc_lx_t++;
                ppv_map_ref_idx_to_poc_lx_b[0] = 0;
                ppv_map_ref_idx_to_poc_lx_b++;
                for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0]; idx++)
                {
                    ps_pic = ps_dec->ps_ref_pic_buf_lx[0][idx];
                    ppv_map_ref_idx_to_poc_lx_t[0] = (ps_pic->pu1_buf1);
                    ppv_map_ref_idx_to_poc_lx_b[1] = (ps_pic->pu1_buf1);

                    ppv_map_ref_idx_to_poc_lx_b[0] = (ps_pic->pu1_buf1) + 1;
                    ppv_map_ref_idx_to_poc_lx_t[1] = (ps_pic->pu1_buf1) + 1;

                    ppv_map_ref_idx_to_poc_lx_t += 2;
                    ppv_map_ref_idx_to_poc_lx_b += 2;
                }

                ppv_map_ref_idx_to_poc_lx_t = ps_dec->ppv_map_ref_idx_to_poc + TOP_LIST_FLD_L1;
                ppv_map_ref_idx_to_poc_lx_b = ps_dec->ppv_map_ref_idx_to_poc + BOT_LIST_FLD_L1;

                ppv_map_ref_idx_to_poc_lx_t[0] = 0;
                ppv_map_ref_idx_to_poc_lx_t++;
                ppv_map_ref_idx_to_poc_lx_b[0] = 0;
                ppv_map_ref_idx_to_poc_lx_b++;
                for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1]; idx++)
                {
                    UWORD8 u1_tmp_idx = idx << 1;
                    ps_pic = ps_dec->ps_ref_pic_buf_lx[1][idx];
                    ppv_map_ref_idx_to_poc_lx_t[u1_tmp_idx] = (ps_pic->pu1_buf1);
                    ppv_map_ref_idx_to_poc_lx_b[u1_tmp_idx + 1] = (ps_pic->pu1_buf1);

                    ppv_map_ref_idx_to_poc_lx_b[u1_tmp_idx] = (ps_pic->pu1_buf1) + 1;
                    ppv_map_ref_idx_to_poc_lx_t[u1_tmp_idx + 1] = (ps_pic->pu1_buf1) + 1;
                }
            }

            /* BS is moved post recon gen in SVC ext*/
            if(ps_dec->u4_num_cores >= 2)
            {
                WORD32 num_entries;
                WORD32 size;
                num_entries = MAX_FRAMES;
                if((1 >= ps_dec->ps_cur_sps->u1_num_ref_frames) && (0 == ps_dec->i4_display_delay))
                {
                    num_entries = 1;
                }

                num_entries = ((2 * num_entries) + 1);
                num_entries *= 2;

                size = num_entries * sizeof(void *);
                size += PAD_MAP_IDX_POC * sizeof(void *);

                memcpy((void *) ps_dec->ps_parse_cur_slice->ppv_map_ref_idx_to_poc,
                       ps_dec->ppv_map_ref_idx_to_poc, size);
            }
        }

        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag &&
           (ps_dec->ps_cur_slice->u1_field_pic_flag == 0))
        {
            ih264d_convert_frm_mbaff_list(ps_dec);
        }

        if(ps_pps->u1_wted_bipred_idc == 1)
        {
            if(!ps_nal_svc_ext->u1_no_inter_layer_pred_flag)
            {
                ps_svc_slice_params->u1_base_pred_weight_table_flag =
                    ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("Slice Header SVC ext: u1_base_pred_weight_table_flag",
                               ps_svc_slice_params->u1_base_pred_weight_table_flag);
            }

            if(ps_nal_svc_ext->u1_no_inter_layer_pred_flag ||
               !ps_svc_slice_params->u1_base_pred_weight_table_flag)
            {
                ih264d_parse_pred_weight_table(ps_slice, ps_bitstrm);

                ih264d_form_pred_weight_matrix(ps_dec);
                ps_dec->pu4_wt_ofsts = ps_dec->pu4_wts_ofsts_mat;
            }
        }
        else if(ps_pps->u1_wted_bipred_idc == 2)
        {
            /* Implicit Weighted prediction */
            ps_slice->u2_log2Y_crwd = 0x0505;
            ps_dec->pu4_wt_ofsts = ps_dec->pu4_wts_ofsts_mat;
            ih264d_get_implicit_weights(ps_dec);
        }
        else
            ps_dec->ps_cur_slice->u2_log2Y_crwd = 0;

        ps_dec->ps_parse_cur_slice->u2_log2Y_crwd = ps_dec->ps_cur_slice->u2_log2Y_crwd;

        /* G050 */
        if(ps_slice->u1_nal_ref_idc != 0)
        {
            if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
            {
                dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
                dec_seq_params_t *ps_sps_tmp = ps_pps->ps_sps;
                UWORD8 u1_nal_unit_type_tmp = ps_dec->u1_nal_unit_type;

                ps_pps->ps_sps = ps_dec->ps_cur_sps;
                if(ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag)
                    ps_dec->u1_nal_unit_type = IDR_SLICE_NAL;

                i_temp = ih264d_read_mmco_commands(ps_dec);

                ps_pps->ps_sps = ps_sps_tmp;
                ps_dec->u1_nal_unit_type = u1_nal_unit_type_tmp;

                if(i_temp < 0)
                {
                    return ERROR_DBP_MANAGER_T;
                }
                ps_dec->u4_bitoffset = i_temp;
            }
            else
                ps_bitstrm->u4_ofst += ps_dec->u4_bitoffset;

            if(!ps_sps_svc_ext->u1_slice_header_restriction_flag)
            {
                ps_svc_slice_params->u1_store_ref_base_pic_flag = ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("SPS_SVC_EXT: u1_store_ref_base_pic_flag",
                               ps_svc_slice_params->u1_store_ref_base_pic_flag);

                if(0 != ps_svc_slice_params->u1_store_ref_base_pic_flag)
                {
                    return NOT_OK;
                }
                if(((1 == ps_nal_svc_ext->u1_use_ref_base_pic_flag) ||
                    (1 == ps_svc_slice_params->u1_store_ref_base_pic_flag)) &&
                   (!ps_nal_svc_ext->u1_idr_flag))
                {
                    i_status = isvcd_dec_ref_base_pic_marking(
                        &ps_svc_slice_params->s_ref_base_pic_marking_svc_ext, ps_bitstrm);
                    if(i_status != OK)
                    {
                        return i_status;
                    }
                }
            }
        }
    }
    /* G050 */
    /*Code is present in standard but omitted in the reference code*/
    if(ps_pps->u1_entropy_coding_mode == CABAC)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > MAX_CABAC_INIT_IDC)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        ps_slice->u1_cabac_init_idc = u4_temp;
        COPYTHECONTEXT("Slice Header SVC ext: cabac_init_idc", ps_slice->u1_cabac_init_idc);
    }

    {
        /* Read slice_qp_delta */
        WORD64 i8_temp =
            (WORD64) ps_pps->u1_pic_init_qp + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP))
        {
            return ERROR_INV_RANGE_QP_T;
        }
        ps_slice->u1_slice_qp = (UWORD8) i8_temp;
        COPYTHECONTEXT("Slice Header SVC ext: slice_qp_delta",
                       (WORD8) (ps_slice->u1_slice_qp - ps_pps->u1_pic_init_qp));
    }
    if(ps_pps->u1_deblocking_filter_parameters_present_flag == 1)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        COPYTHECONTEXT("Slice Header SVC ext: disable_deblocking_filter_idc", u4_temp);
        ps_slice->u1_disable_dblk_filter_idc = u4_temp;
        if(u4_temp != 1)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_alpha_c0_offset = i_temp;
            COPYTHECONTEXT("Slice Header SVC ext: slice_alpha_c0_offset_div2",
                           ps_slice->i1_slice_alpha_c0_offset >> 1);

            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_beta_offset = i_temp;
            COPYTHECONTEXT("Slice Header SVC ext: slice_beta_offset_div2",
                           ps_slice->i1_slice_beta_offset >> 1);
        }
        else
        {
            ps_slice->i1_slice_alpha_c0_offset = 0;
            ps_slice->i1_slice_beta_offset = 0;
        }
    }
    else
    {
        ps_slice->u1_disable_dblk_filter_idc = 0;
        ps_slice->i1_slice_alpha_c0_offset = 0;
        ps_slice->i1_slice_beta_offset = 0;
    }

    /* add the remaining part of the code for svc extension from reference */
    ret = isvcd_set_default_slice_header_ext(ps_svc_lyr_dec);
    if(ret != OK)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    ret = isvcd_parse_slice_header(ps_svc_lyr_dec);
    if(ret != OK)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    ps_dec->u1_slice_header_done = 2;
    if(!ps_svc_slice_params->u1_slice_skip_flag)
    {
        if(ps_pps->u1_entropy_coding_mode)
        {
            SWITCHOFFTRACE;
            SWITCHONTRACECABAC;
            ps_svc_lyr_dec->pf_parse_inter_slice_svc_ext =
                isvcd_parse_inter_slice_data_cabac_enh_lyr;
            ps_svc_lyr_dec->pf_parse_inter_mb_svc_ext = isvcd_parse_bmb_cabac;
            isvcd_init_cabac_contexts(B_SLICE, ps_dec);

            if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
                ps_dec->pf_get_mb_info = ih264d_get_mb_info_cabac_mbaff;
            else
                ps_dec->pf_get_mb_info = isvcd_get_mb_info_cabac_nonmbaff;
        }
        else
        {
            SWITCHONTRACE;
            SWITCHOFFTRACECABAC;
            ps_svc_lyr_dec->pf_parse_inter_slice_svc_ext =
                isvcd_parse_inter_slice_data_cavlc_enh_lyr;
            ps_svc_lyr_dec->pf_parse_inter_mb_svc_ext = isvcd_parse_bmb_cavlc;
            if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
                ps_dec->pf_get_mb_info = ih264d_get_mb_info_cavlc_mbaff;
            else
                ps_dec->pf_get_mb_info = isvcd_get_mb_info_cavlc_nonmbaff;
        }
    }
    else
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ret = ih264d_cal_col_pic(ps_dec);
    if(ret != OK) return ret;
    ps_dec->u1_B = 1;
    ps_dec->pf_mvpred_ref_tfr_nby2mb = isvcd_mv_pred_ref_tfr_nby2_ebmb;
    ret = ps_svc_lyr_dec->pf_parse_inter_slice_svc_ext(ps_svc_lyr_dec, ps_slice,
                                                       u2_first_mb_in_slice);
    if(ret != OK) return ret;

    return OK;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_bslice \endif
*
* \brief
*    Decodes a B Slice
*
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_bslice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_dec->ps_cur_slice;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD8 u1_ref_idx_re_flag_lx;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD64 u8_ref_idx_l0, u8_ref_idx_l1;
    UWORD32 u4_temp;
    WORD32 i_temp;
    WORD32 ret;

    /*--------------------------------------------------------------------*/
    /* Read remaining contents of the slice header                        */
    /*--------------------------------------------------------------------*/
    {
        WORD8 *pi1_buf;
        WORD16 *pi2_mv = ps_dec->s_default_mv_pred.i2_mv;
        WORD32 *pi4_mv = (WORD32 *) pi2_mv;
        WORD16 *pi16_refFrame;
        pi1_buf = ps_dec->s_default_mv_pred.i1_ref_frame;
        pi16_refFrame = (WORD16 *) pi1_buf;
        *pi4_mv = 0;
        *(pi4_mv + 1) = 0;
        *pi16_refFrame = OUT_OF_RANGE_REF;
        ps_dec->s_default_mv_pred.u1_col_ref_pic_idx = (UWORD8) -1;
        ps_dec->s_default_mv_pred.u1_pic_type = (UWORD8) -1;
    }

    ps_slice->u1_num_ref_idx_active_override_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SH: num_ref_idx_override_flag", ps_slice->u1_num_ref_idx_active_override_flag);

    u8_ref_idx_l0 = ps_dec->ps_cur_pps->u1_num_ref_idx_lx_active[0];
    u8_ref_idx_l1 = ps_dec->ps_cur_pps->u1_num_ref_idx_lx_active[1];
    if(ps_slice->u1_num_ref_idx_active_override_flag)
    {
        u8_ref_idx_l0 = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SH: num_ref_idx_l0_active_minus1", u8_ref_idx_l0 - 1);

        u8_ref_idx_l1 = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SH: num_ref_idx_l1_active_minus1", u8_ref_idx_l1 - 1);
    }

    {
        UWORD8 u1_max_ref_idx = H264_MAX_REF_PICS;
        if(ps_slice->u1_field_pic_flag)
        {
            u1_max_ref_idx = H264_MAX_REF_PICS << 1;
        }
        if((u8_ref_idx_l0 >= u1_max_ref_idx) || (u8_ref_idx_l1 >= u1_max_ref_idx))
        {
            return ERROR_NUM_REF;
        }
        ps_slice->u1_num_ref_idx_lx_active[0] = (UWORD8) u8_ref_idx_l0;
        ps_slice->u1_num_ref_idx_lx_active[1] = (UWORD8) u8_ref_idx_l1;
    }

    ih264d_init_ref_idx_lx_b(ps_dec);
    /* Store the value for future slices in the same picture */
    ps_dec->u1_num_ref_idx_lx_active_prev = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0];

    u1_ref_idx_re_flag_lx = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SH: ref_pic_list_reordering_flag_l0", u1_ref_idx_re_flag_lx);

    /* Modified temporarily */
    if(u1_ref_idx_re_flag_lx)
    {
        WORD8 ret;
        ps_dec->ps_ref_pic_buf_lx[0] = ps_dec->ps_dpb_mgr->ps_mod_dpb[0];
        ret = ih264d_ref_idx_reordering(ps_dec, 0);
        if(ret == -1) return ERROR_REFIDX_ORDER_T;
    }
    else
        ps_dec->ps_ref_pic_buf_lx[0] = ps_dec->ps_dpb_mgr->ps_init_dpb[0];

    u1_ref_idx_re_flag_lx = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SH: ref_pic_list_reordering_flag_l1", u1_ref_idx_re_flag_lx);

    /* Modified temporarily */
    if(u1_ref_idx_re_flag_lx)
    {
        WORD8 ret;
        ps_dec->ps_ref_pic_buf_lx[1] = ps_dec->ps_dpb_mgr->ps_mod_dpb[1];
        ret = ih264d_ref_idx_reordering(ps_dec, 1);
        if(ret == -1) return ERROR_REFIDX_ORDER_T;
    }
    else
        ps_dec->ps_ref_pic_buf_lx[1] = ps_dec->ps_dpb_mgr->ps_init_dpb[1];

    /* Create refIdx to POC mapping */
    {
        void **ppv_map_ref_idx_to_poc_lx;
        WORD8 idx;
        struct pic_buffer_t *ps_pic;

        ppv_map_ref_idx_to_poc_lx = ps_dec->ppv_map_ref_idx_to_poc + FRM_LIST_L0;
        ppv_map_ref_idx_to_poc_lx[0] = 0;
        ppv_map_ref_idx_to_poc_lx++;
        for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0]; idx++)
        {
            ps_pic = ps_dec->ps_ref_pic_buf_lx[0][idx];
            ppv_map_ref_idx_to_poc_lx[idx] = (ps_pic->pu1_buf1);
        }

        ppv_map_ref_idx_to_poc_lx = ps_dec->ppv_map_ref_idx_to_poc + FRM_LIST_L1;

        ppv_map_ref_idx_to_poc_lx[0] = 0;
        ppv_map_ref_idx_to_poc_lx++;
        for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1]; idx++)
        {
            ps_pic = ps_dec->ps_ref_pic_buf_lx[1][idx];
            ppv_map_ref_idx_to_poc_lx[idx] = (ps_pic->pu1_buf1);
        }

        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
        {
            void **ppv_map_ref_idx_to_poc_lx_t, **ppv_map_ref_idx_to_poc_lx_b;

            ppv_map_ref_idx_to_poc_lx_t = ps_dec->ppv_map_ref_idx_to_poc + TOP_LIST_FLD_L0;
            ppv_map_ref_idx_to_poc_lx_b = ps_dec->ppv_map_ref_idx_to_poc + BOT_LIST_FLD_L0;

            ppv_map_ref_idx_to_poc_lx_t[0] = 0;
            ppv_map_ref_idx_to_poc_lx_t++;
            ppv_map_ref_idx_to_poc_lx_b[0] = 0;
            ppv_map_ref_idx_to_poc_lx_b++;
            for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0]; idx++)
            {
                ps_pic = ps_dec->ps_ref_pic_buf_lx[0][idx];
                ppv_map_ref_idx_to_poc_lx_t[0] = (ps_pic->pu1_buf1);
                ppv_map_ref_idx_to_poc_lx_b[1] = (ps_pic->pu1_buf1);

                ppv_map_ref_idx_to_poc_lx_b[0] = (ps_pic->pu1_buf1) + 1;
                ppv_map_ref_idx_to_poc_lx_t[1] = (ps_pic->pu1_buf1) + 1;

                ppv_map_ref_idx_to_poc_lx_t += 2;
                ppv_map_ref_idx_to_poc_lx_b += 2;
            }

            ppv_map_ref_idx_to_poc_lx_t = ps_dec->ppv_map_ref_idx_to_poc + TOP_LIST_FLD_L1;
            ppv_map_ref_idx_to_poc_lx_b = ps_dec->ppv_map_ref_idx_to_poc + BOT_LIST_FLD_L1;

            ppv_map_ref_idx_to_poc_lx_t[0] = 0;
            ppv_map_ref_idx_to_poc_lx_t++;
            ppv_map_ref_idx_to_poc_lx_b[0] = 0;
            ppv_map_ref_idx_to_poc_lx_b++;
            for(idx = 0; idx < ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1]; idx++)
            {
                UWORD8 u1_tmp_idx = idx << 1;
                ps_pic = ps_dec->ps_ref_pic_buf_lx[1][idx];
                ppv_map_ref_idx_to_poc_lx_t[u1_tmp_idx] = (ps_pic->pu1_buf1);
                ppv_map_ref_idx_to_poc_lx_b[u1_tmp_idx + 1] = (ps_pic->pu1_buf1);

                ppv_map_ref_idx_to_poc_lx_b[u1_tmp_idx] = (ps_pic->pu1_buf1) + 1;
                ppv_map_ref_idx_to_poc_lx_t[u1_tmp_idx + 1] = (ps_pic->pu1_buf1) + 1;
            }
        }
    }

    if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag && (ps_dec->ps_cur_slice->u1_field_pic_flag == 0))
    {
        ih264d_convert_frm_mbaff_list(ps_dec);
    }

    if(ps_pps->u1_wted_bipred_idc == 1)
    {
        ret = ih264d_parse_pred_weight_table(ps_slice, ps_bitstrm);
        if(ret != OK) return ret;
        ih264d_form_pred_weight_matrix(ps_dec);
        ps_dec->pu4_wt_ofsts = ps_dec->pu4_wts_ofsts_mat;
    }
    else if(ps_pps->u1_wted_bipred_idc == 2)
    {
        /* Implicit Weighted prediction */
        ps_slice->u2_log2Y_crwd = 0x0505;
        ps_dec->pu4_wt_ofsts = ps_dec->pu4_wts_ofsts_mat;
        ih264d_get_implicit_weights(ps_dec);
    }
    else
        ps_dec->ps_cur_slice->u2_log2Y_crwd = 0;

    ps_dec->ps_parse_cur_slice->u2_log2Y_crwd = ps_dec->ps_cur_slice->u2_log2Y_crwd;

    /* G050 */
    if(ps_slice->u1_nal_ref_idc != 0)
    {
        if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
        {
            dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
            dec_seq_params_t *ps_sps_tmp = ps_pps->ps_sps;
            UWORD8 u1_nal_unit_type_tmp = ps_dec->u1_nal_unit_type;

            ps_pps->ps_sps = ps_dec->ps_cur_sps;
            if(ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag)
                ps_dec->u1_nal_unit_type = IDR_SLICE_NAL;

            i_temp = ih264d_read_mmco_commands(ps_dec);

            ps_pps->ps_sps = ps_sps_tmp;
            ps_dec->u1_nal_unit_type = u1_nal_unit_type_tmp;

            if(i_temp < 0)
            {
                return ERROR_DBP_MANAGER_T;
            }
            ps_dec->u4_bitoffset = i_temp;
        }
        else
            ps_bitstrm->u4_ofst += ps_dec->u4_bitoffset;
    }
    /* G050 */

    if(ps_pps->u1_entropy_coding_mode == CABAC)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > MAX_CABAC_INIT_IDC)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        ps_slice->u1_cabac_init_idc = u4_temp;
        COPYTHECONTEXT("SH: cabac_init_idc", ps_slice->u1_cabac_init_idc);
    }
    {
        /* Read slice_qp_delta */
        WORD64 i8_temp =
            (WORD64) ps_pps->u1_pic_init_qp + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP))
        {
            return ERROR_INV_RANGE_QP_T;
        }
        ps_slice->u1_slice_qp = (UWORD8) i8_temp;
        COPYTHECONTEXT("SH: slice_qp_delta",
                       (WORD8) (ps_slice->u1_slice_qp - ps_pps->u1_pic_init_qp));
    }
    if(ps_pps->u1_deblocking_filter_parameters_present_flag == 1)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        COPYTHECONTEXT("SH: disable_deblocking_filter_idc", u4_temp);
        ps_slice->u1_disable_dblk_filter_idc = u4_temp;
        if(u4_temp != 1)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_alpha_c0_offset = i_temp;
            COPYTHECONTEXT("SH: slice_alpha_c0_offset_div2",
                           ps_slice->i1_slice_alpha_c0_offset >> 1);

            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_beta_offset = i_temp;
            COPYTHECONTEXT("SH: slice_beta_offset_div2", ps_slice->i1_slice_beta_offset >> 1);
        }
        else
        {
            ps_slice->i1_slice_alpha_c0_offset = 0;
            ps_slice->i1_slice_beta_offset = 0;
        }
    }
    else
    {
        ps_slice->u1_disable_dblk_filter_idc = 0;
        ps_slice->i1_slice_alpha_c0_offset = 0;
        ps_slice->i1_slice_beta_offset = 0;
    }

    ps_dec->u1_slice_header_done = 2;

    if(ps_pps->u1_entropy_coding_mode)
    {
        SWITCHOFFTRACE;
        SWITCHONTRACECABAC;
        ps_svc_lyr_dec->pf_parse_svc_inter_slice = isvcd_parse_inter_slice_data_cabac;
        ps_dec->pf_parse_inter_mb = ih264d_parse_bmb_cabac;
        ih264d_init_cabac_contexts(B_SLICE, ps_dec);

        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
            ps_dec->pf_get_mb_info = ih264d_get_mb_info_cabac_mbaff;
        else
            ps_dec->pf_get_mb_info = isvcd_get_mb_info_cabac_nonmbaff;
    }
    else
    {
        SWITCHONTRACE;
        SWITCHOFFTRACECABAC;
        ps_svc_lyr_dec->pf_parse_svc_inter_slice = isvcd_parse_inter_slice_data_cavlc;
        ps_dec->pf_parse_inter_mb = ih264d_parse_bmb_cavlc;
        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
            ps_dec->pf_get_mb_info = ih264d_get_mb_info_cavlc_mbaff;
        else
            ps_dec->pf_get_mb_info = isvcd_get_mb_info_cavlc_nonmbaff;
    }

    ret = ih264d_cal_col_pic(ps_dec);
    if(ret != OK) return ret;
    ps_dec->u1_B = 1;
    ps_dec->pf_mvpred_ref_tfr_nby2mb = ih264d_mv_pred_ref_tfr_nby2_bmb;
    ret = ps_svc_lyr_dec->pf_parse_svc_inter_slice(ps_svc_lyr_dec, ps_slice, u2_first_mb_in_slice);
    if(ret != OK) return ret;

    return OK;
}
