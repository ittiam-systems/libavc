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
 *  isvcd_process_epslice.c
 *
 * @brief
 *  Contains routines that decode a I slice type
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <assert.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "isvcd_structs.h"
#include "ih264d_defs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_mvpred.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "isvcd_process_epslice.h"
#include "ih264d_error_handler.h"
#include "ih264d_cabac.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "ih264d_parse_slice.h"
#include "ih264d_utils.h"
#include "ih264d_parse_islice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_intra_mb.h"
#include "isvcd_mode_mv_resamp.h"
#include "ih264_debug.h"

/*!
 **************************************************************************
 * \if Function name : isvcd_retrive_infer_mode_mv \endif
 *
 * \brief
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
void isvcd_retrive_infer_mode_mv(svc_dec_lyr_struct_t *ps_svc_lyr_dec, mv_pred_t *ps_mvpred,
                                 UWORD8 u1_lx, UWORD8 u1_sub_mb_num)
{
    mode_motion_ctxt_t *ps_ctxt;
    mv_pred_t *ps_motion_pred;
    UWORD8 u1_tmp_lx = (u1_lx << 1);

    ps_ctxt = (mode_motion_ctxt_t *) ps_svc_lyr_dec->pv_mode_mv_sample_ctxt;
    ps_motion_pred = ps_ctxt->ps_motion_pred_struct;
    ps_motion_pred += u1_sub_mb_num;
    ps_mvpred->i2_mv[u1_tmp_lx] = ps_motion_pred->i2_mv[u1_tmp_lx];
    ps_mvpred->i2_mv[u1_tmp_lx + 1] = ps_motion_pred->i2_mv[u1_tmp_lx + 1];

    return;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_interlyr_motion_mode_pred \endif
 *
 * \brief
 *
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_interlyr_motion_mode_pred(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                       dec_mb_info_t *ps_cur_mb_info,
                                       dec_svc_mb_info_t *ps_svc_cur_mb_info,
                                       parse_pmbarams_t *ps_mb_part_info,
                                       parse_part_params_t *ps_part)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD32 i4_inter_layer_pred_req_flag;
    WORD32 i4_listx;
    WORD32 i4_mb_mode = -1;
    i4_inter_layer_pred_req_flag = SVCD_FALSE;
    i4_listx = (ps_dec->ps_cur_slice->u1_slice_type == B_SLICE) ? 2 : 1;
    /* check Base mode flag and motion predcition flags */
    if(1 == ps_svc_cur_mb_info->u1_base_mode_flag)
    {
        i4_inter_layer_pred_req_flag = SVCD_TRUE;
    }
    else
    {
        UWORD8 u1_mot_pred_flag;

        /* get the packed the motion pred flag of list 0 */
        u1_mot_pred_flag = ps_svc_cur_mb_info->au1_motion_pred_flag[0];

        /* extract the last 4 bits */
        u1_mot_pred_flag &= 0x0F;

        if(0 != u1_mot_pred_flag)
        {
            i4_inter_layer_pred_req_flag = SVCD_TRUE;
        }

        /* check for list 1 flags if required */
        if((2 == i4_listx) && (SVCD_FALSE == i4_inter_layer_pred_req_flag))
        {
            /* get the packed the motion pred flag of list 1 */
            u1_mot_pred_flag = ps_svc_cur_mb_info->au1_motion_pred_flag[1];

            /* extract the last 4 bits */
            u1_mot_pred_flag &= 0x0F;

            if(0 != u1_mot_pred_flag)
            {
                i4_inter_layer_pred_req_flag = SVCD_TRUE;
            }
        }
    }

    if(SVCD_TRUE == i4_inter_layer_pred_req_flag)
    {
        mode_motion_ctxt_t *ps_ctxt;
        mode_motion_lyr_ctxt *ps_lyr_mem;

        ps_ctxt = (mode_motion_ctxt_t *) ps_svc_lyr_dec->pv_mode_mv_sample_ctxt;
        /* get the current layer ctxt */
        ps_lyr_mem = &ps_ctxt->as_res_lyr_mem[ps_ctxt->i4_res_id];

        {
            ps_ctxt->i4_listx = i4_listx;

            i4_mb_mode =
                ps_lyr_mem->pf_inter_lyr_pred(ps_svc_lyr_dec->pv_mode_mv_sample_ctxt, ps_cur_mb_info,
                                          ps_svc_cur_mb_info, ps_dec, ps_mb_part_info, ps_part);
        }
    }
    return i4_mb_mode;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_mv_pred_ref_tfr_nby2_epmb \endif
 *
 * \brief
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_mv_pred_ref_tfr_nby2_epmb(dec_struct_t *ps_dec, UWORD32 u4_mb_idx, UWORD32 u4_num_mbs)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) ps_dec;
    parse_pmbarams_t *ps_mb_part_info;
    parse_part_params_t *ps_part;
    mv_pred_t *ps_mv_nmb, *ps_mv_nmb_start, *ps_mv_ntop, *ps_mv_ntop_start;
    UWORD32 i, j;
    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    dec_mb_info_t *ps_cur_mb_info;
    dec_svc_mb_info_t *ps_svc_cur_mb_info;
    WORD32 i2_mv_x, i2_mv_y;

    ps_dec->i4_submb_ofst -= (u4_num_mbs - u4_mb_idx) << 4;
    ps_mb_part_info = ps_dec->ps_parse_mb_data;
    ps_part = ps_dec->ps_parse_part_params;

    /* N/2 Mb MvPred and Transfer Setup Loop */
    for(i = u4_mb_idx; i < u4_num_mbs; i++, ps_mb_part_info++)
    {
        UWORD32 u1_colz;
        UWORD32 u1_field;
        mv_pred_t s_mvPred = {0};
        mv_pred_t *ps_mv_pred = &s_mvPred;

        *ps_mv_pred = ps_dec->s_default_mv_pred;

        ps_dec->i4_submb_ofst += SUB_BLK_SIZE;

        /* Restore the slice scratch MbX and MbY context */
        ps_cur_mb_info = ps_dec->ps_nmb_info + i;
        ps_svc_cur_mb_info = ps_svc_lyr_dec->ps_svc_nmb_info + i;
        u1_field = ps_cur_mb_info->u1_mb_field_decodingflag;

        ps_mv_nmb_start = ps_dec->ps_mv_cur + (i << 4);
        ps_dec->u2_mbx = ps_cur_mb_info->u2_mbx;
        ps_dec->u2_mby = ps_cur_mb_info->u2_mby;
        ps_dec->u2_mv_2mb[i & 0x1] = 0;

        /* Look for MV Prediction and Reference Transfer in Non-I Mbs */
        if(!ps_mb_part_info->u4_isI_mb)
        {
            UWORD32 u1_blk_no;
            WORD32 i1_ref_idx, i1_ref_idx1;
            UWORD32 u1_sub_mb_x, u1_sub_mb_y, u1_sub_mb_num;
            UWORD32 u1_num_part, u1_num_ref, u1_wd, u1_ht;
            UWORD32 *pu4_wt_offst, **ppu4_wt_ofst;
            UWORD32 u1_scale_ref, u4_bot_mb;
            WORD8 *pi1_ref_idx = ps_mb_part_info->i1_ref_idx[0];
            pic_buffer_t *ps_ref_frame, **pps_ref_frame;
            deblk_mb_t *ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + i;
            WORD32 i4_mb_mode_svc;
            UWORD8 u1_motion_pred_flag_l0 = ps_svc_cur_mb_info->au1_motion_pred_flag[0];

            /* MB Level initialisations */
            ps_dec->u4_num_pmbair = i >> u1_mbaff;
            ps_dec->u4_mb_idx_mv = i;
            ppu4_wt_ofst = ps_mb_part_info->pu4_wt_offst;
            pps_ref_frame = ps_dec->ps_ref_pic_buf_lx[0];

            i4_mb_mode_svc = isvcd_interlyr_motion_mode_pred(
                ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info, ps_mb_part_info, ps_part);

            if((-1 == i4_mb_mode_svc) || (SVC_INTER_MB == i4_mb_mode_svc))
            {
                ps_mv_ntop_start =
                    ps_mv_nmb_start - (ps_dec->u2_frm_wd_in_mbs << (4 + u1_mbaff)) + 12;

                u1_num_part = ps_mb_part_info->u1_num_part;
                ps_cur_deblk_mb->u1_mb_type |= (u1_num_part > 1) << 1;
                ps_cur_mb_info->u4_pred_info_pkd_idx = ps_dec->u4_pred_info_pkd_idx;
                ps_cur_mb_info->u1_num_pred_parts = 0;

                /****************************************************/
                /* weighted u4_ofst pointer calculations, this loop  */
                /* runs maximum 4 times, even in direct cases       */
                /****************************************************/
                u1_scale_ref = u1_mbaff & u1_field;

                u4_bot_mb = 1 - ps_cur_mb_info->u1_topmb;
                if(ps_dec->ps_cur_pps->u1_wted_pred_flag)
                {
                    u1_num_ref = MIN(u1_num_part, 4);
                    for(u1_blk_no = 0; u1_blk_no < u1_num_ref; u1_blk_no++)
                    {
                        i1_ref_idx = pi1_ref_idx[u1_blk_no];
                        if(u1_scale_ref) i1_ref_idx >>= 1;
                        pu4_wt_offst = (UWORD32 *) &ps_dec->pu4_wt_ofsts[2 * X3(i1_ref_idx)];
                        ppu4_wt_ofst[u1_blk_no] = pu4_wt_offst;
                    }
                }
                else
                {
                    ppu4_wt_ofst[0] = NULL;
                    ppu4_wt_ofst[1] = NULL;
                    ppu4_wt_ofst[2] = NULL;
                    ppu4_wt_ofst[3] = NULL;
                }

                /**************************************************/
                /* Loop on Partitions                             */
                /**************************************************/
                for(j = 0; j < u1_num_part; j++, ps_part++)
                {
                    u1_sub_mb_num = ps_part->u1_sub_mb_num;
                    ps_dec->u1_sub_mb_num = u1_sub_mb_num;

                    if(PART_NOT_DIRECT != ps_part->u1_is_direct)
                    {
                        /* Mb Skip Mode */
                        /* Setting the default and other members of MvPred Structure */
                        s_mvPred.i2_mv[2] = -1;
                        s_mvPred.i2_mv[3] = -1;
                        s_mvPred.i1_ref_frame[0] = 0;
                        i1_ref_idx = (u1_scale_ref && u4_bot_mb) ? MAX_REF_BUFS : 0;
                        ps_ref_frame = pps_ref_frame[i1_ref_idx];
                        s_mvPred.u1_col_ref_pic_idx = ps_ref_frame->u1_mv_buf_id;
                        s_mvPred.u1_pic_type = ps_ref_frame->u1_pic_type;
                        pu4_wt_offst = (UWORD32 *) &ps_dec->pu4_wt_ofsts[0];

                        ps_dec->pf_mvpred(ps_dec, ps_cur_mb_info, ps_mv_nmb_start, ps_mv_ntop_start,
                                          &s_mvPred, 0, 4, 0, 1, MB_SKIP);

                        {
                            pred_info_pkd_t *ps_pred_pkd;
                            ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                            ih264d_fill_pred_info(s_mvPred.i2_mv, 4, 4, 0, PRED_L0, ps_pred_pkd,
                                                  ps_ref_frame->u1_pic_buf_id,
                                                  (i1_ref_idx >> u1_scale_ref), pu4_wt_offst,
                                                  ps_ref_frame->u1_pic_type);

                            ps_dec->u4_pred_info_pkd_idx++;
                            ps_cur_mb_info->u1_num_pred_parts++;
                        }

                        /* Storing colocated zero information */
                        u1_colz = ((ABS(s_mvPred.i2_mv[0]) <= 1) && (ABS(s_mvPred.i2_mv[1]) <= 1)) +
                                  (u1_field << 1);

                        if(ps_mv_nmb_start)
                        {
                            ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0, u1_colz, 4,
                                               4);
                        }
                        else
                        {
                            return NOT_OK;
                        }
                    }
                    else
                    {
                        u1_sub_mb_x = u1_sub_mb_num & 0x03;
                        u1_sub_mb_y = u1_sub_mb_num >> 2;
                        u1_blk_no = (u1_num_part < 4)
                                        ? j
                                        : (((u1_sub_mb_y >> 1) << 1) + (u1_sub_mb_x >> 1));

                        ps_mv_ntop = ps_mv_ntop_start + u1_sub_mb_x;
                        ps_mv_nmb = ps_mv_nmb_start + u1_sub_mb_num;

                        u1_wd = ps_part->u1_partwidth;
                        u1_ht = ps_part->u1_partheight;

                        /* Populate the colpic info and reference frames */
                        i1_ref_idx = pi1_ref_idx[u1_blk_no];
                        /********************************************************************/
                        /* If reference index is inferred from the base layer and it is     */
                        /* exceeding the number of active reference in the current layer.   */
                        /* Then reference index is clipped to the max in the current layer  */
                        /********************************************************************/
                        if(ps_svc_cur_mb_info->u1_base_mode_flag == 1)
                        {
                            if(i1_ref_idx > (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0] - 1))
                            {
                                i1_ref_idx = ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0] - 1;
                            }
                        }
                        s_mvPred.i1_ref_frame[0] = i1_ref_idx;

                        if((1 != ps_svc_cur_mb_info->u1_base_mode_flag) &&
                           (0 == (u1_motion_pred_flag_l0 & (1 << u1_blk_no))))
                        {
                            /********************************************************/
                            /* Predict Mv                                           */
                            /* Add Mv Residuals and store back                      */
                            /********************************************************/
                            ps_dec->pf_mvpred(ps_dec, ps_cur_mb_info, ps_mv_nmb, ps_mv_ntop,
                                              &s_mvPred, u1_sub_mb_num, u1_wd, 0, 1,
                                              ps_cur_mb_info->u1_mb_mc_mode);

                            i2_mv_x = ps_mv_nmb->i2_mv[0];
                            i2_mv_y = ps_mv_nmb->i2_mv[1];
                            i2_mv_x += s_mvPred.i2_mv[0];
                            i2_mv_y += s_mvPred.i2_mv[1];
                            s_mvPred.i2_mv[0] = i2_mv_x;
                            s_mvPred.i2_mv[1] = i2_mv_y;
                        }
                        else
                        {
                            isvcd_retrive_infer_mode_mv(ps_svc_lyr_dec, &s_mvPred, 0,
                                                        u1_sub_mb_num);

                            if(0 != (u1_motion_pred_flag_l0 & (1 << u1_blk_no)))
                            {
                                i2_mv_x = ps_mv_nmb->i2_mv[0];
                                i2_mv_y = ps_mv_nmb->i2_mv[1];
                                i2_mv_x += s_mvPred.i2_mv[0];
                                i2_mv_y += s_mvPred.i2_mv[1];
                                s_mvPred.i2_mv[0] = i2_mv_x;
                                s_mvPred.i2_mv[1] = i2_mv_y;
                            }
                            i2_mv_x = s_mvPred.i2_mv[0];
                            i2_mv_y = s_mvPred.i2_mv[1];
                        }
                        /********************************************************/
                        /* Transfer setup call                                  */
                        /* convert RefIdx if it is MbAff                        */
                        /* Pass Weight Offset and refFrame                      */
                        /********************************************************/
                        i1_ref_idx1 = i1_ref_idx >> u1_scale_ref;
                        if(u1_scale_ref && ((i1_ref_idx & 0x01) != u4_bot_mb))
                            i1_ref_idx1 += MAX_REF_BUFS;
                        if(-1 == i1_ref_idx1) return NOT_OK;
                        ps_ref_frame = pps_ref_frame[i1_ref_idx1];
                        pu4_wt_offst = ppu4_wt_ofst[u1_blk_no];

                        {
                            pred_info_pkd_t *ps_pred_pkd;
                            ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                            ih264d_fill_pred_info(s_mvPred.i2_mv, u1_wd, u1_ht, u1_sub_mb_num,
                                                  PRED_L0, ps_pred_pkd, ps_ref_frame->u1_pic_buf_id,
                                                  (i1_ref_idx >> u1_scale_ref), pu4_wt_offst,
                                                  ps_ref_frame->u1_pic_type);

                            ps_dec->u4_pred_info_pkd_idx++;
                            ps_cur_mb_info->u1_num_pred_parts++;
                        }

                        /* Fill colocated info in MvPred structure */
                        s_mvPred.u1_col_ref_pic_idx = ps_ref_frame->u1_mv_buf_id;
                        s_mvPred.u1_pic_type = ps_ref_frame->u1_pic_type;

                        /* Calculating colocated zero information */
                        u1_colz = (u1_field << 1) |
                                  ((i1_ref_idx == 0) && (ABS(i2_mv_x) <= 1) && (ABS(i2_mv_y) <= 1));
                        u1_colz |= ps_mb_part_info->u1_col_info[u1_blk_no];

                        /* Replicate the motion vectors and colzero u4_flag  */
                        /* for all sub-partitions                         */

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
                ps_cur_deblk_mb->u1_mb_type |= D_INTRA_IBL;
                if((ps_svc_lyr_dec->u1_layer_identifier != TARGET_LAYER) &&
                   (DBLK_ENABLED == ps_dec->ps_cur_slice->u1_disable_dblk_filter_idc))
                {
                    ps_cur_deblk_mb->u1_deblocking_mode = MB_ENABLE_FILTERING;
                }
                /* to take care of 16 parttitions increment for base mode flag case*/
                if(1 != ps_svc_cur_mb_info->u1_base_mode_flag)
                {
                    return NOT_OK;
                }
                {
                    ps_part += (MAX_NUM_MB_PART);
                }
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
 * \if Function name : isvcd_update_intra_mb_inter_layer_info \endif
 *
 * \brief : IT
 *    This function decodes an Inter MB fornfor ot target base layers
 *    Only for Progressive : saves residual for upper enhancement layers
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
void isvcd_update_intra_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                            dec_mb_info_t *ps_cur_mb_info)
{
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
        (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));

    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode = SVC_INTRA_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = 0;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = 0;
}

/*!
**************************************************************************
* \if Function name : isvcd_update_ipcm_mb_inter_layer_info \endif
*
* \brief : IT
*    This function decodes an IPM MB fornfor ot target base layers
*    Only for Progressive : saves residual for upper enhancement layers
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
void isvcd_update_ipcm_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info)
{
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode = SVC_IPCM_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = 0;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = 0;
}

/*!
**************************************************************************
* \if Function name : isvcd_update_ibl_mb_inter_layer_info \endif
*
* \brief : IT
*    This function decodes an IBL MB fornfor ot target base layers
*    Only for Progressive : saves residual for upper enhancement layers
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
void isvcd_update_ibl_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_mb_info_t *ps_cur_mb_info)
{
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode = SVC_IBL_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = 0;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = 0;
}
/*!
**************************************************************************
* \if Function name : isvcd_update_inter_mb_inter_layer_info \endif
*
* \brief : IT
*    This function decodes an IBL MB fornfor ot target base layers
*    Only for Progressive : saves residual for upper enhancement layers
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
void isvcd_update_inter_mb_inter_layer_info(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                            dec_mb_info_t *ps_cur_mb_info, UWORD8 u1_inference_mode)
{
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
        (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode =
        u1_inference_mode ? SVC_IBL_MB : SVC_INTER_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = ps_cur_mb_info->u2_luma_csbp;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz =
        (UWORD8) ps_cur_mb_info->u2_chroma_csbp;
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 1))
    {
        /* Four bits for Cb in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0x0F;
    }
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 2))
    {
        /* Four bits for Cr in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0xF0;
    }
}
/*!
 **************************************************************************
 * \if Function name : isvcd_process_inter_mb_no_rsd_pred_non_target \endif
 *
 * \brief : IT
 *    This function decodes an Inter MB fornfor ot target base layers
 *    Only for Progressive : saves residual for upper enhancement layers
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_process_inter_mb_no_rsd_pred_non_target(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                     dec_mb_info_t *ps_cur_mb_info,
                                                     UWORD8 u1_inference_mode)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD16 u2_luma_stride, u2_chroma_stride;
    WORD16 *pi2_y_coeff, *pi2_luma_res_ptr, *pi2_chroma_res_ptr;
    UWORD32 u4_luma_dc_only_csbp = 0;
    UWORD32 u4_luma_dc_only_cbp = 0;

    if(0 != ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
    {
        return NOT_OK;
    }
    u2_luma_stride = ps_svc_lyr_dec->u2_residual_resample_luma_stride;
    pi2_luma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_luma_frm_start +
                       (ps_cur_mb_info->u2_mbx << 4) +
                       ((ps_cur_mb_info->u2_mby << 4) * u2_luma_stride);

    u2_chroma_stride = ps_svc_lyr_dec->u2_residual_resample_chroma_stride;
    pi2_chroma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_chroma_frm_start +
                         (ps_cur_mb_info->u2_mbx << 4) +
                         ((ps_cur_mb_info->u2_mby << 3) * u2_chroma_stride);

    if(!ps_cur_mb_info->u1_tran_form8x8)
    {
        u4_luma_dc_only_csbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
    }
    else
    {
        if(!ps_dec->ps_cur_pps->u1_entropy_coding_mode)
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
        }
        else
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff8x8_mb(ps_dec, ps_cur_mb_info);
        }
    }

    pi2_y_coeff = ps_dec->pi2_coeff_data;
    /* Inverse Transform and Reconstruction */
    if(ps_cur_mb_info->u1_cbp & 0x0f)
    {
        if(!ps_cur_mb_info->u1_tran_form8x8)
        {
            UWORD32 i;
            WORD16 ai2_tmp[16] = {0};
            for(i = 0; i < 16; i++)
            {
                if(CHECKBIT(ps_cur_mb_info->u2_luma_csbp, i))
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x3) * BLK_SIZE) +
                                      (i >> 2) * (u2_luma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_csbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4_dc(
                                pi2_level, pi2_out, u2_luma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4(
                                pi2_level, pi2_out, u2_luma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                    }
                }
            }
        }
        else
        {
            WORD16 *pi2_scale_matrix_ptr;
            WORD32 i;

            pi2_scale_matrix_ptr = ps_dec->s_high_profile.i2_scalinglist8x8[1];

            for(i = 0; i < 4; i++)
            {
                WORD16 ai2_tmp[64] = {0};
                WORD16 *pi16_levelBlock =
                    pi2_y_coeff + (i << 6); /* move to the next 8x8 adding 64 */

                WORD16 *pi2_out =
                    pi2_luma_res_ptr + ((i & 0x1) * BLK8x8SIZE) + (i >> 1) * (u2_luma_stride << 3);
                if(CHECKBIT(ps_cur_mb_info->u1_cbp, i))
                {
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_cbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8_dc(
                                pi16_levelBlock, pi2_out, u2_luma_stride,
                                gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8(
                                pi16_levelBlock, pi2_out, u2_luma_stride,
                                gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                ai2_tmp, 0, NULL);
                        }
                    }
                }
            }
        }
    }

    /* Decode Chroma Block */
    ih264d_unpack_chroma_coeff4x4_mb(ps_dec, ps_cur_mb_info);
    /*--------------------------------------------------------------------*/
    /* Chroma Blocks decoding                                             */
    /*--------------------------------------------------------------------*/
    {
        UWORD8 u1_chroma_cbp = (UWORD8) (ps_cur_mb_info->u1_cbp >> 4);

        if(u1_chroma_cbp != CBPC_ALLZERO)
        {
            UWORD32 u4_scale_u = ps_cur_mb_info->u1_qpc_div6;
            UWORD32 u4_scale_v = ps_cur_mb_info->u1_qpcr_div6;
            UWORD16 u2_chroma_csbp = ps_cur_mb_info->u2_chroma_csbp;

            pi2_y_coeff = ps_dec->pi2_coeff_data;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_chroma_res_ptr +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4(
                                pi2_level, pi2_out, u2_chroma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4], u4_scale_u,
                                ai2_tmp, pi2_level);
                        }
                        else if(pi2_level[0] != 0)
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4_dc(
                                pi2_level, pi2_out, u2_chroma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4], u4_scale_u,
                                ai2_tmp, pi2_level);
                        }
                    }
                }
            }

            pi2_y_coeff += MB_CHROM_SIZE;
            u2_chroma_csbp >>= 4;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_chroma_res_ptr + 1 +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4(
                                pi2_level, pi2_out, u2_chroma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5], u4_scale_v,
                                ai2_tmp, pi2_level);
                        }
                        else if(pi2_level[0] != 0)
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4_dc(
                                pi2_level, pi2_out, u2_chroma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5], u4_scale_v,
                                ai2_tmp, pi2_level);
                        }
                    }
                }
            }
        }
    }

    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
        (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode =
        u1_inference_mode ? SVC_IBL_MB : SVC_INTER_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = ps_cur_mb_info->u2_luma_csbp;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz =
        (UWORD8) ps_cur_mb_info->u2_chroma_csbp;
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 1))
    {
        /* Four bits for Cb in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0x0F;
    }
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 2))
    {
        /* Four bits for Cr in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0xF0;
    }
    return OK;
}

/*!
 **************************************************************************
 * \if Function name : isvcd_process_inter_mb_rsd_pred_non_target \endif
 *
 * \brief : IT + Residual :
 *    This function decodes an Inter MB for non target layers
 *    Only for Progressive : saves residual + IT for upper enhancement layers
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_process_inter_mb_rsd_pred_non_target(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_mb_info_t *ps_cur_mb_info,
                                                  UWORD8 u1_inference_mode,
                                                  UWORD16 *pu2_res_luma_csbp)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD16 u2_luma_stride, u2_chroma_stride;
    WORD16 *pi2_y_coeff, *pi2_luma_res_ptr, *pi2_chroma_res_ptr;
    UWORD32 u4_luma_dc_only_csbp = 0;
    UWORD32 u4_luma_dc_only_cbp = 0;
    UWORD16 u2_res_luma_csbp = 0;
    UWORD16 u2_res_chroma_csbp = 0, u2_res_chroma_nnz = 0;
    WORD32 ret;

    if(0 != ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
    {
        return NOT_OK;
    }
    u2_luma_stride = ps_svc_lyr_dec->u2_residual_resample_luma_stride;
    pi2_luma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_luma_frm_start +
                       (ps_cur_mb_info->u2_mbx << 4) +
                       ((ps_cur_mb_info->u2_mby << 4) * u2_luma_stride);

    u2_chroma_stride = ps_svc_lyr_dec->u2_residual_resample_chroma_stride;
    pi2_chroma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_chroma_frm_start +
                         (ps_cur_mb_info->u2_mbx << 4) +
                         ((ps_cur_mb_info->u2_mby << 3) * u2_chroma_stride);

    // residual prediction SVC
    ret = isvcd_process_residual_resample_mb(ps_svc_lyr_dec, ps_cur_mb_info);
    if(ret != OK)
    {
        return ret;
    }

    if(!ps_cur_mb_info->u1_tran_form8x8)
    {
        u4_luma_dc_only_csbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
    }
    else
    {
        if(!ps_dec->ps_cur_pps->u1_entropy_coding_mode)
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
        }
        else
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff8x8_mb(ps_dec, ps_cur_mb_info);
        }
    }

    *pu2_res_luma_csbp = 0;
    pi2_y_coeff = ps_dec->pi2_coeff_data;
    /* Inverse Transform and Reconstruction */
    if(ps_cur_mb_info->u1_cbp & 0x0f)
    {
        if(!ps_cur_mb_info->u1_tran_form8x8)
        {
            UWORD32 i;
            WORD16 ai2_tmp[16] = {0};
            for(i = 0; i < 16; i++)
            {
                if(CHECKBIT(ps_cur_mb_info->u2_luma_csbp, i))
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x3) * BLK_SIZE) +
                                      (i >> 2) * (u2_luma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_csbp, i))
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4_dc(
                                    pi2_level, pi2_out, pi2_out, u2_luma_stride, u2_luma_stride,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                    ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            u2_res_luma_csbp = ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4(
                                pi2_level, pi2_out, pi2_out, u2_luma_stride, u2_luma_stride,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                    }
                }
                else
                {
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x3) * BLK_SIZE) +
                                      (i >> 2) * (u2_luma_stride << 2);

                    u2_res_luma_csbp =
                        ps_svc_lyr_dec->pf_residual_luma_4x4(pi2_out, u2_luma_stride);
                }
                *pu2_res_luma_csbp |= (u2_res_luma_csbp << i);
            }
        }
        else
        {
            WORD16 *pi2_scale_matrix_ptr;
            WORD32 i;

            pi2_scale_matrix_ptr = ps_dec->s_high_profile.i2_scalinglist8x8[1];

            for(i = 0; i < 4; i++)
            {
                WORD16 ai2_tmp[64] = {0};
                WORD16 *pi16_levelBlock =
                    pi2_y_coeff + (i << 6); /* move to the next 8x8 adding 64 */

                WORD16 *pi2_out =
                    pi2_luma_res_ptr + ((i & 0x1) * BLK8x8SIZE) + (i >> 1) * (u2_luma_stride << 3);
                if(CHECKBIT(ps_cur_mb_info->u1_cbp, i))
                {
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_cbp, i))
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8_dc(
                                    pi16_levelBlock, pi2_out, pi2_out, u2_luma_stride,
                                    u2_luma_stride,
                                    gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                    ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            u2_res_luma_csbp = ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8(
                                pi16_levelBlock, pi2_out, pi2_out, u2_luma_stride, u2_luma_stride,
                                gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                ai2_tmp, 0, NULL);
                        }
                    }
                }
                else
                {
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x1) * BLK8x8SIZE) +
                                      (i >> 1) * (u2_luma_stride << 3);

                    u2_res_luma_csbp =
                        ps_svc_lyr_dec->pf_residual_luma_8x8(pi2_out, u2_luma_stride);
                }
                *pu2_res_luma_csbp |= (u2_res_luma_csbp << (((i >> 1) << 3) + ((i & 0x01) << 1)));
            }
        }
    }
    else
    {
        WORD16 *pi2_out = pi2_luma_res_ptr;

        *pu2_res_luma_csbp = ps_svc_lyr_dec->pf_residual_luma_16x16(pi2_out, u2_luma_stride);
    }

    /* Decode Chroma Block */
    ih264d_unpack_chroma_coeff4x4_mb(ps_dec, ps_cur_mb_info);
    /*--------------------------------------------------------------------*/
    /* Chroma Blocks decoding                                             */
    /*--------------------------------------------------------------------*/
    {
        UWORD8 u1_chroma_cbp = (UWORD8) (ps_cur_mb_info->u1_cbp >> 4);

        u2_res_chroma_nnz =
            ps_svc_lyr_dec->pf_residual_chroma_cb_cr_8x8(pi2_chroma_res_ptr, u2_chroma_stride);

        if(u1_chroma_cbp != CBPC_ALLZERO)
        {
            UWORD32 u4_scale_u = ps_cur_mb_info->u1_qpc_div6;
            UWORD32 u4_scale_v = ps_cur_mb_info->u1_qpcr_div6;
            UWORD16 u2_chroma_csbp = ps_cur_mb_info->u2_chroma_csbp;

            pi2_y_coeff = ps_dec->pi2_coeff_data;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_chroma_res_ptr +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            u2_res_chroma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4(
                                    pi2_level, pi2_out, pi2_out, u2_chroma_stride, u2_chroma_stride,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4],
                                    u4_scale_u, ai2_tmp, pi2_level);
                            u2_res_chroma_nnz &= (0XFF ^ (1 << i));
                            u2_res_chroma_nnz |= (u2_res_chroma_csbp << i);
                        }
                        else if(pi2_level[0] != 0)
                        {
                            u2_res_chroma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4_dc(
                                    pi2_level, pi2_out, pi2_out, u2_chroma_stride, u2_chroma_stride,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4],
                                    u4_scale_u, ai2_tmp, pi2_level);
                            u2_res_chroma_nnz &= (0XFF ^ (1 << i));
                            u2_res_chroma_nnz |= (u2_res_chroma_csbp << i);
                        }
                    }
                }
            }

            pi2_y_coeff += MB_CHROM_SIZE;
            u2_chroma_csbp >>= 4;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    WORD16 *pi2_out = pi2_chroma_res_ptr + 1 +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            u2_res_chroma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4(
                                    pi2_level, pi2_out, pi2_out, u2_chroma_stride, u2_chroma_stride,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5],
                                    u4_scale_v, ai2_tmp, pi2_level);
                            u2_res_chroma_nnz &= (0XFF ^ (1 << (i + 4)));
                            u2_res_chroma_nnz |= (u2_res_chroma_csbp << (i + 4));
                        }
                        else if(pi2_level[0] != 0)
                        {
                            u2_res_chroma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4_dc(
                                    pi2_level, pi2_out, pi2_out, u2_chroma_stride, u2_chroma_stride,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5],
                                    u4_scale_v, ai2_tmp, pi2_level);
                            u2_res_chroma_nnz &= (0XFF ^ (1 << (i + 4)));
                            u2_res_chroma_nnz |= (u2_res_chroma_csbp << (i + 4));
                        }
                    }
                }
            }
        }
    }

    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
        (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode =
        u1_inference_mode ? SVC_IBL_MB : SVC_INTER_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;

    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = *pu2_res_luma_csbp;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = (UWORD8) u2_res_chroma_nnz;
    return OK;
}
/*!
**************************************************************************
* \if Function name : isvcd_process_ii_mb \endif
*
* \brief
*    This function decodes an intra inter mb
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_process_ii_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                           dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num)
{
    res_prms_t *ps_res_prms;
    WORD32 i4_status;
    UWORD8 u1_ii_mb_mode = 0;
    mb_coord_t s_mb_coord = {0};
    mem_element_t s_ref_mb_mode = {0};
    svc_dec_lyr_struct_t *ps_svc_dec_ref_layer;

    ps_svc_dec_ref_layer = ps_svc_lyr_dec->ps_dec_svc_ref_layer;
    ps_res_prms = &ps_svc_lyr_dec->s_res_prms;
    s_mb_coord.u2_mb_x = ps_cur_mb_info->u2_mbx;
    s_mb_coord.u2_mb_y = ps_cur_mb_info->u2_mby;

    /* Restricted resolution change has significance only */
    /* at resolution change layer                         */
    if(SVCD_FALSE == ps_res_prms->u1_rstrct_res_change_flag)
    {
        s_ref_mb_mode.pv_buffer = ps_svc_dec_ref_layer->ps_inter_lyr_mb_prms_frm_start;
        s_ref_mb_mode.i4_element_size = sizeof(inter_lyr_mb_prms_t);
        s_ref_mb_mode.i4_num_element_stride = ps_svc_dec_ref_layer->u2_inter_lyr_mb_prms_stride;

        i4_status = isvcd_ii_pred_compute_flags_mb(ps_svc_lyr_dec->pv_ii_pred_ctxt, &s_ref_mb_mode,
                                                   &s_mb_coord, ps_cur_mb_info, ps_svc_cur_mb_info,
                                                   &u1_ii_mb_mode);

        if(OK != i4_status)
        {
            return i4_status;
        }
    }

    if(SVC_INTRA_INTER_MB == u1_ii_mb_mode)
    {
        i4_status = isvcd_process_ibl_mb(ps_svc_lyr_dec, ps_cur_mb_info, u1_mb_num, 1);
        if(OK != i4_status)
        {
            return i4_status;
        }
        isvcd_ii_pred_mb(ps_svc_lyr_dec, ps_cur_mb_info);
    }
    return OK;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_decode_recon_tfr_nmb_non_base_lyr \endif
 *
 * \brief
 *
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_decode_recon_tfr_nmb_non_base_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                               UWORD32 u4_mb_idx, UWORD32 u4_num_mbs,
                                               UWORD32 u4_num_mbs_next, UWORD32 u4_tfr_n_mb,
                                               UWORD32 u4_end_of_row)
{
    WORD32 i, j;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD32 u1_end_of_row_next;
    dec_mb_info_t *ps_cur_mb_info;
    dec_svc_mb_info_t *ps_svc_cur_mb_info;
    UWORD16 *pu2_res_luma_csbp;
    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    const UWORD32 u1_slice_type = ps_dec->ps_cur_slice->u1_slice_type;
    const WORD32 u1_skip_th =
        ((u1_slice_type != I_SLICE) ? (ps_dec->u1_B ? B_8x8 : PRED_8x8R0) : -1);
    const UWORD32 u1_ipcm_th = ((u1_slice_type != I_SLICE) ? (ps_dec->u1_B ? 23 : 5) : 0);
    WORD32 ret = OK;

    if(!((0 == ps_svc_lyr_dec->u1_base_res_flag) ||
         ((1 == ps_svc_lyr_dec->u1_base_res_flag) &&
          (1 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_no_inter_layer_pred_flag))))
    {
        return NOT_OK;
    }
    /* N Mb MC Loop */
    for(i = u4_mb_idx; i < u4_num_mbs; i++)
    {
        ps_cur_mb_info = ps_dec->ps_nmb_info + i;
        ps_dec->u4_dma_buf_idx = 0;
        ps_dec->u4_pred_info_idx = 0;

        /*Pointer assignment for Residual NNZ */
        pu2_res_luma_csbp = ps_svc_lyr_dec->pu2_frm_res_luma_csbp + ps_cur_mb_info->u2_mbx;
        pu2_res_luma_csbp += ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride;

        if(ps_cur_mb_info->u1_mb_type <= u1_skip_th)
        {
            {
                WORD32 pred_cnt = 0;
                pred_info_pkd_t *ps_pred_pkd;
                UWORD32 u4_pred_info_pkd_idx;

                u4_pred_info_pkd_idx = ps_cur_mb_info->u4_pred_info_pkd_idx;

                while(pred_cnt < ps_cur_mb_info->u1_num_pred_parts)
                {
                    ps_pred_pkd = ps_dec->ps_pred_pkd + u4_pred_info_pkd_idx;

                    ps_dec->p_form_mb_part_info(ps_pred_pkd, ps_dec, ps_cur_mb_info->u2_mbx,
                                                ps_cur_mb_info->u2_mby, (i >> u1_mbaff),
                                                ps_cur_mb_info);
                    u4_pred_info_pkd_idx++;
                    pred_cnt++;
                }
            }
            if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
            {
                ps_dec->p_motion_compensate(ps_dec, ps_cur_mb_info);
            }
        }
        else if(ps_cur_mb_info->u1_mb_type == MB_SKIP)
        {
            {
                WORD32 pred_cnt = 0;
                pred_info_pkd_t *ps_pred_pkd;
                UWORD32 u4_pred_info_pkd_idx;

                u4_pred_info_pkd_idx = ps_cur_mb_info->u4_pred_info_pkd_idx;

                while(pred_cnt < ps_cur_mb_info->u1_num_pred_parts)
                {
                    ps_pred_pkd = ps_dec->ps_pred_pkd + u4_pred_info_pkd_idx;

                    ps_dec->p_form_mb_part_info(ps_pred_pkd, ps_dec, ps_cur_mb_info->u2_mbx,
                                                ps_cur_mb_info->u2_mby, (i >> u1_mbaff),
                                                ps_cur_mb_info);

                    u4_pred_info_pkd_idx++;
                    pred_cnt++;
                }
            }
            if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
            {
                /* Decode MB skip */
                ps_dec->p_motion_compensate(ps_dec, ps_cur_mb_info);
            }

            *pu2_res_luma_csbp = 0;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
                ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
                (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode = SVC_INTER_MB;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size =
                ps_cur_mb_info->u1_tran_form8x8;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = 0;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = 0;
        }
    }

    /* N Mb IQ IT RECON  Loop */
    for(j = u4_mb_idx; j < i; j++)
    {
        ps_cur_mb_info = ps_dec->ps_nmb_info + j;
        ps_svc_cur_mb_info = ps_svc_lyr_dec->ps_svc_nmb_info + j;

        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
            (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));

        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_slice_id = (WORD8) ps_dec->u2_cur_slice_num;

        /*Pointer assignment for Residual NNZ */
        pu2_res_luma_csbp = ps_svc_lyr_dec->pu2_frm_res_luma_csbp + ps_cur_mb_info->u2_mbx;
        pu2_res_luma_csbp += ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride;

        if(ps_cur_mb_info->u1_mb_type <= u1_skip_th)
        {
            if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
            {
                /* inter intra pred generation */
                if(SVCD_FALSE == ps_svc_lyr_dec->u1_dyadic_flag)
                {
                    ret =
                        isvcd_process_ii_mb(ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info, j);
                    if(ret != OK) return ret;
                }
                if(0 == ps_svc_cur_mb_info->u1_residual_prediction_flag)
                {
                    // IT + Recon
                    ih264d_process_inter_mb(ps_dec, ps_cur_mb_info, j);
                    isvcd_update_inter_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info, 0);
                    *pu2_res_luma_csbp = ps_cur_mb_info->u2_luma_csbp;
                }
                else
                {
                    // IT + Residual + Recon
                    ret = isvcd_process_inter_mb_rsd_pred_target_lyr(ps_svc_lyr_dec, ps_cur_mb_info,
                                                                     j, 0, pu2_res_luma_csbp);
                    if(ret != OK) return ret;
                }
            }
            else if(ps_svc_lyr_dec->u1_layer_identifier == MEDIAL_ENHANCEMENT_LAYER)
            {
                if(0 == ps_svc_cur_mb_info->u1_residual_prediction_flag)
                {
                    // IT : to be consumed by Target
                    ret = isvcd_process_inter_mb_no_rsd_pred_non_target(ps_svc_lyr_dec,
                                                                        ps_cur_mb_info, 0);
                    *pu2_res_luma_csbp = ps_cur_mb_info->u2_luma_csbp;
                    if(ret != OK) return ret;
                }
                else
                {
                    // IT + Residual : to be consumed by target
                    ret = isvcd_process_inter_mb_rsd_pred_non_target(ps_svc_lyr_dec, ps_cur_mb_info,
                                                                     0, pu2_res_luma_csbp);
                    if(ret != OK) return ret;
                }
            }
            else
            {
                return NOT_OK;
            }
        }
        else if((ps_cur_mb_info->u1_mb_type != MB_SKIP) && (ps_cur_mb_info->u1_mb_type != MB_INFER))
        {
            if((u1_ipcm_th + 25) != ps_cur_mb_info->u1_mb_type)
            {
                ps_cur_mb_info->u1_mb_type -= (u1_skip_th + 1);
                ih264d_process_intra_mb(ps_dec, ps_cur_mb_info, j);
                isvcd_update_intra_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info);
            }
            else
            {
                isvcd_update_ipcm_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info);
            }
            *pu2_res_luma_csbp = 0;
        }
        else if(ps_cur_mb_info->u1_mb_type == MB_INFER)
        {
            /* inter layer intra prediction : intra upsample, IQ, IT ,deblock */
            /* Intra resample for IBL mode */
            ret = isvcd_process_ibl_mb(ps_svc_lyr_dec, ps_cur_mb_info, j, 0);
            if(ret != OK) return ret;
            ih264d_process_inter_mb(ps_dec, ps_cur_mb_info, j);
            isvcd_update_inter_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info, 1);
            *pu2_res_luma_csbp = ps_cur_mb_info->u2_luma_csbp;

            ps_dec->pi1_left_pred_mode[0] = DC;
            ps_dec->pi1_left_pred_mode[1] = DC;
            ps_dec->pi1_left_pred_mode[2] = DC;
            ps_dec->pi1_left_pred_mode[3] = DC;

            ps_cur_mb_info->ps_curmb->pi1_intrapredmodes[0] = DC;
            ps_cur_mb_info->ps_curmb->pi1_intrapredmodes[1] = DC;
            ps_cur_mb_info->ps_curmb->pi1_intrapredmodes[2] = DC;
            ps_cur_mb_info->ps_curmb->pi1_intrapredmodes[3] = DC;

            isvcd_update_ibl_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info);
        }

        if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
        {
            if(ps_dec->u4_num_cores < 3)
            {
                if(ps_dec->u4_app_disable_deblk_frm == 0)
                    ps_svc_lyr_dec->pf_svc_compute_bs(ps_svc_lyr_dec, ps_cur_mb_info,
                                                      (UWORD16) (j >> u1_mbaff));
            }
        }
        else if(ps_svc_lyr_dec->u1_layer_identifier == MEDIAL_ENHANCEMENT_LAYER)
        {
            if(ps_dec->u4_num_cores < 3)
            {
                if(ps_dec->u4_app_disable_deblk_frm == 0)
                    ps_svc_lyr_dec->pf_svc_compute_bs(ps_svc_lyr_dec, ps_cur_mb_info,
                                                      (UWORD16) (j >> u1_mbaff));
            }
        }

        if(ps_dec->u4_use_intrapred_line_copy)
        {
            ih264d_copy_intra_pred_line(ps_dec, ps_cur_mb_info, j);
        }
    }

    /*MB deblocking*/
    if(ps_dec->u4_nmb_deblk == 1)
    {
        UWORD32 u4_wd_y, u4_wd_uv;
        tfr_ctxt_t *ps_tfr_cxt = &(ps_dec->s_tran_addrecon);
        UWORD8 u1_field_pic_flag = ps_dec->ps_cur_slice->u1_field_pic_flag;
        const WORD32 i4_cb_qp_idx_ofst = ps_dec->ps_cur_pps->i1_chroma_qp_index_offset;
        const WORD32 i4_cr_qp_idx_ofst = ps_dec->ps_cur_pps->i1_second_chroma_qp_index_offset;

        u4_wd_y = ps_dec->u2_frm_wd_y << u1_field_pic_flag;
        u4_wd_uv = ps_dec->u2_frm_wd_uv << u1_field_pic_flag;

        ps_cur_mb_info = ps_dec->ps_nmb_info + u4_mb_idx;

        ps_dec->u4_deblk_mb_x = ps_cur_mb_info->u2_mbx;
        ps_dec->u4_deblk_mb_y = ps_cur_mb_info->u2_mby;

        for(j = u4_mb_idx; j < i; j++)
        {
            if(ps_dec->u4_cur_deblk_mb_num > ps_dec->ps_cur_sps->u4_max_mb_addr)
            {
                return NOT_OK;
            }
            ih264d_deblock_mb_nonmbaff(ps_dec, ps_tfr_cxt, i4_cb_qp_idx_ofst, i4_cr_qp_idx_ofst,
                                       u4_wd_y, u4_wd_uv);
        }
    }

    if(u4_tfr_n_mb)
    {
        /****************************************************************/
        /* Check for End Of Row in Next iteration                       */
        /****************************************************************/
        u1_end_of_row_next =
            u4_num_mbs_next && (u4_num_mbs_next <= (ps_dec->u4_recon_mb_grp >> u1_mbaff));

        /****************************************************************/
        /* Transfer the Following things                                */
        /* N-Mb DeblkParams Data    ( To Ext DeblkParams Buffer )       */
        /* N-Mb Recon Data          ( To Ext Frame Buffer )             */
        /* N-Mb Intrapredline Data  ( Updated Internally)               */
        /* N-Mb MV Data             ( To Ext MV Buffer )                */
        /* N-Mb MVTop/TopRight Data ( To Int MV Top Scratch Buffers)    */
        /****************************************************************/
        ih264d_transfer_mb_group_data(ps_dec, u4_num_mbs, u4_end_of_row, u1_end_of_row_next);
        ps_dec->u4_num_mbs_prev_nmb = u4_num_mbs;
        ps_dec->u4_pred_info_idx = 0;
        ps_dec->u4_dma_buf_idx = 0;
    }
    return OK;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_decode_recon_tfr_nmb_base_lyr \endif
 *
 * \brief
 *
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_decode_recon_tfr_nmb_base_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD32 u1_mb_idx,
                                           UWORD32 u4_num_mbs, UWORD32 u4_num_mbs_next,
                                           UWORD32 u4_tfr_n_mb, UWORD32 u4_end_of_row)
{
    WORD32 j;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD32 u1_end_of_row_next;
    dec_mb_info_t *ps_cur_mb_info;
    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    const UWORD32 u1_slice_type = ps_dec->ps_cur_slice->u1_slice_type;
    const WORD32 u1_skip_th =
        ((u1_slice_type != I_SLICE) ? (ps_dec->u1_B ? B_8x8 : PRED_8x8R0) : -1);
    const UWORD32 u1_ipcm_th = ((u1_slice_type != I_SLICE) ? (ps_dec->u1_B ? 23 : 5) : 0);
    WORD32 ret = OK;

    if(1 != ps_svc_lyr_dec->u1_base_res_flag)
    {
        return NOT_OK;
    }
    if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
    {
        return NOT_OK;
    }

    /* N Mb IQ IT + Residual Store for Inter / + Recon for Intra Loop */
    for(j = u1_mb_idx; j < u4_num_mbs; j++)
    {
        ps_dec->u4_dma_buf_idx = 0;
        ps_dec->u4_pred_info_idx = 0;
        ps_cur_mb_info = ps_dec->ps_nmb_info + j;

        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
            (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));

        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_slice_id = (WORD8) ps_dec->u2_cur_slice_num;

        if(ps_cur_mb_info->u1_mb_type == MB_SKIP)
        {
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode = SVC_INTER_MB;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size =
                ps_cur_mb_info->u1_tran_form8x8;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = 0;
            ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz = 0;
        }
        else if(ps_cur_mb_info->u1_mb_type <= u1_skip_th)
        {
            /* Only IT : Store Residual (WORD16) for Higher Layers : Base layer*/
            ret = isvcd_process_inter_mb_no_rsd_pred_non_target(ps_svc_lyr_dec, ps_cur_mb_info, 0);
            if(ret != OK) return ret;
        }
        else if(ps_cur_mb_info->u1_mb_type != MB_SKIP)
        {
            if((u1_ipcm_th + 25) != ps_cur_mb_info->u1_mb_type)
            {
                ps_cur_mb_info->u1_mb_type -= (u1_skip_th + 1);
                ih264d_process_intra_mb(ps_dec, ps_cur_mb_info, j);
                isvcd_update_intra_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info);
            }
            else
            {
                isvcd_update_ipcm_mb_inter_layer_info(ps_svc_lyr_dec, ps_cur_mb_info);
            }
        }

        if(ps_dec->u4_use_intrapred_line_copy)
        {
            ih264d_copy_intra_pred_line(ps_dec, ps_cur_mb_info, j);
        }
    }

    /*MB deblocking*/
    if(ps_dec->u4_nmb_deblk == 1)
    {
        UWORD32 u4_wd_y, u4_wd_uv;
        tfr_ctxt_t *ps_tfr_cxt = &(ps_dec->s_tran_addrecon);
        UWORD8 u1_field_pic_flag = ps_dec->ps_cur_slice->u1_field_pic_flag;
        const WORD32 i4_cb_qp_idx_ofst = ps_dec->ps_cur_pps->i1_chroma_qp_index_offset;
        const WORD32 i4_cr_qp_idx_ofst = ps_dec->ps_cur_pps->i1_second_chroma_qp_index_offset;

        u4_wd_y = ps_dec->u2_frm_wd_y << u1_field_pic_flag;
        u4_wd_uv = ps_dec->u2_frm_wd_uv << u1_field_pic_flag;

        ps_cur_mb_info = ps_dec->ps_nmb_info + u1_mb_idx;

        ps_dec->u4_deblk_mb_x = ps_cur_mb_info->u2_mbx;
        ps_dec->u4_deblk_mb_y = ps_cur_mb_info->u2_mby;

        for(j = u1_mb_idx; j < u4_num_mbs; j++)
        {
            /* IN SVC base layers only intra MB's Need to be deblocked*/
            deblk_mb_t *ps_top_mb, *ps_left_mb, *ps_cur_mb;
            ps_cur_mb = ps_dec->ps_cur_deblk_mb;
            if(!(ps_cur_mb->u1_deblocking_mode & MB_DISABLE_FILTERING))
            {
                if(ps_dec->u4_deblk_mb_x)
                {
                    ps_left_mb = ps_cur_mb - 1;
                }
                else
                {
                    ps_left_mb = NULL;
                }
                if(ps_dec->u4_deblk_mb_y != 0)
                {
                    ps_top_mb = ps_cur_mb - (ps_dec->u2_frm_wd_in_mbs);
                }
                else
                {
                    ps_top_mb = NULL;
                }

                if(ps_cur_mb->u1_deblocking_mode & MB_DISABLE_LEFT_EDGE) ps_left_mb = NULL;
                if(ps_cur_mb->u1_deblocking_mode & MB_DISABLE_TOP_EDGE) ps_top_mb = NULL;

                /* Top Horizontal Edge*/
                if(NULL != ps_top_mb)
                {
                    if(!(ps_top_mb->u1_mb_type & D_INTRA_MB))
                    {
                        ps_cur_mb->u4_bs_table[0] = 0;
                    }
                }
                else
                {
                    ps_cur_mb->u4_bs_table[0] = 0;
                }

                /* Left Vertical Edge*/
                if(NULL != ps_left_mb)
                {
                    if(!(ps_left_mb->u1_mb_type & D_INTRA_MB))
                    {
                        ps_cur_mb->u4_bs_table[4] = 0;
                    }
                }
                else
                {
                    ps_cur_mb->u4_bs_table[4] = 0;
                }
            }

            if(ps_dec->u4_cur_deblk_mb_num > ps_dec->ps_cur_sps->u4_max_mb_addr)
            {
                return NOT_OK;
            }

            ih264d_deblock_mb_nonmbaff(ps_dec, ps_tfr_cxt, i4_cb_qp_idx_ofst, i4_cr_qp_idx_ofst,
                                       u4_wd_y, u4_wd_uv);
        }
    }

    if(u4_tfr_n_mb)
    {
        /****************************************************************/
        /* Check for End Of Row in Next iteration                       */
        /****************************************************************/
        u1_end_of_row_next =
            u4_num_mbs_next && (u4_num_mbs_next <= (ps_dec->u4_recon_mb_grp >> u1_mbaff));

        /****************************************************************/
        /* Transfer the Following things                                */
        /* N-Mb DeblkParams Data    ( To Ext DeblkParams Buffer )       */
        /* N-Mb Recon Data          ( To Ext Frame Buffer )             */
        /* N-Mb Intrapredline Data  ( Updated Internally)               */
        /* N-Mb MV Data             ( To Ext MV Buffer )                */
        /* N-Mb MVTop/TopRight Data ( To Int MV Top Scratch Buffers)    */
        /****************************************************************/
        ih264d_transfer_mb_group_data(ps_dec, u4_num_mbs, u4_end_of_row, u1_end_of_row_next);
        ps_dec->u4_num_mbs_prev_nmb = u4_num_mbs;
        ps_dec->u4_pred_info_idx = 0;
        ps_dec->u4_dma_buf_idx = 0;
    }
    return OK;
}
/*!
**************************************************************************
* \if Function name : isvcd_process_ibl_mb \endif
*
* \brief
*    This function decodes an ibl mb
*
*
* \return
*
**************************************************************************
*/
WORD32 isvcd_process_ibl_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                            UWORD32 u4_mb_num, UWORD8 u1_inter_intra_mode)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    intra_sampling_ctxt_t *ps_ctxt;
    svc_dec_lyr_struct_t *ps_svc_dec_ref_layer;
    pic_buffer_t *ps_frame_buf = ps_dec->ps_cur_pic;
    pic_buffer_t *ps_frame_buf_ref_layer;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    mem_element_t s_ref_mb_mode = {0};
    mem_element_t s_inp_luma = {0};
    mem_element_t s_inp_chroma = {0};
    mem_element_t s_out_luma = {0};
    mem_element_t s_out_chroma = {0};
    WORD32 i4_ref_x_luma, i4_ref_y_luma, i4_luma_incr = 0;
    WORD32 i4_ref_x_chroma, i4_ref_y_chroma, i4_chroma_incr = 0;
    UWORD32 u4_cur_y_stride, u4_cur_uv_stride;
    UWORD32 u4_ref_y_stride, u4_ref_uv_stride;
    WORD32 i4_ref_luma_instra_sample_correction_offset = 0;
    WORD32 i4_ref_chroma_instra_sample_correction_offset = 0;
    ref_mb_map_t *ps_x_off_len_luma;
    ref_mb_map_t *ps_y_off_len_luma;
    ref_mb_map_t *ps_x_off_len_chroma;
    ref_mb_map_t *ps_y_off_len_chroma;
    mb_coord_t s_mb_coord = {0};
    WORD32 ret = OK;
    UNUSED(u4_mb_num);

    ps_ctxt = (intra_sampling_ctxt_t *) ps_svc_lyr_dec->pv_intra_sample_ctxt;
    ps_svc_dec_ref_layer = ps_svc_lyr_dec->ps_dec_svc_ref_layer;

    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    u4_cur_y_stride = ps_dec->u2_frm_wd_y;
    u4_cur_uv_stride = ps_dec->u2_frm_wd_uv;
    u4_ref_y_stride = ps_svc_dec_ref_layer->s_dec.u2_frm_wd_y;
    u4_ref_uv_stride = ps_svc_dec_ref_layer->s_dec.u2_frm_wd_uv;

    ps_frame_buf_ref_layer = ps_svc_dec_ref_layer->s_dec.ps_cur_pic;
    if(0 == u1_inter_intra_mode)
    {
        s_out_luma.pv_buffer = ps_frame_buf->pu1_buf1 + (ps_cur_mb_info->u2_mbx << 4) +
                               (u4_cur_y_stride * (ps_cur_mb_info->u2_mby << 4));
        s_out_luma.i4_element_size = 1;
        s_out_luma.i4_num_element_stride = u4_cur_y_stride;

        s_out_chroma.pv_buffer = ps_frame_buf->pu1_buf2 +
                                 (ps_cur_mb_info->u2_mbx << 3) * YUV420SP_FACTOR +
                                 (u4_cur_uv_stride * (ps_cur_mb_info->u2_mby << 3));
        s_out_chroma.i4_element_size = 1;
        s_out_chroma.i4_num_element_stride = u4_cur_uv_stride;
    }
    else
    {
        if(SVCD_TRUE == ps_svc_lyr_dec->s_res_prms.u1_dyadic_flag)
        {
            return NOT_OK;
        }

        s_out_luma.pv_buffer = ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma;
        s_out_luma.i4_element_size = 1;
        s_out_luma.i4_num_element_stride = MB_SIZE;

        s_out_chroma.pv_buffer = ps_svc_lyr_dec->pu1_ii_resamp_buffer_chroma;
        s_out_chroma.i4_element_size = 1;
        s_out_chroma.i4_num_element_stride = MB_SIZE;
    }

    /* get the projected locations buffer pointers */
    {
        intra_samp_map_ctxt_t *ps_luma_map, *ps_chroma_map;

        ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
        ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

        ps_x_off_len_luma = ps_luma_map->ps_x_offset_length;
        ps_y_off_len_luma = ps_luma_map->ps_y_offset_length;
        ps_x_off_len_chroma = ps_chroma_map->ps_x_offset_length;
        ps_y_off_len_chroma = ps_chroma_map->ps_y_offset_length;
    }
    i4_ref_x_luma = ps_svc_lyr_dec->ps_intsam_luma_map_horz[ps_cur_mb_info->u2_mbx].i2_offset;
    i4_ref_y_luma = ps_svc_lyr_dec->ps_intsam_luma_map_vert[ps_cur_mb_info->u2_mby].i2_offset;

    i4_luma_incr = ps_x_off_len_luma[ps_cur_mb_info->u2_mbx].i2_offset - i4_ref_x_luma;
    i4_luma_incr +=
        (ps_y_off_len_luma[ps_cur_mb_info->u2_mby].i2_offset - i4_ref_y_luma) * u4_ref_y_stride;

    i4_ref_x_chroma = ps_svc_lyr_dec->ps_intsam_chroma_map_horz[ps_cur_mb_info->u2_mbx].i2_offset;
    i4_ref_y_chroma = ps_svc_lyr_dec->ps_intsam_chroma_map_vert[ps_cur_mb_info->u2_mby].i2_offset;

    i4_chroma_incr = ps_x_off_len_chroma[ps_cur_mb_info->u2_mbx].i2_offset - i4_ref_x_chroma;
    i4_chroma_incr <<= 1;
    i4_chroma_incr += (ps_y_off_len_chroma[ps_cur_mb_info->u2_mby].i2_offset - i4_ref_y_chroma) *
                      u4_ref_uv_stride;
    if(SVCD_FALSE == ps_svc_lyr_dec->s_res_prms.u1_dyadic_flag)
    {
        i4_ref_x_luma = CLIP3(0, (ps_svc_dec_ref_layer->s_dec.u2_frm_wd_y - 1), i4_ref_x_luma);
        i4_ref_y_luma = CLIP3(0, (ps_svc_dec_ref_layer->s_dec.u2_frm_ht_y - 1), i4_ref_y_luma);
    }

    i4_ref_luma_instra_sample_correction_offset =
        i4_ref_x_luma + (i4_ref_y_luma) * (WORD32) u4_ref_y_stride;

    s_inp_luma.pv_buffer = ps_frame_buf_ref_layer->pu1_buf1 + i4_luma_incr +
                           i4_ref_luma_instra_sample_correction_offset;
    s_inp_luma.i4_element_size = 1;
    s_inp_luma.i4_num_element_stride = u4_ref_y_stride;

    if(SVCD_FALSE == ps_svc_lyr_dec->s_res_prms.u1_dyadic_flag)
    {
        i4_ref_x_chroma = CLIP3(0, (ps_svc_dec_ref_layer->s_dec.u2_frm_wd_uv - 1), i4_ref_x_chroma);
        i4_ref_y_chroma = CLIP3(0, (ps_svc_dec_ref_layer->s_dec.u2_frm_ht_uv - 1), i4_ref_y_chroma);
    }
    i4_ref_chroma_instra_sample_correction_offset =
        (i4_ref_x_chroma << 1) + (i4_ref_y_chroma) * (WORD32) u4_ref_uv_stride;

    s_inp_chroma.pv_buffer = ps_frame_buf_ref_layer->pu1_buf2 + i4_chroma_incr +
                             i4_ref_chroma_instra_sample_correction_offset;
    s_inp_chroma.i4_element_size = 1;
    s_inp_chroma.i4_num_element_stride = u4_ref_uv_stride;

    s_ref_mb_mode.pv_buffer = ps_svc_dec_ref_layer->ps_inter_lyr_mb_prms_frm_start;
    s_ref_mb_mode.i4_element_size = sizeof(inter_lyr_mb_prms_t);
    s_ref_mb_mode.i4_num_element_stride = ps_svc_dec_ref_layer->u2_inter_lyr_mb_prms_stride;

    s_mb_coord.u2_mb_x = ps_cur_mb_info->u2_mbx;
    s_mb_coord.u2_mb_y = ps_cur_mb_info->u2_mby;

    if(SVCD_TRUE == ps_svc_lyr_dec->s_res_prms.u1_dyadic_flag)
    {
        ret = isvcd_intra_resamp_mb_dyadic(ps_ctxt, &s_inp_luma, &s_inp_chroma, &s_ref_mb_mode,
                                           &s_out_luma, &s_out_chroma, &s_mb_coord, ps_svc_lyr_dec);
    }
    else
    {
        ret = isvcd_intra_resamp_mb(ps_ctxt, &s_inp_luma, &s_inp_chroma, &s_ref_mb_mode,
                                    &s_out_luma, &s_out_chroma, &s_mb_coord);
    }
    if(OK != ret) return ret;
    return OK;
}
/*!
**************************************************************************
* \if Function name : isvcd_process_residual_resample_mb \endif
*
* \brief
*    This function decodes a residual resample mb
*
*
* \return
*
**************************************************************************
*/
WORD32 isvcd_process_residual_resample_mb(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                          dec_mb_info_t *ps_cur_mb_info)
{
    residual_sampling_ctxt_t *ps_ctxt;
    svc_dec_lyr_struct_t *ps_svc_dec_ref_layer;
    res_lyr_ctxt *ps_lyr_ctxt;
    mem_element_t s_ref_mb_mode = {0};
    mem_element_t s_inp_luma = {0};
    mem_element_t s_inp_chroma = {0};
    mem_element_t s_out_luma = {0};
    mem_element_t s_out_chroma = {0};

    /* projected locations pointer */
    ref_mb_map_t *ps_x_off_len_luma;
    ref_mb_map_t *ps_y_off_len_luma;
    ref_mb_map_t *ps_x_off_len_chroma;
    ref_mb_map_t *ps_y_off_len_chroma;
    WORD32 i4_ref_x_luma, i4_ref_y_luma;
    WORD32 i4_ref_x_chroma, i4_ref_y_chroma;
    WORD32 i4_ref_luma_ressam_correction_offset = 0;
    WORD32 i4_ref_chroma_ressam_correction_offset = 0;
    WORD32 i4_inp_luma_stride, i4_inp_chroma_stride;
    WORD32 i4_out_luma_stride, i4_out_chroma_stride;
    WORD32 i4_inp_luma_offset = 0, i4_inp_chroma_offset = 0;
    WORD32 ret;

    ps_svc_dec_ref_layer = ps_svc_lyr_dec->ps_dec_svc_ref_layer;
    ps_ctxt = (residual_sampling_ctxt_t *) ps_svc_lyr_dec->pv_residual_sample_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    i4_inp_luma_stride = ps_svc_dec_ref_layer->u2_residual_resample_luma_stride;
    i4_inp_chroma_stride = ps_svc_dec_ref_layer->u2_residual_resample_chroma_stride;
    i4_out_luma_stride = ps_svc_lyr_dec->u2_residual_resample_luma_stride;
    i4_out_chroma_stride = ps_svc_lyr_dec->u2_residual_resample_chroma_stride;

    {
        residual_samp_map_ctxt_t *ps_luma_map, *ps_chroma_map;

        ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
        ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;
        ps_x_off_len_luma = ps_luma_map->ps_x_offset_length;
        ps_y_off_len_luma = ps_luma_map->ps_y_offset_length;
        ps_x_off_len_chroma = ps_chroma_map->ps_x_offset_length;
        ps_y_off_len_chroma = ps_chroma_map->ps_y_offset_length;
    }
    i4_ref_x_luma = ps_svc_lyr_dec->ps_ressam_luma_map_horz[ps_cur_mb_info->u2_mbx].i2_offset;
    i4_ref_y_luma = ps_svc_lyr_dec->ps_ressam_luma_map_vert[ps_cur_mb_info->u2_mby].i2_offset;

    i4_ref_x_chroma = ps_svc_lyr_dec->ps_ressam_chroma_map_horz[ps_cur_mb_info->u2_mbx].i2_offset;
    i4_ref_y_chroma = ps_svc_lyr_dec->ps_ressam_chroma_map_vert[ps_cur_mb_info->u2_mby].i2_offset;

    i4_ref_x_luma = CLIP3(0, (ps_lyr_ctxt->i4_ref_width - 1), i4_ref_x_luma);
    i4_ref_y_luma = CLIP3(0, (ps_lyr_ctxt->i4_ref_height - 1), i4_ref_y_luma);
    i4_ref_x_chroma = CLIP3(0, ((ps_lyr_ctxt->i4_ref_width >> 1) - 1), i4_ref_x_chroma);
    i4_ref_y_chroma = CLIP3(0, ((ps_lyr_ctxt->i4_ref_height >> 1) - 1), i4_ref_y_chroma);

    {
        WORD32 i4_offset_x, i4_offset_y;

        i4_offset_x = ps_x_off_len_luma[ps_cur_mb_info->u2_mbx].i2_offset;
        i4_offset_y = ps_y_off_len_luma[ps_cur_mb_info->u2_mby].i2_offset;

        /* check for offsets inside frame dimensions */
        if(0 <= i4_offset_x)
        {
            /* validity of pointer passed */
            if(!(i4_offset_x >= i4_ref_x_luma))
            {
                return NOT_OK;
            }
            i4_inp_luma_offset += (i4_offset_x - i4_ref_x_luma);
        }

        if(0 <= i4_offset_y)
        {
            /* validity of pointer passed */
            if(!(i4_offset_y >= i4_ref_y_luma))
            {
                return NOT_OK;
            }
            i4_inp_luma_offset += (i4_offset_y - i4_ref_y_luma) * i4_inp_luma_stride;
        }

        i4_offset_x = ps_x_off_len_chroma[ps_cur_mb_info->u2_mbx].i2_offset;
        i4_offset_y = ps_y_off_len_chroma[ps_cur_mb_info->u2_mby].i2_offset;

        /* check for offsets inside frame dimensions */
        if(0 <= i4_offset_x)
        {
            /* validity of pointer passed */
            if(!(i4_offset_x >= i4_ref_x_chroma))
            {
                return NOT_OK;
            }
            i4_inp_chroma_offset += (i4_offset_x - i4_ref_x_chroma) << 1;
        }

        if(0 <= i4_offset_y)
        {
            /* validity of pointer passed */
            if(!(i4_offset_y >= i4_ref_y_chroma))
            {
                return NOT_OK;
            }
            i4_inp_chroma_offset += (i4_offset_y - i4_ref_y_chroma) * (i4_inp_chroma_stride << 1);
        }
    }

    i4_ref_luma_ressam_correction_offset = i4_ref_x_luma + (i4_ref_y_luma) *i4_inp_luma_stride;

    s_inp_luma.pv_buffer = ps_svc_dec_ref_layer->pi2_il_residual_resample_mb_luma_frm_start +
                           i4_inp_luma_offset + i4_ref_luma_ressam_correction_offset;
    s_inp_luma.i4_element_size = 1;
    s_inp_luma.i4_num_element_stride = i4_inp_luma_stride;

    i4_ref_chroma_ressam_correction_offset =
        (i4_ref_x_chroma << 1) + (i4_ref_y_chroma) *i4_inp_chroma_stride;

    s_inp_chroma.pv_buffer = ps_svc_dec_ref_layer->pi2_il_residual_resample_mb_chroma_frm_start +
                             i4_inp_chroma_offset + i4_ref_chroma_ressam_correction_offset;
    s_inp_chroma.i4_element_size = 1;
    s_inp_chroma.i4_num_element_stride = i4_inp_luma_stride;

    s_ref_mb_mode.pv_buffer = ps_svc_dec_ref_layer->ps_inter_lyr_mb_prms_frm_start;
    s_ref_mb_mode.i4_element_size = sizeof(inter_lyr_mb_prms_t);
    s_ref_mb_mode.i4_num_element_stride = ps_svc_dec_ref_layer->u2_inter_lyr_mb_prms_stride;

    s_out_luma.i4_element_size = 1;
    s_out_luma.pv_buffer =
        ps_svc_lyr_dec->pi2_il_residual_resample_mb_luma_frm_start +
        ((ps_cur_mb_info->u2_mbx << 4) +
         (i4_out_luma_stride * (ps_cur_mb_info->u2_mby << 4)) * s_out_luma.i4_element_size);

    s_out_luma.i4_num_element_stride = i4_out_luma_stride;

    s_out_chroma.i4_element_size = 1;
    s_out_chroma.pv_buffer =
        ps_svc_lyr_dec->pi2_il_residual_resample_mb_chroma_frm_start +
        ((ps_cur_mb_info->u2_mbx << 4) +
         (i4_out_chroma_stride * (ps_cur_mb_info->u2_mby << 3)) * s_out_chroma.i4_element_size);
    s_out_chroma.i4_num_element_stride = i4_out_chroma_stride;

    ret = ps_lyr_ctxt->pf_residual_samp_mb(ps_ctxt, &s_inp_luma, &s_inp_chroma, &s_ref_mb_mode,
                                           &s_out_luma, &s_out_chroma, ps_cur_mb_info->u2_mbx,
                                           ps_cur_mb_info->u2_mby);
    if(ret != OK)
    {
        return ret;
    }
    return OK;
}

/*!
 **************************************************************************
 * \if Function name : isvcd_process_inter_mb_rsd_pred_target_lyr \endif
 *
 * \brief IT+ Residual + Recon
 *    This function decodes an Inter MB.
 *
 *
 * \return
 *    0 on Success and Error code otherwise
 **************************************************************************
 */
WORD32 isvcd_process_inter_mb_rsd_pred_target_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                  dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mb_num,
                                                  UWORD8 u1_inference_mode,
                                                  UWORD16 *pu2_res_luma_csbp)
{
    UWORD8 *pu1_rec_y, *pu1_rec_u;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD32 ui_rec_width, u4_recwidth_cr;
    UWORD16 u2_luma_stride, u2_chroma_stride;
    WORD16 *pi2_y_coeff, *pi2_luma_res_ptr, *pi2_chroma_res_ptr;
    UWORD32 u1_mb_field_decoding_flag;
    const UWORD8 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    UWORD32 uc_botMb;
    UWORD32 u4_num_pmbair;
    tfr_ctxt_t *ps_frame_buf = ps_dec->ps_frame_buf_ip_recon;
    UWORD32 u4_luma_dc_only_csbp = 0;
    UWORD32 u4_luma_dc_only_cbp = 0;
    UWORD16 u2_res_luma_csbp = 0;
    WORD32 ret;

    if(0 != ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
    {
        return NOT_OK;
    }
    uc_botMb = 1 - ps_cur_mb_info->u1_topmb;
    u4_num_pmbair = (u4_mb_num >> u1_mbaff);
    u1_mb_field_decoding_flag = ps_cur_mb_info->u1_mb_field_decodingflag;

    pu1_rec_y = ps_frame_buf->pu1_dest_y + (u4_num_pmbair << 4);
    pu1_rec_u = ps_frame_buf->pu1_dest_u + (u4_num_pmbair << 3) * YUV420SP_FACTOR;
    ui_rec_width = ps_dec->u2_frm_wd_y << u1_mb_field_decoding_flag;
    u4_recwidth_cr = ps_dec->u2_frm_wd_uv << u1_mb_field_decoding_flag;

    u2_luma_stride = ps_svc_lyr_dec->u2_residual_resample_luma_stride;
    pi2_luma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_luma_frm_start +
                       (ps_cur_mb_info->u2_mbx << 4) +
                       ((ps_cur_mb_info->u2_mby << 4) * u2_luma_stride);

    u2_chroma_stride = ps_svc_lyr_dec->u2_residual_resample_chroma_stride;
    pi2_chroma_res_ptr = ps_svc_lyr_dec->pi2_il_residual_resample_mb_chroma_frm_start +
                         (ps_cur_mb_info->u2_mbx << 4) +
                         ((ps_cur_mb_info->u2_mby << 3) * u2_chroma_stride);

    ret = isvcd_process_residual_resample_mb(ps_svc_lyr_dec, ps_cur_mb_info);
    if(ret != OK)
    {
        return ret;
    }
    if(u1_mbaff)
    {
        if(uc_botMb)
        {
            pu1_rec_y += (u1_mb_field_decoding_flag ? (ui_rec_width >> 1) : (ui_rec_width << 4));
            pu1_rec_u +=
                (u1_mb_field_decoding_flag ? (u4_recwidth_cr >> 1) : (u4_recwidth_cr << 3));
        }
    }

    if(!ps_cur_mb_info->u1_tran_form8x8)
    {
        u4_luma_dc_only_csbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
    }
    else
    {
        if(!ps_dec->ps_cur_pps->u1_entropy_coding_mode)
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff4x4_mb(ps_dec, ps_cur_mb_info, 0);
        }
        else
        {
            u4_luma_dc_only_cbp = ih264d_unpack_luma_coeff8x8_mb(ps_dec, ps_cur_mb_info);
        }
    }

    *pu2_res_luma_csbp = 0;
    pi2_y_coeff = ps_dec->pi2_coeff_data;

    /* Inverse Transform and Reconstruction */
    if(ps_cur_mb_info->u1_cbp & 0x0f)
    {
        if(!ps_cur_mb_info->u1_tran_form8x8)
        {
            UWORD32 i;
            WORD16 ai2_tmp[16] = {0};
            for(i = 0; i < 16; i++)
            {
                if(CHECKBIT(ps_cur_mb_info->u2_luma_csbp, i))
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    UWORD8 *pu1_pred_sblk =
                        pu1_rec_y + ((i & 0x3) * BLK_SIZE) + (i >> 2) * (ui_rec_width << 2);
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x3) * BLK_SIZE) +
                                      (i >> 2) * (u2_luma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_csbp, i))
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4_dc(
                                    pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, ui_rec_width,
                                    u2_luma_stride, ui_rec_width,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                    ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4(
                                    pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, ui_rec_width,
                                    u2_luma_stride, ui_rec_width,
                                    gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[3],
                                    ps_cur_mb_info->u1_qp_div6, ai2_tmp, 0, NULL);
                        }
                    }
                }
                else
                {
                    UWORD8 *pu1_pred_sblk =
                        pu1_rec_y + ((i & 0x3) * BLK_SIZE) + (i >> 2) * (ui_rec_width << 2);
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x3) * BLK_SIZE) +
                                      (i >> 2) * (u2_luma_stride << 2);

                    u2_res_luma_csbp = ps_svc_lyr_dec->pf_pred_residual_recon_luma_4x4(
                        pu1_pred_sblk, pi2_out, pu1_pred_sblk, ui_rec_width, u2_luma_stride,
                        ui_rec_width);
                }
                *pu2_res_luma_csbp |= (u2_res_luma_csbp << i);
            }
        }
        else
        {
            WORD16 *pi2_scale_matrix_ptr;
            WORD32 i;

            pi2_scale_matrix_ptr = ps_dec->s_high_profile.i2_scalinglist8x8[1];

            for(i = 0; i < 4; i++)
            {
                WORD16 ai2_tmp[64] = {0};
                WORD16 *pi16_levelBlock =
                    pi2_y_coeff + (i << 6); /* move to the next 8x8 adding 64 */

                UWORD8 *pu1_pred_sblk =
                    pu1_rec_y + ((i & 0x1) * BLK8x8SIZE) + (i >> 1) * (ui_rec_width << 3);
                WORD16 *pi2_out =
                    pi2_luma_res_ptr + ((i & 0x1) * BLK8x8SIZE) + (i >> 1) * (u2_luma_stride << 3);
                if(CHECKBIT(ps_cur_mb_info->u1_cbp, i))
                {
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u4_luma_dc_only_cbp, i))
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8_dc(
                                    pi16_levelBlock, pu1_pred_sblk, pi2_out, pu1_pred_sblk,
                                    ui_rec_width, u2_luma_stride, ui_rec_width,
                                    gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                    ai2_tmp, 0, NULL);
                        }
                        else
                        {
                            u2_res_luma_csbp =
                                ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8(
                                    pi16_levelBlock, pu1_pred_sblk, pi2_out, pu1_pred_sblk,
                                    ui_rec_width, u2_luma_stride, ui_rec_width,
                                    gau1_ih264d_dequant8x8_cavlc[ps_cur_mb_info->u1_qp_rem6],
                                    (UWORD16 *) pi2_scale_matrix_ptr, ps_cur_mb_info->u1_qp_div6,
                                    ai2_tmp, 0, NULL);
                        }
                    }
                }
                else
                {
                    UWORD8 *pu1_pred_sblk =
                        pu1_rec_y + ((i & 0x1) * BLK8x8SIZE) + (i >> 1) * (ui_rec_width << 3);
                    WORD16 *pi2_out = pi2_luma_res_ptr + ((i & 0x1) * BLK8x8SIZE) +
                                      (i >> 1) * (u2_luma_stride << 3);

                    u2_res_luma_csbp = ps_svc_lyr_dec->pf_pred_residual_recon_luma_8x8(
                        pu1_pred_sblk, pi2_out, pu1_pred_sblk, ui_rec_width, u2_luma_stride,
                        ui_rec_width);
                }
                *pu2_res_luma_csbp |= (u2_res_luma_csbp << (((i >> 1) << 3) + ((i & 0x01) << 1)));
            }
        }
    }
    else
    {
        UWORD8 *pu1_pred_sblk = pu1_rec_y;
        WORD16 *pi2_out = pi2_luma_res_ptr;

        *pu2_res_luma_csbp = ps_svc_lyr_dec->pf_pred_residual_recon_luma_16x16(
            pu1_pred_sblk, pi2_out, pu1_pred_sblk, ui_rec_width, u2_luma_stride, ui_rec_width);
    }

    /* Decode Chroma Block */
    ih264d_unpack_chroma_coeff4x4_mb(ps_dec, ps_cur_mb_info);
    /*--------------------------------------------------------------------*/
    /* Chroma Blocks decoding                                             */
    /*--------------------------------------------------------------------*/
    {
        UWORD8 u1_chroma_cbp = (UWORD8) (ps_cur_mb_info->u1_cbp >> 4);

        if(u1_chroma_cbp != CBPC_ALLZERO)
        {
            UWORD32 u4_scale_u = ps_cur_mb_info->u1_qpc_div6;
            UWORD32 u4_scale_v = ps_cur_mb_info->u1_qpcr_div6;
            UWORD16 u2_chroma_csbp = ps_cur_mb_info->u2_chroma_csbp;

            pi2_y_coeff = ps_dec->pi2_coeff_data;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    UWORD8 *pu1_pred_sblk = pu1_rec_u + ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                            (i >> 1) * (u4_recwidth_cr << 2);
                    WORD16 *pi2_out = pi2_chroma_res_ptr +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4(
                                pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4], u4_scale_u,
                                ai2_tmp, pi2_level);
                        }
                        else if(pi2_level[0] != 0)
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4_dc(
                                pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpc_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[4], u4_scale_u,
                                ai2_tmp, pi2_level);
                        }
                        else
                        {
                            ps_svc_lyr_dec->pf_pred_residual_recon_chroma_4x4(
                                pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr);
                        }
                    }
                }
            }

            pi2_y_coeff += MB_CHROM_SIZE;
            u2_chroma_csbp >>= 4;

            {
                UWORD32 i;
                WORD16 ai2_tmp[16] = {0};
                for(i = 0; i < 4; i++)
                {
                    WORD16 *pi2_level = pi2_y_coeff + (i << 4);
                    UWORD8 *pu1_pred_sblk = pu1_rec_u + 1 +
                                            ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                            (i >> 1) * (u4_recwidth_cr << 2);
                    WORD16 *pi2_out = pi2_chroma_res_ptr + 1 +
                                      ((i & 0x1) * BLK_SIZE * YUV420SP_FACTOR) +
                                      (i >> 1) * (u2_chroma_stride << 2);
                    PROFILE_DISABLE_IQ_IT_RECON()
                    {
                        if(CHECKBIT(u2_chroma_csbp, i))
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4(
                                pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5], u4_scale_v,
                                ai2_tmp, pi2_level);
                        }
                        else if(pi2_level[0] != 0)
                        {
                            ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4_dc(
                                pi2_level, pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr,
                                gau2_ih264_iquant_scale_4x4[ps_cur_mb_info->u1_qpcr_rem6],
                                (UWORD16 *) ps_dec->s_high_profile.i2_scalinglist4x4[5], u4_scale_v,
                                ai2_tmp, pi2_level);
                        }
                        else
                        {
                            ps_svc_lyr_dec->pf_pred_residual_recon_chroma_4x4(
                                pu1_pred_sblk, pi2_out, pu1_pred_sblk, u4_recwidth_cr,
                                u2_chroma_stride, u4_recwidth_cr);
                        }
                    }
                }
            }
        }
        else
        {
            /* Cr*/
            {
                UWORD8 *pu1_pred_sblk = pu1_rec_u;
                WORD16 *pi2_out = pi2_chroma_res_ptr;

                ps_svc_lyr_dec->pf_pred_residual_recon_chroma_8x8(pu1_pred_sblk, pi2_out,
                                                                  pu1_pred_sblk, u4_recwidth_cr,
                                                                  u2_chroma_stride, u4_recwidth_cr);
            }

            /* Cb*/
            {
                UWORD8 *pu1_pred_sblk = pu1_rec_u + 1;
                WORD16 *pi2_out = pi2_chroma_res_ptr + 1;
                ps_svc_lyr_dec->pf_pred_residual_recon_chroma_8x8(pu1_pred_sblk, pi2_out,
                                                                  pu1_pred_sblk, u4_recwidth_cr,
                                                                  u2_chroma_stride, u4_recwidth_cr);
            }
        }
    }
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start + ps_cur_mb_info->u2_mbx +
        (ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride * (ps_cur_mb_info->u2_mby));
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_mb_mode =
        u1_inference_mode ? SVC_IBL_MB : SVC_INTER_MB;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->i1_tx_size = ps_cur_mb_info->u1_tran_form8x8;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u2_luma_nnz = ps_cur_mb_info->u2_luma_csbp;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz =
        (UWORD8) ps_cur_mb_info->u2_chroma_csbp;
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 1))
    {
        /* Four bits for Cb in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0x0F;
    }
    if(CHECKBIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 2))
    {
        /* Four bits for Cr in DC only cbp */
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_cur_mb->u1_chroma_nnz |= 0xF0;
    }
    return (0);
}
