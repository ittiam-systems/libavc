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
 *  isvcd_process_bslice.c
 *
 * @brief
 *  Contains routines that decode B slice type
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
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"
#include "ih264d_bitstrm.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_mvpred.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "ih264d_error_handler.h"
#include "ih264d_tables.h"
#include "ih264d_parse_slice.h"
#include "ih264d_process_pslice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_tables.h"
#include "ih264d_parse_islice.h"
#include "ih264d_mvpred.h"

/*!
**************************************************************************
* \if Function name : isvcd_one_to_one \endif
*
* \brief
*    Initializes forward and backward refernce lists for B slice decoding.
*
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
void isvcd_one_to_one(svc_dec_lyr_struct_t *ps_svc_lyr_dec, struct pic_buffer_t *ps_col_pic,
                      directmv_t *ps_direct, UWORD8 u1_wd_x, WORD32 u2_sub_mb_ofst,
                      dec_mb_info_t *ps_cur_mb_info)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 *pu1_col_zero_flag_start, u1_col_mb_pred_mode, u1_num_blks, u1_sub_mb_num;
    UWORD8 u1_init_colzero_flag;
    UNUSED(ps_cur_mb_info);
    pu1_col_zero_flag_start = ps_col_pic->pu1_col_zero_flag + u2_sub_mb_ofst;
    u1_col_mb_pred_mode = pu1_col_zero_flag_start[ps_dec->u1_sub_mb_num];
    u1_init_colzero_flag = u1_col_mb_pred_mode & 1;
    u1_col_mb_pred_mode >>= 6;
    ps_direct->u1_vert_mv_scale = ONE_TO_ONE;
    ps_direct->u1_col_zeroflag_change = (ps_svc_lyr_dec->u1_base_res_flag) ? 0 : 1;

    if(u1_wd_x == MB_SIZE)
    {
        ps_dec->u1_currB_type = (!!u1_col_mb_pred_mode);
        if(u1_col_mb_pred_mode == PRED_16x16)
        {
            ps_direct->i1_num_partitions = 1;
            ps_direct->i4_mv_indices[0] = u2_sub_mb_ofst;
            ps_direct->i1_submb_num[0] = 0;
            ps_direct->i1_partitionsize[0] = PRED_16x16;

            return;
        }
        else if(u1_col_mb_pred_mode < PRED_8x8)
        {
            ps_direct->i1_num_partitions = 2;
            ps_direct->i4_mv_indices[0] = u2_sub_mb_ofst;
            ps_direct->i1_submb_num[0] = 0;
            ps_direct->i1_partitionsize[0] = u1_col_mb_pred_mode;
            u1_sub_mb_num = (u1_col_mb_pred_mode == PRED_16x8) ? 8 : 2;
            ps_direct->i1_submb_num[1] = u1_sub_mb_num;
            ps_direct->i4_mv_indices[1] = u2_sub_mb_ofst + ps_direct->i1_submb_num[1];
            ps_direct->i1_partitionsize[1] = u1_col_mb_pred_mode;
            if((pu1_col_zero_flag_start[u1_sub_mb_num] & 1) != u1_init_colzero_flag)
                ps_direct->u1_col_zeroflag_change = 1;
            return;
        }
        else
        {
            u1_num_blks = 4;
        }
    }
    else
    {
        u1_num_blks = 1;
    }

    {
        const UWORD8 *pu1_top_lt_mb_part_idx;
        UWORD8 u1_col_sub_mb_pred_mode, uc_blk, u1_sub_blk, u1_submb_col = 0;
        UWORD8 u1_num_sub_blks, uc_direct8x8inf, *pu1_col_zero_flag, u1_sub_mb_num;
        const UWORD8 *pu1_num_sub_mb_part = (const UWORD8 *) gau1_ih264d_num_submb_part;
        UWORD8 i1_num_partitions = 0, partition_size;
        WORD32 mv_index;
        const UWORD8 *pu1_top_lt_sub_mb_idx = gau1_ih264d_submb_indx_mod_sp_drct;

        u1_sub_mb_num = ps_dec->u1_sub_mb_num;
        uc_direct8x8inf = ps_dec->ps_cur_slice->u1_direct_8x8_inference_flag;
        pu1_top_lt_mb_part_idx = gau1_ih264d_top_left_mb_part_indx_mod + (PRED_8x8 << 1) + 1;

        for(uc_blk = 0; uc_blk < u1_num_blks; uc_blk++)
        {
            partition_size = PRED_8x8;
            pu1_top_lt_sub_mb_idx = gau1_ih264d_submb_indx_mod_sp_drct;
            if(uc_direct8x8inf == 1)
            {
                u1_submb_col = u1_sub_mb_num | (u1_sub_mb_num >> 1);
                mv_index = u2_sub_mb_ofst + u1_submb_col;
                u1_num_sub_blks = 1;
            }
            else
            {
                /* colMbPart is either 8x8, 8x4, 4x8, 4x4 */
                pu1_col_zero_flag = pu1_col_zero_flag_start + u1_sub_mb_num;
                u1_col_sub_mb_pred_mode = *pu1_col_zero_flag;
                u1_col_sub_mb_pred_mode = (u1_col_sub_mb_pred_mode & 0x30) >> 4;
                partition_size = (UWORD8) ((u1_col_sub_mb_pred_mode) | (PRED_8x8 << 2));
                mv_index = u2_sub_mb_ofst + u1_sub_mb_num;
                pu1_top_lt_sub_mb_idx += (u1_col_sub_mb_pred_mode << 1);
                u1_num_sub_blks = pu1_num_sub_mb_part[u1_col_sub_mb_pred_mode];
            }

            for(u1_sub_blk = 0; u1_sub_blk < u1_num_sub_blks; u1_sub_blk++, pu1_top_lt_sub_mb_idx++)
            {
                u1_sub_mb_num += *pu1_top_lt_sub_mb_idx;
                mv_index += *pu1_top_lt_sub_mb_idx;
                ps_direct->i4_mv_indices[i1_num_partitions] = mv_index;
                ps_direct->i1_submb_num[i1_num_partitions] = u1_sub_mb_num;
                ps_direct->i1_partitionsize[i1_num_partitions] = partition_size;
                i1_num_partitions++;
                if(!uc_direct8x8inf) u1_submb_col = u1_sub_mb_num;
                if((pu1_col_zero_flag_start[u1_submb_col] & 1) != u1_init_colzero_flag)
                    ps_direct->u1_col_zeroflag_change = 1;
            }
            u1_sub_mb_num = *pu1_top_lt_mb_part_idx++;
        }
        ps_direct->i1_num_partitions = i1_num_partitions;
    }
}

/*!
 **************************************************************************
 * \if Function name : isvcd_decode_spatial_direct \endif
 *
 * \brief
 *    Decodes spatial direct mode.
 *
 * \return
 *    None.
 *    Vijay
 **************************************************************************
 */
WORD32 isvcd_decode_spatial_direct(dec_struct_t *ps_dec, UWORD8 u1_wd_x,
                                   dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mb_num)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) ps_dec;
    mv_pred_t s_mv_pred = {0};
    mv_pred_t *ps_mv;
    UWORD8 u1_col_zero_flag, u1_direct_zero_pred_flag = 0;
    UWORD32 u4_sub_mb_num;
    UWORD8 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    mv_pred_t *ps_mv_ntop_start;
    mv_pred_t *ps_mv_nmb_start = ps_dec->ps_mv_cur + (u4_mb_num << 4);
    UWORD8 partition_size, sub_partition, u1_mb_partw, u1_mb_parth;
    UWORD8 i;
    WORD8 i1_pred, i1_ref_frame0, i1_ref_frame1;
    struct pic_buffer_t *ps_ref_frame = NULL, *ps_col_pic, *ps_pic_buff0 = NULL,
                        *ps_pic_buff1 = NULL;

    UWORD8 u1_zero_pred_cond_f, u1_zero_pred_cond_b;
    WORD16 i2_spat_pred_mv[4] = {0};
    WORD16 *pi2_final_mv0, *pi2_final_mv1;
    UWORD16 ui2_mask_fwd = 0, ui2_mask_bwd = 0;
    UWORD32 *pui32_weight_ofsts = NULL;
    directmv_t s_mvdirect = {0};
    UWORD8 u1_colz;
    UWORD8 u1_final_ref_idx = 0;
    const UWORD8 *pu1_mb_parth = (const UWORD8 *) gau1_ih264d_mb_parth;
    const UWORD8 *pu1_mb_partw = (const UWORD8 *) gau1_ih264d_mb_partw;

    mv_pred_t s_temp_mv_pred = {0};
    ps_mv_ntop_start =
        ps_dec->ps_mv_cur + (u4_mb_num << 4) - (ps_dec->u2_frm_wd_in_mbs << (4 + u1_mbaff)) + 12;

    u1_direct_zero_pred_flag =
        ps_dec->pf_mvpred(ps_dec, ps_cur_mb_info, (ps_mv_nmb_start + ps_dec->u1_sub_mb_num),
                          ps_mv_ntop_start + (ps_dec->u1_sub_mb_num & 0x03), &s_mv_pred,
                          ps_dec->u1_sub_mb_num, (u1_wd_x >> 2), 0, 1, B_DIRECT_SPATIAL);

    i2_spat_pred_mv[0] = s_mv_pred.i2_mv[0];
    i2_spat_pred_mv[1] = s_mv_pred.i2_mv[1];
    i2_spat_pred_mv[2] = s_mv_pred.i2_mv[2];
    i2_spat_pred_mv[3] = s_mv_pred.i2_mv[3];

    i1_ref_frame0 = s_mv_pred.i1_ref_frame[0];
    i1_ref_frame1 = s_mv_pred.i1_ref_frame[1];

    i1_ref_frame0 = (i1_ref_frame0 < 0) ? -1 : i1_ref_frame0;
    i1_ref_frame1 = (i1_ref_frame1 < 0) ? -1 : i1_ref_frame1;

    i1_pred = 0;

    {
        WORD8 u1_ref_idx, u1_ref_idx1;
        UWORD32 uc_Idx, uc_Idx1;
        UWORD8 u1_scale_ref =
            (ps_dec->ps_cur_slice->u1_mbaff_frame_flag && ps_cur_mb_info->u1_mb_field_decodingflag);
        u1_final_ref_idx = i1_ref_frame0;
        if(i1_ref_frame0 >= 0)
        {
            /* convert RefIdx if it is MbAff */
            u1_ref_idx = i1_ref_frame0;
            u1_ref_idx1 = i1_ref_frame0;
            if(u1_scale_ref)
            {
                u1_ref_idx1 = u1_ref_idx >> 1;
                if((u1_ref_idx & 0x01) != (1 - ps_cur_mb_info->u1_topmb))
                    u1_ref_idx1 += MAX_REF_BUFS;
            }
            /* If i1_ref_frame0 < 0 then refIdxCol is obtained from ps_pic_buff1 */
            ps_pic_buff0 = ps_dec->ps_ref_pic_buf_lx[0][u1_ref_idx1];
            ps_ref_frame = ps_pic_buff0;
            i1_pred = PRED_L0;
        }

        if(i1_ref_frame1 >= 0)
        {
            /* convert RefIdx if it is MbAff */
            u1_ref_idx = i1_ref_frame1;
            u1_ref_idx1 = i1_ref_frame1;
            if(u1_scale_ref)
            {
                u1_ref_idx1 = u1_ref_idx >> 1;
                if((u1_ref_idx & 0x01) != (1 - ps_cur_mb_info->u1_topmb))
                    u1_ref_idx1 += MAX_REF_BUFS;
            }
            ps_pic_buff1 = ps_dec->ps_ref_pic_buf_lx[1][u1_ref_idx1];
            i1_pred = i1_pred | PRED_L1;
        }
        if(i1_ref_frame0 < 0)
        {
            ps_ref_frame = ps_pic_buff1;
            u1_final_ref_idx = i1_ref_frame1;
        }

        u1_zero_pred_cond_f = (u1_direct_zero_pred_flag) || (i1_ref_frame0 < 0);
        u1_zero_pred_cond_b = (u1_direct_zero_pred_flag) || (i1_ref_frame1 < 0);

        if(ps_dec->ps_cur_pps->u1_wted_bipred_idc)
        {
            uc_Idx = ((i1_ref_frame0 < 1) ? 0 : i1_ref_frame0) *
                     ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1];
            if(u1_scale_ref) uc_Idx >>= 1;
            uc_Idx1 = (i1_ref_frame1 < 0) ? 0 : i1_ref_frame1;
            uc_Idx += (u1_scale_ref) ? (uc_Idx1 >> 1) : uc_Idx1;
            pui32_weight_ofsts = (UWORD32 *) &ps_dec->pu4_wt_ofsts[2 * X3(uc_Idx)];

            if(i1_ref_frame0 < 0) pui32_weight_ofsts += 1;

            if(u1_scale_ref && (ps_dec->ps_cur_pps->u1_wted_bipred_idc == 2))
            {
                WORD16 i2_ref_idx;
                i2_ref_idx = MAX(i1_ref_frame0, 0);
                i2_ref_idx *= (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1] << 1);
                i2_ref_idx += MAX(i1_ref_frame1, 0);
                if(!ps_cur_mb_info->u1_topmb)
                    i2_ref_idx += (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[0] << 1) *
                                  (ps_dec->ps_cur_slice->u1_num_ref_idx_lx_active[1] << 1);
                pui32_weight_ofsts = (UWORD32 *) &ps_dec->pu4_mbaff_wt_mat[2 * X3(i2_ref_idx)];
            }
        }
    }

    s_temp_mv_pred.i1_ref_frame[0] = i1_ref_frame0;
    s_temp_mv_pred.i1_ref_frame[1] = i1_ref_frame1;
    s_temp_mv_pred.u1_col_ref_pic_idx = ps_ref_frame->u1_mv_buf_id;
    s_temp_mv_pred.u1_pic_type = ps_ref_frame->u1_pic_type;

    /**********************************************************************/
    /* Call the function which gets the number of partitions and          */
    /* partition info of colocated Mb                                     */
    /**********************************************************************/

    isvcd_one_to_one(ps_svc_lyr_dec, ps_dec->ps_col_pic, &s_mvdirect, u1_wd_x,
                     ps_dec->i4_submb_ofst, ps_cur_mb_info);

    ps_col_pic = ps_dec->ps_col_pic;
    if((s_mvdirect.u1_col_zeroflag_change == 0) || u1_direct_zero_pred_flag)
    {
        WORD16 i2_mv_x, i2_mv_y, i2_mvX1, i2_mvY1;
        /* Most probable case */
        u1_col_zero_flag = *(ps_col_pic->pu1_col_zero_flag + s_mvdirect.i4_mv_indices[0]);
        u1_col_zero_flag = u1_col_zero_flag & 0x01;

        if(u1_zero_pred_cond_f || ((i1_ref_frame0 == 0) && (u1_col_zero_flag == 1)))
        {
            i2_mv_x = 0;
            i2_mv_y = 0;
        }
        else
        {
            i2_mv_x = i2_spat_pred_mv[0];
            i2_mv_y = i2_spat_pred_mv[1];
        }

        if(u1_zero_pred_cond_b || ((i1_ref_frame1 == 0) && (u1_col_zero_flag == 1)))
        {
            i2_mvX1 = 0;
            i2_mvY1 = 0;
        }
        else
        {
            i2_mvX1 = i2_spat_pred_mv[2];
            i2_mvY1 = i2_spat_pred_mv[3];
        }

        u4_sub_mb_num = ps_dec->u1_sub_mb_num;
        u1_mb_partw = (u1_wd_x >> 2);

        if(i1_ref_frame0 >= 0)
        {
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD16 i2_mv[2];
                WORD8 i1_ref_idx = 0;

                i2_mv[0] = i2_mv_x;
                i2_mv[1] = i2_mv_y;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(i2_mv, u1_mb_partw, u1_mb_partw, u4_sub_mb_num, i1_pred,
                                      ps_pred_pkd, ps_pic_buff0->u1_pic_buf_id, i1_ref_idx,
                                      pui32_weight_ofsts, ps_pic_buff0->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }
        }

        if(i1_ref_frame1 >= 0)
        {
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD16 i2_mv[2];
                WORD8 i1_ref_idx = 0;

                i2_mv[0] = i2_mvX1;
                i2_mv[1] = i2_mvY1;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(i2_mv, u1_mb_partw, u1_mb_partw, u4_sub_mb_num, i1_pred,
                                      ps_pred_pkd, ps_pic_buff1->u1_pic_buf_id, i1_ref_idx,
                                      pui32_weight_ofsts, ps_pic_buff1->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }
        }

        /* Replication optimisation */
        s_temp_mv_pred.i2_mv[0] = i2_mv_x;
        s_temp_mv_pred.i2_mv[1] = i2_mv_y;
        s_temp_mv_pred.i2_mv[2] = i2_mvX1;
        s_temp_mv_pred.i2_mv[3] = i2_mvY1;

        /* Calculating colocated zero information */
        {
            /*************************************/
            /* If(bit2 and bit3 set)             */
            /* then                              */
            /*  (bit0 and bit1) => submmbmode    */
            /*  (bit2 and bit3) => mbmode        */
            /* else                              */
            /*  (bit0 and bit1) => mbmode        */
            /*************************************/
            /*UWORD8 u1_packed_mb_sub_mb_mode = sub_partition ?
             (s_mvdirect.i1_partitionsize[0]) : ((s_mvdirect.i1_partitionsize[0]) <<
             2);*/
            UWORD8 u1_packed_mb_sub_mb_mode = (u1_mb_partw == 2) ? 0x03 : 0;

            if(i1_ref_frame0 < 0)
            {
                i2_mv_x = i2_mvX1;
                i2_mv_y = i2_mvY1;
            }

            /* Change from left shift 4 to 6 - Varun */
            u1_colz = (ps_cur_mb_info->u1_mb_field_decodingflag << 1) |
                      ((u1_final_ref_idx == 0) && (ABS(i2_mv_x) <= 1) && (ABS(i2_mv_y) <= 1));
            u1_colz |= (u1_packed_mb_sub_mb_mode << 6);
        }
        ps_mv = ps_mv_nmb_start + u4_sub_mb_num;
        if(ps_mv)
        {
            ih264d_rep_mv_colz(ps_dec, &s_temp_mv_pred, ps_mv, u4_sub_mb_num, u1_colz, u1_mb_partw,
                               u1_mb_partw);
        }
        else
        {
            return NOT_OK;
        }

        if(u1_wd_x == MB_SIZE) ps_dec->u1_currB_type = 0;

        return OK;
    }

    /***************************************************************************/
    /* If present MB is 16x16 and the partition of colocated Mb is >= PRED_8x8 */
    /* i.e 8x8 or less than 8x8 partitions then set up DMA for (0,0) and       */
    /* spatially predicted motion vector and do the multiplexing after         */
    /* motion compensation                                                     */
    /***************************************************************************/

    if((u1_wd_x == MB_SIZE) && (s_mvdirect.i1_num_partitions > 2))
    {
        ps_cur_mb_info->u1_Mux = 1;
        if(i1_ref_frame0 >= 0)
        {
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD8 i1_ref_idx = 0;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(&(i2_spat_pred_mv[0]), 4, 4, 0, i1_pred, ps_pred_pkd,
                                      ps_pic_buff0->u1_pic_buf_id, i1_ref_idx, pui32_weight_ofsts,
                                      ps_pic_buff0->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }

            /******    (0,0) Motion vectors DMA     *****/
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD16 i2_mv[2];
                WORD8 i1_ref_idx = 0;

                i2_mv[0] = 0;
                i2_mv[1] = 0;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(i2_mv, 4, 4, 0, i1_pred, ps_pred_pkd,
                                      ps_pic_buff0->u1_pic_buf_id, i1_ref_idx, pui32_weight_ofsts,
                                      ps_pic_buff0->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }
        }
        if(i1_ref_frame1 >= 0)
        {
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD8 i1_ref_idx = 0;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(&(i2_spat_pred_mv[2]), 4, 4, 0, i1_pred, ps_pred_pkd,
                                      ps_pic_buff1->u1_pic_buf_id, i1_ref_idx, pui32_weight_ofsts,
                                      ps_pic_buff1->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }

            /******    (0,0) Motion vectors DMA     *****/
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD16 i2_mv[2];
                WORD8 i1_ref_idx = 0;

                i2_mv[0] = 0;
                i2_mv[1] = 0;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(i2_mv, 4, 4, 0, i1_pred, ps_pred_pkd,
                                      ps_pic_buff1->u1_pic_buf_id, i1_ref_idx, pui32_weight_ofsts,
                                      ps_pic_buff1->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }
        }
    }

    for(i = 0; i < s_mvdirect.i1_num_partitions; i++)
    {
        partition_size = s_mvdirect.i1_partitionsize[i];
        u4_sub_mb_num = s_mvdirect.i1_submb_num[i];

        sub_partition = partition_size >> 2;
        partition_size &= 0x3;
        u1_mb_partw = pu1_mb_partw[partition_size];
        u1_mb_parth = pu1_mb_parth[partition_size];
        if(sub_partition != 0)
        {
            u1_mb_partw >>= 1;
            u1_mb_parth >>= 1;
        }

        u1_col_zero_flag = *(ps_col_pic->pu1_col_zero_flag + s_mvdirect.i4_mv_indices[i]);
        u1_col_zero_flag = u1_col_zero_flag & 0x01;

        /*if(u1_col != u1_col_zero_flag)
         u1_init = 1;*/

        pi2_final_mv0 = &i2_spat_pred_mv[0];
        pi2_final_mv1 = &i2_spat_pred_mv[2];

        if(ps_cur_mb_info->u1_Mux != 1)
        {
            if(i1_ref_frame0 >= 0)
            {
                {
                    pred_info_pkd_t *ps_pred_pkd;
                    WORD8 i1_ref_idx = 0;

                    ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                    ih264d_fill_pred_info(pi2_final_mv0, u1_mb_partw, u1_mb_parth, u4_sub_mb_num,
                                          i1_pred, ps_pred_pkd, ps_pic_buff0->u1_pic_buf_id,
                                          i1_ref_idx, pui32_weight_ofsts,
                                          ps_pic_buff0->u1_pic_type);
                    ps_dec->u4_pred_info_pkd_idx++;
                    ps_cur_mb_info->u1_num_pred_parts++;
                }
            }

            if(i1_ref_frame1 >= 0)
            {
                pred_info_pkd_t *ps_pred_pkd;
                WORD8 i1_ref_idx = 0;

                ps_pred_pkd = ps_dec->ps_pred_pkd + ps_dec->u4_pred_info_pkd_idx;
                ih264d_fill_pred_info(pi2_final_mv1, u1_mb_partw, u1_mb_parth, u4_sub_mb_num,
                                      i1_pred, ps_pred_pkd, ps_pic_buff1->u1_pic_buf_id, i1_ref_idx,
                                      pui32_weight_ofsts, ps_pic_buff1->u1_pic_type);
                ps_dec->u4_pred_info_pkd_idx++;
                ps_cur_mb_info->u1_num_pred_parts++;
            }
        }

        /* Replication optimisation */
        s_temp_mv_pred.i2_mv[0] = pi2_final_mv0[0];
        s_temp_mv_pred.i2_mv[1] = pi2_final_mv0[1];
        s_temp_mv_pred.i2_mv[2] = pi2_final_mv1[0];
        s_temp_mv_pred.i2_mv[3] = pi2_final_mv1[1];

        /* Calculating colocated zero information */
        {
            WORD16 i2_mv_x = 0, i2_mv_y = 0;
            /*************************************/
            /* If(bit2 and bit3 set)             */
            /* then                              */
            /*  (bit0 and bit1) => submmbmode    */
            /*  (bit2 and bit3) => mbmode        */
            /* else                              */
            /*  (bit0 and bit1) => mbmode        */
            /*************************************/
            UWORD8 u1_packed_mb_sub_mb_mode = sub_partition
                                                  ? (s_mvdirect.i1_partitionsize[i])
                                                  : ((s_mvdirect.i1_partitionsize[i]) << 2);

            if(i1_ref_frame0 >= 0)
            {
                i2_mv_x = pi2_final_mv0[0];
                i2_mv_y = pi2_final_mv0[1];
            }
            else
            {
                i2_mv_x = pi2_final_mv1[0];
                i2_mv_y = pi2_final_mv1[1];
            }

            u1_colz = (ps_cur_mb_info->u1_mb_field_decodingflag << 1) |
                      ((u1_final_ref_idx == 0) && (ABS(i2_mv_x) <= 1) && (ABS(i2_mv_y) <= 1));
            u1_colz |= (u1_packed_mb_sub_mb_mode << 4);
        }
        ps_mv = ps_mv_nmb_start + u4_sub_mb_num;
        if(ps_mv)
        {
            ih264d_rep_mv_colz(ps_dec, &s_temp_mv_pred, ps_mv, u4_sub_mb_num, u1_colz, u1_mb_parth,
                               u1_mb_partw);
        }
        else
        {
            return NOT_OK;
        }
    }
    i = 0;
    if(i1_ref_frame0 >= 0) ps_cur_mb_info->u2_mask[i++] = ui2_mask_fwd;
    if(i1_ref_frame1 >= 0) ps_cur_mb_info->u2_mask[i] = ui2_mask_bwd;

    return OK;
}
