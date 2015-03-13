/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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

/*!
 **************************************************************************
 * \file ih264d_thread_compute_bs.c
 *
 * \brief
 *    Contains routines that for multi-thread decoder
 *
 * Detailed_description
 *
 * \date
 *    20/02/2012
 *
 * \author  ZR
 **************************************************************************
 */
#include "ih264d_error_handler.h"
#include "ih264d_debug.h"
#include <string.h>
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "ih264d_structs.h"
#include "ih264d_defs.h"
#include "ih264d_mb_utils.h"

#include "ih264d_thread_compute_bs.h"
#include "ithread.h"
#include "ih264d_deblocking.h"
#include "ih264d_mb_utils.h"
#include "ih264d_tables.h"
#include "ih264d_format_conv.h"
#include "ih264d_defs.h"
UWORD16 ih264d_update_csbp_8x8(UWORD16 u2_luma_csbp);
void ih264d_fill_bs2_horz_vert(UWORD32 *pu4_bs, /* Base pointer of BS table */
                               WORD32 u4_left_mb_csbp, /* csbp of left mb */
                               WORD32 u4_top_mb_csbp, /* csbp of top mb */
                               WORD32 u4_cur_mb_csbp, /* csbp of current mb */
                               const UWORD32 *pu4_packed_bs2, const UWORD16 *pu2_4x4_v2h_reorder);

#define BS_MB_GROUP 4
#define DEBLK_MB_GROUP 1
#define FORMAT_CONV_MB_GROUP 4

/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264d_compute_bs_non_mbaff_thread                                           */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function for*/
/*                : Boundary Strength Calculation .this function is used     */
/*                : BS being calculated in separate thread                   */
/*  Inputs        : pointer to decoder context,cur_mb_info,u4_mb_num            */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : Produces the Boundary Strength for Current Mb            */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                      ITTIAM                                               */
/*****************************************************************************/

void ih264d_compute_bs_non_mbaff_thread(dec_struct_t * ps_dec,
                                        dec_mb_info_t * ps_cur_mb_info,
                                        UWORD32 u4_mb_num)
{
    /* Mvpred and Nnz for top and Courrent */
    mv_pred_t *ps_cur_mv_pred, *ps_top_mv_pred = NULL, *ps_left_mv_pred;
    /* deblk_mb_t Params */
    deblk_mb_t *ps_cur_mb_params; /*< Parameters of current MacroBlock */
    deblkmb_neighbour_t *ps_deblk_top_mb;

    /* Reference Index to POC mapping*/
    void ** apv_map_ref_idx_to_poc;
    UWORD32 u4_leftmbtype;

    UWORD16 u2_left_csbp, u2_top_csbp, u2_cur_csbp;

    /* Set of flags */
    UWORD32 u4_cur_mb_intra, u1_top_mb_typ, u4_cur_mb_fld;
    UWORD32 u1_cur_mb_type;
    UWORD32 * pu4_bs_table;

    /* Neighbour availability */
    /* Initialization */
    const UWORD32 u2_mbx = ps_cur_mb_info->u2_mbx;
    const UWORD32 u2_mby = ps_cur_mb_info->u2_mby;
    const UWORD32 u1_pingpong = u2_mbx & 0x01;
    ps_deblk_top_mb = ps_dec->ps_deblk_top_mb + u2_mbx;

    /* Pointer assignment for Current DeblkMB, Current Mv Pred  */
    ps_cur_mb_params = ps_dec->ps_deblk_pic + u4_mb_num;
    ps_cur_mv_pred = ps_dec->s_cur_pic.ps_mv + (u4_mb_num << 4);

    apv_map_ref_idx_to_poc =
                    (void **)ps_dec->ps_computebs_cur_slice->ppv_map_ref_idx_to_poc
                                    + 1;
    u1_cur_mb_type = ps_cur_mb_params->u1_mb_type;
    u1_top_mb_typ = ps_deblk_top_mb->u1_mb_type;
    ps_deblk_top_mb->u1_mb_type = u1_cur_mb_type;

    {
        ps_cur_mb_params->u1_topmb_qp = ps_deblk_top_mb->u1_mb_qp;
        ps_deblk_top_mb->u1_mb_qp = ps_cur_mb_params->u1_mb_qp;

        ps_cur_mb_params->u1_left_mb_qp = ps_dec->deblk_left_mb[1].u1_mb_qp;
        ps_dec->deblk_left_mb[1].u1_mb_qp = ps_cur_mb_params->u1_mb_qp;

    }

    /* if no deblocking required for current Mb then continue */
    /* Check next Mbs   in Mb group                           */
    if(ps_cur_mb_params->u1_deblocking_mode & MB_DISABLE_FILTERING)
    {
        void ** pu4_map_ref_idx_to_poc_l1 = apv_map_ref_idx_to_poc +
        POC_LIST_L0_TO_L1_DIFF;
        {
            /* Store Parameter for Top MvPred refernce frame Address */

            void ** ppv_top_mv_pred_addr = ps_cur_mb_info->ps_curmb->u4_pic_addrress;
            WORD8 * p1_refTop0 = (ps_cur_mv_pred + 12)->i1_ref_frame;
            WORD8 * p1_refTop1 = (ps_cur_mv_pred + 14)->i1_ref_frame;

            /* Store Left addresses for Next Mb   */
            void ** ppv_left_mv_pred_addr =
                            ps_dec->ps_left_mvpred_addr[!u1_pingpong][1].u4_add;
            WORD8 * p1_refleft0 = (ps_cur_mv_pred + 3)->i1_ref_frame;


            ppv_top_mv_pred_addr[0] = apv_map_ref_idx_to_poc[p1_refTop0[0]];
            ppv_top_mv_pred_addr[1] = pu4_map_ref_idx_to_poc_l1[p1_refTop0[1]];

            ppv_left_mv_pred_addr[2] = apv_map_ref_idx_to_poc[p1_refTop1[0]];
            ppv_top_mv_pred_addr[2] = apv_map_ref_idx_to_poc[p1_refTop1[0]];
            ppv_left_mv_pred_addr[3] = pu4_map_ref_idx_to_poc_l1[p1_refTop1[1]];
            ppv_top_mv_pred_addr[3] = pu4_map_ref_idx_to_poc_l1[p1_refTop1[1]];

            ppv_left_mv_pred_addr[0] = apv_map_ref_idx_to_poc[p1_refleft0[0]];
            ppv_left_mv_pred_addr[1] = pu4_map_ref_idx_to_poc_l1[p1_refleft0[1]];
            //}
            /* Storing the leftMbtype for next Mb */
            ps_dec->deblk_left_mb[1].u1_mb_type = ps_cur_mb_params->u1_mb_type;
        }

        return;
    }

    /* Flag for extra left Edge */
    ps_cur_mb_params->u1_single_call = 1;

    /* Update the Left deblk_mb_t and Left MvPred Parameters           */
    if(!u2_mbx)
    {
        u4_leftmbtype = 0;

        /* Initialize the ps_left_mv_pred with Junk but Valid Location */
        /* to avoid invalid memory access                           */
        /* this is read only pointer                                */
        ps_left_mv_pred = ps_cur_mv_pred + 3;
    }
    else
    {
        u4_leftmbtype = ps_dec->deblk_left_mb[1].u1_mb_type;

        /* Come to Left Most Edge of the MB */
        ps_left_mv_pred = ps_cur_mv_pred - (1 << 4) + 3;
    }

    if(!u2_mby)
        u1_top_mb_typ = 0;

    /* MvPred Pointer Calculation */
    /* CHANGED CODE */
    ps_top_mv_pred = ps_cur_mv_pred - (ps_dec->u2_frm_wd_in_mbs << 4) + 12;

    u4_cur_mb_intra = u1_cur_mb_type & D_INTRA_MB;
    u4_cur_mb_fld = !!(u1_cur_mb_type & D_FLD_MB);
    /* Compute BS function */
    pu4_bs_table = ps_cur_mb_params->u4_bs_table;

    u2_cur_csbp = ps_cur_mb_info->ps_curmb->u2_luma_csbp;
    u2_left_csbp = ps_cur_mb_info->ps_left_mb->u2_luma_csbp;
    u2_top_csbp = ps_cur_mb_info->ps_top_mb->u2_luma_csbp;

    /* Compute BS function */
    if(ps_dec->ps_cur_sps->u1_profile_idc == HIGH_PROFILE_IDC)
    {
        if(ps_cur_mb_info->u1_tran_form8x8 == 1)
        {
            u2_cur_csbp = ih264d_update_csbp_8x8(
                            ps_cur_mb_info->ps_curmb->u2_luma_csbp);
        }

        if(ps_cur_mb_info->ps_left_mb->u1_tran_form8x8 == 1)
        {
            u2_left_csbp = ih264d_update_csbp_8x8(
                            ps_cur_mb_info->ps_left_mb->u2_luma_csbp);
        }

        if(ps_cur_mb_info->ps_top_mb->u1_tran_form8x8 == 1)
        {
            u2_top_csbp = ih264d_update_csbp_8x8(
                            ps_cur_mb_info->ps_top_mb->u2_luma_csbp);
        }
    }
    if(u4_cur_mb_intra)
    {

        pu4_bs_table[4] = 0x04040404;
        pu4_bs_table[0] = u4_cur_mb_fld ? 0x03030303 : 0x04040404;
        pu4_bs_table[1] = 0x03030303;
        pu4_bs_table[2] = 0x03030303;
        pu4_bs_table[3] = 0x03030303;
        pu4_bs_table[5] = 0x03030303;
        pu4_bs_table[6] = 0x03030303;
        pu4_bs_table[7] = 0x03030303;
    }
    else
    {
        UWORD32 u4_is_non16x16 = !!(u1_cur_mb_type & D_PRED_NON_16x16);
        UWORD32 u4_is_b =
                        (ps_dec->ps_computebs_cur_slice->slice_type == B_SLICE);






        ih264d_fill_bs2_horz_vert(pu4_bs_table, u2_left_csbp, u2_top_csbp,
                                  u2_cur_csbp, gau4_ih264d_packed_bs2,
                                  gau2_ih264d_4x4_v2h_reorder);

        if(u4_leftmbtype & D_INTRA_MB)
            pu4_bs_table[4] = 0x04040404;

        if(u1_top_mb_typ & D_INTRA_MB)
            pu4_bs_table[0] = u4_cur_mb_fld ? 0x03030303 : 0x04040404;

        ps_dec->pf_fill_bs1[u4_is_b][u4_is_non16x16](
                        ps_cur_mv_pred, ps_top_mv_pred, apv_map_ref_idx_to_poc,
                        pu4_bs_table, ps_left_mv_pred,
                        &(ps_dec->ps_left_mvpred_addr[u1_pingpong][1]),
                        ps_cur_mb_info->ps_top_mb->u4_pic_addrress,
                        (4 >> u4_cur_mb_fld));
    }

    {
        void ** pu4_map_ref_idx_to_poc_l1 = apv_map_ref_idx_to_poc +
        POC_LIST_L0_TO_L1_DIFF;
        {
            /* Store Parameter for Top MvPred refernce frame Address */

            void ** ppv_top_mv_pred_addr = ps_cur_mb_info->ps_curmb->u4_pic_addrress;
            WORD8 * p1_refTop0 = (ps_cur_mv_pred + 12)->i1_ref_frame;
            WORD8 * p1_refTop1 = (ps_cur_mv_pred + 14)->i1_ref_frame;

            /* Store Left addresses for Next Mb   */
            void ** ppv_left_mv_pred_addr =
                            ps_dec->ps_left_mvpred_addr[!u1_pingpong][1].u4_add;
            WORD8 * p1_refleft0 = (ps_cur_mv_pred + 3)->i1_ref_frame;

            ppv_top_mv_pred_addr[0] = apv_map_ref_idx_to_poc[p1_refTop0[0]];
            ppv_top_mv_pred_addr[1] = pu4_map_ref_idx_to_poc_l1[p1_refTop0[1]];

            ppv_left_mv_pred_addr[2] = apv_map_ref_idx_to_poc[p1_refTop1[0]];
            ppv_top_mv_pred_addr[2] = apv_map_ref_idx_to_poc[p1_refTop1[0]];
            ppv_left_mv_pred_addr[3] = pu4_map_ref_idx_to_poc_l1[p1_refTop1[1]];
            ppv_top_mv_pred_addr[3] = pu4_map_ref_idx_to_poc_l1[p1_refTop1[1]];

            ppv_left_mv_pred_addr[0] = apv_map_ref_idx_to_poc[p1_refleft0[0]];
            ppv_left_mv_pred_addr[1] = pu4_map_ref_idx_to_poc_l1[p1_refleft0[1]];

            /* Storing the leftMbtype for next Mb */
            ps_dec->deblk_left_mb[1].u1_mb_type = ps_cur_mb_params->u1_mb_type;

        }
    }

    /* For transform 8x8 disable deblocking of the intrernal edges of a 8x8 block */
    if(ps_cur_mb_info->u1_tran_form8x8)
    {
        pu4_bs_table[1] = 0;
        pu4_bs_table[3] = 0;
        pu4_bs_table[5] = 0;
        pu4_bs_table[7] = 0;
    }
}

void ih264d_check_mb_map_deblk(dec_struct_t *ps_dec,
                               UWORD32 deblk_mb_grp,
                               tfr_ctxt_t *ps_tfr_cxt)
{
    UWORD32 i = 0;
    UWORD32 u4_mb_num;
    UWORD32 u4_cur_mb, u4_right_mb;
    volatile UWORD8 *mb_map = ps_dec->pu1_recon_mb_map;
    UWORD32 u4_mb_x, u4_mb_y, u4_image_wd_mb;
    deblk_mb_t *ps_cur_mb = ps_dec->ps_cur_deblk_thrd_mb;
    deblk_mb_t *ps_top_mb;
    deblk_mb_t *ps_left_mb;
    const WORD32 i4_cb_qp_idx_ofst =
                    ps_dec->ps_cur_pps->i1_chroma_qp_index_offset;
    const WORD32 i4_cr_qp_idx_ofst =
                    ps_dec->ps_cur_pps->i1_second_chroma_qp_index_offset;

    UWORD32 u4_wd_y, u4_wd_uv;
    UWORD8 u1_field_pic_flag = ps_dec->ps_cur_slice->u1_field_pic_flag;

    u4_mb_num = ps_dec->u4_cur_deblk_mb_num;
    u4_mb_x = ps_dec->u4_deblk_mb_x;
    u4_mb_y = ps_dec->u4_deblk_mb_y;
    u4_image_wd_mb = ps_dec->u2_frm_wd_in_mbs;
    u4_wd_y = ps_dec->u2_frm_wd_y << u1_field_pic_flag;
    u4_wd_uv = ps_dec->u2_frm_wd_uv << u1_field_pic_flag;
    ps_cur_mb = ps_dec->ps_cur_deblk_thrd_mb;

    for(i = 0; i < deblk_mb_grp; i++)
    {

        //while(1)
        //{
        CHECK_MB_MAP_BYTE(u4_mb_num, mb_map, u4_cur_mb);

        if(ps_dec->u4_cur_bs_mb_num <= u4_mb_num)
            u4_cur_mb = 0;

        if(u4_mb_x < (u4_image_wd_mb - 1))
        {
            CHECK_MB_MAP_BYTE((u4_mb_num + 1), mb_map, u4_right_mb);
        }
        else
            u4_right_mb = 1;

        if((u4_cur_mb && u4_right_mb) == 0)
        {
            break;
        }
        else
        {

        }
        //}

        u4_mb_num++;
        {
            UWORD32 u4_deb_mode, u4_mbs_next;
            u4_deb_mode = ps_cur_mb->u1_deblocking_mode;
            if(!(u4_deb_mode & MB_DISABLE_FILTERING))
            {

                if(u4_mb_x)
                {
                    ps_left_mb = ps_cur_mb - 1;

                }
                else
                {
                    ps_left_mb = NULL;

                }
                if(u4_mb_y != 0)
                {
                    ps_top_mb = ps_cur_mb - (u4_image_wd_mb);
                }
                else
                {
                    ps_top_mb = NULL;
                }

                if(u4_deb_mode & MB_DISABLE_LEFT_EDGE)
                    ps_left_mb = NULL;
                if(u4_deb_mode & MB_DISABLE_TOP_EDGE)
                    ps_top_mb = NULL;

                ih264d_deblock_mb_nonmbaff(ps_dec, ps_tfr_cxt,
                                           i4_cb_qp_idx_ofst, i4_cr_qp_idx_ofst,
                                           ps_cur_mb, u4_wd_y, u4_wd_uv,
                                           ps_top_mb, ps_left_mb);

            }

            ps_cur_mb++;
            u4_mb_x++;
            u4_mbs_next = u4_image_wd_mb - u4_mb_x;

            ps_tfr_cxt->pu1_mb_y += 16;
            ps_tfr_cxt->pu1_mb_u += 8 * YUV420SP_FACTOR;
            ps_tfr_cxt->pu1_mb_v += 8;

            if(!u4_mbs_next)
            {
                ps_tfr_cxt->pu1_mb_y += ps_tfr_cxt->u4_y_inc;
                ps_tfr_cxt->pu1_mb_u += ps_tfr_cxt->u4_uv_inc;
                ps_tfr_cxt->pu1_mb_v += ps_tfr_cxt->u4_uv_inc;
                u4_mb_y++;
                u4_mb_x = 0;
            }
        }

    }

    ps_dec->u4_cur_deblk_mb_num = u4_mb_num;
    ps_dec->u4_deblk_mb_x = u4_mb_x;
    ps_dec->u4_deblk_mb_y = u4_mb_y;
    ps_dec->ps_cur_deblk_thrd_mb = ps_cur_mb;

}

void ih264d_check_mb_map_deblk_wait(dec_struct_t *ps_dec,
                                    UWORD32 deblk_mb_grp,
                                    tfr_ctxt_t *ps_tfr_cxt)
{
    UWORD32 i = 0;
    UWORD32 u4_mb_num;
    UWORD32 u4_cur_mb, u4_right_mb;
    volatile UWORD8 *mb_map = ps_dec->pu1_recon_mb_map;
    UWORD32 u4_mb_x, u4_mb_y, u4_image_wd_mb;
    deblk_mb_t *ps_cur_mb = ps_dec->ps_cur_deblk_thrd_mb;
    deblk_mb_t *ps_top_mb;
    deblk_mb_t *ps_left_mb;
    const WORD32 i4_cb_qp_idx_ofst =
                    ps_dec->ps_cur_pps->i1_chroma_qp_index_offset;
    const WORD32 i4_cr_qp_idx_ofst =
                    ps_dec->ps_cur_pps->i1_second_chroma_qp_index_offset;

    UWORD32 u4_wd_y, u4_wd_uv;
    UWORD8 u1_field_pic_flag = ps_dec->ps_cur_slice->u1_field_pic_flag;

    u4_mb_num = ps_dec->u4_cur_deblk_mb_num;
    u4_mb_x = ps_dec->u4_deblk_mb_x;
    u4_mb_y = ps_dec->u4_deblk_mb_y;
    u4_image_wd_mb = ps_dec->u2_frm_wd_in_mbs;
    u4_wd_y = ps_dec->u2_frm_wd_y << u1_field_pic_flag;
    u4_wd_uv = ps_dec->u2_frm_wd_uv << u1_field_pic_flag;
    ps_cur_mb = ps_dec->ps_cur_deblk_thrd_mb;

    for(i = 0; i < deblk_mb_grp; i++)
    {

        while(1)
        {
            CHECK_MB_MAP_BYTE(u4_mb_num, mb_map, u4_cur_mb);

            if(ps_dec->u4_cur_bs_mb_num <= u4_mb_num)
                u4_cur_mb = 0;

            if(u4_mb_x < (u4_image_wd_mb - 1))
            {
                CHECK_MB_MAP_BYTE((u4_mb_num + 1), mb_map, u4_right_mb);
            }
            else
                u4_right_mb = 1;

            if(ps_dec->u2_mb_skip_error)
            {
                ps_dec->u2_skip_deblock = 1;
                break;
            }


            if(ps_dec->u2_skip_deblock == 1)
            {
                break;
            }
            if((u4_cur_mb && u4_right_mb) == 0)
            {

                if(ps_dec->u4_output_present
                                && ps_dec->u4_fmt_conv_cur_row
                                                < ps_dec->s_disp_frame_info.u4_y_ht)
                {
                    ps_dec->u4_fmt_conv_num_rows =
                                    MIN(ps_dec->u4_fmt_conv_num_rows,
                                        (ps_dec->s_disp_frame_info.u4_y_ht
                                                        - ps_dec->u4_fmt_conv_cur_row));
                    ih264d_format_convert(ps_dec, &(ps_dec->s_disp_op),
                                          ps_dec->u4_fmt_conv_cur_row,
                                          ps_dec->u4_fmt_conv_num_rows);
                    ps_dec->u4_fmt_conv_cur_row += ps_dec->u4_fmt_conv_num_rows;
                }
                else
                    NOP(32);
            }
            else
            {
                break;
            }
        }

        u4_mb_num++;
        {
            UWORD32 u4_deb_mode, u4_mbs_next;
            u4_deb_mode = ps_cur_mb->u1_deblocking_mode;
            if(!(u4_deb_mode & MB_DISABLE_FILTERING))
            {

                if(u4_mb_x)
                {
                    ps_left_mb = ps_cur_mb - 1;

                }
                else
                {
                    ps_left_mb = NULL;

                }
                if(u4_mb_y != 0)
                {
                    ps_top_mb = ps_cur_mb - (u4_image_wd_mb);
                }
                else
                {
                    ps_top_mb = NULL;
                }

                if(u4_deb_mode & MB_DISABLE_LEFT_EDGE)
                    ps_left_mb = NULL;
                if(u4_deb_mode & MB_DISABLE_TOP_EDGE)
                    ps_top_mb = NULL;

                ih264d_deblock_mb_nonmbaff(ps_dec, ps_tfr_cxt,
                                           i4_cb_qp_idx_ofst, i4_cr_qp_idx_ofst,
                                           ps_cur_mb, u4_wd_y, u4_wd_uv,
                                           ps_top_mb, ps_left_mb);
            }

            ps_cur_mb++;
            u4_mb_x++;
            u4_mbs_next = u4_image_wd_mb - u4_mb_x;

            ps_tfr_cxt->pu1_mb_y += 16;
            ps_tfr_cxt->pu1_mb_u += 8 * YUV420SP_FACTOR;
            ps_tfr_cxt->pu1_mb_v += 8;

            if(!u4_mbs_next)
            {
                ps_tfr_cxt->pu1_mb_y += ps_tfr_cxt->u4_y_inc;
                ps_tfr_cxt->pu1_mb_u += ps_tfr_cxt->u4_uv_inc;
                ps_tfr_cxt->pu1_mb_v += ps_tfr_cxt->u4_uv_inc;
                u4_mb_y++;
                u4_mb_x = 0;
            }
        }

    }

    ps_dec->u4_cur_deblk_mb_num = u4_mb_num;
    ps_dec->u4_deblk_mb_x = u4_mb_x;
    ps_dec->u4_deblk_mb_y = u4_mb_y;
    ps_dec->ps_cur_deblk_thrd_mb = ps_cur_mb;

}
void ih264d_computebs_deblk_slice(dec_struct_t *ps_dec, tfr_ctxt_t *ps_tfr_cxt)
{
    dec_mb_info_t *p_cur_mb;
    UWORD32 u4_max_addr = ps_dec->ps_cur_sps->u2_max_mb_addr;
    UWORD32 i;
    UWORD32 u1_mb_aff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    UWORD16 u2_slice_num;
    UWORD32 u4_mb_num;

    ps_dec->u4_cur_slice_bs_done = 0;
    ps_dec->u4_bs_cur_slice_num_mbs = 0;
    ps_dec->u4_cur_bs_mb_num =
                    (ps_dec->ps_computebs_cur_slice->u4_first_mb_in_slice)
                                    << u1_mb_aff;

    while(ps_dec->u4_cur_slice_bs_done != 1)
    {
        UWORD32 bs_mb_grp = BS_MB_GROUP;
        while(1)
        {

            UWORD32 u4_cond = 0;

            u4_mb_num = ps_dec->u4_cur_bs_mb_num;

            /*introducing 1 MB delay*/
            if((u4_mb_num + BS_MB_GROUP) <= u4_max_addr)
                u4_mb_num = u4_mb_num + BS_MB_GROUP;
            else
            {
                bs_mb_grp = u4_max_addr - u4_mb_num + 1;
                u4_mb_num = u4_max_addr;

            }

            CHECK_MB_MAP_BYTE(u4_mb_num, ps_dec->pu1_dec_mb_map, u4_cond);
            if(u4_cond)
            {
                break;
            }

            if(ps_dec->u2_skip_deblock == 0)
            {
                ih264d_check_mb_map_deblk(ps_dec, DEBLK_MB_GROUP, ps_tfr_cxt);
            }
        }

        GET_SLICE_NUM_MAP(ps_dec->pu2_slice_num_map, ps_dec->u4_cur_bs_mb_num,
                          u2_slice_num);

        if(u2_slice_num != ps_dec->u2_cur_slice_num_bs)
        {
            ps_dec->u4_cur_slice_bs_done = 1;
        }

        /* Compute BS for NMB group*/
        for(i = 0; i < bs_mb_grp; i++)
        {
            GET_SLICE_NUM_MAP(ps_dec->pu2_slice_num_map,
                              ps_dec->u4_cur_bs_mb_num, u2_slice_num);

            if(u2_slice_num != ps_dec->u2_cur_slice_num_bs)
            {
                ps_dec->u4_cur_slice_bs_done = 1;
            }

            if(ps_dec->u4_cur_slice_bs_done == 1)
                break;

            p_cur_mb = &ps_dec->ps_frm_mb_info[ps_dec->u4_cur_bs_mb_num
                            & PD_MB_BUF_SIZE_MOD];

            DEBUG_THREADS_PRINTF("ps_dec->u4_cur_bs_mb_num = %d\n",ps_dec->u4_cur_bs_mb_num);
            ih264d_compute_bs_non_mbaff_thread(ps_dec, p_cur_mb,
                                               ps_dec->u4_cur_bs_mb_num);
            ps_dec->u4_cur_bs_mb_num++;
            ps_dec->u4_bs_cur_slice_num_mbs++;

        }

        if(ps_dec->u4_cur_bs_mb_num > u4_max_addr)
        {
            ps_dec->u4_cur_slice_bs_done = 1;
        }

        /*deblock MB group*/
        {
            UWORD32 u4_num_mbs;

            if(ps_dec->u4_cur_bs_mb_num > ps_dec->u4_cur_deblk_mb_num)

                u4_num_mbs = ps_dec->u4_cur_bs_mb_num
                                - ps_dec->u4_cur_deblk_mb_num;
            else
                u4_num_mbs = 0;

            if(u4_num_mbs >= DEBLK_MB_GROUP)
                u4_num_mbs = DEBLK_MB_GROUP;
            if(ps_dec->u2_skip_deblock == 0)
            {
                ih264d_check_mb_map_deblk_wait(ps_dec, u4_num_mbs, ps_tfr_cxt);
            }
        }

    }
}

void ih264d_computebs_deblk_thread(dec_struct_t *ps_dec)
{
    tfr_ctxt_t s_tfr_ctxt;
    tfr_ctxt_t *ps_tfr_cxt = &s_tfr_ctxt; // = &ps_dec->s_tran_addrecon;
    pad_mgr_t *ps_pad_mgr = &ps_dec->s_pad_mgr;

    UWORD32 yield_cnt = 0;

    ithread_set_name("ih264d_computebs_deblk_thread");


    // run the loop till all slices are decoded

    // 0: un-identified state, 1 - bs needed, 2 - bs not needed
    while(1)
    {
        if(ps_dec->u4_start_bs_deblk == 0)
        {
            NOP(128);
            NOP(128);
            NOP(128);
            NOP(128);
        }
        else
        {
            break;
        }
    }

    if(ps_dec->u4_start_bs_deblk == 1)
    {
        ps_dec->u4_cur_deblk_mb_num = 0;
        ps_dec->u4_deblk_mb_x = 0;
        ps_dec->u4_deblk_mb_y = 0;

        ih264d_init_deblk_tfr_ctxt(ps_dec, ps_pad_mgr, ps_tfr_cxt,
                                   ps_dec->u2_frm_wd_in_mbs, 0);

        ps_tfr_cxt->pu1_mb_y = ps_tfr_cxt->pu1_src_y + 4;
        ps_tfr_cxt->pu1_mb_u = ps_tfr_cxt->pu1_src_u + 4;
        ps_tfr_cxt->pu1_mb_v = ps_tfr_cxt->pu1_src_v + 4;

        ps_dec->ps_cur_deblk_thrd_mb = ps_dec->ps_deblk_pic;

        while(1)
        {
            /*Complete all writes before processing next slice*/
            DATA_SYNC();
            /*wait untill all the slice params have been populated*/
            while(ps_dec->ps_computebs_cur_slice->slice_header_done == 0)
            {
                NOP(32); DEBUG_THREADS_PRINTF(" waiting for slice header at compute bs\n");
            }

            DEBUG_THREADS_PRINTF(" Entering compute bs slice\n");
            ih264d_computebs_deblk_slice(ps_dec, ps_tfr_cxt);

            DEBUG_THREADS_PRINTF(" Exit  compute bs slice \n");

            /*Complete all writes before processing next slice*/
            DATA_SYNC();

            while(1)
            {
                volatile void * parse_addr, *computebs_addr;
                volatile UWORD32 last_slice;

                parse_addr = (volatile void *)ps_dec->ps_parse_cur_slice;
                computebs_addr =
                                (volatile void *)ps_dec->ps_computebs_cur_slice;
                last_slice =
                                ps_dec->ps_computebs_cur_slice->last_slice_in_frame;

                if(last_slice == 1)
                    break;

                if(parse_addr != computebs_addr)
                    break;

                DEBUG_THREADS_PRINTF("Waiting at compute bs for next slice  or end of frame\n");

                NOP(32);

            }

            DEBUG_THREADS_PRINTF("CBS thread:Got next slice/end of frame signal \n ");

            if((void *)ps_dec->ps_parse_cur_slice
                            > (void *)ps_dec->ps_computebs_cur_slice)
            {
                ps_dec->ps_computebs_cur_slice++;
                ps_dec->u2_cur_slice_num_bs++;
            }
            else
            {
                /*Last slice in frame*/
                break;
            }

        }

        /*deblock remaining MBs*/
        {
            UWORD32 u4_num_mbs;

            u4_num_mbs = ps_dec->ps_cur_sps->u2_max_mb_addr
                            - ps_dec->u4_cur_deblk_mb_num + 1;

            DEBUG_PERF_PRINTF("mbs left for deblocking= %d \n",u4_num_mbs);

            if(u4_num_mbs != 0)
                if(ps_dec->u2_skip_deblock == 0)
                    ih264d_check_mb_map_deblk_wait(ps_dec, u4_num_mbs,
                                                   ps_tfr_cxt);
        }
    }

    ps_dec->u4_start_bs_deblk = 0;
    ithread_exit(0);
}


