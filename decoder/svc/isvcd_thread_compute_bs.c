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
 *  isvcd_thread_compute_bs.c
 *
 * @brief
 *  Contains routines that for multi-thread decoder
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
#include "ih264d_error_handler.h"
#include "ih264d_debug.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "isvcd_structs.h"
#include "ih264d_defs.h"
#include "ih264d_mb_utils.h"
#include "ih264d_thread_compute_bs.h"
#include "isvcd_thread_compute_bs.h"
#include "ithread.h"
#include "ih264d_deblocking.h"
#include "ih264d_process_pslice.h"
#include "isvcd_process_epslice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264d_mb_utils.h"
#include "ih264d_tables.h"
#include "ih264d_format_conv.h"
#include "ih264d_defs.h"

UWORD16 ih264d_update_csbp_8x8(UWORD16 u2_luma_csbp);
void ih264d_fill_bs2_horz_vert(UWORD32 *pu4_bs,        /* Base pointer of BS table */
                               WORD32 u4_left_mb_csbp, /* csbp of left mb */
                               WORD32 u4_top_mb_csbp,  /* csbp of top mb */
                               WORD32 u4_cur_mb_csbp,  /* csbp of current mb */
                               const UWORD32 *pu4_packed_bs2, const UWORD16 *pu2_4x4_v2h_reorder);
void isvcd_fill_bs_ibl(deblk_mb_t *ps_deblk_mb, UWORD8 u1_top_mb_type, UWORD8 u1_left_mb_type,
                       dec_mb_info_t *ps_cur_mb_info, UWORD16 *pu2_curr_res_luma_csbp,
                       UWORD16 *pu2_left_res_luma_csbp, UWORD16 *pu2_top_res_luma_csbp);

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compute_bs_non_mbaff_thread */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function
 * for*/
/*                : Boundary Strength Calculation .this function is used     */
/*                : BS being calculated in separate thread                   */
/*  Inputs        : pointer to decoder context,cur_mb_info,u4_mb_num */
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

void isvcd_compute_bs_non_mbaff_thread(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                       dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mb_num)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    /* Mvpred and Nnz for top and Courrent */
    mv_pred_t *ps_cur_mv_pred, *ps_top_mv_pred = NULL, *ps_left_mv_pred;
    /* deblk_mb_t Params */
    deblk_mb_t *ps_cur_mb_params; /*< Parameters of current MacroBlock */
    deblkmb_neighbour_t *ps_deblk_top_mb;

    /* Reference Index to POC mapping*/
    void **apv_map_ref_idx_to_poc;
    UWORD32 u4_leftmbtype;

    UWORD16 u2_left_csbp, u2_top_csbp, u2_cur_csbp;

    /* Set of flags */
    UWORD32 u4_cur_mb_intra, u1_top_mb_typ, u4_cur_mb_fld;
    UWORD32 u4_cur_mb_ibl;
    UWORD32 u1_cur_mb_type;
    UWORD32 *pu4_bs_table;

    UWORD16 *pu2_curr_res_luma_csbp;
    UWORD16 *pu2_left_res_luma_csbp;
    UWORD16 *pu2_top_res_luma_csbp;

    /* Neighbour availability */
    /* Initialization */
    const UWORD32 u2_mbx = ps_cur_mb_info->u2_mbx;
    const UWORD32 u2_mby = ps_cur_mb_info->u2_mby;
    const UWORD32 u1_pingpong = u2_mbx & 0x01;

    PROFILE_DISABLE_BOUNDARY_STRENGTH()
    ps_deblk_top_mb = ps_dec->ps_deblk_top_mb + u2_mbx;

    /* Pointer assignment for Current DeblkMB, Current Mv Pred  */
    ps_cur_mb_params = ps_dec->ps_deblk_pic + u4_mb_num;
    ps_cur_mv_pred = ps_dec->s_cur_pic.ps_mv + (u4_mb_num << 4);

    /*Pointer assignment for Residual NNZ */
    pu2_curr_res_luma_csbp = ps_svc_lyr_dec->pu2_frm_res_luma_csbp + ps_cur_mb_info->u2_mbx;
    pu2_curr_res_luma_csbp += ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride;

    pu2_left_res_luma_csbp = pu2_curr_res_luma_csbp - (ps_cur_mb_info->u2_mbx != 0);
    pu2_top_res_luma_csbp = pu2_curr_res_luma_csbp - ((ps_cur_mb_info->u2_mby != 0) *
                                                      ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride);

    apv_map_ref_idx_to_poc = (void **) ps_dec->ps_computebs_cur_slice->ppv_map_ref_idx_to_poc + 1;
    u1_cur_mb_type = ps_cur_mb_params->u1_mb_type;
    u1_top_mb_typ = ps_deblk_top_mb->u1_mb_type;
    ps_deblk_top_mb->u1_mb_type = u1_cur_mb_type;

    ps_cur_mb_params->u1_topmb_qp = ps_deblk_top_mb->u1_mb_qp;
    ps_deblk_top_mb->u1_mb_qp = ps_cur_mb_params->u1_mb_qp;
    ps_cur_mb_params->u1_left_mb_qp = ps_dec->deblk_left_mb[1].u1_mb_qp;
    ps_dec->deblk_left_mb[1].u1_mb_qp = ps_cur_mb_params->u1_mb_qp;

    /* if no deblocking required for current Mb then continue */
    /* Check next Mbs   in Mb group                           */
    if(ps_cur_mb_params->u1_deblocking_mode & MB_DISABLE_FILTERING)
    {
        void **pu4_map_ref_idx_to_poc_l1 = apv_map_ref_idx_to_poc + POC_LIST_L0_TO_L1_DIFF;
        {
            /* Store Parameter for Top MvPred refernce frame Address */
            void **ppv_top_mv_pred_addr = ps_cur_mb_info->ps_curmb->u4_pic_addrress;
            WORD8 *p1_refTop0 = (ps_cur_mv_pred + 12)->i1_ref_frame;
            WORD8 *p1_refTop1 = (ps_cur_mv_pred + 14)->i1_ref_frame;

            /* Store Left addresses for Next Mb   */
            void **ppv_left_mv_pred_addr = ps_dec->ps_left_mvpred_addr[!u1_pingpong][1].u4_add;
            WORD8 *p1_refleft0 = (ps_cur_mv_pred + 3)->i1_ref_frame;

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

        return;
    }

    /* Flag for extra left Edge */
    ps_cur_mb_params->u1_single_call = 1;

    /* Update the Left deblk_mb_t and Left MvPred Parameters */
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

    if(!u2_mby) u1_top_mb_typ = 0;

    /* MvPred Pointer Calculation */
    ps_top_mv_pred = ps_cur_mv_pred - (ps_dec->u2_frm_wd_in_mbs << 4) + 12;
    u4_cur_mb_intra = u1_cur_mb_type & D_INTRA_MB;
    u4_cur_mb_ibl = u1_cur_mb_type & D_INTRA_IBL;
    u4_cur_mb_fld = !!(u1_cur_mb_type & D_FLD_MB);
    /* Compute BS function */
    pu4_bs_table = ps_cur_mb_params->u4_bs_table;

    u2_cur_csbp = ps_cur_mb_info->ps_curmb->u2_luma_csbp;
    u2_left_csbp = ps_cur_mb_info->ps_left_mb->u2_luma_csbp;
    u2_top_csbp = ps_cur_mb_info->ps_top_mb->u2_luma_csbp;

    /* Compute BS function */
    if((ps_dec->ps_cur_sps->u1_profile_idc == HIGH_PROFILE_IDC) ||
       (ps_dec->ps_cur_sps->u1_profile_idc == HIGH_PROFILE_IDC) ||
       (ps_dec->ps_cur_sps->u1_profile_idc == SCALABLE_HIGH_PROFILE_IDC) ||
       (ps_dec->ps_cur_sps->u1_profile_idc == SCALABLE_BASELINE_PROFILE_IDC))
    {
        if(ps_cur_mb_info->u1_tran_form8x8 == 1)
        {
            u2_cur_csbp = ih264d_update_csbp_8x8(ps_cur_mb_info->ps_curmb->u2_luma_csbp);
            ps_cur_mb_info->ps_curmb->u2_luma_csbp = u2_cur_csbp;
        }
    }
    u2_cur_csbp |= *pu2_curr_res_luma_csbp;
    u2_left_csbp |= *pu2_left_res_luma_csbp;
    u2_top_csbp |= *pu2_top_res_luma_csbp;

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
        isvcd_fill_bs_ibl(ps_cur_mb_params, u1_top_mb_typ, u4_leftmbtype, ps_cur_mb_info,
                          pu2_curr_res_luma_csbp, pu2_left_res_luma_csbp, pu2_top_res_luma_csbp);

        if(!u4_cur_mb_ibl)
        {
            UWORD32 u4_is_non16x16 = !!(u1_cur_mb_type & D_PRED_NON_16x16);
            UWORD32 u4_is_b = (ps_dec->ps_computebs_cur_slice->slice_type == B_SLICE);
            UWORD32 u4_bs_0, u4_bs_4;

            u4_bs_0 = pu4_bs_table[0];
            u4_bs_4 = pu4_bs_table[4];

            ih264d_fill_bs2_horz_vert(pu4_bs_table, u2_left_csbp, u2_top_csbp, u2_cur_csbp,
                                      gau4_ih264d_packed_bs2, gau2_ih264d_4x4_v2h_reorder);

            if(u4_leftmbtype & D_INTRA_MB)
            {
                pu4_bs_table[4] = 0x04040404;
            }
            else if(u4_leftmbtype & D_INTRA_IBL)
            {
                pu4_bs_table[4] = u4_bs_4;
            }

            if(u1_top_mb_typ & D_INTRA_MB)
            {
                pu4_bs_table[0] = u4_cur_mb_fld ? 0x03030303 : 0x04040404;
            }
            else if(u1_top_mb_typ & D_INTRA_IBL)
            {
                pu4_bs_table[0] = u4_bs_0;
            }

            ps_dec->pf_fill_bs1[u4_is_b][u4_is_non16x16](
                ps_cur_mv_pred, ps_top_mv_pred, apv_map_ref_idx_to_poc, pu4_bs_table,
                ps_left_mv_pred, &(ps_dec->ps_left_mvpred_addr[u1_pingpong][1]),
                ps_cur_mb_info->ps_top_mb->u4_pic_addrress, (4 >> u4_cur_mb_fld));
        }
    }

    {
        void **pu4_map_ref_idx_to_poc_l1 = apv_map_ref_idx_to_poc + POC_LIST_L0_TO_L1_DIFF;
        {
            /* Store Parameter for Top MvPred refernce frame Address */
            void **ppv_top_mv_pred_addr = ps_cur_mb_info->ps_curmb->u4_pic_addrress;
            WORD8 *p1_refTop0 = (ps_cur_mv_pred + 12)->i1_ref_frame;
            WORD8 *p1_refTop1 = (ps_cur_mv_pred + 14)->i1_ref_frame;

            /* Store Left addresses for Next Mb   */
            void **ppv_left_mv_pred_addr = ps_dec->ps_left_mvpred_addr[!u1_pingpong][1].u4_add;
            WORD8 *p1_refleft0 = (ps_cur_mv_pred + 3)->i1_ref_frame;

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
