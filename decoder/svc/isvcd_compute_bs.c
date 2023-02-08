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
 *  isvcd_compute_bs.c
 *
 * @brief
 *  This file contains bit flags and info on target layer.
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_deblk_extract_bit_flags()
 *  - isvcd_fill_bs_ibl()
 *  - isvcd_compute_bs_non_mbaff_target_lyr_no_inter_layer()
 *  - isvcd_compute_bs_non_mbaff()
 *  - isvcd_compute_bs_non_mbaff_target_lyr()
 *  - isvcd_compute_bs_non_mbaff_medial_lyr()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"
#include "ih264d_defs.h"
#include "ih264d_deblocking.h"
#include "string.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
UWORD16 ih264d_update_csbp_8x8(UWORD16 u2_luma_csbp);
void ih264d_fill_bs2_horz_vert(UWORD32 *pu4_bs,        /* Base pointer of BS table */
                               WORD32 u4_left_mb_csbp, /* csbp of left mb */
                               WORD32 u4_top_mb_csbp,  /* csbp of top mb */
                               WORD32 u4_cur_mb_csbp,  /* csbp of current mb */
                               const UWORD32 *pu4_packed_bs2, const UWORD16 *pu2_4x4_v2h_reorder);
const UWORD32 g_au4_extract_set[NUM_SUB_MB_PARTS] = {0x00000000F, 0x0000000F0, 0x000000F00,
                                                     0x00000F000};

/** \brief extracts the bits from given offset and packs it to last 4 bits  */
UWORD16 isvcd_deblk_extract_bit_flags(UWORD16 u2_bit_field, WORD32 i4_initial_bit_mask)
{
    WORD32 i4_i;
    WORD32 i4_bit_mask;
    UWORD16 u2_result = 0;

    i4_bit_mask = i4_initial_bit_mask;

    for(i4_i = 0; i4_i < NUM_SUB_MB_PARTS; i4_i++)
    {
        WORD32 i4_bit;
        /* extract the bits of the last column 4x4 blocks */
        if(0 == (i4_bit_mask & u2_bit_field))
        {
            i4_bit = 0;
        }
        else
        {
            i4_bit = 1;
        }
        /* store the result */
        u2_result |= i4_bit << i4_i;
        i4_bit_mask <<= 4;

    } /* end of loop over num sub Mb parts */
    return (u2_result);
}

/** \brief Fills the BS for edges falling on a IBL boundary  */
void isvcd_fill_bs_ibl(deblk_mb_t *ps_deblk_mb, UWORD8 u1_top_mb_type, UWORD8 u1_left_mb_type,
                       dec_mb_info_t *ps_cur_mb_info, UWORD16 *pu2_curr_res_luma_csbp,
                       UWORD16 *pu2_left_res_luma_csbp, UWORD16 *pu2_top_res_luma_csbp)
{
    /*! Flow of the module is as follows                                  */
    /*! 1. checks if MB edge is falling on IBL boundary                   */
    /*! 2. if only Mb edge then it fills the BS based on INTRA or INTER
           stauts                                                         */
    /*! 3. if the current MB is IBL and neighbours are also neighbours
           then it uses the current layer t_coeff flag to decide the
           BS of a particular edge                                        */
    /*!4. fills the BS for all the edges in curretn MB if IBL             */

    UWORD16 u2_top_horz_nnz;
    UWORD8 u1_top_mb_ibl, u1_left_mb_ibl;
    UWORD32 i4_i, i4_edge;
    UWORD8 u1_bs;
    UWORD8 u1_cnd;
    UWORD8 u1_top_intra;
    UWORD8 u1_left_intra;
    UWORD8 u1_p_nnz, u1_q_nnz;
    UWORD8 u1_curr_mb_ibl;
    UWORD32 *pu4_bs_table;
    UWORD16 u2_curr_nnz;
    UWORD8 u1_left_mb_nnz = 0, u1_left_nnz;
    WORD32 i4_horz_start = 0;
    WORD32 i4_vertical_start = 0;

    pu4_bs_table = &(ps_deblk_mb->u4_bs_table[0]);

    u1_top_mb_ibl = u1_top_mb_type & D_INTRA_IBL;
    u1_left_mb_ibl = u1_left_mb_type & D_INTRA_IBL;

    u1_curr_mb_ibl = ps_deblk_mb->u1_mb_type & D_INTRA_IBL;

    u1_top_intra = u1_top_mb_type & D_INTRA_MB;
    u1_left_intra = u1_left_mb_type & D_INTRA_MB;

    /* return if none of the current top and left is IBL */
    if((0 == u1_curr_mb_ibl) && (0 == u1_top_mb_ibl) && (0 == u1_left_mb_ibl))
    {
        return;
    }

    /* set up the vertical and horz MB edge skip flags */
    if(0 != u1_curr_mb_ibl)
    {
        /* if top is not IBL */
        if(0 == u1_top_mb_ibl)
        {
            i4_horz_start = 1;
        }

        /* if left in not IBL */
        if(0 == u1_left_mb_ibl)
        {
            i4_vertical_start = 1;
        }
    }

    /*******************************************************/
    /* Fill BS for mb egdex assuming non IBL case          */
    /*******************************************************/

    /* only the  MB edges fall across IBL boundary */
    if((0 != u1_curr_mb_ibl) || (0 != u1_top_mb_ibl) || (0 != u1_left_mb_ibl))
    {
        UWORD16 u2_temp, u2_i;
        u2_temp = *pu2_left_res_luma_csbp;
        for(u2_i = 0; u2_i < 4; u2_i++)
        {
            u1_left_mb_nnz |= ((u2_temp & 0x08) >> (3 - u2_i));
            u2_temp >>= 4;
        }
        u2_curr_nnz = *pu2_curr_res_luma_csbp;
        u2_top_horz_nnz = *pu2_top_res_luma_csbp >> 12;

        /* top is intra and not ibl */
        if(0 != u1_top_intra)
        {
            pu4_bs_table[0] = 0x04040404;
        }
        /* left is intra and not ibl */
        if(0 != u1_left_intra)
        {
            pu4_bs_table[4] = 0x04040404;
        }

        /* assume neighbours are inter and update bs */

        /* Edge = 0 means Vert Edges and Edge = 1 means Horz edges */
        for(i4_edge = 0; i4_edge < 2; i4_edge++)
        {
            UWORD8 u1_p_nnz, u1_q_nnz;
            UWORD32 u4_bs_edge = 0;
            WORD32 i4_bit_mask;
            WORD32 i4_curr_intra_flag;
            WORD32 i4_neibor_intra_flag;

            i4_curr_intra_flag = (0 != u1_curr_mb_ibl);

            if(0 != i4_edge)
            {
                /* Initialize for the TOP edge */
                u1_p_nnz = (UWORD8) u2_top_horz_nnz;
                u1_q_nnz = (UWORD8) (u2_curr_nnz & g_au4_extract_set[0]);
                i4_neibor_intra_flag = (u1_top_mb_ibl || u1_top_intra);
            }
            else
            {
                u1_p_nnz = u1_left_mb_nnz;
                u1_q_nnz = (UWORD8) isvcd_deblk_extract_bit_flags(u2_curr_nnz, 0x01);
                i4_neibor_intra_flag = (u1_left_mb_ibl || u1_left_intra);
            }

            i4_bit_mask = 1;
            /* find bs of 4 edges */
            for(i4_i = 0; i4_i < 4; i4_i++)
            {
                UWORD8 u1_p_nnz_temp, u1_q_nnz_temp;

                u1_p_nnz_temp = (u1_p_nnz & i4_bit_mask);
                u1_q_nnz_temp = (u1_q_nnz & i4_bit_mask);

                u1_cnd = ((u1_p_nnz_temp && (!i4_neibor_intra_flag)) ||
                          (u1_q_nnz_temp && (!i4_curr_intra_flag)));

                u1_bs = u1_cnd ? 2 : 1;

                /* update the bs of the edge */
                u4_bs_edge = (u4_bs_edge << 8) + u1_bs;
                i4_bit_mask <<= 1;

            } /* end of loop over blk edges */

            /* update the bs of edges */
            if(i4_edge && !u1_top_intra)
            {
                pu4_bs_table[0] = u4_bs_edge;
            }
            else if(!i4_edge && !u1_left_intra)
            {
                pu4_bs_table[4] = u4_bs_edge;
            }
        } /* end of loop over v1 vetical and horizontal edge */
    }
    /* current MB is IBL */
    if(0 != u1_curr_mb_ibl)
    {
        UWORD16 u2_temp, u2_i;
        WORD32 i4_bit_mask_edge = 1;

        u1_left_mb_nnz = 0;
        u2_temp = ps_cur_mb_info->ps_left_mb->u2_luma_csbp;
        for(u2_i = 0; u2_i < 4; u2_i++)
        {
            u1_left_mb_nnz |= ((u2_temp & 0x08) >> (3 - u2_i));
            u2_temp >>= 4;
        }
        u2_curr_nnz = ps_cur_mb_info->ps_curmb->u2_luma_csbp;
        u2_top_horz_nnz = ps_cur_mb_info->ps_top_mb->u2_luma_csbp >> 12;
        /* all are IBL edges then use only t_coeff of current layer*/
        /* loop over all edges */
        for(i4_edge = 0; i4_edge < 4; i4_edge++)
        {
            UWORD16 u2_curr_horz_nnz = 0;
            WORD32 i4_bit_mask = 1;

            u2_curr_horz_nnz = u2_curr_nnz & g_au4_extract_set[i4_edge];

            u2_curr_horz_nnz = (u2_curr_horz_nnz >> (i4_edge * 4));

            u1_left_nnz = (u1_left_mb_nnz & i4_bit_mask_edge);

            for(i4_i = 0; i4_i < 4; i4_i++)
            {
                UWORD8 u1_curr_nnz, u1_top_nnz;

                u1_curr_nnz = (u2_curr_horz_nnz & i4_bit_mask);
                u1_top_nnz = (u2_top_horz_nnz & i4_bit_mask);
                /* update bs horizontal */

                if(!((1 == i4_horz_start) && (0 == i4_edge)))
                {
                    u1_p_nnz = u1_top_nnz;
                    u1_q_nnz = u1_curr_nnz;
                    u1_cnd = !(u1_p_nnz || u1_q_nnz);
                    u1_bs = u1_cnd ? 0 : 1;
                    pu4_bs_table[i4_edge] = (pu4_bs_table[i4_edge] << 8) + u1_bs;
                }

                /* update bs vertical */
                if(!((1 == i4_vertical_start) && (0 == i4_i)))
                {
                    u1_p_nnz = u1_left_nnz;
                    u1_q_nnz = u1_curr_nnz;
                    u1_cnd = !(u1_p_nnz || u1_q_nnz);
                    u1_bs = u1_cnd ? 0 : 1;
                    pu4_bs_table[i4_i + 4] = (pu4_bs_table[i4_i + 4] << 8) + u1_bs;
                }
                /* store the current nnz to left nnz */
                u1_left_nnz = u1_curr_nnz;
                i4_bit_mask <<= 1;
            }
            /* store the current row nnz to top row nnz */
            u2_top_horz_nnz = u2_curr_horz_nnz;
            i4_bit_mask_edge <<= 1;
        }
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compute_bs_non_mbaff_target_lyr_no_inter_layer     */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function */
/*                : for Boundary Strength Calculation                        */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Processing    : This functions calls deblock MB in the MB increment order*/
/*                                                                           */
/*  Outputs       : Produces the Boundary Strength for Current Mb            */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore                                              */
/*****************************************************************************/

void isvcd_compute_bs_non_mbaff_target_lyr_no_inter_layer(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                                          dec_mb_info_t *ps_cur_mb_info,
                                                          const UWORD16 u2_mbxn_mb)
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
    UWORD32 u1_cur_mb_type;
    UWORD32 *pu4_bs_table;

    /* Neighbour availability */
    /* Initialization */
    const UWORD32 u2_mbx = ps_cur_mb_info->u2_mbx;
    const UWORD32 u2_mby = ps_cur_mb_info->u2_mby;
    const UWORD32 u1_pingpong = u2_mbx & 0x01;

    PROFILE_DISABLE_BOUNDARY_STRENGTH()

    ps_deblk_top_mb = ps_dec->ps_deblk_top_mb + u2_mbx;

    /* Pointer assignment for Current DeblkMB, Current Mv Pred  */
    ps_cur_mb_params = ps_dec->ps_deblk_mbn + u2_mbxn_mb;
    ps_cur_mv_pred = ps_dec->ps_mv_cur + (u2_mbxn_mb << 4);

    apv_map_ref_idx_to_poc = ps_dec->ppv_map_ref_idx_to_poc + 1;
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

    /* Update the Left deblk_mb_t and Left MvPred Parameters           */
    if(!u2_mbx)
    {
        u4_leftmbtype = 0;

        /* Initialize the ps_left_mv_pred with Junk but Valid Location */
        /* to avoid invalid memory access                           */
        /* this is read only pointer                                */
        ps_left_mv_pred = ps_dec->ps_mv_cur + 3;
    }
    else
    {
        u4_leftmbtype = ps_dec->deblk_left_mb[1].u1_mb_type;

        /* Come to Left Most Edge of the MB */
        ps_left_mv_pred =
            (u2_mbxn_mb) ? ps_dec->ps_mv_cur + ((u2_mbxn_mb - 1) << 4) + 3 : ps_dec->ps_mv_left + 3;
    }

    if(!u2_mby) u1_top_mb_typ = 0;

    /* MvPred Pointer Calculation */
    ps_top_mv_pred = ps_cur_mv_pred - (ps_dec->u2_frm_wd_in_mbs << 4) + 12;

    u4_cur_mb_intra = u1_cur_mb_type & D_INTRA_MB;
    u4_cur_mb_fld = !!(u1_cur_mb_type & D_FLD_MB);
    /* Compute BS function */
    pu4_bs_table = ps_cur_mb_params->u4_bs_table;

    u2_cur_csbp = ps_cur_mb_info->ps_curmb->u2_luma_csbp;
    u2_left_csbp = ps_cur_mb_info->ps_left_mb->u2_luma_csbp;
    u2_top_csbp = ps_cur_mb_info->ps_top_mb->u2_luma_csbp;
    /* Compute BS function */
    if((ps_dec->ps_cur_sps->u1_profile_idc == HIGH_PROFILE_IDC) ||
       (ps_dec->ps_cur_sps->u1_profile_idc == SCALABLE_HIGH_PROFILE_IDC) ||
       (ps_dec->ps_cur_sps->u1_profile_idc == SCALABLE_BASELINE_PROFILE_IDC))
    {
        if(ps_cur_mb_info->u1_tran_form8x8 == 1)
        {
            u2_cur_csbp = ih264d_update_csbp_8x8(ps_cur_mb_info->ps_curmb->u2_luma_csbp);
        }

        if(ps_cur_mb_info->ps_left_mb->u1_tran_form8x8 == 1)
        {
            u2_left_csbp = ih264d_update_csbp_8x8(ps_cur_mb_info->ps_left_mb->u2_luma_csbp);
        }

        if(ps_cur_mb_info->ps_top_mb->u1_tran_form8x8 == 1)
        {
            u2_top_csbp = ih264d_update_csbp_8x8(ps_cur_mb_info->ps_top_mb->u2_luma_csbp);
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
        UWORD32 u4_is_b = ps_dec->u1_B;

        ih264d_fill_bs2_horz_vert(pu4_bs_table, u2_left_csbp, u2_top_csbp, u2_cur_csbp,
                                  (const UWORD32 *) (gau4_ih264d_packed_bs2),
                                  (const UWORD16 *) (gau2_ih264d_4x4_v2h_reorder));

        if(u4_leftmbtype & D_INTRA_MB) pu4_bs_table[4] = 0x04040404;

        if(u1_top_mb_typ & D_INTRA_MB) pu4_bs_table[0] = u4_cur_mb_fld ? 0x03030303 : 0x04040404;

        ps_dec->pf_fill_bs1[u4_is_b][u4_is_non16x16](
            ps_cur_mv_pred, ps_top_mv_pred, apv_map_ref_idx_to_poc, pu4_bs_table, ps_left_mv_pred,
            &(ps_dec->ps_left_mvpred_addr[u1_pingpong][1]),
            ps_cur_mb_info->ps_top_mb->u4_pic_addrress, (4 >> u4_cur_mb_fld));
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
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compute_bs_non_mbaff                               */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function */
/*                : for Boundary Strength Calculation                        */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Processing    : This functions calls deblock MB in the MB increment order*/
/*                                                                           */
/*  Outputs       : Produces the Boundary Strength for Current Mb            */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore                                              */
/*****************************************************************************/

void isvcd_compute_bs_non_mbaff(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                                const UWORD16 u2_mbxn_mb)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    ps_dec->pf_compute_bs(ps_dec, ps_cur_mb_info, u2_mbxn_mb);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compute_bs_non_mbaff_target_lyr                    */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function */
/*                : for Boundary Strength Calculation                        */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Processing    : This functions calls deblock MB in the MB increment order*/
/*                                                                           */
/*  Outputs       : Produces the Boundary Strength for Current Mb            */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore                                              */
/*****************************************************************************/

void isvcd_compute_bs_non_mbaff_target_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info, const UWORD16 u2_mbxn_mb)
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
    ps_cur_mb_params = ps_dec->ps_deblk_mbn + u2_mbxn_mb;
    ps_cur_mv_pred = ps_dec->ps_mv_cur + (u2_mbxn_mb << 4);

    /*Pointer assignment for Residual NNZ */
    pu2_curr_res_luma_csbp = ps_svc_lyr_dec->pu2_frm_res_luma_csbp + ps_cur_mb_info->u2_mbx;
    pu2_curr_res_luma_csbp += ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride;

    pu2_left_res_luma_csbp = pu2_curr_res_luma_csbp - (ps_cur_mb_info->u2_mbx != 0);
    pu2_top_res_luma_csbp = pu2_curr_res_luma_csbp - ((ps_cur_mb_info->u2_mby != 0) *
                                                      ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride);

    apv_map_ref_idx_to_poc = ps_dec->ppv_map_ref_idx_to_poc + 1;
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
        ps_left_mv_pred = ps_dec->ps_mv_cur + 3;
    }
    else
    {
        u4_leftmbtype = ps_dec->deblk_left_mb[1].u1_mb_type;

        /* Come to Left Most Edge of the MB */
        ps_left_mv_pred =
            (u2_mbxn_mb) ? ps_dec->ps_mv_cur + ((u2_mbxn_mb - 1) << 4) + 3 : ps_dec->ps_mv_left + 3;
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
            UWORD32 u4_is_b = ps_dec->u1_B;
            UWORD32 u4_bs_0, u4_bs_4;

            u4_bs_0 = pu4_bs_table[0];
            u4_bs_4 = pu4_bs_table[4];

            ih264d_fill_bs2_horz_vert(pu4_bs_table, u2_left_csbp, u2_top_csbp, u2_cur_csbp,
                                      (const UWORD32 *) (gau4_ih264d_packed_bs2),
                                      (const UWORD16 *) (gau2_ih264d_4x4_v2h_reorder));

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

    /* For transform 8x8 disable deblocking of the intrernal edges of a 8x8 block */
    if(ps_cur_mb_info->u1_tran_form8x8)
    {
        pu4_bs_table[1] = 0;
        pu4_bs_table[3] = 0;
        pu4_bs_table[5] = 0;
        pu4_bs_table[7] = 0;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compute_bs_non_mbaff_medial_lyr                    */
/*                                                                           */
/*  Description   : This function computes the pointers of left,top & current*/
/*                : Nnz, MvPred & deblk_mb_t and supplies to FillBs function */
/*                : for Boundary Strength Calculation                        */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Processing    : This functions calls deblock MB in the MB increment order*/
/*                                                                           */
/*  Outputs       : Produces the Boundary Strength for Current Mb            */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore                                              */
/*****************************************************************************/

void isvcd_compute_bs_non_mbaff_medial_lyr(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           dec_mb_info_t *ps_cur_mb_info, const UWORD16 u2_mbxn_mb)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    /* deblk_mb_t Params */
    deblk_mb_t *ps_cur_mb_params; /*< Parameters of current MacroBlock */
    deblkmb_neighbour_t *ps_deblk_top_mb;
    UWORD32 u4_leftmbtype;
    UWORD16 u2_cur_csbp;

    /* Set of flags */
    UWORD32 u4_cur_mb_intra, u1_top_mb_typ, u4_cur_mb_fld;
    UWORD32 u1_cur_mb_type;
    UWORD32 *pu4_bs_table;
    UWORD32 mb_type_intra = 0;

    UWORD16 *pu2_curr_res_luma_csbp;
    UWORD16 *pu2_left_res_luma_csbp;
    UWORD16 *pu2_top_res_luma_csbp;

    /* Neighbour availability */
    const UWORD32 u2_mbx = ps_cur_mb_info->u2_mbx;
    const UWORD32 u2_mby = ps_cur_mb_info->u2_mby;

    PROFILE_DISABLE_BOUNDARY_STRENGTH()

    ps_deblk_top_mb = ps_dec->ps_deblk_top_mb + u2_mbx;

    /* Pointer assignment for Current DeblkMB, Current Mv Pred  */
    ps_cur_mb_params = ps_dec->ps_deblk_mbn + u2_mbxn_mb;

    /*Pointer assignment for Residual NNZ */
    pu2_curr_res_luma_csbp = ps_svc_lyr_dec->pu2_frm_res_luma_csbp + ps_cur_mb_info->u2_mbx;
    pu2_curr_res_luma_csbp += ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride;

    pu2_left_res_luma_csbp = pu2_curr_res_luma_csbp - (ps_cur_mb_info->u2_mbx != 0);
    pu2_top_res_luma_csbp = pu2_curr_res_luma_csbp - ((ps_cur_mb_info->u2_mby != 0) *
                                                      ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride);

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
        /* Storing the leftMbtype for next Mb */
        ps_dec->deblk_left_mb[1].u1_mb_type = ps_cur_mb_params->u1_mb_type;
        return;
    }

    /* Flag for extra left Edge */
    ps_cur_mb_params->u1_single_call = 1;

    /* Update the Left deblk_mb_t */
    if(!u2_mbx)
    {
        u4_leftmbtype = 0;
    }
    else
    {
        u4_leftmbtype = ps_dec->deblk_left_mb[1].u1_mb_type;
    }

    if(!u2_mby) u1_top_mb_typ = 0;

    u4_cur_mb_intra = u1_cur_mb_type & D_INTRA_MB;
    u4_cur_mb_fld = !!(u1_cur_mb_type & D_FLD_MB);
    /* Compute BS function */
    pu4_bs_table = ps_cur_mb_params->u4_bs_table;

    u2_cur_csbp = ps_cur_mb_info->ps_curmb->u2_luma_csbp;
    /* Compute BS function */
    if((ps_dec->ps_cur_sps->u1_profile_idc == HIGH_PROFILE_IDC) ||
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
    }

    mb_type_intra = (u1_top_mb_typ & D_INTRA_MB) || (u1_top_mb_typ & D_INTRA_IBL);

    /* if Top MB or current MB is INTER */
    if(!mb_type_intra)
    {
        pu4_bs_table[0] = 0;
        /* disable the processing of top edge  */
        ps_cur_mb_params->u1_deblocking_mode |= MB_DISABLE_TOP_EDGE;
    }

    mb_type_intra = (u4_leftmbtype & D_INTRA_MB) || (u4_leftmbtype & D_INTRA_IBL);
    /* if Left MB current MB is INTER */
    if(!mb_type_intra)
    {
        pu4_bs_table[4] = 0;
        /* disable the processing of left edge  */
        ps_cur_mb_params->u1_deblocking_mode |= MB_DISABLE_LEFT_EDGE;
    }

    /* overwrite the BS 0 values for corner cases */
    if(0 == u2_mbx)
    {
        pu4_bs_table[4] = 0;
    }
    if(0 == u2_mby)
    {
        pu4_bs_table[0] = 0;
    }

    /* Storing the leftMbtype for next Mb */
    ps_dec->deblk_left_mb[1].u1_mb_type = ps_cur_mb_params->u1_mb_type;

    /* For transform 8x8 disable deblocking of the intrernal edges of a 8x8 block */
    if(ps_cur_mb_info->u1_tran_form8x8)
    {
        pu4_bs_table[1] = 0;
        pu4_bs_table[3] = 0;
        pu4_bs_table[5] = 0;
        pu4_bs_table[7] = 0;
    }
}
