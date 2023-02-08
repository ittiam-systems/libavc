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
 *  isvcd_mb_utils.c
 *
 * @brief
 *  Contains utitlity functions needed for Macroblock decoding
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_get_mb_info_cabac_nonmbaff()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include <stdlib.h>
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "isvcd_structs.h"
#include "ih264d_defs.h"
#include "ih264d_mb_utils.h"
#include "ih264d_parse_slice.h"
#include "ih264d_error_handler.h"
#include "ih264d_parse_mb_header.h"
#include "ih264d_cabac.h"
#include "ih264d_defs.h"
#include "ih264d_tables.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : get_mb_info_cabac                                        */
/*                                                                           */
/*  Description   : This function sets the following information of cur MB   */
/*                  (a) mb_x and mb_y                                        */
/*                  (b) Neighbour availablity                                */
/*                  (c) Macroblock location in the frame buffer              */
/*                  (e) leftMb parama and TopMb params of curMB              */
/*                  (f) For Mbaff case leftMb params and TopMb params of     */
/*                      bottomMb are also set if curMB is top                */
/*                  (g) For mbaff predicts field/frame u4_flag for topMb     */
/*                      and sets the field/frame for botMb. This is          */
/*                      written in ps_dec->u1_cur_mb_fld_dec_flag            */
/*                                                                           */
/*  Inputs        : pointer to decstruct                                     */
/*                  pointer to current mb info                               */
/*                  currentMbaddress                                         */
/*                                                                           */
/*  Processing    : leftMb and TopMb params are used by DecMbskip and        */
/*                  DecCtxMbfield  modules so that these modules do not      */
/*                  check for neigbour availability and then find the        */
/*                  neigbours for context increments                         */
/*                                                                           */
/*  Returns       : OK                                                       */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
UWORD32 isvcd_get_mb_info_cabac_nonmbaff(dec_struct_t *ps_dec, const UWORD16 u2_cur_mb_address,
                                         dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mbskip)
{
    WORD32 mb_x;
    WORD32 mb_y;
    UWORD32 u1_mb_ngbr_avail = 0;
    UWORD32 u2_frm_width_in_mb = ps_dec->u2_frm_wd_in_mbs;
    UWORD32 u1_top_mb = 1;
    WORD32 i2_prev_slice_mbx = ps_dec->i2_prev_slice_mbx;
    UWORD32 u2_top_right_mask = TOP_RIGHT_DEFAULT_AVAILABLE;
    UWORD32 u2_top_left_mask = TOP_LEFT_DEFAULT_AVAILABLE;
    ctxt_inc_mb_info_t *const p_ctx_inc_mb_map = ps_dec->p_ctxt_inc_mb_map;

    /*--------------------------------------------------------------------*/
    /* Calculate values of mb_x and mb_y                                  */
    /*--------------------------------------------------------------------*/
    mb_x = (WORD16) ps_dec->u2_mbx;
    mb_y = (WORD16) ps_dec->u2_mby;
    ps_dec->u2_cur_mb_addr = u2_cur_mb_address;

    mb_x++;
    if((UWORD32) mb_x == u2_frm_width_in_mb)
    {
        mb_x = 0;
        mb_y++;
        if(mb_y >= ps_dec->u2_frm_ht_in_mbs)
        {
            mb_y = ps_dec->u2_frm_ht_in_mbs - 1;
        }
    }
    /*********************************************************************/
    /* Cabac Context Initialisations                                     */
    /*********************************************************************/
    ps_dec->ps_curr_ctxt_mb_info = p_ctx_inc_mb_map + mb_x;
    ps_dec->p_left_ctxt_mb_info = p_ctx_inc_mb_map - 1;
    ps_dec->p_top_ctxt_mb_info = p_ctx_inc_mb_map - 1;

    /********************************************************************/
    /* neighbour availablility                                          */
    /********************************************************************/
    if(mb_y > ps_dec->i2_prev_slice_mby)
    {
        /* if not in the immemdiate row of prev slice end then top will be available */
        if(mb_y > (ps_dec->i2_prev_slice_mby + 1)) i2_prev_slice_mbx = -1;

        if(mb_x > i2_prev_slice_mbx)
        {
            u1_mb_ngbr_avail |= TOP_MB_AVAILABLE_MASK;
            u2_top_right_mask |= TOP_RIGHT_TOP_AVAILABLE;
            u2_top_left_mask |= TOP_LEFT_TOP_AVAILABLE;
            ps_dec->p_top_ctxt_mb_info = ps_dec->ps_curr_ctxt_mb_info;
        }
        if((mb_x > (i2_prev_slice_mbx - 1)) && ((UWORD32) mb_x != (u2_frm_width_in_mb - 1)))
        {
            u1_mb_ngbr_avail |= TOP_RIGHT_MB_AVAILABLE_MASK;
            u2_top_right_mask |= TOP_RIGHT_TOPR_AVAILABLE;
        }

        if(mb_x > (i2_prev_slice_mbx + 1))
        {
            u1_mb_ngbr_avail |= TOP_LEFT_MB_AVAILABLE_MASK;
            u2_top_left_mask |= TOP_LEFT_TOPL_AVAILABLE;
        }
        /* Next row */
        i2_prev_slice_mbx = -1;
    }
    /* Same row */
    if(mb_x > (i2_prev_slice_mbx + 1))
    {
        u1_mb_ngbr_avail |= LEFT_MB_AVAILABLE_MASK;
        u2_top_left_mask |= TOP_LEFT_LEFT_AVAILABLE;
        ps_dec->p_left_ctxt_mb_info = ps_dec->ps_curr_ctxt_mb_info - 1;
    }
    {
        mb_neigbour_params_t *ps_cur_mb_row = ps_dec->ps_cur_mb_row;
        mb_neigbour_params_t *ps_top_mb_row = ps_dec->ps_top_mb_row;
        /* copy the parameters of topleft Mb */
        ps_cur_mb_info->u1_topleft_mbtype = ps_dec->u1_topleft_mbtype;
        /* Neighbour pointer assignments*/
        ps_cur_mb_info->ps_curmb = ps_cur_mb_row + mb_x;
        ps_cur_mb_info->ps_left_mb = ps_cur_mb_row + mb_x - 1;
        ps_cur_mb_info->ps_top_mb = ps_top_mb_row + mb_x;
        ps_cur_mb_info->ps_top_right_mb = ps_top_mb_row + mb_x + 1;

        /* Update the parameters of topleftmb*/
        ps_dec->u1_topleft_mbtype = ps_cur_mb_info->ps_top_mb->u1_mb_type;
    }

    ps_dec->u2_mby = mb_y;
    ps_dec->u2_mbx = mb_x;
    ps_cur_mb_info->u2_mbx = mb_x;
    ps_cur_mb_info->u2_mby = mb_y;
    ps_cur_mb_info->u1_topmb = u1_top_mb;
    ps_dec->i4_submb_ofst += SUB_BLK_SIZE;
    ps_dec->u1_mb_ngbr_availablity = u1_mb_ngbr_avail;
    ps_cur_mb_info->u1_mb_ngbr_availablity = u1_mb_ngbr_avail;
    ps_cur_mb_info->ps_curmb->u1_mb_fld = ps_dec->u1_cur_mb_fld_dec_flag;
    ps_cur_mb_info->u1_mb_field_decodingflag = ps_dec->u1_cur_mb_fld_dec_flag;
    ps_cur_mb_info->u2_top_left_avail_mask = u2_top_left_mask;
    ps_cur_mb_info->u2_top_right_avail_mask = u2_top_right_mask;

    /*********************************************************************/
    /*                  Assign the neigbours                             */
    /*********************************************************************/
    if(u4_mbskip)
    {
        UWORD8 u1_a, u1_b;
        UWORD32 u4_ctx_inc;

        u1_a = (ps_dec->p_top_ctxt_mb_info->u1_mb_type != CAB_INFERRED)
                   ? (!!(ps_dec->p_top_ctxt_mb_info->u1_mb_type & CAB_SKIP_MASK))
                   : 0;
        u1_b = (ps_dec->p_left_ctxt_mb_info->u1_mb_type != CAB_INFERRED)
                   ? (!!(ps_dec->p_left_ctxt_mb_info->u1_mb_type & CAB_SKIP_MASK))
                   : 0;
        u4_ctx_inc = 2 - (u1_a + u1_b);

        u4_mbskip = ih264d_decode_bin(u4_ctx_inc, ps_dec->p_mb_skip_flag_t, ps_dec->ps_bitstrm,
                                      &ps_dec->s_cab_dec_env);

        if(!u4_mbskip)
        {
            if(!(u1_mb_ngbr_avail & LEFT_MB_AVAILABLE_MASK))
            {
                UWORD32 *pu4_buf;
                UWORD8 *pu1_buf;

                pu1_buf = ps_dec->pu1_left_nnz_y;
                pu4_buf = (UWORD32 *) pu1_buf;
                *pu4_buf = 0;
                pu1_buf = ps_dec->pu1_left_nnz_uv;
                pu4_buf = (UWORD32 *) pu1_buf;
                *pu4_buf = 0;

                *(ps_dec->pu1_left_yuv_dc_csbp) = 0;
                MEMSET_16BYTES(&ps_dec->pu1_left_mv_ctxt_inc[0][0], 0);
                *(UWORD32 *) ps_dec->pi1_left_ref_idx_ctxt_inc = 0;
            }
            if(!(u1_mb_ngbr_avail & TOP_MB_AVAILABLE_MASK))
            {
                MEMSET_16BYTES(ps_dec->ps_curr_ctxt_mb_info->u1_mv, 0);
                memset(ps_dec->ps_curr_ctxt_mb_info->i1_ref_idx, 0, 4);
            }
        }
    }
    return (u4_mbskip);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_mb_info_cavlc_nonmbaff                         */
/*                                                                           */
/*  Description   : This function sets the following information of cur MB   */
/*                  (a) mb_x and mb_y                                        */
/*                  (b) Neighbour availablity                                */
/*                  (c) Macroblock location in the frame buffer              */
/*                  (e) For mbaff predicts field/frame u4_flag for topMb     */
/*                      and sets the field/frame for botMb. This is          */
/*                      written in ps_dec->u1_cur_mb_fld_dec_flag            */
/*                                                                           */
/*  Inputs        : pointer to decstruct                                     */
/*                  pointer to current mb info                               */
/*                  currentMbaddress                                         */
/*                                                                           */
/*  Processing    : leftMb and TopMb params are used by DecMbskip and        */
/*                  DecCtxMbfield  modules so that these modules do not      */
/*                  check for neigbour availability and then find the        */
/*                  neigbours for context increments                         */
/*                                                                           */
/*  Returns       : OK                                                       */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 01 2023   Kishore             Draft                            */
/*                                                                           */
/*****************************************************************************/

UWORD32 isvcd_get_mb_info_cavlc_nonmbaff(dec_struct_t *ps_dec, const UWORD16 u2_cur_mb_address,
                                         dec_mb_info_t *ps_cur_mb_info, UWORD32 u4_mbskip_run)
{
    WORD32 mb_x;
    WORD32 mb_y;
    UWORD8 u1_mb_ngbr_avail = 0;
    UWORD16 u2_frm_width_in_mb = ps_dec->u2_frm_wd_in_mbs;
    WORD16 i2_prev_slice_mbx = ps_dec->i2_prev_slice_mbx;
    UWORD16 u2_top_right_mask = TOP_RIGHT_DEFAULT_AVAILABLE;
    UWORD16 u2_top_left_mask = TOP_LEFT_DEFAULT_AVAILABLE;
    UNUSED(u4_mbskip_run);
    /*--------------------------------------------------------------------*/
    /* Calculate values of mb_x and mb_y                                  */
    /*--------------------------------------------------------------------*/
    mb_x = (WORD16) ps_dec->u2_mbx;
    mb_y = (WORD16) ps_dec->u2_mby;

    ps_dec->u2_cur_mb_addr = u2_cur_mb_address;

    mb_x++;

    if(mb_x == u2_frm_width_in_mb)
    {
        mb_x = 0;
        mb_y++;
        if(mb_y >= ps_dec->u2_frm_ht_in_mbs)
        {
            mb_y = ps_dec->u2_frm_ht_in_mbs - 1;
        }
    }
    if(mb_y > ps_dec->i2_prev_slice_mby)
    {
        /* if not in the immemdiate row of prev slice end then top
         will be available */
        if(mb_y > (ps_dec->i2_prev_slice_mby + 1)) i2_prev_slice_mbx = -1;

        if(mb_x > i2_prev_slice_mbx)
        {
            u1_mb_ngbr_avail |= TOP_MB_AVAILABLE_MASK;
            u2_top_right_mask |= TOP_RIGHT_TOP_AVAILABLE;
            u2_top_left_mask |= TOP_LEFT_TOP_AVAILABLE;
        }

        if((mb_x > (i2_prev_slice_mbx - 1)) && (mb_x != (u2_frm_width_in_mb - 1)))
        {
            u1_mb_ngbr_avail |= TOP_RIGHT_MB_AVAILABLE_MASK;
            u2_top_right_mask |= TOP_RIGHT_TOPR_AVAILABLE;
        }

        if(mb_x > (i2_prev_slice_mbx + 1))
        {
            u1_mb_ngbr_avail |= TOP_LEFT_MB_AVAILABLE_MASK;
            u2_top_left_mask |= TOP_LEFT_TOPL_AVAILABLE;
        }

        /* Next row  Left will be available*/
        i2_prev_slice_mbx = -1;
    }

    /* Same row */
    if(mb_x > (i2_prev_slice_mbx + 1))
    {
        u1_mb_ngbr_avail |= LEFT_MB_AVAILABLE_MASK;
        u2_top_left_mask |= TOP_LEFT_LEFT_AVAILABLE;
    }

    {
        mb_neigbour_params_t *ps_cur_mb_row = ps_dec->ps_cur_mb_row;
        mb_neigbour_params_t *ps_top_mb_row = ps_dec->ps_top_mb_row;

        /* copy the parameters of topleft Mb */
        ps_cur_mb_info->u1_topleft_mbtype = ps_dec->u1_topleft_mbtype;
        /* Neighbour pointer assignments*/
        ps_cur_mb_info->ps_curmb = ps_cur_mb_row + mb_x;
        ps_cur_mb_info->ps_left_mb = ps_cur_mb_row + mb_x - 1;
        ps_cur_mb_info->ps_top_mb = ps_top_mb_row + mb_x;
        ps_cur_mb_info->ps_top_right_mb = ps_top_mb_row + mb_x + 1;

        /* Update the parameters of topleftmb*/
        ps_dec->u1_topleft_mbtype = ps_cur_mb_info->ps_top_mb->u1_mb_type;
    }

    ps_dec->u2_mby = mb_y;
    ps_dec->u2_mbx = mb_x;
    ps_cur_mb_info->u2_mbx = mb_x;
    ps_cur_mb_info->u2_mby = mb_y;
    ps_cur_mb_info->u1_topmb = 1;
    ps_dec->i4_submb_ofst += SUB_BLK_SIZE;
    ps_dec->u1_mb_ngbr_availablity = u1_mb_ngbr_avail;
    ps_cur_mb_info->u1_mb_ngbr_availablity = u1_mb_ngbr_avail;
    ps_cur_mb_info->ps_curmb->u1_mb_fld = ps_dec->u1_cur_mb_fld_dec_flag;
    ps_cur_mb_info->u1_mb_field_decodingflag = ps_dec->u1_cur_mb_fld_dec_flag;
    ps_cur_mb_info->u2_top_left_avail_mask = u2_top_left_mask;
    ps_cur_mb_info->u2_top_right_avail_mask = u2_top_right_mask;
    return (OK);
}
