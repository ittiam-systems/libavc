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
 *  isvce_deblk.c
 *
 * @brief
 *  This file contains functions that are associated with deblocking
 *
 * @author
 *  ittiam
 *
 * @par List of Functions:
 *  - isvce_fill_bs_1mv_1ref_non_mbaff
 *  - isvce_compute_bs
 *  - isvce_filter_top_edge
 *  - isvce_filter_left_edge
 *  - isvce_deblock_mb
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
#include <stdio.h>
#include <string.h>

#include "ih264e_config.h"
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "isvc_macros.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "ih264_trans_data.h"
#include "isvc_trans_quant_itrans_iquant.h"
#include "isvc_inter_pred_filters.h"
#include "isvc_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "isvc_cabac_tables.h"
#include "ih264_deblk_tables.h"
#include "isvce_defs.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "ime_distortion_metrics.h"
#include "ime_defs.h"
#include "ime_structs.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "isvce_rate_control.h"
#include "isvce_cabac_structs.h"
#include "isvce_structs.h"
#include "isvce_deblk.h"
#include "isvce_globals.h"

static const UWORD32 gau4_isvce_packed_bs2[(1 << MAX_TU_IN_MB_COL) * 2] = {
    /* BS TABLES FOR NORMAL EDGES */
    0x00000000, 0x02000000, 0x00020000, 0x02020000, 0x00000200, 0x02000200, 0x00020200, 0x02020200,
    0x00000002, 0x02000002, 0x00020002, 0x02020002, 0x00000202, 0x02000202, 0x00020202, 0x02020202,

    /* BS TABLES FOR XTRA LEFT MB EDGES IN MBAFF CASE */
    0x01010101, 0x02010101, 0x01020101, 0x02020101, 0x01010201, 0x02010201, 0x01020201, 0x02020201,
    0x01010102, 0x02010102, 0x01020102, 0x02020102, 0x01010202, 0x02010202, 0x01020202, 0x02020202};

static const UWORD16 gau2_isvce_4x4_v2h_reorder[(1 << MAX_TU_IN_MB_COL)] = {
    0x0000, 0x0001, 0x0010, 0x0011, 0x0100, 0x0101, 0x0110, 0x0111,
    0x1000, 0x1001, 0x1010, 0x1011, 0x1100, 0x1101, 0x1110, 0x1111};

static void isvce_fill_bs1_16x16mb_pslice(isvce_mb_info_t *ps_cur_mb, isvce_mb_info_t *ps_top_mb,
                                          isvce_mb_info_t *ps_left_mb, UWORD32 *pu4_bs_table,
                                          coordinates_t *ps_mb_pos)
{
    WORD16 i2_q_mv0, i2_q_mv1;
    WORD16 i2_p_mv0, i2_p_mv1;
    UWORD32 i;
    UWORD32 u4_bs_horz = pu4_bs_table[0];
    UWORD32 u4_bs_vert = pu4_bs_table[4];

    i2_q_mv0 = ps_cur_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
    i2_q_mv1 = ps_cur_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;

    if(ps_mb_pos->i4_ordinate)
    {
        /* Computing Bs for the top edge */
        for(i = 0; i < 4; i++)
        {
            UWORD32 u4_idx = 24 - (i << 3);

            /* check if Bs is already set */
            if(!((u4_bs_horz >> u4_idx) & 0xf))
            {
                /************************************************************/
                /* If Bs is not set, use left edge and current edge mvs and */
                /* reference pictures addresses to evaluate Bs==1           */
                /************************************************************/
                UWORD32 u4_bs_temp1;
                UWORD32 u4_bs;

                /*********************************************************/
                /* If any motion vector component differs by more than 1 */
                /* integer pel or if reference pictures are different Bs */
                /* is set to 1. Note that this condition shall be met for*/
                /* both (fwd-fwd,bwd-bwd) and (fwd-bwd,bwd-fwd) direction*/
                /*********************************************************/
                i2_p_mv0 = ps_top_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
                i2_p_mv1 = ps_top_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;

                u4_bs_temp1 =
                    ((ABS((i2_p_mv0 - i2_q_mv0)) >= 4) || (ABS((i2_p_mv1 - i2_q_mv1)) >= 4));

                u4_bs = ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_top_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         u4_bs_temp1);

                u4_bs_horz |= (u4_bs << u4_idx);
            }
        }

        pu4_bs_table[0] = u4_bs_horz;
    }

    if(ps_mb_pos->i4_abscissa)
    {
        /* Computing Bs for the left edge */
        for(i = 0; i < 4; i++)
        {
            UWORD32 u4_idx = 24 - (i << 3);

            /* check if Bs is already set */
            if(!((u4_bs_vert >> u4_idx) & 0xf))
            {
                /* If Bs is not set, evalaute conditions for Bs=1 */
                UWORD32 u4_bs_temp1;
                UWORD32 u4_bs;
                /*********************************************************/
                /* If any motion vector component differs by more than 1 */
                /* integer pel or if reference pictures are different Bs */
                /* is set to 1. Note that this condition shall be met for*/
                /* both (fwd-fwd,bwd-bwd) and (fwd-bwd,bwd-fwd) direction*/
                /*********************************************************/

                i2_p_mv0 = ps_left_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
                i2_p_mv1 = ps_left_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;

                u4_bs_temp1 =
                    ((ABS((i2_p_mv0 - i2_q_mv0)) >= 4) | (ABS((i2_p_mv1 - i2_q_mv1)) >= 4));

                u4_bs = ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_left_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         u4_bs_temp1);

                u4_bs_vert |= (u4_bs << u4_idx);
            }
        }

        pu4_bs_table[4] = u4_bs_vert;
    }
}

static void isvce_fill_bs1_16x16mb_bslice(isvce_mb_info_t *ps_cur_mb, isvce_mb_info_t *ps_top_mb,
                                          isvce_mb_info_t *ps_left_mb, UWORD32 *pu4_bs_table,
                                          coordinates_t *ps_mb_pos)
{
    WORD16 i2_q_mv0, i2_q_mv1, i2_q_mv2, i2_q_mv3;
    WORD16 i2_p_mv0, i2_p_mv1, i2_p_mv2, i2_p_mv3;
    UWORD32 i;
    UWORD32 u4_bs_horz = pu4_bs_table[0];
    UWORD32 u4_bs_vert = pu4_bs_table[4];

    i2_q_mv0 = ps_cur_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
    i2_q_mv1 = ps_cur_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;
    i2_q_mv2 = ps_cur_mb->as_pu->as_me_info[L1].s_mv.i2_mvx;
    i2_q_mv3 = ps_cur_mb->as_pu->as_me_info[L1].s_mv.i2_mvy;

    /* Computing Bs for the top edge */
    if(ps_mb_pos->i4_ordinate)
    {
        for(i = 0; i < 4; i++)
        {
            UWORD32 u4_idx = 24 - (i << 3);

            /* check if Bs is already set */
            if(!((u4_bs_horz >> u4_idx) & 0xf))
            {
                /************************************************************/
                /* If Bs is not set, use left edge and current edge mvs and */
                /* reference pictures addresses to evaluate Bs==1           */
                /************************************************************/
                UWORD32 u4_bs_temp1, u4_bs_temp2;
                UWORD32 u4_bs;

                /*********************************************************/
                /* If any motion vector component differs by more than 1 */
                /* integer pel or if reference pictures are different Bs */
                /* is set to 1. Note that this condition shall be met for*/
                /* both (fwd-fwd,bwd-bwd) and (fwd-bwd,bwd-fwd) direction*/
                /*********************************************************/
                i2_p_mv0 = ps_top_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
                i2_p_mv1 = ps_top_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;
                i2_p_mv2 = ps_top_mb->as_pu->as_me_info[L1].s_mv.i2_mvx;
                i2_p_mv3 = ps_top_mb->as_pu->as_me_info[L1].s_mv.i2_mvy;

                u4_bs_temp1 =
                    ((ABS((i2_p_mv0 - i2_q_mv0)) >= 4) | (ABS((i2_p_mv1 - i2_q_mv1)) >= 4) |
                     (ABS((i2_p_mv2 - i2_q_mv2)) >= 4) | (ABS((i2_p_mv3 - i2_q_mv3)) >= 4));

                u4_bs_temp2 =
                    ((ABS((i2_p_mv0 - i2_q_mv2)) >= 4) | (ABS((i2_p_mv1 - i2_q_mv3)) >= 4) |
                     (ABS((i2_p_mv2 - i2_q_mv0)) >= 4) | (ABS((i2_p_mv3 - i2_q_mv1)) >= 4));

                u4_bs = ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_top_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         (ps_cur_mb->as_pu->as_me_info[L1].i1_ref_idx !=
                          ps_top_mb->as_pu->as_me_info[L1].i1_ref_idx) ||
                         u4_bs_temp1) &&
                        ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_top_mb->as_pu->as_me_info[L1].i1_ref_idx) ||
                         (ps_cur_mb->as_pu->as_me_info[L1].i1_ref_idx !=
                          ps_top_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         u4_bs_temp2);

                u4_bs_horz |= (u4_bs << u4_idx);
            }
        }

        pu4_bs_table[0] = u4_bs_horz;
    }

    /* Computing Bs for the left edge */
    if(ps_mb_pos->i4_abscissa)
    {
        for(i = 0; i < 4; i++)
        {
            UWORD32 u4_idx = 24 - (i << 3);

            /* check if Bs is already set */
            if(!((u4_bs_vert >> u4_idx) & 0xf))
            {
                /* If Bs is not set, evalaute conditions for Bs=1 */
                UWORD32 u4_bs_temp1, u4_bs_temp2;
                UWORD32 u4_bs;
                /*********************************************************/
                /* If any motion vector component differs by more than 1 */
                /* integer pel or if reference pictures are different Bs */
                /* is set to 1. Note that this condition shall be met for*/
                /* both (fwd-fwd,bwd-bwd) and (fwd-bwd,bwd-fwd) direction*/
                /*********************************************************/

                i2_p_mv0 = ps_left_mb->as_pu->as_me_info[L0].s_mv.i2_mvx;
                i2_p_mv1 = ps_left_mb->as_pu->as_me_info[L0].s_mv.i2_mvy;
                i2_p_mv2 = ps_left_mb->as_pu->as_me_info[L1].s_mv.i2_mvx;
                i2_p_mv3 = ps_left_mb->as_pu->as_me_info[L1].s_mv.i2_mvy;

                u4_bs_temp1 =
                    ((ABS((i2_p_mv0 - i2_q_mv0)) >= 4) | (ABS((i2_p_mv1 - i2_q_mv1)) >= 4) |
                     (ABS((i2_p_mv2 - i2_q_mv2)) >= 4) | (ABS((i2_p_mv3 - i2_q_mv3)) >= 4));

                u4_bs_temp2 =
                    ((ABS((i2_p_mv0 - i2_q_mv2)) >= 4) | (ABS((i2_p_mv1 - i2_q_mv3)) >= 4) |
                     (ABS((i2_p_mv2 - i2_q_mv0)) >= 4) | (ABS((i2_p_mv3 - i2_q_mv1)) >= 4));

                u4_bs = ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_left_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         (ps_cur_mb->as_pu->as_me_info[L1].i1_ref_idx !=
                          ps_left_mb->as_pu->as_me_info[L1].i1_ref_idx) ||
                         u4_bs_temp1) &&
                        ((ps_cur_mb->as_pu->as_me_info[L0].i1_ref_idx !=
                          ps_left_mb->as_pu->as_me_info[L1].i1_ref_idx) ||
                         (ps_cur_mb->as_pu->as_me_info[L1].i1_ref_idx !=
                          ps_left_mb->as_pu->as_me_info[L0].i1_ref_idx) ||
                         u4_bs_temp2);

                u4_bs_vert |= (u4_bs << u4_idx);
            }
        }

        pu4_bs_table[4] = u4_bs_vert;
    }
}

static void isvce_fill_bs2_horz_vert(UWORD32 *pu4_bs, WORD32 u4_left_mb_csbp, WORD32 u4_top_mb_csbp,
                                     WORD32 u4_cur_mb_csbp, coordinates_t *ps_mb_pos,
                                     const UWORD32 *pu4_packed_bs2,
                                     const UWORD16 *pu2_4x4_v2h_reorder)
{
    UWORD32 u4_nbr_horz_csbp, u4_nbr_vert_csbp;
    UWORD32 u4_horz_bs2_dec, u4_vert_bs2_dec;
    UWORD32 u4_left_mb_masked_csbp, u4_cur_mb_masked_csbp;

    UWORD32 u4_reordered_vert_bs2_dec, u4_temp;

    WORD32 u4_cur_mb_csbp_seq = 0;
    WORD32 u4_top_mb_csbp_seq = 0;
    WORD32 u4_left_mb_csbp_seq = 0;

    /* Convert the csbp packed data in sequential pattern from raster order */
    u4_cur_mb_csbp_seq |= u4_cur_mb_csbp & 3;  // 0 1
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 4;  // 4 5
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 2;  // 2 3
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 6;  // 6 7
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 8;  // 8 9
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 12;  // 12 13
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 10;  // 10 11
    u4_cur_mb_csbp >>= 2;
    u4_cur_mb_csbp_seq |= (u4_cur_mb_csbp & 3) << 14;  // 14 15

    u4_left_mb_csbp_seq |= u4_left_mb_csbp & 3;  // 0 1
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 4;  // 4 5
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 2;  // 2 3
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 6;  // 6 7
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 8;  // 8 9
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 12;  // 12 13
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 10;  // 10 11
    u4_left_mb_csbp >>= 2;
    u4_left_mb_csbp_seq |= (u4_left_mb_csbp & 3) << 14;  // 14 15

    /* Required only the last row of top MB */
    u4_top_mb_csbp = u4_top_mb_csbp >> 10;  // 12 13
    u4_top_mb_csbp_seq |= (u4_top_mb_csbp & 3);
    u4_top_mb_csbp = u4_top_mb_csbp >> 4;  // 14 15
    u4_top_mb_csbp_seq |= ((u4_top_mb_csbp & 3) << 2);

    /* u4_nbr_horz_csbp=11C|10C|9C|8C|7C|6C|5C|4C|3C|2C|1C|0C|15T|14T|13T|12T */
    u4_nbr_horz_csbp = (u4_cur_mb_csbp_seq << 4) | u4_top_mb_csbp_seq;
    u4_horz_bs2_dec = u4_cur_mb_csbp_seq | u4_nbr_horz_csbp;

    /* u4_left_mb_masked_csbp = 15L|0|0|0|11L|0|0|0|7L|0|0|0|3L|0|0|0 */
    u4_left_mb_masked_csbp = u4_left_mb_csbp_seq & CSBP_RIGHT_BLOCK_MASK;

    /* u4_cur_mb_masked_csbp =14C|13C|12C|x|10C|9C|8C|x|6C|5C|4C|x|2C|1C|0C|x */
    u4_cur_mb_masked_csbp = (u4_cur_mb_csbp_seq << 1) & (~CSBP_LEFT_BLOCK_MASK);

    /* u4_nbr_vert_csbp=14C|13C|12C|15L|10C|9C|8C|11L|6C|5C|4C|7L|2C|1C|0C|3L */
    u4_nbr_vert_csbp = (u4_cur_mb_masked_csbp) | (u4_left_mb_masked_csbp >> 3);

    u4_vert_bs2_dec = u4_cur_mb_csbp_seq | u4_nbr_vert_csbp;

    /* Fill horz edges (0,1,2,3) boundary strengths 2 using look up table */
    if(ps_mb_pos->i4_ordinate)
    {
        pu4_bs[0] = pu4_packed_bs2[u4_horz_bs2_dec & 0xF];
    }

    pu4_bs[1] = pu4_packed_bs2[(u4_horz_bs2_dec >> 4) & 0xF];
    pu4_bs[2] = pu4_packed_bs2[(u4_horz_bs2_dec >> 8) & 0xF];
    pu4_bs[3] = pu4_packed_bs2[(u4_horz_bs2_dec >> 12) & 0xF];

    /* Do 4x4 tranpose of u4_vert_bs2_dec by using look up table for reorder */
    u4_reordered_vert_bs2_dec = pu2_4x4_v2h_reorder[u4_vert_bs2_dec & 0xF];
    u4_temp = pu2_4x4_v2h_reorder[(u4_vert_bs2_dec >> 4) & 0xF];
    u4_reordered_vert_bs2_dec |= (u4_temp << 1);
    u4_temp = pu2_4x4_v2h_reorder[(u4_vert_bs2_dec >> 8) & 0xF];
    u4_reordered_vert_bs2_dec |= (u4_temp << 2);
    u4_temp = pu2_4x4_v2h_reorder[(u4_vert_bs2_dec >> 12) & 0xF];
    u4_reordered_vert_bs2_dec |= (u4_temp << 3);

    /* Fill vert edges (4,5,6,7) boundary strengths 2 using look up table */
    if(ps_mb_pos->i4_abscissa)
    {
        pu4_bs[4] = pu4_packed_bs2[u4_reordered_vert_bs2_dec & 0xF];
    }

    pu4_bs[5] = pu4_packed_bs2[(u4_reordered_vert_bs2_dec >> 4) & 0xF];
    pu4_bs[6] = pu4_packed_bs2[(u4_reordered_vert_bs2_dec >> 8) & 0xF];
    pu4_bs[7] = pu4_packed_bs2[(u4_reordered_vert_bs2_dec >> 12) & 0xF];
}

/* brief Fills the BS for edges falling on a IBL boundary */
static void isvce_fill_bs_ibl(isvce_mb_info_t *ps_cur_mb, isvce_mb_info_t *ps_top_mb,
                              isvce_mb_info_t *ps_left_mb, UWORD32 *pu4_bs_table)
{
    /*! Flow of the module is as follows                                  */
    /*! 1. checks if MB edge is falling on IBL boundary                   */
    /*! 2. if only Mb edge then it fills the BS based on INTRA or INTER
           stauts                                                         */
    /*! 3. if the current MB is IBL and neighbours are also neighbours
           then it uses the current layer t_coeff flag to decide the
           BS of a particular edge                                        */
    /*! 4. fills the BS for all the edges in curretn MB if IBL            */

    UWORD16 u2_top_horz_nnz;
    UWORD8 u1_top_mb_ibl, u1_left_mb_ibl;
    UWORD32 i4_i, i4_edge;
    UWORD8 u1_bs;
    UWORD8 u1_cnd;
    UWORD8 u1_top_intra;
    UWORD8 u1_left_intra;
    UWORD8 u1_p_nnz, u1_q_nnz;
    UWORD8 u1_curr_mb_ibl;
    UWORD16 u2_curr_nnz;
    UWORD8 u1_left_mb_nnz = 0, u1_left_nnz;
    WORD32 i4_horz_start = 0;
    WORD32 i4_vertical_start = 0;

    u1_top_mb_ibl = ps_top_mb ? (ps_top_mb->u1_base_mode_flag && ps_top_mb->u1_is_intra) : 0;
    u1_left_mb_ibl = ps_left_mb ? (ps_left_mb->u1_base_mode_flag && ps_left_mb->u1_is_intra) : 0;

    u1_curr_mb_ibl = ps_cur_mb ? (ps_cur_mb->u1_base_mode_flag && ps_cur_mb->u1_is_intra) : 0;

    u1_top_intra = ps_top_mb ? ps_top_mb->u1_is_intra : 0;
    u1_left_intra = ps_left_mb ? ps_left_mb->u1_is_intra : 0;

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

    /* Fill BS for mb egdex assuming non IBL case */

    /* only the  MB edges fall across IBL boundary */
    if((0 != u1_curr_mb_ibl) || (0 != u1_top_mb_ibl) || (0 != u1_left_mb_ibl))
    {
        UWORD16 u2_temp, u2_i, u1_i;
        u2_temp = ps_left_mb ? ps_left_mb->u4_res_csbp : 0;
        for(u2_i = 0; u2_i < MAX_TU_IN_MB_COL; u2_i++)
        {
            UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[u2_i * 4 + MAX_TU_IN_MB_ROW - 1];
            u1_left_mb_nnz |= ((u2_temp & (1 << u1_zscan_idx)) ? 1 << u2_i : 0);
        }

        u2_curr_nnz = ps_cur_mb->u4_res_csbp;

        u2_top_horz_nnz = 0;
        if(ps_top_mb)
        {
            /* last row of top MB */
            for(u1_i = 12; u1_i < 16; u1_i++)
            {
                UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[u1_i];
                u2_top_horz_nnz |=
                    ((ps_top_mb->u4_res_csbp & (1 << u1_zscan_idx)) ? 1 << (u1_i - 12) : 0);
            }
        }
        else
        {
            u2_top_horz_nnz = 0;
        }

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
            UWORD8 u1_p_nnz = 0, u1_q_nnz = 0;
            UWORD32 u4_bs_edge = 0;
            WORD32 i4_bit_mask;
            WORD32 i4_curr_intra_flag;
            WORD32 i4_neibor_intra_flag;

            if(((1 == i4_horz_start) && (i4_edge == 1))) continue;
            if(((1 == i4_vertical_start) && (i4_edge == 0))) continue;

            i4_curr_intra_flag = (0 != u1_curr_mb_ibl);

            if(0 != i4_edge)
            {
                /* initialize for the TOP edge */
                u1_p_nnz = (UWORD8) u2_top_horz_nnz;
                for(i4_i = 0; i4_i < MAX_TU_IN_MB_ROW; i4_i++)
                {
                    UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[i4_i];
                    u1_q_nnz |= ((u2_curr_nnz & (1 << u1_zscan_idx)) ? (1 << i4_i) : 0);
                }

                i4_neibor_intra_flag = (u1_top_mb_ibl || u1_top_intra);
            }
            else
            {
                u1_p_nnz = u1_left_mb_nnz;
                for(u2_i = 0; u2_i < MAX_TU_IN_MB_COL; u2_i++)
                {
                    UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[u2_i * 4];
                    u1_q_nnz |= ((u2_curr_nnz & (1 << u1_zscan_idx)) ? 1 << u2_i : 0);
                }

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
        WORD32 i4_bit_mask_edge = 1;
        UWORD16 u2_temp, u2_i, u1_i;

        u1_left_mb_nnz = 0;
        u2_temp = ps_left_mb ? ps_left_mb->u4_csbp : 0;
        for(u2_i = 0; u2_i < MAX_TU_IN_MB_COL; u2_i++)
        {
            UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[u2_i * 4 + MAX_TU_IN_MB_ROW - 1];
            u1_left_mb_nnz |= ((u2_temp & (1 << u1_zscan_idx)) ? 1 << u2_i : 0);
        }

        u2_curr_nnz = ps_cur_mb->u4_csbp;

        u2_top_horz_nnz = 0;
        if(ps_top_mb)
        {
            for(u1_i = 12; u1_i < 16; u1_i++)
            {
                UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[u1_i];
                u2_top_horz_nnz |=
                    ((ps_top_mb->u4_csbp & (1 << u1_zscan_idx)) ? 1 << (u1_i - 12) : 0);
            }
        }
        else
        {
            u2_top_horz_nnz = 0;
        }

        /* all are IBL edges then use only t_coeff of current layer */
        /* loop over all edges */
        for(i4_edge = 0; i4_edge < 4; i4_edge++)
        {
            UWORD16 u2_curr_horz_nnz = 0;
            WORD32 i4_bit_mask = 1;

            u1_left_nnz = (u1_left_mb_nnz & i4_bit_mask_edge);

            for(i4_i = 0; i4_i < 4; i4_i++)
            {
                UWORD8 u1_curr_nnz, u1_top_nnz;
                UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[(4 * i4_edge) + i4_i];

                u2_curr_horz_nnz |= ((ps_cur_mb->u4_csbp & (1 << u1_zscan_idx)) ? (1 << i4_i) : 0);
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
}

void isvce_compute_bs(isvce_process_ctxt_t *ps_proc, UWORD8 u1_inter_layer_deblk_flag)
{
    coordinates_t s_mb_pos;

    UWORD32 *pu4_pic_vert_bs;
    UWORD32 *pu4_pic_horz_bs;

    isvce_bs_ctxt_t *ps_bs = &(ps_proc->s_deblk_ctxt.s_bs_ctxt);
    block_neighbors_t *ps_ngbr_avbl = ps_proc->ps_ngbr_avbl;
    nbr_info_t *ps_nbr_info = &ps_proc->s_nbr_info;
    isvce_mb_info_t *ps_left_mb = ps_ngbr_avbl->u1_mb_a ? ps_nbr_info->ps_left_mb_info : NULL;
    isvce_mb_info_t *ps_top_mb =
        ps_ngbr_avbl->u1_mb_b ? &ps_nbr_info->ps_top_row_mb_info[ps_bs->i4_mb_x] : NULL;
    isvce_mb_info_t *ps_cur_mb = ps_proc->ps_mb_info;

    UWORD32 u1_left_mb_intra, u1_left_mb_ibl;

    UWORD16 u2_left_csbp, u2_top_csbp, u2_cur_csbp;

    UWORD32 u4_cur_mb_intra, u1_top_mb_intra, u4_cur_mb_fld;
    UWORD32 u4_cur_mb_ibl, u1_top_mb_ibl;
    UWORD32 au4_bs_table[8];
    UWORD32 *pu4_bs_table;

    u4_cur_mb_intra = ps_cur_mb->u1_is_intra;
    u4_cur_mb_ibl = ps_cur_mb->u1_base_mode_flag && ps_cur_mb->u1_is_intra;
    u4_cur_mb_fld = 0;

    u1_top_mb_intra = ps_top_mb ? ps_top_mb->u1_is_intra : 0;
    u1_top_mb_ibl = ps_top_mb ? (ps_top_mb->u1_base_mode_flag && ps_top_mb->u1_is_intra) : 0;

    u1_left_mb_intra = ps_left_mb ? ps_left_mb->u1_is_intra : 0;
    u1_left_mb_ibl = ps_left_mb ? (ps_left_mb->u1_base_mode_flag && ps_left_mb->u1_is_intra) : 0;

    pu4_bs_table = au4_bs_table;
    memset(pu4_bs_table, 0, sizeof(pu4_bs_table[0]) * NUM_EDGES_IN_MB * 2);

    s_mb_pos.i4_abscissa = ps_bs->i4_mb_x;
    s_mb_pos.i4_ordinate = ps_bs->i4_mb_y;

    if(!u1_inter_layer_deblk_flag)
    {
        pu4_pic_vert_bs =
            ps_bs->pu4_pic_vert_bs +
            ((s_mb_pos.i4_ordinate * ps_proc->i4_wd_mbs) + s_mb_pos.i4_abscissa) * NUM_EDGES_IN_MB;
        pu4_pic_horz_bs =
            ps_bs->pu4_pic_horz_bs +
            ((s_mb_pos.i4_ordinate * ps_proc->i4_wd_mbs) + s_mb_pos.i4_abscissa) * NUM_EDGES_IN_MB;
    }
    else
    {
        pu4_pic_vert_bs =
            ps_bs->pu4_intra_base_vert_bs +
            ((s_mb_pos.i4_ordinate * ps_proc->i4_wd_mbs) + s_mb_pos.i4_abscissa) * NUM_EDGES_IN_MB;
        pu4_pic_horz_bs =
            ps_bs->pu4_intra_base_horz_bs +
            ((s_mb_pos.i4_ordinate * ps_proc->i4_wd_mbs) + s_mb_pos.i4_abscissa) * NUM_EDGES_IN_MB;
    }

    if(u4_cur_mb_intra && !(u4_cur_mb_ibl))
    {
        pu4_bs_table[4] = ps_bs->i4_mb_x ? 0x04040404 : 0;
        pu4_bs_table[0] = ps_bs->i4_mb_y ? 0x04040404 : 0;
        pu4_bs_table[1] = 0x03030303;
        pu4_bs_table[2] = 0x03030303;
        pu4_bs_table[3] = 0x03030303;
        pu4_bs_table[5] = 0x03030303;
        pu4_bs_table[6] = 0x03030303;
        pu4_bs_table[7] = 0x03030303;
    }
    else
    {
        isvce_fill_bs_ibl(ps_cur_mb, ps_top_mb, ps_left_mb, pu4_bs_table);

        if(!u4_cur_mb_ibl)
        {
            UWORD32 u4_bs_0, u4_bs_4;

            UWORD32 u4_is_b = (ps_proc->i4_slice_type == BSLICE);

            u2_cur_csbp = ps_cur_mb->u4_csbp;
            u2_left_csbp = ps_left_mb ? ps_left_mb->u4_csbp : 0;
            u2_top_csbp = ps_top_mb ? ps_top_mb->u4_csbp : 0;

            u2_cur_csbp |= (ps_cur_mb->u4_res_csbp);
            u2_left_csbp |= ps_left_mb ? ps_left_mb->u4_res_csbp : 0;
            u2_top_csbp |= ps_top_mb ? ps_top_mb->u4_res_csbp : 0;

            u4_bs_0 = pu4_bs_table[0];
            u4_bs_4 = pu4_bs_table[4];

            isvce_fill_bs2_horz_vert(pu4_bs_table, u2_left_csbp, u2_top_csbp, u2_cur_csbp,
                                     &s_mb_pos, (gau4_isvce_packed_bs2),
                                     (gau2_isvce_4x4_v2h_reorder));

            if(u1_left_mb_intra)
            {
                pu4_bs_table[4] = 0x04040404;
            }
            else if(u1_left_mb_ibl)
            {
                pu4_bs_table[4] = u4_bs_4;
            }

            if(u1_top_mb_intra)
            {
                pu4_bs_table[0] = u4_cur_mb_fld ? 0x03030303 : 0x04040404;
            }
            else if(u1_top_mb_ibl)
            {
                pu4_bs_table[0] = u4_bs_0;
            }

            if(!u4_is_b)
            {
                isvce_fill_bs1_16x16mb_pslice(ps_cur_mb, ps_top_mb, ps_left_mb, pu4_bs_table,
                                              &s_mb_pos);
            }
            else
            {
                isvce_fill_bs1_16x16mb_bslice(ps_cur_mb, ps_top_mb, ps_left_mb, pu4_bs_table,
                                              &s_mb_pos);
            }
        }
    }

    pu4_pic_horz_bs[0] = pu4_bs_table[0];
    pu4_pic_horz_bs[1] = pu4_bs_table[1];
    pu4_pic_horz_bs[2] = pu4_bs_table[2];
    pu4_pic_horz_bs[3] = pu4_bs_table[3];

    pu4_pic_vert_bs[0] = pu4_bs_table[4];
    pu4_pic_vert_bs[1] = pu4_bs_table[5];
    pu4_pic_vert_bs[2] = pu4_bs_table[6];
    pu4_pic_vert_bs[3] = pu4_bs_table[7];
}

/**
*******************************************************************************
*
* @brief This function performs deblocking of top horizontal edge
*
* @par Description:
*  This function performs deblocking of top horizontal edge
*
* @param[in] ps_codec
*  pointer to codec context
*
* @param[in] ps_proc
*  pointer to proc context
*
* @param[in] pu1_mb_qp
*  pointer to mb quantization param
*
* @param[in] pu1_cur_pic_luma
*  pointer to recon buffer luma
*
* @param[in] pu1_cur_pic_chroma
*  pointer to recon buffer chroma
*
* @param[in] pu4_pic_horz_bs
*  pointer to horizontal blocking strength
*
* @returns  none
*
* @remarks none
*
*******************************************************************************
*/
static void isvce_filter_top_edge(isvce_codec_t *ps_codec, UWORD8 u1_qp_p, UWORD8 u1_qp_q,
                                  UWORD8 *pu1_cur_pic_luma, WORD32 i4_luma_stride,
                                  UWORD8 *pu1_cur_pic_chroma, WORD32 i4_chroma_stride,
                                  UWORD32 *pu4_pic_horz_bs)
{
    UWORD32 u4_alpha_luma, u4_beta_luma, u4_qp_luma, u4_idx_A_luma, u4_idx_B_luma;
    UWORD32 u4_alpha_chroma, u4_beta_chroma, u4_qp_chroma, u4_idx_A_chroma, u4_idx_B_chroma;

    /********/
    /* luma */
    /********/
    u4_qp_luma = (u1_qp_p + u1_qp_q + 1) >> 1;

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_luma = MIN(51, u4_qp_luma + 0);
    u4_idx_B_luma = MIN(51, u4_qp_luma + 0);

    /* alpha, beta computation */
    u4_alpha_luma = gu1_ih264_alpha_table[u4_idx_A_luma];
    u4_beta_luma = gu1_ih264_beta_table[u4_idx_B_luma];

    /**********/
    /* chroma */
    /**********/
    u4_qp_chroma = (gu1_qpc_fqpi[u1_qp_p] + gu1_qpc_fqpi[u1_qp_q] + 1) >> 1;

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_chroma = MIN(51, u4_qp_chroma + 0);
    u4_idx_B_chroma = MIN(51, u4_qp_chroma + 0);

    /* alpha, beta computation */
    u4_alpha_chroma = gu1_ih264_alpha_table[u4_idx_A_chroma];
    u4_beta_chroma = gu1_ih264_beta_table[u4_idx_B_chroma];

    /* deblk edge */
    /* top Horizontal edge - allowed to be deblocked ? */
    if(pu4_pic_horz_bs[0] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_horz_bs4(pu1_cur_pic_luma, i4_luma_stride, u4_alpha_luma,
                                         u4_beta_luma);
        ps_codec->pf_deblk_chroma_horz_bs4(pu1_cur_pic_chroma, i4_chroma_stride, u4_alpha_chroma,
                                           u4_beta_chroma, u4_alpha_chroma, u4_beta_chroma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_horz_bslt4(pu1_cur_pic_luma, i4_luma_stride, u4_alpha_luma,
                                           u4_beta_luma, pu4_pic_horz_bs[0],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);

        ps_codec->pf_deblk_chroma_horz_bslt4(
            pu1_cur_pic_chroma, i4_chroma_stride, u4_alpha_chroma, u4_beta_chroma, u4_alpha_chroma,
            u4_beta_chroma, pu4_pic_horz_bs[0], gu1_ih264_clip_table[u4_idx_A_chroma],
            gu1_ih264_clip_table[u4_idx_A_chroma]);
    }
}

/**
*******************************************************************************
*
* @brief This function performs deblocking of left vertical edge
*
* @par Description:
*  This function performs deblocking of top horizontal edge
*
* @param[in] ps_codec
*  pointer to codec context
*
* @param[in] ps_proc
*  pointer to proc context
*
* @param[in] pu1_mb_qp
*  pointer to mb quantization param
*
* @param[in] pu1_cur_pic_luma
*  pointer to recon buffer luma
*
* @param[in] pu1_cur_pic_chroma
*  pointer to recon buffer chroma
*
* @param[in] pu4_pic_vert_bs
*  pointer to vertical blocking strength
*
* @returns  none
*
* @remarks none
*
*******************************************************************************
*/
static void isvce_filter_left_edge(isvce_codec_t *ps_codec, UWORD8 u1_qp_p, UWORD8 u1_qp_q,
                                   UWORD8 *pu1_cur_pic_luma, WORD32 i4_luma_stride,
                                   UWORD8 *pu1_cur_pic_chroma, WORD32 i4_chroma_stride,
                                   UWORD32 *pu4_pic_vert_bs)
{
    UWORD32 u4_alpha_luma, u4_beta_luma, u4_qp_luma, u4_idx_A_luma, u4_idx_B_luma;
    UWORD32 u4_alpha_chroma, u4_beta_chroma, u4_qp_chroma, u4_idx_A_chroma, u4_idx_B_chroma;

    /********/
    /* luma */
    /********/
    u4_qp_luma = (u1_qp_p + u1_qp_q + 1) >> 1;

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_luma = MIN(51, u4_qp_luma + 0);
    u4_idx_B_luma = MIN(51, u4_qp_luma + 0);

    /* alpha, beta computation */
    u4_alpha_luma = gu1_ih264_alpha_table[u4_idx_A_luma];
    u4_beta_luma = gu1_ih264_beta_table[u4_idx_B_luma];

    /**********/
    /* chroma */
    /**********/
    u4_qp_chroma = (gu1_qpc_fqpi[u1_qp_p] + gu1_qpc_fqpi[u1_qp_q] + 1) >> 1;

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_chroma = MIN(51, u4_qp_chroma + 0);
    u4_idx_B_chroma = MIN(51, u4_qp_chroma + 0);

    /* alpha, beta computation */
    u4_alpha_chroma = gu1_ih264_alpha_table[u4_idx_A_chroma];
    u4_beta_chroma = gu1_ih264_beta_table[u4_idx_B_chroma];

    /* deblk edge */
    if(pu4_pic_vert_bs[0] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_vert_bs4(pu1_cur_pic_luma, i4_luma_stride, u4_alpha_luma,
                                         u4_beta_luma);
        ps_codec->pf_deblk_chroma_vert_bs4(pu1_cur_pic_chroma, i4_chroma_stride, u4_alpha_chroma,
                                           u4_beta_chroma, u4_alpha_chroma, u4_beta_chroma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_vert_bslt4(pu1_cur_pic_luma, i4_luma_stride, u4_alpha_luma,
                                           u4_beta_luma, pu4_pic_vert_bs[0],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);

        ps_codec->pf_deblk_chroma_vert_bslt4(
            pu1_cur_pic_chroma, i4_chroma_stride, u4_alpha_chroma, u4_beta_chroma, u4_alpha_chroma,
            u4_beta_chroma, pu4_pic_vert_bs[0], gu1_ih264_clip_table[u4_idx_A_chroma],
            gu1_ih264_clip_table[u4_idx_A_chroma]);
    }
}

static UWORD8 isvce_get_deblk_mb_qp(isvce_process_ctxt_t *ps_proc, coordinates_t *ps_mb_pos)
{
    UWORD8 u1_mb_qp;

    isvce_deblk_ctxt_t *ps_deblk = &ps_proc->s_deblk_ctxt;
    isvce_bs_ctxt_t *ps_bs_ctxt = &ps_deblk->s_bs_ctxt;
    coordinates_t s_cur_mb_pos = {ps_deblk->i4_mb_x, ps_deblk->i4_mb_y};

    UWORD32 u4_mb_idx = ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * ps_proc->i4_wd_mbs;

    if((s_cur_mb_pos.i4_abscissa != ps_mb_pos->i4_abscissa) ||
       (s_cur_mb_pos.i4_ordinate != ps_mb_pos->i4_ordinate))
    {
        u1_mb_qp = ps_bs_ctxt->pu1_pic_qp[u4_mb_idx];
    }
    else
    {
        isvce_mb_info_t *ps_mb_info =
            ps_proc->ps_cur_mv_buf->ps_svc_layer_data[ps_proc->u1_spatial_layer_id].ps_mb_info +
            u4_mb_idx;

        if((0 == ps_mb_pos->i4_abscissa) && (0 == ps_mb_pos->i4_ordinate))
        {
            u1_mb_qp = ps_mb_info->u1_mb_qp;
        }
        else
        {
            if((ps_mb_info->u4_cbp > 0) || (I16x16 == ps_mb_info->u2_mb_type))
            {
                u1_mb_qp = ps_mb_info->u1_mb_qp;
            }
            else
            {
                u1_mb_qp = ps_bs_ctxt->pu1_pic_qp[u4_mb_idx - 1];
            }
        }
    }

    return u1_mb_qp;
}

/**
*******************************************************************************
*
* @brief This function performs deblocking on an mb
*
* @par Description:
*  This function performs deblocking on an mb
*
* @param[in] ps_proc
*  process context corresponding to the job
*
* @param[in] ps_deblk
*  pointer to deblock context
*
* @returns  none
*
* @remarks none
*
*******************************************************************************
*/
void isvce_deblock_mb(isvce_process_ctxt_t *ps_proc, isvce_deblk_ctxt_t *ps_deblk,
                      UWORD8 u1_inter_layer_deblk_flag)
{
    UWORD8 u1_mb_a, u1_mb_b;
    UWORD32 *pu4_pic_vert_bs;
    UWORD32 *pu4_pic_horz_bs;
    UWORD8 u1_cur_mb_qp;
    UWORD8 u1_left_mb_qp;
    UWORD8 u1_top_mb_qp;
    UWORD32 u4_alpha_luma, u4_beta_luma, u4_idx_A_luma, u4_idx_B_luma;
    UWORD32 u4_alpha_chroma, u4_beta_chroma, u4_qp_chroma, u4_idx_A_chroma, u4_idx_B_chroma;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    coordinates_t s_cur_mb_pos = {ps_deblk->i4_mb_x, ps_deblk->i4_mb_y};
    coordinates_t s_left_mb_pos = {ps_deblk->i4_mb_x - 1, ps_deblk->i4_mb_y};
    coordinates_t s_top_mb_pos = {ps_deblk->i4_mb_x, ps_deblk->i4_mb_y - 1};

    WORD32 i4_mb_x = ps_deblk->i4_mb_x, i4_mb_y = ps_deblk->i4_mb_y;
    WORD32 i4_luma_stride = ps_deblk->s_rec_pic_buf_props.as_component_bufs[0].i4_data_stride;
    UWORD8 *pu1_cur_pic_luma =
        (UWORD8 *) (ps_deblk->s_rec_pic_buf_props.as_component_bufs[0].pv_data) +
        (i4_mb_x * MB_SIZE) + ((i4_mb_y * MB_SIZE) * i4_luma_stride);
    WORD32 i4_chroma_stride = ps_deblk->s_rec_pic_buf_props.as_component_bufs[1].i4_data_stride;
    UWORD8 *pu1_cur_pic_chroma =
        (UWORD8 *) (ps_deblk->s_rec_pic_buf_props.as_component_bufs[1].pv_data) +
        (i4_mb_x * MB_SIZE) + (i4_mb_y * (MB_SIZE / 2) * i4_chroma_stride);
    UWORD32 push_ptr = (i4_mb_y * ps_proc->i4_wd_mbs) + i4_mb_x;

    if(!u1_inter_layer_deblk_flag)
    {
        pu4_pic_vert_bs = ps_deblk->s_bs_ctxt.pu4_pic_vert_bs;
        pu4_pic_horz_bs = ps_deblk->s_bs_ctxt.pu4_pic_horz_bs;
    }
    else
    {
        pu4_pic_vert_bs = ps_deblk->s_bs_ctxt.pu4_intra_base_vert_bs;
        pu4_pic_horz_bs = ps_deblk->s_bs_ctxt.pu4_intra_base_horz_bs;
    }

    /* derive neighbor availability */
    /* In slice mode the edges of mbs that lie on the slice boundary are not
     * deblocked */
    /* deblocking filter idc '2' */
    if(ps_codec->s_cfg.e_slice_mode != IVE_SLICE_MODE_NONE)
    {
        /* slice index */
        UWORD8 *pu1_slice_idx = ps_deblk->pu1_slice_idx;

        pu1_slice_idx += (i4_mb_y * ps_proc->i4_wd_mbs);
        /* left macroblock availability */
        u1_mb_a = (i4_mb_x == 0 || (pu1_slice_idx[i4_mb_x - 1] != pu1_slice_idx[i4_mb_x])) ? 0 : 1;
        /* top macroblock availability */
        u1_mb_b = (i4_mb_y == 0 ||
                   (pu1_slice_idx[i4_mb_x - ps_proc->i4_wd_mbs] != pu1_slice_idx[i4_mb_x]))
                      ? 0
                      : 1;
    }
    else
    {
        /* left macroblock availability */
        u1_mb_a = (i4_mb_x == 0) ? 0 : 1;
        /* top macroblock availability */
        u1_mb_b = (i4_mb_y == 0) ? 0 : 1;
    }

    pu4_pic_vert_bs += push_ptr * NUM_EDGES_IN_MB;
    pu4_pic_horz_bs += push_ptr * NUM_EDGES_IN_MB;

    /********/
    /* luma */
    /********/
    u1_cur_mb_qp = isvce_get_deblk_mb_qp(ps_proc, &s_cur_mb_pos);
    ps_deblk->s_bs_ctxt.pu1_pic_qp[push_ptr] = u1_cur_mb_qp;

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_luma = MIN(51, u1_cur_mb_qp + 0);
    u4_idx_B_luma = MIN(51, u1_cur_mb_qp + 0);

    /* alpha, beta computation */
    u4_alpha_luma = gu1_ih264_alpha_table[u4_idx_A_luma];
    u4_beta_luma = gu1_ih264_beta_table[u4_idx_B_luma];

    /**********/
    /* chroma */
    /**********/
    u4_qp_chroma = gu1_qpc_fqpi[u1_cur_mb_qp];

    /* filter offset A and filter offset B have to be received from slice header
     */
    /* TODO : for now lets set these offsets as zero */

    u4_idx_A_chroma = MIN(51, u4_qp_chroma + 0);
    u4_idx_B_chroma = MIN(51, u4_qp_chroma + 0);

    /* alpha, beta computation */
    u4_alpha_chroma = gu1_ih264_alpha_table[u4_idx_A_chroma];
    u4_beta_chroma = gu1_ih264_beta_table[u4_idx_B_chroma];

    /* Deblock vertical edges */
    /* left vertical edge 0 - allowed to be deblocked ? */
    if(u1_mb_a)
    {
        u1_left_mb_qp = isvce_get_deblk_mb_qp(ps_proc, &s_left_mb_pos);

        isvce_filter_left_edge(ps_codec, u1_left_mb_qp, u1_cur_mb_qp, pu1_cur_pic_luma,
                               i4_luma_stride, pu1_cur_pic_chroma, i4_chroma_stride,
                               pu4_pic_vert_bs);
    }

    /* vertical edge 1 */
    if(pu4_pic_vert_bs[1] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_vert_bs4(pu1_cur_pic_luma + 4, i4_luma_stride, u4_alpha_luma,
                                         u4_beta_luma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_vert_bslt4(pu1_cur_pic_luma + 4, i4_luma_stride, u4_alpha_luma,
                                           u4_beta_luma, pu4_pic_vert_bs[1],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);
    }

    /* vertical edge 2 */
    if(pu4_pic_vert_bs[2] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_vert_bs4(pu1_cur_pic_luma + 8, i4_luma_stride, u4_alpha_luma,
                                         u4_beta_luma);
        ps_codec->pf_deblk_chroma_vert_bs4(pu1_cur_pic_chroma + 8, i4_chroma_stride,
                                           u4_alpha_chroma, u4_beta_chroma, u4_alpha_chroma,
                                           u4_beta_chroma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_vert_bslt4(pu1_cur_pic_luma + 8, i4_luma_stride, u4_alpha_luma,
                                           u4_beta_luma, pu4_pic_vert_bs[2],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);

        ps_codec->pf_deblk_chroma_vert_bslt4(
            pu1_cur_pic_chroma + 8, i4_chroma_stride, u4_alpha_chroma, u4_beta_chroma,
            u4_alpha_chroma, u4_beta_chroma, pu4_pic_vert_bs[2],
            gu1_ih264_clip_table[u4_idx_A_chroma], gu1_ih264_clip_table[u4_idx_A_chroma]);
    }

    /* vertical edge 3 */
    if(pu4_pic_vert_bs[3] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_vert_bs4(pu1_cur_pic_luma + 12, i4_luma_stride, u4_alpha_luma,
                                         u4_beta_luma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_vert_bslt4(pu1_cur_pic_luma + 12, i4_luma_stride, u4_alpha_luma,
                                           u4_beta_luma, pu4_pic_vert_bs[3],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);
    }

    /* Deblock Horizontal edges */
    /* Horizontal edge 0 */
    if(u1_mb_b)
    {
        u1_top_mb_qp = isvce_get_deblk_mb_qp(ps_proc, &s_top_mb_pos);

        isvce_filter_top_edge(ps_codec, u1_top_mb_qp, u1_cur_mb_qp, pu1_cur_pic_luma,
                              i4_luma_stride, pu1_cur_pic_chroma, i4_chroma_stride,
                              pu4_pic_horz_bs);
    }

    /* horizontal edge 1 */
    if(pu4_pic_horz_bs[1] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_horz_bs4(pu1_cur_pic_luma + 4 * i4_luma_stride, i4_luma_stride,
                                         u4_alpha_luma, u4_beta_luma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_horz_bslt4(pu1_cur_pic_luma + 4 * i4_luma_stride, i4_luma_stride,
                                           u4_alpha_luma, u4_beta_luma, pu4_pic_horz_bs[1],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);
    }

    /* horizontal edge 2 */
    if(pu4_pic_horz_bs[2] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_horz_bs4(pu1_cur_pic_luma + 8 * i4_luma_stride, i4_luma_stride,
                                         u4_alpha_luma, u4_beta_luma);
        ps_codec->pf_deblk_chroma_horz_bs4(pu1_cur_pic_chroma + 4 * i4_chroma_stride,
                                           i4_chroma_stride, u4_alpha_chroma, u4_beta_chroma,
                                           u4_alpha_chroma, u4_beta_chroma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_horz_bslt4(pu1_cur_pic_luma + 8 * i4_luma_stride, i4_luma_stride,
                                           u4_alpha_luma, u4_beta_luma, pu4_pic_horz_bs[2],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);

        ps_codec->pf_deblk_chroma_horz_bslt4(
            pu1_cur_pic_chroma + 4 * i4_chroma_stride, i4_chroma_stride, u4_alpha_chroma,
            u4_beta_chroma, u4_alpha_chroma, u4_beta_chroma, pu4_pic_horz_bs[2],
            gu1_ih264_clip_table[u4_idx_A_chroma], gu1_ih264_clip_table[u4_idx_A_chroma]);
    }

    /* horizontal edge 3 */
    if(pu4_pic_horz_bs[3] == 0x04040404)
    {
        /* strong filter */
        ps_codec->pf_deblk_luma_horz_bs4(pu1_cur_pic_luma + 12 * i4_luma_stride, i4_luma_stride,
                                         u4_alpha_luma, u4_beta_luma);
    }
    else
    {
        /* normal filter */
        ps_codec->pf_deblk_luma_horz_bslt4(pu1_cur_pic_luma + 12 * i4_luma_stride, i4_luma_stride,
                                           u4_alpha_luma, u4_beta_luma, pu4_pic_horz_bs[3],
                                           gu1_ih264_clip_table[u4_idx_A_luma]);
    }
}
