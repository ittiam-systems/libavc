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

/**
 *******************************************************************************
 * @file
 *  ih264e_me.c
 *
 * @brief
 *  Contains definition of functions for motion estimation
 *
 * @author
 *  ittiam
 *
 * @par List of Functions:
 *  - ih264e_init_mv_bits()
 *  - ih264e_skip_analysis_chroma()
 *  - ih264e_skip_analysis_luma()
 *  - ih264e_analyse_skip()
 *  - ih264e_get_search_candidates()
 *  - ih264e_find_skip_motion_vector()
 *  - ih264e_get_mv_predictor()
 *  - ih264e_mv_pred()
 *  - ih264e_mv_pred_me()
 *  - ih264e_init_me()
 *  - ih264e_compute_me()
 *  - ih264e_compute_me_nmb()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <assert.h>
#include <limits.h>

/* User include files */
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ithread.h"
#include "ih264_platform_macros.h"
#include "ih264_defs.h"
#include "ime_defs.h"
#include "ime_distortion_metrics.h"
#include "ime_structs.h"
#include "ih264_structs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "ih264_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "ih264e_defs.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "ih264e_rate_control.h"
#include "ih264e_structs.h"
#include "ih264e_globals.h"
#include "ih264_macros.h"
#include "ih264e_me.h"
#include "ime.h"
#include "ime_distortion_metrics.h"
#include "ih264_debug.h"
#include "ithread.h"
#include "ih264e_intra_modes_eval.h"
#include "ih264e_core_coding.h"
#include "ih264e_mc.h"
#include "ih264e_debug.h"
#include "ih264e_half_pel.h"
#include "ime_statistics.h"
#include "ih264e_platform_macros.h"


/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
*******************************************************************************
*
* @brief
*  This function populates the length of the codewords for motion vectors in the
*  range (-search range, search range) in pixels
*
* @param[in] ps_me
*  Pointer to me ctxt
*
* @param[out] pu1_mv_bits
*  length of the codeword for all mv's
*
* @remarks The length of the code words are derived from signed exponential
* goloumb codes.
*
*******************************************************************************
*/
void ih264e_init_mv_bits(me_ctxt_t *ps_me_ctxt)
{
    /* temp var */
    WORD32 i, codesize = 3, diff, limit;
    UWORD32 u4_code_num, u4_range;
    UWORD32 u4_uev_min, u4_uev_max, u4_sev_min, u4_sev_max;

    /* max srch range */
    diff = MAX(DEFAULT_MAX_SRCH_RANGE_X, DEFAULT_MAX_SRCH_RANGE_Y);
    /* sub pel */
    diff <<= 2;
    /* delta mv */
    diff <<= 1;

    /* codeNum for positive integer     =  2x-1     : Table9-3  */
    u4_code_num = (diff << 1);

    /* get range of the bit string and put using put_bits()                 */
    GETRANGE(u4_range, u4_code_num);

    limit = 2*u4_range - 1;

    /* init mv bits */
    ps_me_ctxt->pu1_mv_bits[0] = 1;

    while (codesize < limit)
    {
        u4_uev_min = (1 << (codesize >> 1));
        u4_uev_max = 2*u4_uev_min - 1;

        u4_sev_min = u4_uev_min >> 1;
        u4_sev_max = u4_uev_max >> 1;

        DEBUG("\n%d min, %d max %d codesize", u4_sev_min, u4_sev_max, codesize);

        for (i = u4_sev_min; i <= (WORD32)u4_sev_max; i++)
        {
            ps_me_ctxt->pu1_mv_bits[-i] = ps_me_ctxt->pu1_mv_bits[i] = codesize;
        }

        codesize += 2;
    }
}

/**
*******************************************************************************
*
* @brief Determines the valid candidates for which the initial search shall happen.
* The best of these candidates is used to center the diamond pixel search.
*
* @par Description: The function sends the skip, (0,0), left, top and top-right
* neighbouring MBs MVs. The left, top and top-right MBs MVs are used because
* these are the same MVs that are used to form the MV predictor. This initial MV
* search candidates need not take care of slice boundaries and hence neighbor
* availability checks are not made here.
*
* @param[in] ps_left_mb_pu
*  pointer to left mb motion vector info
*
* @param[in] ps_top_mb_pu
*  pointer to top & top right mb motion vector info
*
* @param[in] ps_top_left_mb_pu
*  pointer to top left mb motion vector info
*
* @param[out] ps_skip_mv
*  pointer to skip motion vectors for the curr mb
*
* @param[in] i4_mb_x
*  mb index x
*
* @param[in] i4_mb_y
*  mb index y
*
* @param[in] i4_wd_mbs
*  pic width in mbs
*
* @param[in] ps_motionEst
*  pointer to me context
*
* @returns  The list of MVs to be used of priming the full pel search and the
* number of such MVs
*
* @remarks
*   Assumptions : 1. Assumes Single reference frame
*                 2. Assumes Only partition of size 16x16
*
*******************************************************************************
*/
static void ih264e_get_search_candidates(process_ctxt_t *ps_proc,
                                         me_ctxt_t *ps_me_ctxt)
{
    /* curr mb indices */
    WORD32 i4_mb_x = ps_proc->i4_mb_x;

    /* left mb motion vector */
    mv_t *ps_left_mv;

    /* top left mb motion vector */
    mv_t *ps_top_mv;

    /* top left mb motion vector */
    mv_t *ps_top_left_mv;

    /* top left mb motion vector */
    mv_t *ps_top_right_mv;

    /* skip mv */
    mv_t *ps_skip_mv = ps_proc->ps_skip_mv;

    /* mb part info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;

    /* num of candidate search candidates */
    UWORD32 u4_num_candidates = 0;

    /* mvs */
    WORD32 mvx, mvy;

    /* ngbr availability */
    block_neighbors_t *ps_ngbr_avbl = ps_proc->ps_ngbr_avbl;

    /* srch range*/
    WORD32 i4_srch_range_n = ps_me_ctxt->i4_srch_range_n;
    WORD32 i4_srch_range_s = ps_me_ctxt->i4_srch_range_s;
    WORD32 i4_srch_range_e = ps_me_ctxt->i4_srch_range_e;
    WORD32 i4_srch_range_w = ps_me_ctxt->i4_srch_range_w;

    ps_left_mv = &ps_proc->s_left_mb_pu_ME.s_l0_mv;
    ps_top_mv = &(ps_proc->ps_top_row_pu_ME + i4_mb_x)->s_l0_mv;
    ps_top_left_mv = &ps_proc->s_top_left_mb_pu_ME.s_l0_mv;
    ps_top_right_mv = &(ps_proc->ps_top_row_pu_ME + i4_mb_x + 1)->s_l0_mv;

    /************************************************************/
    /* Taking the Zero motion vector as one of the candidates   */
    /************************************************************/
    ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = 0;
    ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = 0;

    u4_num_candidates++;

    /************************************************************/
    /* Taking the Left MV Predictor as one of the candidates    */
    /************************************************************/
    if (ps_ngbr_avbl->u1_mb_a)
    {
        mvx      = (ps_left_mv->i2_mvx + 2) >> 2;
        mvy      = (ps_left_mv->i2_mvy + 2) >> 2;

        mvx = CLIP3(i4_srch_range_w, i4_srch_range_e, mvx);
        mvy = CLIP3(i4_srch_range_n, i4_srch_range_s, mvy);

        ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = mvx;
        ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = mvy;

        u4_num_candidates ++;
    }
    /*else
    {
        ps_me_ctxt->as_mv_init_search[LEFT_CAND].i2_mvx = 0;
        ps_me_ctxt->as_mv_init_search[LEFT_CAND].i2_mvy = 0;
    }*/

    /************************************************************/
    /* Taking the Top MV Predictor as one of the candidates     */
    /************************************************************/
    if (ps_ngbr_avbl->u1_mb_b)
    {
        mvx      = (ps_top_mv->i2_mvx + 2) >> 2;
        mvy      = (ps_top_mv->i2_mvy + 2) >> 2;

        mvx = CLIP3(i4_srch_range_w, i4_srch_range_e, mvx);
        mvy = CLIP3(i4_srch_range_n, i4_srch_range_s, mvy);

        ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = mvx;
        ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = mvy;

        u4_num_candidates ++;

        /************************************************************/
        /* Taking the TopRt MV Predictor as one of the candidates   */
        /************************************************************/
        if (ps_ngbr_avbl->u1_mb_c)
        {
            mvx      = (ps_top_right_mv->i2_mvx + 2) >> 2;
            mvy      = (ps_top_right_mv->i2_mvy + 2)>> 2;

            mvx = CLIP3(i4_srch_range_w, i4_srch_range_e, mvx);
            mvy = CLIP3(i4_srch_range_n, i4_srch_range_s, mvy);

            ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = mvx;
            ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = mvy;

            u4_num_candidates ++;
        }
        /************************************************************/
        /* Taking the TopLt MV Predictor as one of the candidates   */
        /************************************************************/
        else if (ps_ngbr_avbl->u1_mb_d)
        {
            mvx      = (ps_top_left_mv->i2_mvx + 2) >> 2;
            mvy      = (ps_top_left_mv->i2_mvy + 2) >> 2;

            mvx = CLIP3(i4_srch_range_w, i4_srch_range_e, mvx);
            mvy = CLIP3(i4_srch_range_n, i4_srch_range_s, mvy);

            ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = mvx;
            ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = mvy;

            u4_num_candidates ++;
        }
        /*else
        {
            ps_me_ctxt->as_mv_init_search[TOPR_CAND].i2_mvx = 0;
            ps_me_ctxt->as_mv_init_search[TOPR_CAND].i2_mvy = 0;
        }*/
    }
    /*else
    {
        ps_me_ctxt->as_mv_init_search[TOP_CAND].i2_mvx = 0;
        ps_me_ctxt->as_mv_init_search[TOP_CAND].i2_mvy = 0;

        ps_me_ctxt->as_mv_init_search[TOPR_CAND].i2_mvx = 0;
        ps_me_ctxt->as_mv_init_search[TOPR_CAND].i2_mvy = 0;
    }*/


    /********************************************************************/
    /*                            MV Prediction                         */
    /********************************************************************/
    ih264e_mv_pred_me(ps_proc);

    ps_mb_part->s_mv_pred.i2_mvx = ps_proc->ps_pred_mv->i2_mvx;
    ps_mb_part->s_mv_pred.i2_mvy = ps_proc->ps_pred_mv->i2_mvy;

    /************************************************************/
    /* Get the skip motion vector                               */
    /************************************************************/
    ih264e_find_skip_motion_vector(ps_proc, 1);

    /************************************************************/
    /* Taking the Skip motion vector as one of the candidates   */
    /************************************************************/
    mvx = (ps_skip_mv->i2_mvx + 2) >> 2;
    mvy = (ps_skip_mv->i2_mvy + 2) >> 2;

    mvx = CLIP3(i4_srch_range_w, i4_srch_range_e, mvx);
    mvy = CLIP3(i4_srch_range_n, i4_srch_range_s, mvy);

    ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvx = mvx;
    ps_me_ctxt->as_mv_init_search[u4_num_candidates].i2_mvy = mvy;

    u4_num_candidates++;

    ASSERT(u4_num_candidates <= 5);

    ps_me_ctxt->u4_num_candidates = u4_num_candidates;
}

/**
*******************************************************************************
*
* @brief The function gives the skip motion vector
*
* @par Description:
*  The function gives the skip motion vector
*
* @param[in] ps_left_mb_pu
*  pointer to left mb motion vector info
*
* @param[in] ps_top_row_pu
*  pointer to top & top right mb motion vector info
*
* @param[out] ps_pred_mv
*  pointer to candidate predictors for the current block
*
* @returns The x & y components of the MV predictor.
*
* @remarks The code implements the logic as described in sec 8.4.1.1 in H264
*   specification.
*
*******************************************************************************
*/
void ih264e_find_skip_motion_vector(process_ctxt_t *ps_proc, UWORD32 u4_for_me)
{
    /* left mb motion vector */
    enc_pu_t *ps_left_mb_pu ;

    /* top mb motion vector */
    enc_pu_t *ps_top_mb_pu ;

    /* skip mv */
    mv_t *ps_skip_mv = ps_proc->ps_skip_mv;

    if (u4_for_me == 1)
    {
        ps_left_mb_pu = &ps_proc->s_left_mb_pu_ME;
        ps_top_mb_pu = ps_proc->ps_top_row_pu_ME + ps_proc->i4_mb_x;
    }
    else
    {
        ps_left_mb_pu = &ps_proc->s_left_mb_pu ;
        ps_top_mb_pu = ps_proc->ps_top_row_pu + ps_proc->i4_mb_x;
    }

    if (  (!ps_proc->ps_ngbr_avbl->u1_mb_a) ||
          (!ps_proc->ps_ngbr_avbl->u1_mb_b) ||
          ((ps_left_mb_pu->i1_l0_ref_idx | ps_left_mb_pu->s_l0_mv.i2_mvx | ps_left_mb_pu->s_l0_mv.i2_mvy) == 0) ||
          ((ps_top_mb_pu->i1_l0_ref_idx | ps_top_mb_pu->s_l0_mv.i2_mvx | ps_top_mb_pu->s_l0_mv.i2_mvy) == 0) )
    {
        ps_skip_mv->i2_mvx = 0;
        ps_skip_mv->i2_mvy = 0;
    }
    else
    {
        ps_skip_mv->i2_mvx = ps_proc->ps_pred_mv->i2_mvx;
        ps_skip_mv->i2_mvy = ps_proc->ps_pred_mv->i2_mvy;
    }
}

/**
*******************************************************************************
*
* @brief motion vector predictor
*
* @par Description:
*  The routine calculates the motion vector predictor for a given block,
*  given the candidate MV predictors.
*
* @param[in] ps_left_mb_pu
*  pointer to left mb motion vector info
*
* @param[in] ps_top_row_pu
*  pointer to top & top right mb motion vector info
*
* @param[out] ps_pred_mv
*  pointer to candidate predictors for the current block
*
* @returns  The x & y components of the MV predictor.
*
* @remarks The code implements the logic as described in sec 8.4.1.3 in H264
*   specification.
*   Assumptions : 1. Assumes Single reference frame
*                 2. Assumes Only partition of size 16x16
*
*******************************************************************************
*/
void ih264e_get_mv_predictor(enc_pu_t *ps_left_mb_pu,
                             enc_pu_t *ps_top_row_pu,
                             mv_t *ps_pred_mv)
{
    /* curr frame ref idx */
    /* we are assuming that we are operating on single reference frame
     * hence the ref idx is insignificant during mv prediction.
     */
    WORD32 u4_ref_idx = 0;

    /* temp var */
    WORD32 pred_algo = 3, a, b, c;

    /* If only one of the candidate blocks has a reference frame equal to
     * the current block then use the same block as the final predictor */
    a = (ps_left_mb_pu->i1_l0_ref_idx == u4_ref_idx)? 0:-1;
    b = (ps_top_row_pu[0].i1_l0_ref_idx == u4_ref_idx)? 0:-1;
    c = (ps_top_row_pu[1].i1_l0_ref_idx == u4_ref_idx)? 0:-1;

    if (a == 0 && b == -1 && c == -1)
        pred_algo = 0; /* LEFT */
    else if (a == -1 && b == 0 && c == -1)
        pred_algo = 1; /* TOP */
    else if (a == -1 && b == -1 && c == 0)
        pred_algo = 2; /* TOP RIGHT */

    switch (pred_algo)
    {
        case 0:
            /* left */
            ps_pred_mv->i2_mvx = ps_left_mb_pu->s_l0_mv.i2_mvx;
            ps_pred_mv->i2_mvy = ps_left_mb_pu->s_l0_mv.i2_mvy;
            break;
        case 1:
            /* top */
            ps_pred_mv->i2_mvx = ps_top_row_pu[0].s_l0_mv.i2_mvx;
            ps_pred_mv->i2_mvy = ps_top_row_pu[0].s_l0_mv.i2_mvy;
            break;
        case 2:
            /* top right */
            ps_pred_mv->i2_mvx = ps_top_row_pu[1].s_l0_mv.i2_mvx;
            ps_pred_mv->i2_mvy = ps_top_row_pu[1].s_l0_mv.i2_mvy;
            break;
        case 3:
            /* median */
            MEDIAN(ps_left_mb_pu->s_l0_mv.i2_mvx,
                   ps_top_row_pu[0].s_l0_mv.i2_mvx,
                   ps_top_row_pu[1].s_l0_mv.i2_mvx,
                   ps_pred_mv->i2_mvx);
            MEDIAN(ps_left_mb_pu->s_l0_mv.i2_mvy,
                   ps_top_row_pu[0].s_l0_mv.i2_mvy,
                   ps_top_row_pu[1].s_l0_mv.i2_mvy,
                   ps_pred_mv->i2_mvy);

            break;
        default:
            break;
    }
}

/**
*******************************************************************************
*
* @brief This function performs MV prediction
*
* @par Description:
*
* @param[in] ps_proc
*  Process context corresponding to the job
*
* @returns  none
*
* @remarks none
*  This function will update the MB availability since intra inter decision
*  should be done before the call
*
*******************************************************************************
*/
void ih264e_mv_pred(process_ctxt_t *ps_proc)
{

    /* left mb motion vector */
    enc_pu_t *ps_left_mb_pu ;

    /* top left mb motion vector */
    enc_pu_t *ps_top_left_mb_pu ;

    /* top row motion vector info */
    enc_pu_t *ps_top_row_pu;

    /* predicted motion vector */
    mv_t *ps_pred_mv = ps_proc->ps_pred_mv;

    /* zero mv */
    mv_t zero_mv = {0, 0};

    /*  mb neighbor availability */
    block_neighbors_t *ps_ngbr_avbl = ps_proc->ps_ngbr_avbl;

    /* mb syntax elements of neighbors */
    mb_info_t   *ps_top_syn = ps_proc->ps_top_row_mb_syntax_ele + ps_proc->i4_mb_x;
    mb_info_t   *ps_top_left_syn;
    UWORD32     u4_left_is_intra;

    ps_top_left_syn = &(ps_proc->s_top_left_mb_syntax_ele);
    u4_left_is_intra = ps_proc->s_left_mb_syntax_ele.u2_is_intra;
    ps_left_mb_pu = &ps_proc->s_left_mb_pu;
    ps_top_left_mb_pu = &ps_proc->s_top_left_mb_pu;
    ps_top_row_pu = (ps_proc->ps_top_row_pu + ps_proc->i4_mb_x);

    /* Before performing mv prediction prepare the ngbr information and
     * reset motion vectors basing on their availability */
    if (!ps_ngbr_avbl->u1_mb_a || (u4_left_is_intra == 1) )
    {
        /* left mv */
        ps_left_mb_pu->i1_l0_ref_idx = -1;
        ps_left_mb_pu->s_l0_mv = zero_mv;
    }
    if (!ps_ngbr_avbl->u1_mb_b || ps_top_syn->u2_is_intra)
    {
        /* top mv */
        ps_top_row_pu[0].i1_l0_ref_idx = -1;
        ps_top_row_pu[0].s_l0_mv = zero_mv;
    }
    if (!ps_ngbr_avbl->u1_mb_c)
    {
        /* top right mv - When top right partition is not available for
         * prediction if top left is available use it for prediction else
         * set the mv information to -1 and (0, 0)
         * */
        if (!ps_ngbr_avbl->u1_mb_d || ps_top_left_syn->u2_is_intra)
        {
            ps_top_row_pu[1].i1_l0_ref_idx = -1;
            ps_top_row_pu[1].s_l0_mv = zero_mv;
        }
        else
        {
            ps_top_row_pu[1].i1_l0_ref_idx = ps_top_left_mb_pu->i1_l0_ref_idx;
            ps_top_row_pu[1].s_l0_mv = ps_top_left_mb_pu->s_l0_mv;
        }
    }
    else if (ps_top_syn[1].u2_is_intra)
    {
        ps_top_row_pu[1].i1_l0_ref_idx = -1;
        ps_top_row_pu[1].s_l0_mv = zero_mv;
    }

    ih264e_get_mv_predictor(ps_left_mb_pu, ps_top_row_pu, ps_pred_mv);
}

/**
*******************************************************************************
*
* @brief This function approximates Pred. MV
*
* @par Description:
*
* @param[in] ps_proc
*  Process context corresponding to the job
*
* @returns  none
*
* @remarks none
*  Motion estimation happens at nmb level. For cost calculations, mv is appro
*  ximated using this function
*
*******************************************************************************
*/
void ih264e_mv_pred_me(process_ctxt_t *ps_proc)
{
    /* left mb motion vector */
    enc_pu_t *ps_left_mb_pu ;

    /* top left mb motion vector */
    enc_pu_t *ps_top_left_mb_pu ;

    /* top row motion vector info */
    enc_pu_t *ps_top_row_pu;

    enc_pu_t s_top_row_pu[2];

    /* predicted motion vector */
    mv_t *ps_pred_mv = ps_proc->ps_pred_mv;

    /* zero mv */
    mv_t zero_mv = {0, 0};

    /*  mb neighbor availability */
    block_neighbors_t *ps_ngbr_avbl = ps_proc->ps_ngbr_avbl;

    ps_left_mb_pu = &ps_proc->s_left_mb_pu_ME;
    ps_top_left_mb_pu = &ps_proc->s_top_left_mb_pu_ME;
    ps_top_row_pu = (ps_proc->ps_top_row_pu_ME + ps_proc->i4_mb_x);

    s_top_row_pu[0] = ps_top_row_pu[0];
    s_top_row_pu[1] = ps_top_row_pu[1];

    /* Before performing mv prediction prepare the ngbr information and
     * reset motion vectors basing on their availability */
    if (!ps_ngbr_avbl->u1_mb_a  )
    {
        /* left mv */
        ps_left_mb_pu->i1_l0_ref_idx = -1;
        ps_left_mb_pu->s_l0_mv = zero_mv;
    }
    if (!ps_ngbr_avbl->u1_mb_b )
    {
        /* top mv */
        s_top_row_pu[0].i1_l0_ref_idx = -1;
        s_top_row_pu[0].s_l0_mv = zero_mv;
    }
    if (!ps_ngbr_avbl->u1_mb_c)
    {
        /* top right mv - When top right partition is not available for
         * prediction if top left is available use it for prediction else
         * set the mv information to -1 and (0, 0)
         * */
        if (!ps_ngbr_avbl->u1_mb_d)
        {
            s_top_row_pu[1].i1_l0_ref_idx = -1;
            s_top_row_pu[1].s_l0_mv = zero_mv;
        }
        else
        {
            s_top_row_pu[1].i1_l0_ref_idx = ps_top_left_mb_pu->i1_l0_ref_idx;
            s_top_row_pu[1].s_l0_mv = ps_top_left_mb_pu->s_l0_mv;
        }
    }

    ih264e_get_mv_predictor(ps_left_mb_pu, &(s_top_row_pu[0]), ps_pred_mv);
}

/**
*******************************************************************************
*
* @brief This function initializes me ctxt
*
* @par Description:
*  Before dispatching the current job to me thread, the me context associated
*  with the job is initialized.
*
* @param[in] ps_proc
*  Process context corresponding to the job
*
* @returns  none
*
* @remarks none
*
*******************************************************************************
*/
void ih264e_init_me(process_ctxt_t *ps_proc)
{
    /* me ctxt */
    me_ctxt_t *ps_me_ctxt = &ps_proc->s_me_ctxt;

    /* src ptr */
    ps_me_ctxt->pu1_src_buf_luma = ps_proc->pu1_src_buf_luma;

    /* ref ptr */
    ps_me_ctxt->pu1_ref_buf_luma = ps_proc->pu1_ref_buf_luma;

    /* lagrange param */
    ps_me_ctxt->u4_lambda_motion = gu1_qp0[ps_me_ctxt->u1_mb_qp];
}

/**
*******************************************************************************
*
* @brief This function performs motion estimation for the current mb
*
* @par Description:
*  The current mb is compared with a list of mb's in the reference frame for
*  least cost. The mb that offers least cost is chosen as predicted mb and the
*  displacement of the predicted mb from index location of the current mb is
*  signaled as mv. The list of the mb's that are chosen in the reference frame
*  are dependent on the speed of the ME configured.
*
* @param[in] ps_proc
*  Process context corresponding to the job
*
* @returns  motion vector of the pred mb, sad, cost.
*
* @remarks none
*
*******************************************************************************
*/
void ih264e_compute_me(process_ctxt_t *ps_proc)
{
    /* me ctxt */
    me_ctxt_t *ps_me_ctxt = &ps_proc->s_me_ctxt;

    /* codec context */
    codec_t *ps_codec = ps_proc->ps_codec;

//    /* mb syntax elements of neighbors */
//    mb_info_t *ps_top_syn = ps_proc->ps_top_row_mb_syntax_ele + ps_proc->i4_mb_x;
//    mb_info_t *ps_top_left_syn = &(ps_proc->s_top_left_mb_syntax_ME);

    /* mb part info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;
    mb_part_ctxt skip_mb_part_info;

    /* temp var */
    WORD32 rows_above, rows_below, columns_left, columns_right,u4_use_stat_sad;

    /* Motion vectors in full-pel units */
    WORD16 mv_x, mv_y;

    /* recon stride */
    WORD32 i4_rec_strd = ps_proc->i4_rec_strd;

    /* source buffer for halp pel generation functions */
    UWORD8 *pu1_hpel_src;

    /* quantization parameters */
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[0];

    /* Sad therholds */
    ps_me_ctxt->pu2_sad_thrsh = ps_qp_params->pu2_sad_thrsh;

    /*Best half pel buffer*/
    UWORD8 *pu1_best_subpel_buf = ps_proc->pu1_best_subpel_buf;
    UWORD32 u4_bst_spel_strd = ps_proc->u4_bst_spel_buf_strd;

    /* During evaluation for motion vectors do not search through padded regions */
    /* Obtain number of rows and columns that are effective for computing for me evaluation */
    rows_above = MB_SIZE + ps_proc->i4_mb_y * MB_SIZE;
    rows_below = (ps_proc->i4_ht_mbs - ps_proc->i4_mb_y) * MB_SIZE;
    columns_left = MB_SIZE + ps_proc->i4_mb_x * MB_SIZE;
    columns_right = (ps_proc->i4_wd_mbs - ps_proc->i4_mb_x) * MB_SIZE;

    /* init srch range */
    /* NOTE : For now, lets limit the search range by DEFAULT_MAX_SRCH_RANGE_X / 2
     * on all sides.
     */
//    ps_me_ctxt->i4_srch_range_w = -MIN(columns_left, ps_me_ctxt->ai2_srch_boundaries[0]);
//    ps_me_ctxt->i4_srch_range_e = MIN(columns_right, ps_me_ctxt->ai2_srch_boundaries[0]);
//    ps_me_ctxt->i4_srch_range_n = -MIN(rows_above, ps_me_ctxt->ai2_srch_boundaries[1]);
//    ps_me_ctxt->i4_srch_range_s = MIN(rows_below, ps_me_ctxt->ai2_srch_boundaries[1]);

    ps_me_ctxt->i4_srch_range_w = -MIN(columns_left, DEFAULT_MAX_SRCH_RANGE_X >> 1);
    ps_me_ctxt->i4_srch_range_e = MIN(columns_right, DEFAULT_MAX_SRCH_RANGE_X >> 1);
    ps_me_ctxt->i4_srch_range_n = -MIN(rows_above, DEFAULT_MAX_SRCH_RANGE_Y >> 1);
    ps_me_ctxt->i4_srch_range_s = MIN(rows_below, DEFAULT_MAX_SRCH_RANGE_Y >> 1);

    /* this is to facilitate fast sub pel computation with minimal loads */
    if (ps_me_ctxt->u4_enable_hpel)
    {
        ps_me_ctxt->i4_srch_range_w += 1;
        ps_me_ctxt->i4_srch_range_e -= 1;
        ps_me_ctxt->i4_srch_range_n += 1;
        ps_me_ctxt->i4_srch_range_s -= 1;
    }

    /*Initialize the min sad option*/
    ps_me_ctxt->u4_min_sad_reached  = 0;    /*Not yet found min sad*/
    ps_me_ctxt->i4_min_sad          = ps_proc->ps_cur_mb->u4_min_sad;

    /************************************************************/
    /* Get the seed motion vector candidates                    */
    /************************************************************/
    ih264e_get_search_candidates(ps_proc, ps_me_ctxt);

    /************************************************************/
    /* Init the MB part ctxt structure                          */
    /************************************************************/
    ps_mb_part->s_mv_curr.i2_mvx = 0;
    ps_mb_part->s_mv_curr.i2_mvy = 0;
    ps_mb_part->i4_mb_cost = INT_MAX;
    ps_mb_part->i4_mb_distortion = INT_MAX;

    /* With NMB changes this logic will not work as we cannot exit NME in between*/
    /********************************************************************/
    /*                  Analyse skip                                    */
    /********************************************************************/
//    if (ps_proc->ps_codec->s_cfg.u4_enable_satqd == 0
//                    && u4_frame_level_me == 0)
//    {
//        if ( (ps_proc->ps_ngbr_avbl->u1_mb_a && (ps_me_ctxt->u4_left_is_skip == 1)) ||
//                        (ps_proc->ps_ngbr_avbl->u1_mb_b && ps_top_syn->u2_mb_type == PSKIP) ||
//                        (ps_proc->ps_ngbr_avbl->u1_mb_d && ps_top_left_syn->u2_mb_type == PSKIP) )
//        {
//            if ( 0 == ih264e_analyse_skip(ps_proc, ps_me_ctxt) )
//            {
//                return;
//            }
//        }
//    }

    /********************************************************************/
    /*                  compute skip cost                               */
    /********************************************************************/
    /* See if we need to use modified sad */
    u4_use_stat_sad = (ps_proc->ps_codec->s_cfg.u4_enable_satqd == 1);

    /* init the cost of skip MB */
    skip_mb_part_info.i4_mb_cost = INT_MAX;
    ime_compute_skip_cost(ps_me_ctxt, ps_proc->ps_skip_mv, &skip_mb_part_info, u4_use_stat_sad);


    if (ps_me_ctxt->u4_min_sad_reached == 0)
    {
        /************************************************************/
        /* Evaluate search candidates for initial mv pt.            */
        /************************************************************/
        ime_evaluate_init_srchposn_16x16(ps_me_ctxt);

        /********************************************************************/
        /*                  full pel motion estimation                      */
        /********************************************************************/
        ime_full_pel_motion_estimation_16x16(ps_me_ctxt);

        DEBUG_MV_HISTOGRAM_ADD((ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvx >> 2),
                               (ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvy >> 2));

        DEBUG_SAD_HISTOGRAM_ADD(ps_me_ctxt->s_mb_part.i4_mb_distortion, 1);
        /********************************************************************/
        /*                   sub pel motion estimation                      */
        /********************************************************************/
        if (ps_me_ctxt->u4_enable_hpel)
        {
            /* motion vectors in terms of full pel values */
            mv_x = ps_mb_part->s_mv_curr.i2_mvx >> 2;
            mv_y = ps_mb_part->s_mv_curr.i2_mvy >> 2;

            /* moving src pointer to the converged motion vector location*/
            pu1_hpel_src = ps_me_ctxt->pu1_ref_buf_luma + mv_x + (mv_y * i4_rec_strd);

            ps_me_ctxt->pu1_half_x = ps_proc->pu1_half_x;
            ps_me_ctxt->pu1_half_y = ps_proc->pu1_half_y;
            ps_me_ctxt->pu1_half_xy = ps_proc->pu1_half_xy;
            ps_me_ctxt->u4_hp_buf_strd = HP_BUFF_WD;

            /* half  pel search is done for both sides of full pel,
             * hence half_x of width x height = 17x16 is created
             * starting from left half_x of converged full pel */
            pu1_hpel_src -= 1;

            /* computing half_x */
            ps_codec->pf_ih264e_sixtapfilter_horz(pu1_hpel_src,
                                                  ps_proc->pu1_half_x,
                                                  i4_rec_strd,
                                                  ps_me_ctxt->u4_hp_buf_strd);

            /*
             * Halfpel search is done for both sides of full pel,
             * hence half_y of width x height = 16x17 is created
             * starting from top half_y of converged full pel
             * for half_xy top_left is required
             * hence it starts from pu1_hpel_src = full_pel_converged_point - i4_rec_strd - 1
             */

            pu1_hpel_src -= i4_rec_strd;

            /* computing half_y , and half_xy*/
            ps_codec->pf_ih264e_sixtap_filter_2dvh_vert(
                            pu1_hpel_src, ps_proc->pu1_half_y,
                            ps_proc->pu1_half_xy, i4_rec_strd,
                            ps_me_ctxt->u4_hp_buf_strd, ps_proc->ai16_pred1 + 3,
                            ps_me_ctxt->u4_hp_buf_strd);

            ime_sub_pel_motion_estimation_16x16(ps_me_ctxt);
        }
    }

    {

        /* if skip gives a better cost than other search, copy the cost accordingly*/
        if (skip_mb_part_info.i4_mb_cost < ps_mb_part->i4_mb_cost)
        {
            ps_mb_part->i4_mb_cost = skip_mb_part_info.i4_mb_cost;
            ps_mb_part->i4_mb_distortion = skip_mb_part_info.i4_mb_distortion;
            ps_mb_part->s_mv_curr.i2_mvx = skip_mb_part_info.s_mv_curr.i2_mvx;
            ps_mb_part->s_mv_curr.i2_mvy = skip_mb_part_info.s_mv_curr.i2_mvy;
        }
        else
        {
            /*
             * If the current MB has a sub pel component,
             * we need to copy that to the best subpel buffer
             */
            if (ps_me_ctxt->u4_enable_hpel && ps_mb_part->pu1_best_hpel_buf)
            {
                ps_codec->pf_inter_pred_luma_copy(ps_mb_part->pu1_best_hpel_buf,
                                                  pu1_best_subpel_buf,
                                                  ps_me_ctxt->u4_hp_buf_strd,
                                                  u4_bst_spel_strd, MB_SIZE,
                                                  MB_SIZE, NULL, 0);
            }
        }
    }

    DEBUG_SAD_HISTOGRAM_ADD(ps_me_ctxt->s_mb_part.i4_mb_distortion, 0);

    /* update the type of the mb if necessary */
    if (ps_me_ctxt->s_mb_part.i4_mb_cost < ps_proc->ps_cur_mb->i4_mb_cost)
    {
        /* mb cost */
        ps_proc->ps_cur_mb->i4_mb_cost = ps_me_ctxt->s_mb_part.i4_mb_cost;

        /* mb distortion */
        ps_proc->ps_cur_mb->i4_mb_distortion = ps_me_ctxt->s_mb_part.i4_mb_distortion;

        /* mb type */
        ps_proc->ps_cur_mb->u4_mb_type  = P16x16;
    }

    /* number of partitions */
    ps_proc->u4_num_sub_partitions = 1;
    *(ps_proc->pu4_mb_pu_cnt) = 1;

    /* position in-terms of PU */
    ps_proc->ps_pu->b4_pos_x = 0;
    ps_proc->ps_pu->b4_pos_y = 0;

    /* PU size */
    ps_proc->ps_pu->b4_wd = 3;
    ps_proc->ps_pu->b4_ht = 3;

    /* ref idx */
    ps_proc->ps_pu->i1_l0_ref_idx = 0;

    /* motion vector L0 */
    ps_proc->ps_pu->s_l0_mv.i2_mvx = ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvx;
    ps_proc->ps_pu->s_l0_mv.i2_mvy = ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvy;

    /* Update min sad conditions */
    if (ps_me_ctxt->u4_min_sad_reached == 1)
    {
        ps_proc->ps_cur_mb->u4_min_sad_reached = 1;
        ps_proc->ps_cur_mb->u4_min_sad = ps_me_ctxt->i4_min_sad;
    }
}

/**
*******************************************************************************
*
* @brief This function performs motion estimation for the current NMB
*
* @par Description:
* Intializes input and output pointers required by the function ih264e_compute_me
* and calls the function ih264e_compute_me in a loop to process NMBs.
*
* @param[in] ps_proc
*  Process context corresponding to the job
*
* @returns
*
* @remarks none
*
*******************************************************************************
*/
void ih264e_compute_me_nmb(process_ctxt_t *ps_proc, UWORD32 u4_nmb_count)
{
    /* pic pu */
    enc_pu_t *ps_pu_begin = ps_proc->ps_pu;

    /* ME map */
    UWORD8 *pu1_me_map = ps_proc->pu1_me_map + (ps_proc->i4_mb_y * ps_proc->i4_wd_mbs);

    /* temp var */
    UWORD32 u4_i;

    ps_proc->s_me_ctxt.u4_left_is_intra = ps_proc->s_left_mb_syntax_ele.u2_is_intra;
    ps_proc->s_me_ctxt.u4_left_is_skip = (ps_proc->s_left_mb_syntax_ele.u2_mb_type == PSKIP);

    for (u4_i = 0; u4_i < u4_nmb_count; u4_i++)
    {
        /* Wait for ME map */
        if (ps_proc->i4_mb_y > 0)
        {
            /* Wait for top right ME to be done */
            UWORD8 *pu1_me_map_tp_rw = ps_proc->pu1_me_map + (ps_proc->i4_mb_y - 1) * ps_proc->i4_wd_mbs;

            while (1)
            {
                volatile UWORD8 *pu1_buf;
                WORD32 idx = ps_proc->i4_mb_x + u4_i + 1;

                idx = MIN(idx, (ps_proc->i4_wd_mbs - 1));
                pu1_buf =  pu1_me_map_tp_rw + idx;
                if(*pu1_buf)
                    break;
                ithread_yield();
            }
        }

        ps_proc->ps_skip_mv = &(ps_proc->ps_nmb_info[u4_i].s_skip_mv);
        ps_proc->ps_ngbr_avbl = &(ps_proc->ps_nmb_info[u4_i].s_ngbr_avbl);
        ps_proc->ps_pred_mv = &(ps_proc->ps_nmb_info[u4_i].s_pred_mv);

        ps_proc->ps_cur_mb = &(ps_proc->ps_nmb_info[u4_i]);

        ps_proc->ps_cur_mb->u4_min_sad = ps_proc->u4_min_sad;
        ps_proc->ps_cur_mb->u4_min_sad_reached = 0;

        ps_proc->ps_cur_mb->i4_mb_cost = INT_MAX;
        ps_proc->ps_cur_mb->i4_mb_distortion = SHRT_MAX;

        /* Set the best subpel buf to the correct mb so that the buffer can be copied */
        ps_proc->pu1_best_subpel_buf = ps_proc->ps_nmb_info[u4_i].pu1_best_sub_pel_buf;
        ps_proc->u4_bst_spel_buf_strd = ps_proc->ps_nmb_info[u4_i].u4_bst_spel_buf_strd;

        /* Set the min sad conditions */
        ps_proc->ps_cur_mb->u4_min_sad = ps_proc->ps_codec->u4_min_sad;
        ps_proc->ps_cur_mb->u4_min_sad_reached = 0;

        /* Derive neighbor availability for the current macroblock */
        ih264e_derive_nghbr_avbl_of_mbs(ps_proc);

        /* init me */
        ih264e_init_me(ps_proc);

        ih264e_compute_me(ps_proc);

        /* update top and left structs */
        {
            mb_info_t *ps_top_syn = ps_proc->ps_top_row_mb_syntax_ele + ps_proc->i4_mb_x;
            mb_info_t *ps_top_left_syn = &(ps_proc->s_top_left_mb_syntax_ME);
            enc_pu_t *ps_left_mb_pu = &ps_proc->s_left_mb_pu_ME;
            enc_pu_t *ps_top_left_mb_pu = &ps_proc->s_top_left_mb_pu_ME;
            enc_pu_t *ps_top_mv = ps_proc->ps_top_row_pu_ME + ps_proc->i4_mb_x;

            *ps_top_left_syn = *ps_top_syn;

            *ps_top_left_mb_pu = *ps_top_mv;
            *ps_left_mb_pu = *ps_proc->ps_pu;
        }

        ps_proc->ps_pu += *ps_proc->pu4_mb_pu_cnt;

        /* Copy the min sad reached info */
        ps_proc->ps_nmb_info[u4_i].u4_min_sad_reached = ps_proc->ps_cur_mb->u4_min_sad_reached;
        ps_proc->ps_nmb_info[u4_i].u4_min_sad   = ps_proc->ps_cur_mb->u4_min_sad;

        /*
         * To make sure that the MV map is properly sync to the
         * cache we need to do a DDB
         */
        {
            DATA_SYNC();

            pu1_me_map[ps_proc->i4_mb_x] = 1;
        }
        ps_proc->i4_mb_x++;

        ps_proc->s_me_ctxt.u4_left_is_intra = 0;
        ps_proc->s_me_ctxt.u4_left_is_skip = (ps_proc->ps_cur_mb->u4_mb_type  == PSKIP);

        /* update buffers pointers */
        ps_proc->pu1_src_buf_luma += MB_SIZE;
        ps_proc->pu1_rec_buf_luma += MB_SIZE;
        ps_proc->pu1_ref_buf_luma += MB_SIZE;

        /*
         * Note: Although chroma mb size is 8, as the chroma buffers are interleaved,
         * the stride per MB is MB_SIZE
         */
        ps_proc->pu1_src_buf_chroma += MB_SIZE;
        ps_proc->pu1_rec_buf_chroma += MB_SIZE;
        ps_proc->pu1_ref_buf_chroma += MB_SIZE;

        ps_proc->pu4_mb_pu_cnt += 1;
    }


    ps_proc->ps_pu = ps_pu_begin;
    ps_proc->i4_mb_x = ps_proc->i4_mb_x - u4_nmb_count;

    /* update buffers pointers */
    ps_proc->pu1_src_buf_luma -= MB_SIZE * u4_nmb_count;
    ps_proc->pu1_rec_buf_luma -= MB_SIZE * u4_nmb_count;
    ps_proc->pu1_ref_buf_luma -= MB_SIZE * u4_nmb_count;

    /*
     * Note: Although chroma mb size is 8, as the chroma buffers are interleaved,
     * the stride per MB is MB_SIZE
     */
    ps_proc->pu1_src_buf_chroma -= MB_SIZE * u4_nmb_count;
    ps_proc->pu1_rec_buf_chroma -= MB_SIZE * u4_nmb_count;
    ps_proc->pu1_ref_buf_chroma -= MB_SIZE * u4_nmb_count;

    ps_proc->pu4_mb_pu_cnt -= u4_nmb_count;
}
