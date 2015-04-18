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
 *
 *
 * @author
 *  Ittiam
 *
 * @par List of Functions:
 *  -
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
#include <string.h>

/* User include files */
#include "ime_typedefs.h"
#include "ime_distortion_metrics.h"
#include "ime_structs.h"
#include "ime_defs.h"
#include "ime_macros.h"
#include "ime.h"
#include "ime_statistics.h"

/**
*******************************************************************************
*
* @brief Diamond Search
*
* @par Description:
*  This function computes the sad at vertices of several layers of diamond grid
*  at a time. The number of layers of diamond grid that would be evaluated is
*  configurable.The function computes the sad at vertices of a diamond grid. If
*  the sad at the center of the diamond grid is lesser than the sad at any other
*  point of the diamond grid, the function marks the candidate Mb partition as
*  mv.
*
* @param[in] ps_mb_part
*  pointer to current mb partition ctxt with respect to ME
*
* @param[in] ps_me_ctxt
*  pointer to me context
*
* @param[in] u4_lambda_motion
*  lambda motion
*
* @param[in] u4_enable_fast_sad
*  enable/disable fast sad computation
*
* @returns  mv pair & corresponding distortion and cost
*
* @remarks Diamond Srch, radius is 1
*
*******************************************************************************
*/
void ime_diamond_search_16x16(me_ctxt_t *ps_me_ctxt)
{
    /* MB partition info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;

    /* lagrange parameter */
    UWORD32 u4_lambda_motion = ps_me_ctxt->u4_lambda_motion;

    /* srch range*/
    WORD32 i4_srch_range_n = ps_me_ctxt->i4_srch_range_n;
    WORD32 i4_srch_range_s = ps_me_ctxt->i4_srch_range_s;
    WORD32 i4_srch_range_e = ps_me_ctxt->i4_srch_range_e;
    WORD32 i4_srch_range_w = ps_me_ctxt->i4_srch_range_w;

    /* enabled fast sad computation */
//    UWORD32 u4_enable_fast_sad = ps_me_ctxt->u4_enable_fast_sad;

    /* pointer to src macro block */
    UWORD8 *pu1_curr_mb = ps_me_ctxt->pu1_src_buf_luma;
    UWORD8 *pu1_ref_mb = ps_me_ctxt->pu1_ref_buf_luma;

    /* strides */
    WORD32 i4_src_strd = ps_me_ctxt->i4_src_strd;
    WORD32 i4_ref_strd = ps_me_ctxt->i4_rec_strd;

    /* least cost */
    WORD32 i4_cost_least = ps_mb_part->i4_mb_cost;

    /* least sad */
    WORD32 i4_distortion_least = ps_mb_part->i4_mb_distortion;

    /* mv pair */
    WORD16 i2_mvx, i2_mvy;

    /* mv bits */
    UWORD8 *pu1_mv_bits = ps_me_ctxt->pu1_mv_bits;

    /* temp var */
    WORD32 i4_cost[4];
    WORD32 i4_sad[4];
    UWORD8 *pu1_ref;
    WORD16 i2_mv_u_x, i2_mv_u_y;

    /* Diamond search Iteration Max Cnt */
    UWORD32 u4_num_layers = ps_me_ctxt->u4_num_layers;

    /* temp var */
//    UWORD8 u1_prev_jump = NONE;
//    UWORD8 u1_curr_jump = NONE;
//    UWORD8 u1_next_jump;
//    WORD32 mask_arr[5] = {15, 13, 14, 7, 11};
//    WORD32 mask;
//    UWORD8 *apu1_ref[4];
//    WORD32 i, cnt;
//    WORD32 dia[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    /* mv with best sad during initial evaluation */
    i2_mvx = ps_mb_part->s_mv_curr.i2_mvx;
    i2_mvy = ps_mb_part->s_mv_curr.i2_mvy;

    i2_mv_u_x = i2_mvx;
    i2_mv_u_y = i2_mvy;

    while (u4_num_layers--)
    {
        /* FIXME : is this the write way to check for out of bounds ? */
        if ( (i2_mvx - 1 < i4_srch_range_w) ||
                        (i2_mvx + 1 > i4_srch_range_e) ||
                        (i2_mvy - 1 < i4_srch_range_n) ||
                        (i2_mvy + 1 > i4_srch_range_s) )
        {
            break;
        }

        pu1_ref = pu1_ref_mb + i2_mvx + (i2_mvy * i4_ref_strd);

        ps_me_ctxt->pf_ime_compute_sad4_diamond(pu1_ref,
                                                pu1_curr_mb,
                                                i4_ref_strd,
                                                i4_src_strd,
                                                i4_sad);

        DEBUG_SAD_HISTOGRAM_ADD(i4_sad[0], 2);
        DEBUG_SAD_HISTOGRAM_ADD(i4_sad[1], 2);
        DEBUG_SAD_HISTOGRAM_ADD(i4_sad[2], 2);
        DEBUG_SAD_HISTOGRAM_ADD(i4_sad[3], 2);

        /* compute cost */
        i4_cost[0] = i4_sad[0] + u4_lambda_motion * ( pu1_mv_bits[ ((i2_mvx - 1) << 2) - ps_mb_part->s_mv_pred.i2_mvx]
                                                                   + pu1_mv_bits[(i2_mvy << 2) - ps_mb_part->s_mv_pred.i2_mvy] );
        i4_cost[1] = i4_sad[1] + u4_lambda_motion * ( pu1_mv_bits[ ((i2_mvx + 1) << 2) - ps_mb_part->s_mv_pred.i2_mvx]
                                                                   + pu1_mv_bits[(i2_mvy << 2) - ps_mb_part->s_mv_pred.i2_mvy] );
        i4_cost[2] = i4_sad[2] + u4_lambda_motion * ( pu1_mv_bits[ (i2_mvx << 2) - ps_mb_part->s_mv_pred.i2_mvx]
                                                                   + pu1_mv_bits[((i2_mvy - 1) << 2) - ps_mb_part->s_mv_pred.i2_mvy] );
        i4_cost[3] = i4_sad[3] + u4_lambda_motion * ( pu1_mv_bits[ (i2_mvx << 2) - ps_mb_part->s_mv_pred.i2_mvx]
                                                                   + pu1_mv_bits[((i2_mvy + 1) << 2) - ps_mb_part->s_mv_pred.i2_mvy] );


        if (i4_cost_least > i4_cost[0])
        {
            i4_cost_least = i4_cost[0];
            i4_distortion_least = i4_sad[0];

            i2_mv_u_x = (i2_mvx - 1);
            i2_mv_u_y = i2_mvy;
        }

        if (i4_cost_least > i4_cost[1])
        {
            i4_cost_least = i4_cost[1];
            i4_distortion_least = i4_sad[1];

            i2_mv_u_x = (i2_mvx + 1);
            i2_mv_u_y = i2_mvy;
        }

        if (i4_cost_least > i4_cost[2])
        {
            i4_cost_least = i4_cost[2];
            i4_distortion_least = i4_sad[2];

            i2_mv_u_x = i2_mvx;
            i2_mv_u_y = i2_mvy - 1;
        }

        if (i4_cost_least > i4_cost[3])
        {
            i4_cost_least = i4_cost[3];
            i4_distortion_least = i4_sad[3];

            i2_mv_u_x = i2_mvx;
            i2_mv_u_y = i2_mvy + 1;
        }

        if( (i2_mv_u_x == i2_mvx) && (i2_mv_u_y == i2_mvy))
        {
            ps_mb_part->u4_exit = 1;
            break;
        }
        else
        {
            i2_mvx = i2_mv_u_x;
            i2_mvy = i2_mv_u_y;
        }


    }

    if (i4_cost_least < ps_mb_part->i4_mb_cost)
    {
        ps_mb_part->i4_mb_cost = i4_cost_least;
        ps_mb_part->i4_mb_distortion = i4_distortion_least;
        ps_mb_part->s_mv_curr.i2_mvx = i2_mvx;
        ps_mb_part->s_mv_curr.i2_mvy = i2_mvy;
    }

}


/**
*******************************************************************************
*
* @brief This function computes the best motion vector among the tentative mv
* candidates chosen.
*
* @par Description:
*  This function determines the position in the search window at which the motion
*  estimation should begin in order to minimise the number of search iterations.
*
* @param[in] ps_mb_part
*  pointer to current mb partition ctxt with respect to ME
*
* @param[in] u4_lambda_motion
*  lambda motion
*
* @param[in] u4_fast_flag
*  enable/disable fast sad computation
*
* @returns  mv pair & corresponding distortion and cost
*
* @remarks none
*
*******************************************************************************
*/
void ime_evaluate_init_srchposn_16x16
        (
            me_ctxt_t *ps_me_ctxt
        )
{
    UWORD32 u4_lambda_motion = ps_me_ctxt->u4_lambda_motion;

    /* candidate mv cnt */
    UWORD32 u4_num_candidates = ps_me_ctxt->u4_num_candidates;

    /* list of candidate mvs */
    ime_mv_t *ps_mv_list = ps_me_ctxt->as_mv_init_search;

    /* pointer to src macro block */
    UWORD8 *pu1_curr_mb = ps_me_ctxt->pu1_src_buf_luma;
    UWORD8 *pu1_ref_mb = ps_me_ctxt->pu1_ref_buf_luma;

    /* strides */
    WORD32 i4_src_strd = ps_me_ctxt->i4_src_strd;
    WORD32 i4_ref_strd = ps_me_ctxt->i4_rec_strd;

    /* enabled fast sad computation */
    UWORD32 u4_enable_fast_sad = ps_me_ctxt->u4_enable_fast_sad;

    /* SAD(distortion metric) of an 8x8 block */
    WORD32 i4_mb_distortion;

    /* cost = distortion + u4_lambda_motion * rate */
    WORD32 i4_mb_cost, i4_mb_cost_least = INT_MAX, i4_distortion_least = INT_MAX;

    /* mb partitions info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;

    /* mv bits */
    UWORD8 *pu1_mv_bits = ps_me_ctxt->pu1_mv_bits;

    /* temp var */
    UWORD32  i, j, u4_srch_pos_idx = 0;
    UWORD8 *pu1_ref = NULL;
    WORD16 mv_x, mv_y;

    if (0)
    {
        /************************************************************/
        /* Compute SKIP Cost                                        */
        /************************************************************/
        mv_x = ps_mv_list[SKIP_CAND].i2_mvx;
        mv_y = ps_mv_list[SKIP_CAND].i2_mvy;

        /* adjust ref pointer */
        pu1_ref = pu1_ref_mb + mv_x + (mv_y * i4_ref_strd);

        /* compute distortion */
        ps_me_ctxt->pf_ime_compute_sad_16x16[u4_enable_fast_sad](pu1_curr_mb, pu1_ref, i4_src_strd, i4_ref_strd, i4_mb_cost_least, &i4_mb_distortion);

        /* for skip mode cost & distortion are identical
         * But we shall add a bias to favor skip mode.
         * Doc. JVT B118 Suggests SKIP_BIAS as 16.
         * TODO : Empirical analysis of SKIP_BIAS is necessary */

        i4_distortion_least = i4_mb_distortion;

        u4_srch_pos_idx = 0;

#define SKIP_BIAS 8

        i4_mb_cost_least = i4_mb_distortion - (u4_lambda_motion * SKIP_BIAS);

#undef SKIP_BIAS
    }


    /* Carry out a search using each of the motion vector pairs identified above as predictors. */
    /* TODO : Just like Skip, Do we need to add any bias to zero mv as well */
    for(i = 0; i < u4_num_candidates; i++)
    {
        /* compute sad */
        WORD32 c_sad = 1;

        for(j = 0; j < i; j++ )
        {
            if ( (ps_mv_list[i].i2_mvx == ps_mv_list[j].i2_mvx) &&
                            (ps_mv_list[i].i2_mvy == ps_mv_list[j].i2_mvy) )
            {
                c_sad = 0;
                break;
            }
        }
        if(c_sad)
        {
            /* adjust ref pointer */
            pu1_ref = pu1_ref_mb + ps_mv_list[i].i2_mvx + (ps_mv_list[i].i2_mvy * i4_ref_strd);

            /* compute distortion */
            ps_me_ctxt->pf_ime_compute_sad_16x16[u4_enable_fast_sad](pu1_curr_mb, pu1_ref, i4_src_strd, i4_ref_strd, i4_mb_cost_least, &i4_mb_distortion);
            DEBUG_SAD_HISTOGRAM_ADD(i4_mb_distortion, 3);
            /* compute cost */
            i4_mb_cost = i4_mb_distortion + u4_lambda_motion * ( pu1_mv_bits[ (ps_mv_list[i].i2_mvx << 2) - ps_mb_part->s_mv_pred.i2_mvx]
                            + pu1_mv_bits[(ps_mv_list[i].i2_mvy << 2) - ps_mb_part->s_mv_pred.i2_mvy] );

            if (i4_mb_cost < i4_mb_cost_least)
            {
                i4_mb_cost_least = i4_mb_cost;

                i4_distortion_least = i4_mb_distortion;

                u4_srch_pos_idx = i;
            }
        }
    }

    if (i4_mb_cost_least < ps_mb_part->i4_mb_cost)
    {
        ps_mb_part->u4_srch_pos_idx = u4_srch_pos_idx;
        ps_mb_part->i4_mb_cost = i4_mb_cost_least;
        ps_mb_part->i4_mb_distortion = i4_distortion_least;
        ps_mb_part->s_mv_curr.i2_mvx = ps_mv_list[u4_srch_pos_idx].i2_mvx;
        ps_mb_part->s_mv_curr.i2_mvy = ps_mv_list[u4_srch_pos_idx].i2_mvy;
    }
}


/**
*******************************************************************************
*
* @brief Searches for the best matching full pixel predictor within the search
* range
*
* @par Description:
*  This function begins by computing the mv predict vector for the current mb.
*  This is used for cost computations. Further basing on the algo. chosen, it
*  looks through a set of candidate vectors that best represent the mb a least
*  cost and returns this information.
*
* @param[in] ps_proc
*  pointer to current proc ctxt
*
* @param[in] ps_me_ctxt
*  pointer to me context
*
* @returns  mv pair & corresponding distortion and cost
*
* @remarks none
*
*******************************************************************************
*/
void ime_full_pel_motion_estimation_16x16
    (
        me_ctxt_t *ps_me_ctxt
    )
{
    /* mb part info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;

    /******************************************************************/
    /* Modify Search range about initial candidate instead of zero mv */
    /******************************************************************/
    /*
     * FIXME: The motion vectors in a way can become unbounded. It may so happen that
     * MV might exceed the limit of the profile configured.
     */
    ps_me_ctxt->i4_srch_range_w = MAX(ps_me_ctxt->i4_srch_range_w,
                                      -ps_me_ctxt->ai2_srch_boundaries[0] + ps_mb_part->s_mv_curr.i2_mvx);
    ps_me_ctxt->i4_srch_range_e = MIN(ps_me_ctxt->i4_srch_range_e,
                                       ps_me_ctxt->ai2_srch_boundaries[0] + ps_mb_part->s_mv_curr.i2_mvx);
    ps_me_ctxt->i4_srch_range_n = MAX(ps_me_ctxt->i4_srch_range_n,
                                      -ps_me_ctxt->ai2_srch_boundaries[1] + ps_mb_part->s_mv_curr.i2_mvy);
    ps_me_ctxt->i4_srch_range_s = MIN(ps_me_ctxt->i4_srch_range_s,
                                       ps_me_ctxt->ai2_srch_boundaries[1] + ps_mb_part->s_mv_curr.i2_mvy);

    /************************************************************/
    /* Traverse about best initial candidate for mv             */
    /************************************************************/

    switch (ps_me_ctxt->u4_me_speed_preset)
    {
        case DMND_SRCH:
            ime_diamond_search_16x16(ps_me_ctxt);
            break;
        default:
            assert(0);
            break;
    }

    ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvx = ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvx << 2;
    ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvy = ps_me_ctxt->s_mb_part.s_mv_curr.i2_mvy << 2;

}


/**
*******************************************************************************
*
* @brief Searches for the best matching sub pixel predictor within the search
* range
*
* @par Description:
*  This function begins by searching across all sub pixel sample points
*  around the full pel motion vector. The vector with least cost is chosen as
*  the mv for the current mb. If the skip mode is not evaluated while analysing
*  the initial search candidates then analyse it here and update the mv.
*
* @param[in] ps_proc
*  pointer to current proc ctxt
*
* @param[in] ps_me_ctxt
*  pointer to me context
*
* @returns none
*
* @remarks none
*
*******************************************************************************
*/
void ime_sub_pel_motion_estimation_16x16
    (
        me_ctxt_t *ps_me_ctxt
    )
{
    /* pointers to src & ref macro block */
    UWORD8 *pu1_curr_mb = ps_me_ctxt->pu1_src_buf_luma;


    /* pointers to ref. half pel planes */
    UWORD8 *pu1_ref_mb_half_x;
    UWORD8 *pu1_ref_mb_half_y;
    UWORD8 *pu1_ref_mb_half_xy;

    /* pointers to ref. half pel planes */
    UWORD8 *pu1_ref_mb_half_x_temp;
    UWORD8 *pu1_ref_mb_half_y_temp;
    UWORD8 *pu1_ref_mb_half_xy_temp;

    /* strides */
    WORD32 i4_src_strd = ps_me_ctxt->i4_src_strd;

    WORD32 i4_ref_strd = ps_me_ctxt->u4_hp_buf_strd;

    /* mb partitions info */
    mb_part_ctxt *ps_mb_part = &ps_me_ctxt->s_mb_part;

    /* SAD(distortion metric) of an mb */
    WORD32 i4_mb_distortion;
    WORD32 i4_distortion_least = ps_mb_part->i4_mb_distortion;

    /* cost = distortion + u4_lambda_motion * rate */
    WORD32 i4_mb_cost;
    WORD32 i4_mb_cost_least = ps_mb_part->i4_mb_cost;

    /*Best half pel buffer*/
    UWORD8 *pu1_best_hpel_buf = NULL;


    /* mv bits */
    UWORD8 *pu1_mv_bits = ps_me_ctxt->pu1_mv_bits;

    /* Motion vectors in full-pel units */
    WORD16 mv_x, mv_y;

    /* lambda - lagrange constant */
    UWORD32 u4_lambda_motion = ps_me_ctxt->u4_lambda_motion;

    /* Flags to check if half pel points needs to be evaluated */
    /**************************************/
    /* 1 bit for each half pel candidate  */
    /* bit 0 - half x = 1, half y = 0     */
    /* bit 1 - half x = -1, half y = 0    */
    /* bit 2 - half x = 0, half y = 1     */
    /* bit 3 - half x = 0, half y = -1    */
    /* bit 4 - half x = 1, half y = 1     */
    /* bit 5 - half x = -1, half y = 1    */
    /* bit 6 - half x = 1, half y = -1    */
    /* bit 7 - half x = -1, half y = -1   */
    /**************************************/
    /* temp var */
    WORD16 i2_mv_u_x, i2_mv_u_y;
    WORD32 i, j;
    WORD32 ai4_sad[8];

    i2_mv_u_x = ps_mb_part->s_mv_curr.i2_mvx;
    i2_mv_u_y = ps_mb_part->s_mv_curr.i2_mvy;

    /************************************************************/
    /* Evaluate half pel                                        */
    /************************************************************/
    mv_x = ps_mb_part->s_mv_curr.i2_mvx >> 2;
    mv_y = ps_mb_part->s_mv_curr.i2_mvy >> 2;


    /**************************************************************/
    /* ps_me_ctxt->pu1_half_x points to the half pel pixel on the */
    /* left side of full pel                                      */
    /* ps_me_ctxt->pu1_half_y points to the half pel pixel on the */
    /* top  side of full pel                                      */
    /* ps_me_ctxt->pu1_half_xy points to the half pel pixel       */
    /* on the top left side of full pel                           */
    /* for the function pf_ime_sub_pel_compute_sad_16x16 the      */
    /* default postions are                                       */
    /* ps_me_ctxt->pu1_half_x = right halp_pel                    */
    /*  ps_me_ctxt->pu1_half_y = bottom halp_pel                  */
    /*  ps_me_ctxt->pu1_half_xy = bottom right halp_pel           */
    /* Hence corresponding adjustments made here                  */
    /**************************************************************/

    pu1_ref_mb_half_x_temp = pu1_ref_mb_half_x = ps_me_ctxt->pu1_half_x + 1;
    pu1_ref_mb_half_y_temp = pu1_ref_mb_half_y = ps_me_ctxt->pu1_half_y + 1 + i4_ref_strd;
    pu1_ref_mb_half_xy_temp = pu1_ref_mb_half_xy = ps_me_ctxt->pu1_half_xy + 1 + i4_ref_strd;


    ps_me_ctxt->pf_ime_sub_pel_compute_sad_16x16(pu1_curr_mb, pu1_ref_mb_half_x,
                                                 pu1_ref_mb_half_y,
                                                 pu1_ref_mb_half_xy,
                                                 i4_src_strd, i4_ref_strd,
                                                 ai4_sad);

    /* Half x plane */
    for(i = 0; i < 2; i++)
    {
        WORD32 mv_x_tmp = (mv_x << 2) + 2;
        WORD32 mv_y_tmp = (mv_y << 2);

        mv_x_tmp -= (i * 4);

        i4_mb_distortion = ai4_sad[i];

        /* compute cost */
        i4_mb_cost = i4_mb_distortion + u4_lambda_motion * ( pu1_mv_bits[ mv_x_tmp - ps_mb_part->s_mv_pred.i2_mvx]
                        + pu1_mv_bits[mv_y_tmp - ps_mb_part->s_mv_pred.i2_mvy] );

        if (i4_mb_cost < i4_mb_cost_least)
        {
            i4_mb_cost_least = i4_mb_cost;

            i4_distortion_least = i4_mb_distortion;

            i2_mv_u_x = mv_x_tmp;

            i2_mv_u_y = mv_y_tmp;

#ifndef HP_PL /*choosing whether left or right half_x*/
            ps_me_ctxt->pu1_half_x = pu1_ref_mb_half_x_temp - i;
            pu1_best_hpel_buf = pu1_ref_mb_half_x_temp - i;
#endif
        }

    }

    /* Half y plane */
    for(i = 0; i < 2; i++)
    {
        WORD32 mv_x_tmp = (mv_x << 2);
        WORD32 mv_y_tmp = (mv_y << 2) + 2;

        mv_y_tmp -= (i * 4);

        i4_mb_distortion = ai4_sad[2 + i];

        /* compute cost */
        i4_mb_cost = i4_mb_distortion + u4_lambda_motion * ( pu1_mv_bits[ mv_x_tmp - ps_mb_part->s_mv_pred.i2_mvx]
                        + pu1_mv_bits[mv_y_tmp - ps_mb_part->s_mv_pred.i2_mvy] );

        if (i4_mb_cost < i4_mb_cost_least)
        {
            i4_mb_cost_least = i4_mb_cost;

            i4_distortion_least = i4_mb_distortion;

            i2_mv_u_x = mv_x_tmp;

            i2_mv_u_y = mv_y_tmp;

#ifndef HP_PL/*choosing whether top or bottom half_y*/
            ps_me_ctxt->pu1_half_y = pu1_ref_mb_half_y_temp  - i*(i4_ref_strd);
            pu1_best_hpel_buf = pu1_ref_mb_half_y_temp  - i*(i4_ref_strd);
#endif
        }

    }

    /* Half xy plane */
    for(j = 0; j < 2; j++)
    {
        for(i = 0; i < 2; i++)
        {
            WORD32 mv_x_tmp = (mv_x << 2) + 2;
            WORD32 mv_y_tmp = (mv_y << 2) + 2;

            mv_x_tmp -= (i * 4);
            mv_y_tmp -= (j * 4);

            i4_mb_distortion = ai4_sad[4 + i + 2 * j];

            /* compute cost */
            i4_mb_cost = i4_mb_distortion + u4_lambda_motion * ( pu1_mv_bits[ mv_x_tmp - ps_mb_part->s_mv_pred.i2_mvx]
                            + pu1_mv_bits[mv_y_tmp - ps_mb_part->s_mv_pred.i2_mvy] );

            if (i4_mb_cost < i4_mb_cost_least)
            {
                i4_mb_cost_least = i4_mb_cost;

                i4_distortion_least = i4_mb_distortion;

                i2_mv_u_x = mv_x_tmp;

                i2_mv_u_y = mv_y_tmp;

#ifndef HP_PL /*choosing between four half_xy */
                ps_me_ctxt->pu1_half_xy = pu1_ref_mb_half_xy_temp  - j*(i4_ref_strd) - i;
                pu1_best_hpel_buf =  pu1_ref_mb_half_xy_temp  - j*(i4_ref_strd) - i;
#endif
            }

        }
    }

    ps_mb_part->i4_mb_cost = i4_mb_cost_least;
    ps_mb_part->i4_mb_distortion = i4_distortion_least;
    ps_mb_part->s_mv_curr.i2_mvx = i2_mv_u_x;
    ps_mb_part->s_mv_curr.i2_mvy = i2_mv_u_y;
    ps_mb_part->pu1_best_hpel_buf = pu1_best_hpel_buf;

}


/**
*******************************************************************************
*
* @brief This function computes cost of skip macroblocks
*
* @par Description:
*
* @param[in] ps_me_ctxt
*  pointer to me ctxt
*
* @param[in] ps_skip_mv
*  pointer to skip mv
*
* @returns  none
*
* @remarks
* NOTE: while computing the skip cost, do not enable early exit from compute
* sad function because, a negative bias gets added later
*
*******************************************************************************
*/
void ime_compute_skip_cost
    (
         me_ctxt_t *ps_me_ctxt,
         void *pv_skip_mv,
         mb_part_ctxt *ps_smb_part_info,
         UWORD32 u4_use_stat_sad
    )
{

    /* pointers to src & ref macro block */
    UWORD8 *pu1_curr_mb = ps_me_ctxt->pu1_src_buf_luma;
    UWORD8 *pu1_ref_mb = ps_me_ctxt->pu1_ref_buf_luma;

    /* strides */
    WORD32 i4_src_strd = ps_me_ctxt->i4_src_strd;
    WORD32 i4_ref_strd = ps_me_ctxt->i4_rec_strd;

    /* enabled fast sad computation */
    UWORD32 u4_enable_fast_sad = ps_me_ctxt->u4_enable_fast_sad;

    /* SAD(distortion metric) of an mb */
    WORD32 i4_mb_distortion;

    /* cost = distortion + u4_lambda_motion * rate */
    WORD32 i4_mb_cost;

    /* Motion vectors in full-pel units */
    WORD16 mv_x, mv_y;

    /* lambda - lagrange constant */
    UWORD32 u4_lambda_motion = ps_me_ctxt->u4_lambda_motion;

    /* skip mv */
    ime_mv_t *ps_skip_mv = pv_skip_mv, s_clip_skip_mv;

    /* temp var */
    UWORD8 *pu1_ref = NULL;
    UWORD32 u4_is_nonzero;

    s_clip_skip_mv.i2_mvx = CLIP3(ps_me_ctxt->i4_srch_range_w, ps_me_ctxt->i4_srch_range_e, ps_skip_mv->i2_mvx);
    s_clip_skip_mv.i2_mvy = CLIP3(ps_me_ctxt->i4_srch_range_n, ps_me_ctxt->i4_srch_range_s, ps_skip_mv->i2_mvy);

    if ((s_clip_skip_mv.i2_mvx != ps_skip_mv->i2_mvx) ||
                    (s_clip_skip_mv.i2_mvy != ps_skip_mv->i2_mvy))
    {
        /* skip motion vector not with in bounds */
        /* it is possible that mv is already evaluated */
        return ;
    }

    mv_x = (ps_skip_mv->i2_mvx + 2) >> 2;
    mv_y = (ps_skip_mv->i2_mvy + 2) >> 2;

    if ((mv_x << 2) != ps_skip_mv->i2_mvx || (mv_y << 2) != ps_skip_mv->i2_mvy)
    {


        return ;


    }
    else
    {
        /* adjust ref pointer */
        pu1_ref = pu1_ref_mb + mv_x + (mv_y * i4_ref_strd);
    }

    if(u4_use_stat_sad == 1)
    {
        ps_me_ctxt->pf_ime_compute_sad_stat_luma_16x16(pu1_curr_mb, pu1_ref, i4_src_strd, i4_ref_strd,
                ps_me_ctxt->pu2_sad_thrsh, &i4_mb_distortion,&u4_is_nonzero);

        /*
         *NOTE The check here is two fold
         * One is checking if the sad has been reached, ie min sad, which a configurable parameter
         * If that is reached,we need not do any mode evaluation
         * Similary if we find a distortion of zero there is no point of doing any further mode evaluation
         * as sad is a non negative quantity
         * hence in this case too, no further evaluation is necessary
         */
        /*
         *NOTE in case we need to disable the zero check using satdq,
         *  we need only to set the u4_is_zero to a non zero value
         */
        if(u4_is_nonzero==0 || i4_mb_distortion <= ps_me_ctxt->i4_min_sad)
        {
            ps_me_ctxt->u4_min_sad_reached = 1;    /* found min sad*/
            ps_me_ctxt->i4_min_sad =  (u4_is_nonzero == 0)?0:i4_mb_distortion;
        }
    }
    else
    {
        ps_me_ctxt->pf_ime_compute_sad_16x16[u4_enable_fast_sad](pu1_curr_mb, pu1_ref, i4_src_strd, i4_ref_strd, INT_MAX, &i4_mb_distortion);

        if(i4_mb_distortion <= ps_me_ctxt->i4_min_sad)
        {
            ps_me_ctxt->i4_min_sad = i4_mb_distortion;
            ps_me_ctxt->u4_min_sad_reached = 1;    /* found min sad*/
        }
    }

    /* for skip mode cost & distortion are identical
     * But we shall add a bias to favor skip mode.
     * Doc. JVT B118 Suggests SKIP_BIAS as 16.
     * TODO : Empirical analysis of SKIP_BIAS is necessary */
#define SKIP_BIAS 8
    i4_mb_cost = i4_mb_distortion - (u4_lambda_motion * SKIP_BIAS);
#undef SKIP_BIAS

    if (i4_mb_cost <= ps_smb_part_info->i4_mb_cost)
    {
        ps_smb_part_info->i4_mb_cost = i4_mb_cost;
        ps_smb_part_info->i4_mb_distortion = i4_mb_distortion;
        ps_smb_part_info->s_mv_curr.i2_mvx = ps_skip_mv->i2_mvx;
        ps_smb_part_info->s_mv_curr.i2_mvy = ps_skip_mv->i2_mvy;
    }
}

