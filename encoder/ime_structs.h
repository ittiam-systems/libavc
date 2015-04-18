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
 *  ih264e_me.h
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

#ifndef _IME_STRUCTS_H_
#define _IME_STRUCTS_H_

/**
 * Motion vector
 */
typedef struct
{
    /**
     * Horizontal Motion Vector
     */
    WORD16 i2_mvx;

    /**
     * Vertical Motion Vector
     */
    WORD16 i2_mvy;
} ime_mv_t;


/**
**************************************************************************
*   @brief   mb_part_ctxt
*
*   Structure that would hold the information for individual MB partitions
*   gathered during the full pel ME stage
**************************************************************************
*/
typedef struct
{
    /**
     * best mvs
     */
    ime_mv_t  s_mv_curr;

    /**
     * mv predictor
     */
    ime_mv_t  s_mv_pred;

    /**
     * SAD associated with the MB partition
     */
    WORD32 i4_mb_distortion;

    /**
     * cost for the MB partition
     */
    WORD32 i4_mb_cost;

    /**
     * Search position for least cost among the list of candidates
     */
    UWORD32 u4_srch_pos_idx;

    /**
     * Search position for least cost among the list of candidates
     */
    UWORD32 u4_exit;

    /*
     * Buffer corresponding to best half pel cost
     */
    UWORD8 *pu1_best_hpel_buf;

} mb_part_ctxt;


/**
**************************************************************************
*   @brief   me_ctxt_t
*
*   Structure encapsulating the parameters used in the motion estimation
*   context
**************************************************************************
*/
typedef struct
{
    /**
     * Ref pointer to current MB luma
     */
    UWORD8 *pu1_ref_buf_luma;

    /**
     * Src pointer to current MB luma
     */
    UWORD8 *pu1_src_buf_luma;

    /**
     * source stride
     * (strides for luma and chroma are the same)
     */
    WORD32 i4_src_strd;

    /**
     * recon stride
     * (strides for luma and chroma are the same)
     */
    WORD32 i4_rec_strd;

    /**
     * Offset for half pel x plane from the pic buf
     */
    UWORD32 u4_half_x_offset;

    /**
     * Offset for half pel y plane from half x plane
     */
    UWORD32 u4_half_y_offset;

    /**
     * Offset for half pel xy plane from half y plane
     */
    UWORD32 u4_half_xy_offset;

    /**
     *  Search range in the X, Y axis in terms of pixels
     */
    WORD32 ai2_srch_boundaries[2];

    /**
     *  Search range in the north direction in terms of pixels
     */
    WORD32 i4_srch_range_n;

    /**
     *  Search range in the south direction in terms of pixels
     */
    WORD32 i4_srch_range_s;

    /**
     *  Search range in the east direction in terms of pixels
     */
    WORD32 i4_srch_range_e;

    /**
     *  Search range in the west direction in terms of pixels
     */
    WORD32 i4_srch_range_w;

    /**
     * left mb motion vector
     */
    ime_mv_t s_left_mv;

    /**
     * top left mb motion vector
     */
    ime_mv_t s_top_left_mv;

    /**
     * Number of valid candidates for the Initial search position
     */
    UWORD32 u4_num_candidates;

    /**
     * Motion vector predictors derived from neighbouring
     * blocks for each of the six block partitions
     */
    ime_mv_t as_mv_init_search[5];

    /**
     * mv bits
     */
    UWORD8 *pu1_mv_bits;

    /**
     * lambda (lagrange multiplier for cost computation)
     */
    UWORD32 u4_lambda_motion;

    /**
     * enabled fast sad computation
     */
    UWORD32 u4_enable_fast_sad;

    /*
     * Enable SKIP block prediction based on SATQD
     */
    UWORD32 u4_enable_stat_sad;

    /*
     * Minimum distortion to search for
     * */
    WORD32 i4_min_sad;

    /*
     * Signal that minimum sad has been reached in ME
     * */
    UWORD32 u4_min_sad_reached;

    /**
     * Flag to enable/disbale half pel motion estimation
     */
    UWORD32 u4_enable_hpel;

    /**
     * Diamond search Iteration Max Cnt
     */
    UWORD32 u4_num_layers;

    /**
     * encoder me speed
     */
    UWORD32 u4_me_speed_preset;

    UWORD32 u4_left_is_intra;

    UWORD32 u4_left_is_skip;

    /**
     * Structure to store the MB partition info
     */
    mb_part_ctxt s_mb_part;
    /*
     * Threshold to compare the sad with
     */
    UWORD16 *pu2_sad_thrsh;

    /**
     * fn ptrs for compute sad routines
     */
    ime_compute_sad_ft *pf_ime_compute_sad_16x16[2];
    ime_compute_sad_ft *pf_ime_compute_sad_16x8;
    ime_compute_sad4_diamond *pf_ime_compute_sad4_diamond;
    ime_compute_sad3_diamond *pf_ime_compute_sad3_diamond;
    ime_compute_sad2_diamond *pf_ime_compute_sad2_diamond;
    ime_sub_pel_compute_sad_16x16_ft *pf_ime_sub_pel_compute_sad_16x16;

    /*
     * Function poitners for SATQD
     */
    ime_compute_sad_stat *pf_ime_compute_sad_stat_luma_16x16;

    /**
     * Qp
     */
    UWORD8 u1_mb_qp;

    /*
     * Buffers for holding half_x , half_y and half_xy
     * values when halfpel generation
     *  for the entire plane is not enabled
     */
    UWORD8 *pu1_half_x;
    UWORD8 *pu1_half_y;
    UWORD8 *pu1_half_xy;


    /*
     * Buffers to store the best halfpel plane*
     */
    UWORD8 *pu1_hpel_buf;

    /*
     * Stride for hpel buffer
     */
    UWORD32 u4_hpel_buf_strd;

    WORD32 u4_hp_buf_strd;

} me_ctxt_t;


#endif  // _IME_STRUCTS_H_

