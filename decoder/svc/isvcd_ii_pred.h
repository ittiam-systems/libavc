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
 *  isvcd_ii_pred.h
 *
 * @brief
 *  Contains structures and function definitions required for SVC resampling
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_II_PRED_H_
#define _ISVCD_II_PRED_H_

typedef struct
{
    /*  Figure shows the projection of current MB onto reference layer
        mapping to a certain case

                  MB_WIDTH
       <--------------------------->
        ---------------------------   ^
       |                |          |  |
       |  TOP_L         |  TOP_R   |  |
       |                |          |  |
       |                |          |  |
       |                |          |  |
       |          (x,y) |          |  | MB_HEIGHT
       |---------------------------|  |
       |                |          |  |
       |                |          |  |
       |   BOT_L        |   BOT_R  |  |
       |                |          |  |
        ---------------------------   |
                                      ^
    */

    UWORD8 u1_top_left_intra_flag; /*!< flag to  inidicate the TOP_L
                                       partition is falling in to a INTRA region
                                       in base layer
                                    */
    UWORD8 u1_top_rt_intra_flag;   /*!< flag to  inidicate the TOP_R
                                      partition is falling in to a INTRA region
                                      in base layer
                                   */
    UWORD8 u1_bot_rt_intra_flag;   /*!< flag to  inidicate the BOT_R
                                      partition is falling in to a INTRA region
                                      in base layer
                                   */
    UWORD8 u1_bot_left_intra_flag; /*!< flag to  inidicate the BOT_L
                                      partition is falling in to a INTRA region
                                      in base layer
                                   */
    UWORD8 u1_intersection_x;      /*!< Horizontal point where the projection
                                        falls into a different MB in reference
                                        layer
                                    */
    UWORD8 u1_intersection_y;      /*!< Vertical point where the projection
                                        falls into a different MB in reference
                                        layer
                                    */
} intra_inter_mb_t;

typedef struct
{
    WORD16 *pi2_ref_loc_x;                  /*!< buffer pointer which holds
                                               the projected location on reference
                                               layer in horizontal direction
                                               for each pixel in current layer
                                             */
    WORD16 *pi2_ref_loc_y;                  /*!< buffer pointer which holds
                                              the projected location on reference
                                              layer in vertical direction
                                              for each pixel in current layer
                                            */
    intra_inter_mb_t s_intra_inter_mb_prms; /*!< array structure
                                       to hold the intersection points
                                       and the intra status for the
                                       all 4 parts around the
                                       intersection point
                                   */
    WORD32 i4_ref_res_lyr_wd;               /*!< Width of reference layer */
    WORD32 i4_ref_res_lyr_ht;               /*!< Height of reference layer */
    WORD32 i4_cur_res_lyr_wd;               /*!< Width of reference layer */
    WORD32 i4_cur_res_lyr_ht;               /*!< Height of reference layer */

} intra_inter_pred_ctxt_t;

typedef struct
{
    UWORD8 *pu1_recon_luma;
    UWORD8 *pu1_mc_pred_luma;
    UWORD8 *pu1_intra_pred_luma;
    WORD16 *pi2_res_luma;
    WORD32 i4_recon_luma_stride;
    WORD32 i4_mc_pred_luma_stride;
    WORD32 i4_intra_pred_luma_stride;
    WORD32 i4_res_luma_stride;
    UWORD8 *pu1_recon_chroma;
    UWORD8 *pu1_mc_pred_chroma;
    UWORD8 *pu1_intra_pred_chroma;
    WORD16 *pi2_res_chroma;
    WORD32 i4_recon_chroma_stride;
    WORD32 i4_mc_pred_chroma_stride;
    WORD32 i4_intra_pred_chroma_stride;
    WORD32 i4_res_chroma_stride;
} intra_inter_mb_buff_t;

WORD32 isvcd_ii_pred_compute_flags_mb(void *pv_ii_pred_ctxt, mem_element_t *ps_ref_mb_mode,
                                      mb_coord_t *ps_coord, void *pv_mb_prms, void *pv_svc_mb_prms,
                                      UWORD8 *pu1_ii_mb_mode);

WORD32 isvcd_ii_pred_res_init(void *pv_svc_dec);
#endif /* _ISVCD_II_PRED_H_ */
