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
#ifndef _ISVCD_MODE_MV_RESAMPLE_H_
#define _ISVCD_MODE_MV_RESAMPLE_H_

/**
 *******************************************************************************
 * @file
 *  isvcd_mode_mv_resamp.h
 *
 * @brief
 *  Contains routines that resample for SVC resampling
 *
 * @author
 *  Kishore
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

#define MAX_NEIGHBORS 3 /*!< left top and top right neighbours */

#define LEFT                                     \
    0 /*!<Index for accessing the left MB in the \
         MV predictor array */

#define TOP                                      \
    1 /*!< Index for accessing the top MB in the \
         MV predictor array */

#define TOP_R                                   \
    2 /*!< Index for accessing the top right MB \
           in the MV predictor array */

#define MB_SKIP 255

typedef void ftype_comp_mode_mv_mb(void *pv_comp_mode_mv_hdl, void *pv_mb_params,
                                   void *pv_curr_mb_motion_mem_elements,
                                   void *pv_top_mb_motion_mem_elements, WORD32 i4_lyr_id);

typedef void(ftype_pred_direct)(void *pv_comp_mode_mv_ctxt, void *pv_top_mb_motion_l0,
                                void *pv_top_mb_motion_l1, void *pv_curr_part_motion_pred_l0,
                                void *pv_curr_part_motion_pred_l1,
                                void *pv_curr_part_motion_params_l0,
                                void *pv_curr_part_motion_params_l1,
                                WORD32 i4_curr_mb_motion_stride, WORD32 i4_availability,
                                WORD32 i4_part_width, WORD32 i4_num_mb_parts,
                                UWORD16 u2_col_zero_flag);

typedef WORD32 ftype_inter_lyr_motion_pred(void *pv_comp_mode_mv_ctxt, void *pv_mb_params,
                                           void *pv_svc_mb_params, void *ps_dec,
                                           void *ps_mb_part_info, void *pv_part);

WORD32 isvcd_interlyr_motion_mode_pred_dyadic(void *pv_comp_mode_mv_ctxt, void *pv_mb_params,
                                              void *pv_svc_mb_params, void *ps_dec,
                                              void *ps_mb_part_info, void *pv_part);

WORD32 isvcd_compute_interlyr_motion_mode(void *pv_comp_mode_mv_ctxt, void *pv_mb_params,
                                          void *pv_svc_mb_params, void *ps_dec,
                                          void *ps_mb_part_info, void *pv_part);

WORD32 isvcd_comp_mode_mv_res_init(void *pv_svc_dec);
typedef struct
{
    res_prms_t *ps_curr_lyr_res_prms; /*!< pointer to current layer
    resolution    parameters */

    WORD32 i4_offset_x;               /*!< scaled ref layer left offset. will
                  be used during mv calculation for crop window change case */

    WORD32 i4_offset_y;               /*!< scaled ref layer top offset. will
                  be used during mv calculation for crop window change case */

    WORD32 i4_scale_mv_x;             /*!< scale factor x for motion upscaling*/

    WORD32 i4_scale_mv_y;             /*!< scale factor y for motion upscaling*/

    WORD32 i4_dyadic_flag;            /*!< flag inidcates the
               dyadic upscaling or non - dyadic upscaling  */

    /* projected location pointers */
    WORD16 *pi2_ref_loc_x;                           /*!< pointer to the bufer having the
                              projected locations on horizontal direction.This used as look up
                              during calculation of reference locations */

    WORD16 *pi2_ref_loc_y;                           /*!< pointer to the bufer having the
                              projected locations on vertical direction.This used as look up
                              during calculation of reference locations */

    ref_lyr_scaled_offset_t *ps_ref_pic_lyr_offsets; /*!<
    array of structure pointer for MAX_NUM_RES - 1 to hold the scaled ref
    offset of reference picture for each resolution */

    WORD32 i4_ref_width;                             /*!< reference layer width in
                                terms of luma samples */

    WORD32 i4_ref_height;                            /*!< reference layer height in
                               terms of luma samples */

    ftype_inter_lyr_motion_pred *pf_inter_lyr_pred;  /*!< function pointer
     for dyadic optimization*/

    void *pv_ref_mv_bank_l0;     /*!< pointer to reference layer MV bank List 0*/

    mem_element_t s_ref_mb_mode; /*!< pointer to reference layer mb mode */

} mode_motion_lyr_ctxt;

typedef struct
{
    mv_pred_t *ps_motion_pred_struct; /*!< temporary motion structure
     array for 16 sub partitions used in both inter layer and spatial
     MV prediction */

    UWORD8 u1_direct_8x8_inference_flag;
    /*!< flag to indicate
    the corner mv have to be inherited in B slice */

    WORD32 i4_listx; /*!< number of lists to be procesed. this
    is set on B_SLICE or PSLICE */

    /* array to store the ref layer part idc */
    WORD32 ai4_ref_part_idc[4][4]; /*!< reference layer partition
    identifications stored of all 16 sub partitions */

    /*!< default mv pred used to store
    default motion params for intra cases in motion map */

    mode_motion_lyr_ctxt as_res_lyr_mem[MAX_NUM_RES_LYRS]; /*!< array of
    strcuture each structure for a pair of resolution except the base layer */

    WORD32 i4_res_id;                                      /*!< resolution ID used to access the
                                         Layer presistant memory */

    WORD32 i4_ref_width;                                   /*!< reference layer width in
                                      terms luma samples */

    WORD32 i4_ref_height;                                  /*!< reference layer height in
                                     terms luma samples */

    WORD32 ai4_dqid[MAX_NUM_LYRS_IN_RES];                  /*! Array of DQID of all the layers in
                                             the resolution */

} mode_motion_ctxt_t;

#endif /* _ISVCD_MODE_MV_RESAMPLE_H_ */
