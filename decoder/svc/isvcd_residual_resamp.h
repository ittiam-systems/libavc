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
#ifndef _ISVCD_RESIDUAL_RESAMP_H_
#define _ISVCD_RESIDUAL_RESAMP_H_

/**
 *******************************************************************************
 * @file
 *  isvcd_residual_resamp.h
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

#define REF_ARRAY_WIDTH_RES_SAMP (MB_WIDTH + 6)
#define REF_ARRAY_HEIGHT_RES_SAMP (MB_HEIGHT + 6)

typedef void i264_residual_reflayer_const_non_boundary_mb(
    WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride, WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
    WORD32 i4_refarray_ht, WORD32 i4_ref_mb_type_q0, WORD32 i4_ref_mb_type_q1,
    WORD32 i4_ref_mb_type_q2, WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
    WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag);

typedef void i264_residual_reflayer_const_boundary_mb(
    WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride, WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
    WORD32 i4_refarray_ht, WORD32 i4_ref_wd, WORD32 i4_ref_ht, WORD32 i4_x_offset,
    WORD32 i4_y_offset, WORD32 i4_ref_mb_type_q0, WORD32 i4_ref_mb_type_q1,
    WORD32 i4_ref_mb_type_q2, WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
    WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag);

typedef void i264_interpolate_residual(void *pv_residual_samp_ctxt, WORD16 *pi2_out,
                                       WORD32 i4_out_stride, WORD32 i4_refarray_wd, UWORD16 u2_mb_x,
                                       UWORD16 u2_mb_y, WORD32 i4_chroma_flag);

typedef void i264_residual_luma_dyadic(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                       WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                       WORD32 i4_out_res_stride, mem_element_t *ps_ref_mb_mode,
                                       UWORD16 u2_mb_x, UWORD16 u2_mb_y, WORD32 i4_ref_nnz,
                                       WORD32 i4_ref_tx_size);

typedef void i264_residual_chroma_dyadic(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                         WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                         WORD32 i4_out_res_stride);

typedef void i264_residual_chroma_dyadic_alt(void *pv_residual_samp_ctxt, UWORD16 u2_mb_x,
                                             UWORD16 u2_mb_y, mem_element_t *ps_ref_mb_mode,
                                             WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride,
                                             WORD16 *pi2_out_res, WORD32 i4_out_res_stride,
                                             WORD32 i4_cr_flag);

/*C Declarations*/
i264_residual_luma_dyadic isvcd_residual_luma_dyadic;
i264_residual_chroma_dyadic isvcd_residual_chroma_dyadic;
i264_residual_chroma_dyadic_alt isvcd_residual_chroma_dyadic_alt;

i264_interpolate_residual isvcd_interpolate_residual;
i264_residual_reflayer_const_non_boundary_mb isvcd_residual_reflayer_const_non_boundary_mb;
i264_residual_reflayer_const_boundary_mb isvcd_residual_reflayer_const_boundary_mb;

/*ARM Declarations*/
i264_residual_luma_dyadic isvcd_residual_luma_dyadic_neonintr;
i264_interpolate_residual isvcd_interpolate_residual_neonintr;
i264_residual_reflayer_const_non_boundary_mb isvcd_residual_reflayer_const_non_boundary_mb_neonintr;

/*x86 Declarations*/
i264_residual_luma_dyadic isvcd_residual_luma_dyadic_sse42;
i264_interpolate_residual isvcd_interpolate_residual_sse42;
i264_residual_reflayer_const_non_boundary_mb isvcd_residual_reflayer_const_non_boundary_mb_sse42;

typedef WORD32 ftype_residual_samp_mb(void *pv_residual_samp_ctxt, mem_element_t *ps_ref_luma,
                                      mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode,
                                      mem_element_t *ps_out_luma, mem_element_t *ps_out_chroma,
                                      UWORD16 u2_mb_x, UWORD16 u2_mb_y);

WORD32 isvcd_residual_samp_mb_dyadic(void *pv_residual_samp_ctxt, mem_element_t *ps_ref_luma,
                                     mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode,
                                     mem_element_t *ps_out_luma, mem_element_t *ps_out_chroma,
                                     UWORD16 u2_mb_x, UWORD16 u2_mb_y);

WORD32 isvcd_residual_samp_mb(void *pv_residual_samp_ctxt, mem_element_t *ps_ref_luma,
                              mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode,
                              mem_element_t *ps_out_luma, mem_element_t *ps_out_chroma,
                              UWORD16 u2_mb_x, UWORD16 u2_mb_y);

typedef struct
{
    /* used for mapping purpose */
    ref_pixel_map_t *ps_x_pos_phase;  /*!< buffers to store the projected
                                         referecne X and phase X for each
                                         pixel in current layer in
                                         horizontal direction
                                      */
    ref_pixel_map_t *ps_y_pos_phase;  /*!< buffers to store the projected
                                         referecne Y and phase Y for each
                                         pixel in current layer in
                                         vertical direction
                                      */
    ref_mb_map_t *ps_x_offset_length; /*!< buffers to store the projected
                                      start point of reference window and
                                      reference array  width in
                                      horizontal direction for each MB in
                                      current layer
                                  */
    ref_mb_map_t *ps_y_offset_length; /*!< buffers to store the projected
                                      start point of reference window and
                                      reference array  height in
                                      vertical direction for each MB in
                                      current layer
                                  */
} residual_samp_map_ctxt_t;

typedef struct
{
    residual_samp_map_ctxt_t s_luma_map_ctxt;    /*!< map structure for luma
                                                     projected locations
                                                     for curr resolution layer
                                                   */
    residual_samp_map_ctxt_t s_chroma_map_ctxt;  /*!< map structure for chroma
                                                    projected locations
                                                    for curr resolution layer
                                                  */
    WORD32 i4_ref_width;                         /*!< reference layer width in
                                                   terms luma samples
                                                 */
    WORD32 i4_ref_height;                        /*!< reference layer height in
                                                    terms luma samples
                                                  */
    WORD32 i4_curr_width;                        /*!< current layer width in
                                                   terms luma samples
                                                 */
    WORD32 i4_curr_height;                       /*!< current layer height in
                                                   terms luma samples
                                                 */
    WORD32 i4_dyadic_flag;                       /*!< flag to indicate whether
                                                   the upscaling factor is 2
                                                   in both directions
                                                 */
    ftype_residual_samp_mb *pf_residual_samp_mb; /*!< function pointer
                                                 for dyadic optimization*/

    /* following variables are for Dyadic cases only */
    WORD32 i4_chrm_alt_proc;      /*!< Alternate processing
                                  for chroma for specific
                                  values of phases
                                  */

    WORD32 i4_chrm_vert_int_mode; /*!< Chroma horizontal
                                  interpolation alternate
                                  mode
                                 */

    WORD32 i4_chrm_horz_int_mode; /*!<Chroma vertical
                                  interpolation alternate
                                  modes
                                  */
} res_lyr_ctxt;

typedef struct
{
    res_lyr_ctxt as_res_lyrs[MAX_NUM_RES_LYRS]; /*!< Array of resolutoin layer
                                                  ctxt.The first strcuture in the
                                                  array will be used for storing
                                                  the "second resolution" map in
                                                  an access unit w.r.t to its
                                                  base resolution, and for base
                                                  resolution nothing will be
                                                  computed or stored
                                              */

    WORD16 *pi2_refarray_buffer;                /*!< buffer to store the reference
                                                      layer data before residual
                                                      sampling
                                                 */

    UWORD8 *pu1_ref_x_ptr_incr;                 /*!< buffer to store the reference
                                                    array ptr increments for
                                                    operand 2 of interpolation
                                              */
    UWORD8 *pu1_ref_y_ptr_incr;                 /*!< buffer to store the reference
                                                    array ptr increments for
                                                    operand 2 of interpolation
                                              */

    WORD32 i4_res_lyr_id;                       /*!< resolution id of the layer
                                                     which is to be processed
                                                 */
    WORD32 i4_ref_width;                        /*!< reference layer width in
                                                  terms luma samples
                                                 */

    WORD32 i4_ref_height;                       /*!< reference layer height in
                                                 terms luma samples
                                                  */

    /*Dyadic Residual Resamp*/
    i264_residual_luma_dyadic *pf_residual_luma_dyadic;
    i264_residual_chroma_dyadic *pf_residual_chroma_dyadic;
    i264_residual_chroma_dyadic_alt *pf_residual_chroma_dyadic_alt;

    /*Non-dyadic Residual Resamp*/
    i264_interpolate_residual *pf_interpolate_residual;
    i264_residual_reflayer_const_non_boundary_mb *pf_residual_reflayer_const_non_boundary_mb;
    i264_residual_reflayer_const_boundary_mb *pf_residual_reflayer_const_boundary_mb;
} residual_sampling_ctxt_t;

WORD32 isvcd_residual_samp_res_init(void *pv_dec);

#endif /* _ISVCD_RESIDUAL_RESAMP_H_ */
