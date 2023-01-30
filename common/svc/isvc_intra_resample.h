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

#ifndef _ISVC_INTRA_RESAMPLE_H_
#define _ISVC_INTRA_RESAMPLE_H_

#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_platform_macros.h"
#include "isvc_structs.h"

#define DYADIC_REF_W_Y 20
#define DYADIC_REF_H_Y 20
#define DYADIC_REF_W_C 10
#define DYADIC_REF_H_C 10

#define MAX_NUM_RES_LYRS 4

#define MAX_PIX_FILL_LUMA 4
#define MAX_PIX_FILL_CHROMA 2

#define MAX_REF_ARR_WD_HT 48
#define MAX_REF_IDX_ARRAY (MAX_REF_ARR_WD_HT + MB_SIZE)

#define CLIPUCHAR(x) CLIP3(0, 255, (x))

#define REF_ARRAY_WIDTH 48
#define REF_ARRAY_HEIGHT 48

typedef void FT_INTERPOLATE_LUMA_2X(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                    UWORD8 *pu1_out_buf, WORD32 i4_out_stride);

typedef void FT_VERT_INTERPOLATE_CHROMA_2X(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                           WORD32 i4_phase_0, WORD32 i4_phase_1);

typedef void FT_HORZ_INTERPOLATE_CHROMA_2X(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                           WORD32 i4_out_stride, WORD32 i4_phase_0,
                                           WORD32 i4_phase_1);

typedef struct mem_element_t
{
    /* Buffer pointer */
    void *pv_buffer;

    /* size of the structure or unit */
    WORD32 i4_element_size;

    /* Stride of buffer in terms of number of elements.*/
    WORD32 i4_num_element_stride;
} mem_element_t;

typedef struct seg_description_t
{
    /* describes segment dimension */
    UWORD8 u1_seg_dim;

    /* describes offset from start */
    UWORD8 u1_seg_off;

    /* describes whether mb is adjoining the segment
       0 => not adjoining 1 => adjoining */
    UWORD8 u1_mb_adjoin;

    /* distance to nearest MB */
    WORD8 i1_dist_idx;

    /* describes the nearest mb boundary
       +1 => rightMB/bottomMB
       -1 => leftMB/topMB	*/
    WORD8 i1_nearst_mb_bdry;
} seg_description_t;

typedef struct seg_lookup_desc_t
{
    /* place holder to store the number of segments */
    UWORD8 u1_num_segments;

    /* this variable indicates where is start locatiion of the segment with
       respect to less the block_width or greater than block width*/
    UWORD8 u4_start_pos;

    /* place holder to store per segment description */
    seg_description_t s_segments[4];
} seg_lookup_desc_t;

typedef struct intra_samp_lyr_ctxt
{
    /* mb position */
    coordinates_t *ps_mb_pos;

    /* reference layer width in terms luma samples */
    WORD32 i4_ref_width;

    /* reference layer height in terms luma samples */
    WORD32 i4_ref_height;

    /* Constrained intra resampling flag. Range is [0,1]. */
    WORD8 i1_constrained_intra_rsmpl_flag;

    /* Chroma xPhase for even values of x for dyadic cases */
    WORD32 i4_x_phase_0;

    /* Chroma xPhase for odd values of x for dyadic cases */
    WORD32 i4_x_phase_1;

    /* Chroma yPhase for even values of y for dyadic cases */
    WORD32 i4_y_phase_0;

    /* Chroma yPhase for odd values of y for dyadic cases */
    WORD32 i4_y_phase_1;

    FT_INTERPOLATE_LUMA_2X *pf_interpolate_luma;

    FT_VERT_INTERPOLATE_CHROMA_2X *pf_vert_interpol_chroma;

    FT_HORZ_INTERPOLATE_CHROMA_2X *pf_horz_interpol_chroma;

    WORD16 i2_x_min_pos;

    WORD16 i2_x_max_pos;

    WORD16 i2_y_min_pos;

    WORD16 i2_y_max_pos;

    coordinates_t *ps_phase;

    WORD32 *pi4_ref_array_positions_x;

    WORD32 *pi4_ref_array_positions_y;

    coordinates_t *ps_offsets;

    coordinates_t *ps_ref_array_dims;

    /* buffers to store lookup for horizontal segment description  */
    seg_lookup_desc_t as_seg_lookup_horz[MB_SIZE];

    /* buffers to store lookup for vertical segment description  */
    seg_lookup_desc_t as_seg_lookup_vert[MB_SIZE];

    /* buffers to store lookup for x indexes to get
       availability from 4x4 availability grid */
    UWORD8 au1_refarray_x_idx[MAX_REF_IDX_ARRAY];

    /* buffers to store lookup for y indexes to get
       availability from 4x4 availability grid */
    UWORD8 au1_refarray_y_idx[MAX_REF_IDX_ARRAY];
} intra_samp_lyr_ctxt;

typedef struct intra_sampling_ctxt_t
{
    /* Array of resolution layer ctxt. */
    intra_samp_lyr_ctxt as_res_lyrs[MAX_NUM_RES_LYRS];

    /* pointer to array of SPS */
    void *ps_sps;

    /* buffer to store the reference layer data before intra sampling */
    UWORD8 *pu1_refarray_buffer;

    /* buffer to hold the reference layer Cb data before intra
       resampling (used for dyadic cases only) */
    UWORD8 *pu1_refarray_cb;

    /* buffer to hold the reference layer Cr data before intra
       resampling (used for dyadic cases only) */
    UWORD8 *pu1_refarray_cr;

    /* intermideate buffer for interpolation */
    WORD32 *pi4_temp_interpolation_buffer;

    /* resolution id of the layer which is to be processed */
    WORD32 i4_res_lyr_id;

    /* reference layer width in terms luma samples */
    WORD32 i4_ref_width;

    /* reference layer width in terms luma samples */
    WORD32 i4_refarray_stride;

    /* reference layer height in terms luma samples */
    WORD32 i4_ref_height;
} intra_sampling_ctxt_t;

typedef struct inter_lyr_mb_prms_t
{
    /* NNZs of Chroma. Here each bit corresonds
       to a NNZs of 4x4 sub block. Lower 4 bits are
       used for Cb and upper are used for Cr */
    UWORD8 u1_chroma_nnz;

    /* NNZs of Luma. Here each bit corresonds
       to a NNZs of 4x4 sub block in raster scan order. */
    UWORD16 u2_luma_nnz;

    /* Packed MB mode transform size of an MB */
    WORD8 i1_mb_mode;
} inter_lyr_mb_prms_t;

/* Function declarations */
extern void isvc_intra_samp_mb_dyadic(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                                      mem_element_t *ps_ref_chroma,
                                      mem_element_t *ps_ref_mb_mode_map,
                                      mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                                      UWORD16 u2_mb_x, UWORD16 u2_mb_y,
                                      WORD32 i4_scaled_ref_layer_left_offset,
                                      WORD32 i4_scaled_ref_layer_top_offset);

extern void isvc_intra_samp_mb(void *pv_intra_samp_ctxt_luma, void *pv_intra_samp_ctxt_chroma,
                               mem_element_t *ps_ref_luma, mem_element_t *ps_ref_chroma,
                               mem_element_t *ps_ref_mb_mode_map, mem_element_t *ps_curr_luma,
                               mem_element_t *ps_curr_chroma);

extern void isvc_intra_resamp_generate_segment_lookup(seg_lookup_desc_t *ps_seg_lookup_table,
                                                      WORD32 i4_dimension, WORD32 i4_mb_size,
                                                      WORD32 i4_shift_val);

/* C Declarations */
extern FT_INTERPOLATE_LUMA_2X isvc_interpolate_base_luma_dyadic;
extern FT_VERT_INTERPOLATE_CHROMA_2X isvc_vert_interpol_chroma_dyadic;
extern FT_HORZ_INTERPOLATE_CHROMA_2X isvc_horz_interpol_chroma_dyadic;

/* SSE42 Declarations */
extern FT_INTERPOLATE_LUMA_2X isvc_interpolate_base_luma_dyadic_sse42;
extern FT_VERT_INTERPOLATE_CHROMA_2X isvc_vert_interpol_chroma_dyadic_sse42;
extern FT_HORZ_INTERPOLATE_CHROMA_2X isvc_horz_interpol_chroma_dyadic_sse42;

/* NEON Declarations */
extern FT_INTERPOLATE_LUMA_2X isvc_interpolate_base_luma_dyadic_neon;
extern FT_VERT_INTERPOLATE_CHROMA_2X isvc_vert_interpol_chroma_dyadic_neon;
extern FT_HORZ_INTERPOLATE_CHROMA_2X isvc_horz_interpol_chroma_dyadic_neon;

#endif
