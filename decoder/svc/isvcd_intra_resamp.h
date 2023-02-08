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
 *  isvcd_intra_resamp.h
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

#ifndef _ISVCD_INTRA_RESAMPLE_H_
#define _ISVCD_INTRA_RESAMPLE_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"

#define SVCD_FALSE 0
#define SVCD_TRUE 1

#define MAP_BUFF_WIDTH 48
#define MAP_BUFF_HEIGHT 48
#define INTERMEDIATE_BUFF_WIDTH 48
#define INTERMEDIATE_BUFF_HEIGHT (MB_HEIGHT + 4)

#define MAX_REF_ARR_WD_HT 48

#define MAX_REF_IDX_ARRAY (MAX_REF_ARR_WD_HT + MB_WIDTH)
#define DYADIC_REF_W_Y 20
#define DYADIC_REF_H_Y 20
#define DYADIC_REF_W_C 10
#define DYADIC_REF_H_C 10

#define SUB_BLOCK_WIDTH 4
#define SUB_BLOCK_HEIGHT 4
#define SUB_BLOCK_SIZE (SUB_BLOCK_WIDTH * SUB_BLOCK_HEIGHT)
#define BLOCK_WIDTH 8
#define BLOCK_HEIGHT 8
#define BLOCK_SIZE (BLOCK_WIDTH * BLOCK_HEIGHT)
#define MB_WIDTH 16
#define MB_HEIGHT 16
#define CLIPUCHAR(x) CLIP3(0, 255, (x))
#define MB_WIDTH_SHIFT 4
#define MB_HEIGHT_SHIFT 4

#define REF_ARRAY_WIDTH 48
#define REF_ARRAY_HEIGHT 48
#define MAX_PIX_FILL_LUMA 4
#define MAX_PIX_FILL_CHROMA 2

extern const WORD8 g_ai1_interp_filter_luma[64];
extern const UWORD8 g_au1_interp_filter_chroma[32];
extern WORD32 ref_pos_luma[4][16];
extern WORD32 ref_pos_chroma[4][8];
extern WORD32 phase_luma[3][16];
extern UWORD8 phase_luma_u8[3][16];
extern WORD8 phase_luma_x86[6][16];
extern WORD32 phase_chroma[3][8];
extern UWORD8 phase_chroma_u8[3][8];
extern UWORD8 ref_pos_luma_mask_m48[8][16];
extern UWORD8 ref_pos_luma_mask_m16[8][16];
extern UWORD8 ref_pos_luma_mask_m32[8][16];
extern UWORD8 ref_pos_chroma_mask_m24[2][16];
extern UWORD8 ref_pos_chroma_mask_m8[2][16];
extern UWORD8 ref_pos_chroma_mask_m16[2][16];

void isvcd_copy_data(UWORD8 *pu1_src, WORD32 i4_src_stride, UWORD8 *pu1_dst, WORD32 i4_dst_stride,
                     WORD32 i4_num_bytes, WORD32 i4_num_lines);

static __inline int isvcd_left_most_bit_detect(UWORD32 u4_num)
{
    WORD32 i4_number = 0;
    if(0xff == u4_num)
    {
        return 32;
    }

    do
    {
        if(0 == (u4_num & 0x80000000))
        {
            return i4_number;
        }
        u4_num <<= 1;
        i4_number++;
    } while(1);
}

typedef void i264_vert_interpol_chroma_dyadic(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                              WORD32 i4_phase_0, WORD32 i4_phase_1);

typedef void i264_horz_interpol_chroma_dyadic(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                              WORD32 i4_out_stride, WORD32 i4_phase_0,
                                              WORD32 i4_phase_1);

typedef void (*pf_vert_interpol_chroma_dyadic)(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                               WORD32 i4_phase_0, WORD32 i4_phase_1);

typedef void (*pf_horz_interpol_chroma_dyadic)(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                               WORD32 i4_out_stride, WORD32 i4_phase_0,
                                               WORD32 i4_phase_1);

typedef void(ftype_intra_samp_padding)(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                       WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                       UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                       WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                       WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available);
void isvcd_left_right_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                              UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                              UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                              WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                              WORD32 i4_corner_pixel_available);

void isvcd_left_right_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                     UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                     UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                     WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                     WORD32 i4_corner_pixel_available);

void isvcd_top_bot_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                           UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                           UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                           WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available);

void isvcd_top_bot_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                  UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                  UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                  WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                  WORD32 i4_corner_pixel_available);

void isvcd_diag_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                        UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                        UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                        WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available);

void isvcd_diag_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                               UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                               UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                               WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                               WORD32 i4_corner_pixel_available);

void isvcd_diag_reconstruction(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                               UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                               UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                               WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                               WORD32 i4_corner_pixel_available);

void isvcd_diag_reconstruction_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                      WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                      UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                      WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                      WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available);

typedef void i264_interpolate_base_luma_dyadic(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                               UWORD8 *pu1_out_buf, WORD32 i4_out_stride);

typedef void i264_interpolate_intra_base(void *pv_intra_samp_ctxt, UWORD8 *pu1_out,
                                         WORD32 i4_out_stride, WORD32 i4_refarray_wd,
                                         WORD32 u2_mb_x, WORD32 u2_mb_y, WORD32 i4_chroma_flag,
                                         WORD32 i4_refarray_flag);

/*C Declarations*/
i264_interpolate_base_luma_dyadic isvcd_interpolate_base_luma_dyadic;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_1;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_2;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_3;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_1;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_2;

/*ARM Declarations*/
i264_interpolate_base_luma_dyadic isvcd_interpolate_base_luma_dyadic_neonintr;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_1_neonintr;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_2_neonintr;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_3_neonintr;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_1_neonintr;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_2_neonintr;

i264_interpolate_intra_base isvcd_interpolate_intra_base;
i264_interpolate_intra_base isvcd_interpolate_intra_base_sse42;
i264_interpolate_intra_base isvcd_interpolate_intra_base_neonintr;

/*x86 Declarations*/
i264_interpolate_base_luma_dyadic isvcd_interpolate_base_luma_dyadic_sse42;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_1_sse42;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_2_sse42;
i264_vert_interpol_chroma_dyadic isvcd_vert_interpol_chroma_dyadic_3_sse42;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_1_sse42;
i264_horz_interpol_chroma_dyadic isvcd_horz_interpol_chroma_dyadic_2_sse42;

typedef struct
{
    UWORD16 u2_mb_x; /*!< MB X of the MB which has to
                          be processed
                      */
    UWORD16 u2_mb_y; /*!< MB Y of the MB which has to
                          be processed
                      */
} mb_coord_t;

typedef struct
{
    void *pv_buffer;              /*!< Buffer pointer */
    WORD32 i4_element_size;       /*!< size of the structure or unit */
    WORD32 i4_num_element_stride; /*!< Stride of buffer in terms of number
                                   of elements.                                 */
} mem_element_t;

typedef struct
{
    void *pv_buffer1;             /*!< Buffer pointer */
    void *pv_buffer2;             /*!< Buffer pointer */
    WORD32 i4_element_size;       /*!< size of the structure or unit */
    WORD32 i4_num_element_stride; /*!< Stride of buffer in terms of number
                                   of elements.
                               */
} mem_element_2_t;
typedef struct
{
    UWORD8 u1_seg_dim;       /*!< describes segment dimension */
    UWORD8 u1_seg_off;       /*!< describes offset from start */
    UWORD8 u1_mb_adjoin;     /*!< describes whether mb is adjoining
                                 the segment 0 => not adjoining
                                             1 => adjoining    */

    WORD8 i1_dist_idx;       /*!< distance to nearest MB */

    WORD8 i1_nearst_mb_bdry; /*!< describes the nearest mb boundary
                                 +1 => rightMB/bottomMB
                                 -1 => leftMB/topMB    */
} seg_description_t;

typedef struct
{
    UWORD8 u1_num_segments;          /*!< place holder to store the number of
                                     segments */

    UWORD8 u4_start_pos;             /*!< this variable indicates where is
                                      start locatiion of the segment with
                                      respect to less the block_width or
                                      greater than block width            */

    seg_description_t s_segments[4]; /*!< place holder to store per segment
                                           description    */
} seg_lookup_desc_t;
typedef struct
{
    WORD16 i2_min_pos; /*!< place holder to store the projected
                           MIN referecne position for a MB in
                           current layer. can be used to store
                           either horizontal or vertical positions
                       */
    WORD16 i2_max_pos; /*!< place holder to store the projected
                           MAX referecne position for a MB in
                           current layer. can be used to store
                           either horizontal or vertical positions
                                  */
} ref_min_max_map_t;

typedef struct
{
    WORD16 i2_left; /*!< Horizontal offset of upper left luma sample
                       after resampling process on reference
                       layer with respect to upper left luma
                       sample of current layer.
                   */
    WORD16 i2_top;  /*!< Vertical offset of upper left luma pixel
                       after resampling process on reference
                       layer
                   */
    WORD16 i2_rt;   /*!< Horizontal offset of bottom right luma
                       sample after resampling process on
                       reference layer with respect to bottom
                       right luma sample.
                   */
    WORD16 i2_bot;  /*!< Vertical offset of bottom right luma
                       pixel after resampling process on
                       reference layer
                   */
} ref_lyr_scaled_offset_t;

typedef struct
{
    UWORD8 i2_ref_pos; /*!<  place holder to store the projected
                             referecne position for a pixel in
                             current layer. can be used to store
                             either horizontal or vertical positions
                        */
    UWORD8 i2_phase;   /*!<  place holder to store the projected
                             phase for a pixel in current layer.
                             can be used to store either
                             horizontal or vertical phase
                        */
} ref_pixel_map_t;

typedef struct
{
    WORD16 i2_offset; /*!<  place holder to store the projected
                            start point of reference window
                            for each MB in current layer.can be
                            used to store either horizontal or
                            vertical offset
                       */
    WORD16 i2_length; /*!<  place holder to store reference array
                            length of the reference window
                            for each MB in current layer.can be
                            used to store either horizontal width
                            or vertical height
                       */
} ref_mb_map_t;

typedef struct
{
    WORD32 i4_res_width;                             /*!< Frame width of top most layer in the
                                                          resolution. It's expressed in terms
                                                          of luma samples.
                                                      */
    WORD32 i4_res_height;                            /*!< Frame height of top most layer in the
                                                          resolution. It's expressed in terms
                                                          of luma samples.
                                                      */
    ref_lyr_scaled_offset_t s_ref_lyr_scaled_offset; /*!< Scaled offset
                                      parameters of reference layer considering
                                      bottom most layer of the resolution as
                                      current layer. Offsets are with respect
                                      to upper left luma samples in top most
                                      layer in the resolution.
                                  */
    UWORD16 u2_scaled_ref_width;                     /*!< Considering bottom most layer of the
                                                       resolution as current layer, scaled
                                                       width of reference layer in terms of
                                                       luma pixels. It's inferred parameter
                                                       based on reference layer offsets.
                                                   */
    UWORD16 u2_scaled_ref_height;                    /*!< Considering bottom most layer of the
                                                      resolution as current layer, scaled
                                                      height of reference layer in terms of
                                                      luma pixels. It's inferred parameter
                                                      based on reference layer offsets.
                                                   */
    UWORD8 u1_rstrct_res_change_flag;                /*!< restrictedResolutionChangeFlag for
                                                  bottom most layer of the resolution. It's
                                                  a inferred parameter.
                                               */
    UWORD8 u1_cropping_change_flag;                  /*!< croppingChangeFlag for bottom most
                                                    layer of the resolution. It's a inferred
                                                    parameter.
                                                 */
    UWORD8 u1_disable_inter_lyr_dblk_filter_idc;     /*!< Mode of operation for
                                         inter layer de-blocking. It shall be
                                         set for bottom most layer of the top
                                         resolution. By top resolution, it
                                         referes to the resolution just above
                                         the current spatial resolution. This
                                         shall be valid for all resolutions
                                         except target resolution.
                                     */
    WORD8 i1_inter_lyr_alpha_c0_offset;              /*!< specifies the offset used in
                                                accessing the alpha and tC0 deblocking
                                                filter tables for filtering operations
                                                in inter layer de-blocking. Applicable
                                                for bottom most layer of the top
                                                resolution. By top resolution, it referes
                                                to the resolution just above the current
                                                spatial resolution. This shall be valid
                                                for resolutions except target resolution.
                                            */
    WORD8 i1_inter_lyr_beta_offset;                  /*!< specifies the offset used in
                                                    accessing the beta deblocking filter table
                                                    for filtering operations in inter layer
                                                    de-blocking. Applicable for bottom most
                                                    layer of the top resolution. By top
                                                    resolution, it referes to the resolution
                                                    just above the current spatial resolution.
                                                    This shall be valid for resolutions
                                                    except target resolution.
                                                */
    WORD8 i1_constrained_intra_rsmpl_flag;           /*!< Constrained intra resampling
                                             flag. Range is [0,1].
                                         */
    WORD8 i1_ref_lyr_chroma_phase_x_plus1_flag;      /*!< specifies the horizontal
                                         phase shift of the chroma components in
                                        units of half luma samples of a layer
                                        frame for the layer pictures that may be
                                        used for inter-layer prediction
                                    */

    WORD8 i1_ref_lyr_chroma_phase_y_plus1;           /*!< specifies the vertical phase
                                              shift of the chroma components in units of
                                              half luma samples of a layer frame for the
                                              layer pictures that may be used for
                                             inter-layer prediction
                                         */
    UWORD8 u1_direct_8x8_inference_flag;             /*!< Direct 8x8 inference flag
                                                      . Range is [0,1].
                                            */

    UWORD8 u1_remap_req_flag;                        /*!< this flag signifies to the
                                                          upsampling modules whether the Map
                                                          buffers have to recomputed for current
                                                          access unit or to retain the previous
                                                          access units values
                                                          */
    UWORD8 u1_dyadic_flag;                           /*!< this flag signifies the scaling
                                                          ratios are 2 in both directions
                                                          and the cropping is MB aligned
                                                        */
} res_prms_t;

typedef struct
{
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
    ref_min_max_map_t *ps_x_min_max;
    /*!< Buffer to store the projected
                                           MIN and MAX referecne position for a
                                           MB in current layer in
                                           horizontal direction
                                       */

    ref_min_max_map_t *ps_y_min_max;
    /*!< Buffer to store the projected
                                           MIN and MAX referecne position for a
                                           MB in current layer in
                                           Vertical direction
                                       */

    WORD16 *pi2_xd_index;                  /*!< buffers to store the projected
                                                 XD for each pixel in an MB
                                            */
    WORD16 *pi2_yd_index;                  /*!< buffers to store the projected
                                                YD for each pixel in an MB
                                            */
    WORD16 *pi2_ya_index;                  /*!< buffers to store the projected
                                                YA for each pixel in an MB
                                          */

    seg_lookup_desc_t *ps_seg_lookup_horz; /*!< buffers to store lookup for
                                          horizontal segment description  */

    seg_lookup_desc_t *ps_seg_lookup_vert; /*!< buffers to store lookup for
                                          vertical segment description  */

    UWORD8 *pu1_refarray_x_idx;            /*!< buffers to store lookup for
                                          x indexes to get availability
                                           from 4x4 availability grid */

    UWORD8 *pu1_refarray_y_idx;            /*!< buffers to store lookup for
                                          y indexes to get availability
                                           from 4x4 availability grid */
} intra_samp_map_ctxt_t;

typedef struct
{
    intra_samp_map_ctxt_t s_luma_map_ctxt;   /*!< map structure for luma
                                                projected locations
                                                for curr resolution layer
                                              */
    intra_samp_map_ctxt_t s_chroma_map_ctxt; /*!< map structure for chroma
                                                projected locations
                                                for curr resolution layer
                                              */
    WORD32 i4_ref_width;                     /*!< reference layer width in
                                                  terms luma samples
                                              */
    WORD32 i4_ref_height;                    /*!< reference layer height in
                                                  terms luma samples
                                              */
    WORD32 i4_curr_width;                    /*!< current layer width in
                                                 terms luma samples
                                             */
    WORD32 i4_curr_height;                   /*!< current layer height in
                                                 terms luma samples
                                             */
    WORD8 i1_constrained_intra_rsmpl_flag;   /*!< Constrained intra resampling
                                                  flag. Range is [0,1].
                                              */
    WORD32 i4_x_phase_0;                     /*!< Chroma xPhase for even
                                          values of x for dyadic cases
                                          */
    WORD32 i4_x_phase_1;                     /*!< Chroma xPhase for odd
                                          values of x for dyadic cases
                                          */
    WORD32 i4_y_phase_0;                     /*!< Chroma yPhase for even
                                          values of y for dyadic cases
                                          */
    WORD32 i4_y_phase_1;                     /*!< Chroma yPhase for odd
                                          values of y for dyadic cases
                                          */

    pf_vert_interpol_chroma_dyadic pf_vert_chroma_interpol;
    /*!< Vertical interpolation
    function for chroma
    */
    pf_horz_interpol_chroma_dyadic pf_horz_chroma_interpol;
    /*!< Horizontal interpolation
    function for chroma
    */
} intra_samp_lyr_ctxt;
typedef struct
{
    intra_samp_lyr_ctxt as_res_lyrs[MAX_NUM_RES_LYRS]; /*!< Array of resolution
                                                  layer ctxt.
                                                  The first strcuture in the
                                                  array will be used for storing
                                                  the "second resolution" map in
                                                  an access unit w.r.t to its
                                                  base resolution, and for base
                                                  resolution nothing will be
                                                  computed or stored
                                              */
    void *ps_sps;                                      /*!< pointer to array of SPS
                                                        */

    UWORD8 *pu1_refarray_buffer;                       /*!< buffer to store the reference
                                                           layer data before intra
                                                           sampling
                                                      */
    UWORD8 *pu1_refarray_cb;                           /*!< buffer to hold the reference
                                                         layer Cb data before intra
                                                         resampling (used for dyadic
                                                         cases only)
                                                    */
    UWORD8 *pu1_refarray_cr;                           /*!< buffer to hold the reference
                                                         layer Cr data before intra
                                                         resampling (used for dyadic
                                                         cases only)
                                                    */
    WORD32 *pi4_temp_interpolation_buffer;             /*!< intermideate buffer
                                                for interpolation
                                           */

    WORD32 i4_res_lyr_id;                              /*!< resolution id of the layer
                                                            which is to be processed
                                                        */
    WORD32 i4_ref_width;                               /*!< reference layer width in
                                                        terms luma samples
                                                       */

    WORD32 i4_refarray_stride;                         /*!< reference layer width in
                                                       terms luma samples
                                                       */

    WORD32 i4_ref_height;                              /*!< reference layer height in
                                                       terms luma samples
                                                        */
    res_prms_t *ps_res_prms;                           /*!< Current resolution params
                                                        */

    i264_interpolate_base_luma_dyadic *pf_interpolate_base_luma_dyadic;
    i264_interpolate_intra_base *pf_interpolate_intra_base;

    /*!< Vertical interpolation
    function for chroma
    */
    i264_vert_interpol_chroma_dyadic *pf_vert_chroma_interpol[3];

    /*!< Horizontal interpolation
    function for chroma
    */
    i264_horz_interpol_chroma_dyadic *pf_horz_chroma_interpol[2];

} intra_sampling_ctxt_t;

WORD32 isvcd_get_ceil_log2(WORD32 i4_x);

void isvcd_2d_memset(void *pv_buf, WORD32 i4_width, WORD32 i4_ht, WORD32 i4_stride, WORD32 i4_val);

WORD32 isvcd_intra_resamp_mb(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                             mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode_map,
                             mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                             mb_coord_t *ps_mb_coord);

WORD32 isvcd_intra_resamp_mb_dyadic(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                                    mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode_map,
                                    mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                                    mb_coord_t *ps_mb_coord, void *ps_svc_dec);

WORD32 isvcd_populate_res_prms(void *ps_svc_dec);

void isvcd_crop_wnd_flag_res_int(void *ps_svc_dec);

WORD32 isvcd_intra_resamp_res_init(void *ps_svc_dec);

#endif /* _ISVCD_INTRA_RESAMPLE_H_ */
