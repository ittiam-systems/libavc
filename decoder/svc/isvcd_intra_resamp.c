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
 *  isvcd_intra_resamp.c
 *
 * @brief
 *  Contains routines that resample for SVC resampling
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_get_ceil_log2()
 *  - isvcd_2d_memset()
 *  - isvcd_copy_data()
 *  - isvcd_copy_data_semiplanr()
 *  - isvcd_get_ref_layer_avlblty_dyadic()
 *  - isvcd_diagonal_construct_dyadic()
 *  - isvcd_left_right_padding()
 *  - isvcd_left_right_padding_chroma()
 *  - isvcd_top_bot_padding()
 *  - isvcd_top_bot_padding_chroma()
 *  - isvcd_diag_reconstruction()
 *  - isvcd_diag_reconstruction_chroma()
 *  - isvcd_diag_padding()
 *  - isvcd_diag_padding_chroma()
 *  - isvcd_corner_samp_dyadic()
 *  - isvcd_fill_non_avail_pixel()
 *  - isvcd_get_ref_layer_mbtype()
 *  - isvcd_reflayer_construction()
 *  - isvcd_reflayer_construction_dyadic()
 *  - isvcd_interpolate_base_luma_dyadic()
 *  - isvcd_vert_interpol_chroma_dyadic_1()
 *  - isvcd_vert_interpol_chroma_dyadic_2()
 *  - isvcd_vert_interpol_chroma_dyadic_3()
 *  - isvcd_horz_interpol_chroma_dyadic_1()
 *  - isvcd_horz_interpol_chroma_dyadic_2()
 *  - isvcd_intra_resamp_mb_dyadic()
 *  - isvcd_interpolate_intra_base()
 *  - isvcd_intra_resamp_mb()
 *  - isvcd_intra_resamp_generate_segment_lookup()
 *  - isvcd_intra_resamp_populate_list()
 *  - isvcd_populate_res_prms()
 *  - isvcd_crop_wnd_flag_res_int()
 *  - isvcd_intra_resamp_res_init()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <assert.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_mvpred.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "ih264d_error_handler.h"
#include "ih264d_cabac.h"
#include "ih264d_tables.h"
#include "ih264d_parse_slice.h"
#include "ih264d_utils.h"
#include "ih264d_parse_islice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264_debug.h"

ftype_intra_samp_padding *gpf_lookup_fxns_luma[32] = {
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &isvcd_left_right_padding,
    NULL,
    &isvcd_diag_reconstruction,
    NULL,
    &isvcd_left_right_padding,
    NULL,
    &isvcd_diag_reconstruction,
    NULL,
    NULL,
    &isvcd_top_bot_padding,
    &isvcd_diag_reconstruction,
    NULL,
    NULL,
    &isvcd_top_bot_padding,
    &isvcd_diag_reconstruction,
    NULL,
    &isvcd_left_right_padding,
    &isvcd_top_bot_padding,
    &isvcd_diag_reconstruction,
    &isvcd_diag_padding,
    &isvcd_left_right_padding,
    &isvcd_top_bot_padding,
    &isvcd_diag_reconstruction,
};

ftype_intra_samp_padding *gpf_lookup_fxns_chroma[32] = {
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &isvcd_left_right_padding_chroma,
    NULL,
    &isvcd_diag_reconstruction_chroma,
    NULL,
    &isvcd_left_right_padding_chroma,
    NULL,
    &isvcd_diag_reconstruction_chroma,
    NULL,
    NULL,
    &isvcd_top_bot_padding_chroma,
    &isvcd_diag_reconstruction_chroma,
    NULL,
    NULL,
    &isvcd_top_bot_padding_chroma,
    &isvcd_diag_reconstruction_chroma,
    NULL,
    &isvcd_left_right_padding_chroma,
    &isvcd_top_bot_padding_chroma,
    &isvcd_diag_reconstruction_chroma,
    &isvcd_diag_padding_chroma,
    &isvcd_left_right_padding_chroma,
    &isvcd_top_bot_padding_chroma,
    &isvcd_diag_reconstruction_chroma,
};

const UWORD32 gu4_valid_segs_lookup[16] = {
    0x0F000000, 0xCF000000, 0x3F000000, 0xFF000000, 0x0F000000, 0xCF000000, 0x3F000000, 0xFF000000,
    0x0F000000, 0x8F000000, 0x6F000000, 0xEF000000, 0x1F000000, 0x9F000000, 0x7F000000, 0xFF000000};

const WORD8 g_ai1_interp_filter_luma[64] = {
    0,  -1, -2, -3, -3, -4, -4, -3, -3, -3, -2, -1, -1, -1, -1, -1, 32, 32, 31, 30, 28, 26,
    24, 22, 19, 16, 14, 11, 8,  6,  4,  2,  0,  2,  4,  6,  8,  11, 14, 16, 19, 22, 24, 26,
    28, 30, 31, 32, 0,  -1, -1, -1, -1, -1, -2, -3, -3, -3, -4, -4, -3, -3, -2, -1};

const UWORD8 g_au1_interp_filter_chroma[32] = {32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12,
                                               10, 8,  6,  4,  2,  0,  2,  4,  6,  8,  10,
                                               12, 14, 16, 18, 20, 22, 24, 26, 28, 30};

WORD32 ref_pos_luma[4][16] = {{10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18, 18, 19, 20, 20},
                              {10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18, 18, 19, 20},
                              {2, 3, 3, 4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12},
                              {2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12}};
WORD32 ref_pos_chroma[4][8] = {{6, 7, 8, 8, 9, 10, 10, 11},
                               {6, 7, 7, 8, 9, 9, 10, 11},
                               {6, 6, 7, 8, 8, 9, 10, 10},
                               {2, 3, 4, 4, 5, 6, 6, 7}};

WORD32 phase_luma[3][16] = {{13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13},
                            {8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8},
                            {3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3}};

UWORD8 phase_luma_u8[3][16] = {{13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13},
                               {8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8},
                               {3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3, 13, 8, 3}};

WORD8 phase_luma_x86[6][16] = {{13, 8, 3, 13, 8, 3, 13, 8, 0, 0, 0, 0, 0, 0, 0, 0},
                               {3, 13, 8, 3, 13, 8, 3, 13, 0, 0, 0, 0, 0, 0, 0, 0},
                               {8, 3, 13, 8, 3, 13, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0},
                               {13, 8, 3, 13, 8, 3, 13, 8, 0, 0, 0, 0, 0, 0, 0, 0},
                               {3, 13, 8, 3, 13, 8, 3, 13, 0, 0, 0, 0, 0, 0, 0, 0},
                               {8, 3, 13, 8, 3, 13, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0}};

WORD32 phase_chroma[3][8] = {
    {13, 8, 3, 13, 8, 3, 13, 8}, {3, 13, 8, 3, 13, 8, 3, 13}, {8, 3, 13, 8, 3, 13, 8, 3}};

UWORD8 phase_chroma_u8[3][8] = {
    {13, 8, 3, 13, 8, 3, 13, 8}, {3, 13, 8, 3, 13, 8, 3, 13}, {8, 3, 13, 8, 3, 13, 8, 3}};

UWORD8 ref_pos_luma_mask_m48[8][16] = {{0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11},
                                       {2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11},
                                       {2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13},
                                       {4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13, 12, 13},
                                       {4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13, 12, 13, 14, 15},
                                       {6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13, 14, 15, 14, 15},
                                       {0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11},
                                       {4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13, 12, 13}};

UWORD8 ref_pos_luma_mask_m16[8][16] = {{0, 1, 2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11},
                                       {0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11},
                                       {2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13},
                                       {2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13},
                                       {4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13, 14, 15},
                                       {4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13, 12, 13, 14, 15},
                                       {0, 1, 2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11},
                                       {2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13}};

UWORD8 ref_pos_luma_mask_m32[8][16] = {{0, 1, 0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9},
                                       {0, 1, 2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11},
                                       {2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11},
                                       {2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13},
                                       {4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13, 12, 13},
                                       {4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13, 14, 15},
                                       {0, 1, 0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9},
                                       {2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13}};

UWORD8 ref_pos_chroma_mask_m24[2][16] = {{0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11},
                                         {2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11, 12, 13}};
UWORD8 ref_pos_chroma_mask_m8[2][16] = {{0, 1, 0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9},
                                        {2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11, 10, 11}};
UWORD8 ref_pos_chroma_mask_m16[2][16] = {{0, 1, 2, 3, 2, 3, 4, 5, 6, 7, 6, 7, 8, 9, 10, 11},
                                         {2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 8, 9, 10, 11, 12, 13}};

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_get_ceil_log2                                      */
/*                                                                           */
/*  Description   : this function returns the CeilLog2 of the given number   */
/*                                                                           */
/*                                                                           */
/*  Inputs        : i4_input : input number                                  */
/*  Globals       : none                                                     */
/*  Processing    : it calculate the bits and returns it                     */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : ceil of log to base 2                                    */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_ceil_log2(WORD32 i4_input)
{
    WORD32 i4_bits = 0;

    i4_input--;
    while(i4_input > 0)
    {
        i4_bits++;
        i4_input >>= 1;
    }
    return (i4_bits);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_2d_memset                                           */
/*                                                                           */
/*  Description   : Function performs 2D memset operation                    */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Buffer pointer                                        */
/*                  2. width                                                 */
/*                  3. Height                                                */
/*                  4. Stride                                                */
/*                  5. value                                                 */
/*  Globals       : None                                                     */
/*  Processing    : calls memset fucntion                                    */
/*                                                                           */
/*  Outputs       : Updates the buffer                                       */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 09 2021   Vijay           Draft                                */
/*                                                                           */
/*****************************************************************************/
void isvcd_2d_memset(void *pv_buf, WORD32 i4_width, WORD32 i4_ht, WORD32 i4_stride, WORD32 i4_val)
{
    WORD32 i4_y;
    UWORD8 *pu1_buf;

    pu1_buf = (UWORD8 *) pv_buf;

    for(i4_y = 0; i4_y < i4_ht; i4_y++)
    {
        memset(pu1_buf, i4_val, i4_width);
        /* Increment the pointer */
        pu1_buf += i4_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_copy_data                                           */
/*                                                                           */
/*  Description   : this module copies the data from source to destination   */
/*                  the amount of data to be copied is passed as input       */
/*                                                                           */
/*  Inputs        : pu1_src : pointer to the source buffer                   */
/*                  u2_src_stride : source buffer stride                     */
/*                  pu1_dst : pointer to the destination buffer              */
/*                  u2_dst_stride : destination buffer stride                */
/*                  u4_num_bytes : number of bytes to be copied              */
/*                  u4_num_lines : number of lines to be copied              */
/*  Globals       : none                                                     */
/*  Processing    : it does a memcpy from source to destination              */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : both buffers are assumed to be 2-D buffers               */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         29 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_copy_data(UWORD8 *pu1_src, WORD32 i4_src_stride, UWORD8 *pu1_dst, WORD32 i4_dst_stride,
                     WORD32 i4_num_bytes, WORD32 i4_num_lines)
{
    WORD32 i4_vert_lines;

    /* loop for copy all the lines requried */
    for(i4_vert_lines = 0; i4_vert_lines < i4_num_lines; i4_vert_lines++)
    {
        memcpy(pu1_dst, pu1_src, i4_num_bytes);
        pu1_src += i4_src_stride;
        pu1_dst += i4_dst_stride;
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_copy_data_semiplanr                                 */
/*                                                                           */
/*  Description   : this module copies the data from source to destination   */
/*                  the amount of data to be copied is passed as input       */
/*                                                                           */
/*  Inputs        : pu1_src : pointer to the source buffer                   */
/*                  i4_src_stride : source buffer stride                     */
/*                  pu1_dst1 : pointer to the destination buffer 1           */
/*                  pu1_dst2 : pointer to the destination buffer 2           */
/*                  i4_dst_stride : destination buffer stride                */
/*                  i4_num_bytes : number of bytes to be copied              */
/*                  i4_num_lines : number of lines to be copied              */
/*  Globals       : none                                                     */
/*  Processing    : it does a memcpy from source to destination              */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : both buffers are assumed to be 2-D buffers               */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         29 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_copy_data_semiplanr(UWORD8 *pu1_src, WORD32 i4_src_stride, UWORD8 *pu1_dst1,
                               UWORD8 *pu1_dst2, WORD32 i4_dst_stride, WORD32 i4_num_bytes,
                               WORD32 i4_num_lines)
{
    WORD32 i4_vert_lines, u4_i;

    /* loop for copy all the lines requried */
    for(i4_vert_lines = 0; i4_vert_lines < i4_num_lines; i4_vert_lines++)
    {
        for(u4_i = 0; u4_i < i4_num_bytes; u4_i++)
        {
            *(pu1_dst1 + u4_i) = *(pu1_src + (2 * u4_i));
            *(pu1_dst2 + u4_i) = *(pu1_src + (2 * u4_i) + 1);
        }

        pu1_src += i4_src_stride;
        pu1_dst1 += i4_dst_stride;
        pu1_dst2 += i4_dst_stride;
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_ref_layer_avlblty_dyadic                        */
/*                                                                           */
/*  Description   : This function is used to find the mb type of the         */
/*                  corresponding MB in the reference layer for dyadic cases */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra samp context                  */
/*                  pi1_ref_mb_modes : ref mb modes buffer pointer           */
/*                  i4_ref_mode_stride : mb mode buffer stride               */
/*                  i4_ref_mb_x : reference MB location X                    */
/*                  i4_ref_mb_y : reference MB location Y                    */
/*                  pi4_mb_type : pointer to store the mb type               */
/*                  i1_curr_slice_id : slice id of current MB                */
/*                  i1_cons_intr_samp_flag :constrained intra resampling flag*/
/*  Globals       : none                                                     */
/*  Processing    : it derives the bit corresponding to reference MB and     */
/*                  stores the mbtype as INTRA if the bit is set             */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_get_ref_layer_avlblty_dyadic(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                                        WORD32 i4_element_size, WORD32 i4_ref_mb_x,
                                        WORD32 i4_ref_mb_y, WORD32 *pi4_avlblty,
                                        WORD8 i1_curr_slice_id, WORD8 i1_cons_intr_samp_flag)
{
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
    WORD8 i1_mb_mode;

    /* get the location of the byte which has the current mb mode */
    pi1_ref_mb_modes += (i4_ref_mb_y * i4_ref_mode_stride * i4_element_size);
    pi1_ref_mb_modes += (i4_ref_mb_x * i4_element_size);
    ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes;
    i1_mb_mode = ps_inter_lyr_mb_prms->i1_mb_mode;

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        /* INTER */
        *pi4_avlblty = 0;
    }
    else
    {
        /* INTRA */
        *pi4_avlblty = 1;
    }

    /* if constrained intra flag is 1 then check for same slice id */
    if(1 == i1_cons_intr_samp_flag)
    {
        if(1 == *pi4_avlblty)
        {
            /* check for different slice idc */
            if(ps_inter_lyr_mb_prms->i1_slice_id != i1_curr_slice_id)
            {
                /* store the mode as not available for upsampling */
                *pi4_avlblty = 0;
            }
        }
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_diagonal_construct_dyadic                           */
/*                                                                           */
/*  Description   : This function fills the unavaible pixels in the reference*/
/*                    array with diagonally constructed samples              */
/*  Inputs        : i4_x :current position in reference array X to be filled */
/*                  i4_y :current position in reference array Y to be filled */
/*                  i4_xd_index : diagonal index in horizontal direction     */
/*                  i4_yd_index : diagonal index in vertical direction       */
/*                  pu1_refarray : popinter to reference array               */
/*                  i4_refarray_wd: width of the reference array             */
/*  Globals       : none                                                     */
/*  Processing    : Fills the sample which is unavailable with filtered      */
/*                    diagonal samples                                       */
/*  Outputs       : pixel filled                                             */
/*  Returns       : constructed pixel                                        */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         creation                             */
/*                                                                           */
/*****************************************************************************/
UWORD8 isvcd_diagonal_construct_dyadic(WORD32 i4_x, WORD32 i4_y, WORD32 i4_xd_index,
                                       WORD32 i4_yd_index, UWORD8 *pu1_refarray,
                                       WORD32 i4_refarray_wd)
{
    WORD32 i4_diff_hor_ver, i4_sgn_xy;
    WORD32 i4_xc, i4_yc;
    WORD32 i4_samp1, i4_samp2, i4_samp3;
    WORD32 i4_result;
    UWORD8 *pu1_tmp;

    i4_diff_hor_ver = ABS(i4_xd_index) - ABS(i4_yd_index);
    i4_sgn_xy = SIGN(i4_xd_index * i4_yd_index);

    if(i4_diff_hor_ver > 0)
    {
        i4_xc = i4_x - (i4_sgn_xy * i4_yd_index);
        i4_yc = i4_y - i4_yd_index;
        pu1_tmp = pu1_refarray + (i4_yc * i4_refarray_wd);
        i4_samp1 = pu1_tmp[i4_xc - 1];
        i4_samp2 = pu1_tmp[i4_xc];
        i4_samp3 = pu1_tmp[i4_xc + 1];
    }
    else if(i4_diff_hor_ver < 0)
    {
        i4_xc = i4_x - i4_xd_index;
        i4_yc = i4_y - (i4_sgn_xy * i4_xd_index);
        pu1_tmp = pu1_refarray + ((i4_yc - 1) * i4_refarray_wd);
        i4_samp1 = pu1_tmp[i4_xc];
        pu1_tmp += i4_refarray_wd;
        i4_samp2 = pu1_tmp[i4_xc];
        pu1_tmp += i4_refarray_wd;
        i4_samp3 = pu1_tmp[i4_xc];
    }
    else
    {
        WORD32 i4_ref_xd, i4_ref_yd;

        i4_ref_xd = i4_x - i4_xd_index;
        i4_ref_yd = i4_y - i4_yd_index;
        i4_xc = i4_ref_xd + SIGN(i4_xd_index);
        i4_yc = i4_ref_yd + SIGN(i4_yd_index);
        pu1_tmp = pu1_refarray + (i4_ref_yd * i4_refarray_wd);
        i4_samp1 = pu1_tmp[i4_xc];
        i4_samp2 = pu1_tmp[i4_ref_xd];
        pu1_tmp = pu1_refarray + (i4_yc * i4_refarray_wd);
        i4_samp3 = pu1_tmp[i4_ref_xd];
    }

    i4_result = (i4_samp1 + (i4_samp2 << 1) + i4_samp3 + 2) >> 2;
    pu1_tmp = pu1_refarray + (i4_y * i4_refarray_wd);
    /* Store the filled sample */
    pu1_tmp[i4_x] = i4_result;

    return (i4_result);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_left_right_padding                                  */
/*  Description   :  This function does the left/right padding for intra     */
/*                      upsampling                                           */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the left/right padding for intra      */
/*                      upsampling                                           */
/*  Outputs       : none                                                     */
/*  Returns       : number of leading zeroes                                 */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_left_right_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                              UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                              UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                              WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                              WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(i1_yd_index);
    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + i1_xd_index;
    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - i1_xd_index;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst, *pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
        pu1_src += i4_refarray_stride;
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_left_right_padding_chroma                           */
/*  Description   :  This function does the left/right padding for intra     */
/*                      upsampling for chroma                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the left/right padding for intra      */
/*                      upsampling for chroma                                */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2011   Vijay                creation                        */
/*****************************************************************************/
void isvcd_left_right_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                     UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                     UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                     WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                     WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i1_yd_index);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_src_cb = pu1_dst_cb + i1_xd_index;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i1_xd_index;

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_CHROMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_CHROMA);
    pu1_dst_cb = pu1_src_cb - i1_xd_index;
    pu1_dst_cr = pu1_src_cr - i1_xd_index;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst_cb, *pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;
        pu1_src_cb += i4_refarray_stride;
        memset(pu1_dst_cr, *pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
        pu1_src_cr += i4_refarray_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_top_bot_padding                                     */
/*  Description   :  This function does the top/bottom padding for intra     */
/*                      upsampling for chroma                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the top/bottom padding for intra      */
/*                      upsampling for chroma                                */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_top_bot_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                           UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                           UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                           WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(i1_xd_index);
    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + (i1_yd_index * i4_refarray_stride);
    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - (i1_yd_index * i4_refarray_stride);

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memcpy(pu1_dst, pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_top_bot_padding                                     */
/*  Description   :  This function does the top/bottom padding for intra     */
/*                      upsampling for chroma                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the top/bottom padding for intra      */
/*                      upsampling for chroma                                */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_top_bot_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                                  UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                                  UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                                  WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                                  WORD32 i4_corner_pixel_available)
{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i1_xd_index);
    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;

    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_src_cb = pu1_dst_cb + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i4_tmp;
    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_CHROMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_CHROMA);

    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_dst_cb = pu1_src_cb - i4_tmp;
    pu1_dst_cr = pu1_src_cr - i4_tmp;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memcpy(pu1_dst_cb, pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;
        memcpy(pu1_dst_cr, pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_diag_reconstruction                                 */
/*  Description   :  This function does the diagonal reconstuction for intra */
/*                      upsampling for luma                                  */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the diagonal reconstruction for intra */
/*                      upsampling for luma                                  */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_diag_reconstruction(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                               UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                               UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                               WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                               WORD32 i4_corner_pixel_available)
{
    WORD32 i4_i;
    UWORD8 *pu1_src_1, *pu1_src_2, *pu1_dst;
    UWORD8 u1_filter_delay_buf[18] = {0};
    UWORD8 u1_out_buf[16] = {0};
    WORD32 i4_width, i4_height;
    WORD32 i4_x_off, i4_y_off;
    WORD32 i4_block_size = BLOCK_WIDTH;

    UNUSED(pu1_refarray_2);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src_1 = pu1_dst + i1_xd_index;
    pu1_src_2 = pu1_dst + (i1_yd_index * i4_refarray_stride);

    i4_width = MAX(u1_seg_wd, (((i4_mb_adjoin_x >> 3) ^ 1) * i4_block_size));
    i4_height = MAX(u1_seg_ht, (((i4_mb_adjoin_y >> 4) ^ 1) * i4_block_size));
    i4_x_off = (i4_width - u1_seg_wd);
    i4_y_off = (i4_height - u1_seg_ht);

    if(i1_xd_index < 0 && i1_yd_index > 0)
    {
        /* Quadrant 1 Processing load the pixel in the filter delay buffer */
        for(i4_i = 0; i4_i < (i4_height + 1); i4_i++)
        {
            u1_filter_delay_buf[i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        pu1_src_2 -= i4_x_off;
        memcpy(&u1_filter_delay_buf[i4_i], pu1_src_2, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_i - 1] =
                (u1_filter_delay_buf[i4_i] + u1_filter_delay_buf[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index > 0)
    {
        /* Quadrant 2 Processing */
        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf[0], pu1_src_2, (i4_width + 1));
        for(i4_i = i4_height; i4_i > 0; i4_i--)
        {
            u1_filter_delay_buf[i4_width + i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_width] =
                (u1_filter_delay_buf[i4_width - 1] + u1_filter_delay_buf[i4_width + 1] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_height - i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index < 0)
    {
        /* Quadrant 3 Processing */
        /* load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf[0], pu1_src_2, (i4_width + 1));

        pu1_src_1 -= (i4_y_off * i4_refarray_stride);
        for(i4_i = 1; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf[i4_width + i4_i] = *pu1_src_1;
            pu1_src_1 += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_width] =
                (u1_filter_delay_buf[i4_width - 1] + u1_filter_delay_buf[i4_width + 1] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
    else
    {
        /* Quadrant 4 Processing */
        /* load the pixel in the filter delay buffer  */
        pu1_src_1 += ((u1_seg_ht - 1) * i4_refarray_stride);
        for(i4_i = 0; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf[i4_i] = *pu1_src_1;
            pu1_src_1 -= i4_refarray_stride;
        }

        pu1_src_2 -= i4_x_off;
        memcpy(&u1_filter_delay_buf[i4_i], pu1_src_2, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf[i4_i - 1] =
                (u1_filter_delay_buf[i4_i] + u1_filter_delay_buf[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf[i4_i] = ((u1_filter_delay_buf[i4_i]) + (u1_filter_delay_buf[i4_i + 1] * 2) +
                                (u1_filter_delay_buf[i4_i + 2]) + 2) >>
                               2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst, &u1_out_buf[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst += i4_refarray_stride;
        }
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_diag_reconstruction_chroma                          */
/*  Description   :  This function does the diagonal reconstuction for intra */
/*                      upsampling for chroma                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the diagonal reconstruction for intra */
/*                      upsampling for chroma                                */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_diag_reconstruction_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index,
                                      WORD8 i1_yd_index, UWORD8 u1_seg_wd, UWORD8 u1_seg_ht,
                                      UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                      WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                                      WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)
{
    WORD32 i4_i;
    UWORD8 u1_filter_delay_buf_cb[18] = {0};
    UWORD8 u1_filter_delay_buf_cr[18] = {0};
    UWORD8 u1_out_buf_cb[16] = {0};
    UWORD8 u1_out_buf_cr[16] = {0};
    WORD32 i4_width, i4_height;
    WORD32 i4_x_off, i4_y_off;
    WORD32 i4_block_size = BLOCK_WIDTH >> 1;
    UWORD8 *pu1_src_1_cb, *pu1_src_2_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_1_cr, *pu1_src_2_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;

    pu1_src_1_cb = pu1_dst_cb + i1_xd_index;
    pu1_src_1_cr = pu1_dst_cr + i1_xd_index;
    i4_tmp = (i1_yd_index * i4_refarray_stride);
    pu1_src_2_cb = pu1_dst_cb + i4_tmp;
    pu1_src_2_cr = pu1_dst_cr + i4_tmp;

    i4_width = MAX(u1_seg_wd, (((i4_mb_adjoin_x >> 3) ^ 1) * i4_block_size));
    i4_height = MAX(u1_seg_ht, (((i4_mb_adjoin_y >> 4) ^ 1) * i4_block_size));
    i4_x_off = (i4_width - u1_seg_wd);
    i4_y_off = (i4_height - u1_seg_ht);

    if(i1_xd_index < 0 && i1_yd_index > 0)
    {
        /* Quadrant 1 Processing load the pixel in the filter delay buffer */
        for(i4_i = 0; i4_i < (i4_height + 1); i4_i++)
        {
            u1_filter_delay_buf_cb[i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;
            u1_filter_delay_buf_cr[i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        pu1_src_2_cb -= i4_x_off;
        pu1_src_2_cr -= i4_x_off;
        memcpy(&u1_filter_delay_buf_cb[i4_i], pu1_src_2_cb, i4_width);
        memcpy(&u1_filter_delay_buf_cr[i4_i], pu1_src_2_cr, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_i - 1] =
                (u1_filter_delay_buf_cb[i4_i] + u1_filter_delay_buf_cb[i4_i - 2] + 1) >> 1;

            u1_filter_delay_buf_cr[i4_i - 1] =
                (u1_filter_delay_buf_cr[i4_i] + u1_filter_delay_buf_cr[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;
            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_x_off + i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index > 0)
    {
        /* Quadrant 2 Processing load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf_cb[0], pu1_src_2_cb, (i4_width + 1));
        memcpy(&u1_filter_delay_buf_cr[0], pu1_src_2_cr, (i4_width + 1));

        for(i4_i = i4_height; i4_i > 0; i4_i--)
        {
            u1_filter_delay_buf_cb[i4_width + i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_width + i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_width] =
                (u1_filter_delay_buf_cb[i4_width - 1] + u1_filter_delay_buf_cb[i4_width + 1] + 1) >>
                1;

            u1_filter_delay_buf_cr[i4_width] =
                (u1_filter_delay_buf_cr[i4_width - 1] + u1_filter_delay_buf_cr[i4_width + 1] + 1) >>
                1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_height - i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;

            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_height - i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else if(i1_xd_index > 0 && i1_yd_index < 0)
    {
        /* Quadrant 3 Processing load the pixel in the filter delay buffer  */
        memcpy(&u1_filter_delay_buf_cb[0], pu1_src_2_cb, (i4_width + 1));
        memcpy(&u1_filter_delay_buf_cr[0], pu1_src_2_cr, (i4_width + 1));

        i4_tmp = (i4_y_off * i4_refarray_stride);
        pu1_src_1_cb -= i4_tmp;
        pu1_src_1_cr -= i4_tmp;
        for(i4_i = 1; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf_cb[i4_width + i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb += i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_width + i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr += i4_refarray_stride;
        }

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_width] =
                (u1_filter_delay_buf_cb[i4_width - 1] + u1_filter_delay_buf_cb[i4_width + 1] + 1) >>
                1;

            u1_filter_delay_buf_cr[i4_width] =
                (u1_filter_delay_buf_cr[i4_width - 1] + u1_filter_delay_buf_cr[i4_width + 1] + 1) >>
                1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 0; i4_i < u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;
            memcpy(pu1_dst_cr, &u1_out_buf_cr[i4_y_off + i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
    else
    {
        /* Quadrant 4 Processing load the pixel in the filter delay buffer  */
        i4_tmp = ((u1_seg_ht - 1) * i4_refarray_stride);
        pu1_src_1_cb += i4_tmp;
        pu1_src_1_cr += i4_tmp;

        for(i4_i = 0; i4_i <= i4_height; i4_i++)
        {
            u1_filter_delay_buf_cb[i4_i] = *pu1_src_1_cb;
            pu1_src_1_cb -= i4_refarray_stride;

            u1_filter_delay_buf_cr[i4_i] = *pu1_src_1_cr;
            pu1_src_1_cr -= i4_refarray_stride;
        }

        pu1_src_2_cb -= i4_x_off;
        pu1_src_2_cr -= i4_x_off;
        memcpy(&u1_filter_delay_buf_cb[i4_i], pu1_src_2_cb, i4_width);
        memcpy(&u1_filter_delay_buf_cr[i4_i], pu1_src_2_cr, i4_width);

        if(0 == i4_corner_pixel_available)
        {
            /* interpolate the unavailable corner pixel */
            u1_filter_delay_buf_cb[i4_i - 1] =
                (u1_filter_delay_buf_cb[i4_i] + u1_filter_delay_buf_cb[i4_i - 2] + 1) >> 1;

            u1_filter_delay_buf_cr[i4_i - 1] =
                (u1_filter_delay_buf_cr[i4_i] + u1_filter_delay_buf_cr[i4_i - 2] + 1) >> 1;
        }

        for(i4_i = 0; i4_i < (i4_width + i4_height - 1); i4_i++)
        {
            /* get the filtered output */
            u1_out_buf_cb[i4_i] =
                ((u1_filter_delay_buf_cb[i4_i]) + (u1_filter_delay_buf_cb[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cb[i4_i + 2]) + 2) >>
                2;

            u1_out_buf_cr[i4_i] =
                ((u1_filter_delay_buf_cr[i4_i]) + (u1_filter_delay_buf_cr[i4_i + 1] * 2) +
                 (u1_filter_delay_buf_cr[i4_i + 2]) + 2) >>
                2;
        }

        /* fill the segment with diagonal reconstructed output */
        for(i4_i = 1; i4_i <= u1_seg_ht; i4_i++)
        {
            memcpy(pu1_dst_cb, &u1_out_buf_cb[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst_cb += i4_refarray_stride;
            memcpy(pu1_dst_cr, &u1_out_buf_cr[(u1_seg_ht + i4_x_off) - i4_i], u1_seg_wd);
            pu1_dst_cr += i4_refarray_stride;
        }
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_diag_padding                                        */
/*  Description   :  This function does the diagonal padding for intra       */
/*                      upsampling for luma                                  */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the diagonal padding for intra        */
/*                      upsampling for luma                                  */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_diag_padding(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                        UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                        UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride, WORD32 i4_mb_adjoin_x,
                        WORD32 i4_mb_adjoin_y, WORD32 i4_corner_pixel_available)

{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src, *pu1_dst;

    UNUSED(pu1_refarray_2);
    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    pu1_dst = pu1_refarray_1 + i4_x + (i4_y * i4_refarray_stride);
    pu1_src = pu1_dst + i1_xd_index + (i1_yd_index * i4_refarray_stride);
    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);
    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);
    pu1_dst = pu1_src - i1_xd_index - (i1_yd_index * i4_refarray_stride);

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst, *pu1_src, u1_seg_wd);
        pu1_dst += i4_refarray_stride;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_diag_padding_chroma                                 */
/*  Description   :  This function does the diagonal padding for intra       */
/*                      upsampling for chroma                                */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    : This function does the diagonal padding for intra        */
/*                      upsampling for chroma                                */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*  Issues        : none                                                     */
/*  Revision History:                                                        */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*****************************************************************************/
void isvcd_diag_padding_chroma(WORD32 i4_x, WORD32 i4_y, WORD8 i1_xd_index, WORD8 i1_yd_index,
                               UWORD8 u1_seg_wd, UWORD8 u1_seg_ht, UWORD8 *pu1_refarray_1,
                               UWORD8 *pu1_refarray_2, WORD32 i4_refarray_stride,
                               WORD32 i4_mb_adjoin_x, WORD32 i4_mb_adjoin_y,
                               WORD32 i4_corner_pixel_available)

{
    WORD32 i4_idx_i;
    UWORD8 *pu1_src_cb, *pu1_dst_cb;
    UWORD8 *pu1_src_cr, *pu1_dst_cr;
    WORD32 i4_tmp;

    UNUSED(i4_mb_adjoin_x);
    UNUSED(i4_mb_adjoin_y);
    UNUSED(i4_corner_pixel_available);

    i4_tmp = i4_x + (i4_y * i4_refarray_stride);
    pu1_dst_cb = pu1_refarray_1 + i4_tmp;
    pu1_dst_cr = pu1_refarray_2 + i4_tmp;
    i4_tmp = i1_xd_index + (i1_yd_index * i4_refarray_stride);
    pu1_src_cb = pu1_dst_cb + i4_tmp;
    pu1_src_cr = pu1_dst_cr + i4_tmp;

    i1_xd_index = MIN(i1_xd_index, MAX_PIX_FILL_LUMA);
    u1_seg_wd = MIN(u1_seg_wd, MAX_PIX_FILL_LUMA);
    i1_yd_index = MIN(i1_yd_index, MAX_PIX_FILL_LUMA);
    u1_seg_ht = MIN(u1_seg_ht, MAX_PIX_FILL_LUMA);

    i4_tmp = (i1_xd_index + (i1_yd_index * i4_refarray_stride));
    pu1_dst_cb = pu1_src_cb - i4_tmp;
    pu1_dst_cr = pu1_src_cr - i4_tmp;

    for(i4_idx_i = 0; i4_idx_i < u1_seg_ht; i4_idx_i++)
    {
        memset(pu1_dst_cb, *pu1_src_cb, u1_seg_wd);
        pu1_dst_cb += i4_refarray_stride;
        memset(pu1_dst_cr, *pu1_src_cr, u1_seg_wd);
        pu1_dst_cr += i4_refarray_stride;
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_corner_samp_dyadic                                  */
/*                                                                           */
/*  Description   : This function fills the corner sample in the reference   */
/*                    array with diagonally constructed samples              */
/*  Inputs        : i4_x :current position in reference array X to be filled */
/*                  i4_y :current position in reference array Y to be filled */
/*                  i4_xd_index : diagonal index in horizontal direction     */
/*                  i4_yd_index : diagonal index in vertical direction       */
/*                  pu1_refarray_y : pointer to luma reference array         */
/*                    pu1_refarray_cb : pointer to Cb reference array        */
/*                    pu1_refarray_cr : pointer to Cr reference array        */
/*  Globals       : none                                                     */
/*  Processing    : Fills the sample which is unavailable with filtered      */
/*                    diagonal samples                                       */
/*  Outputs       : pixel filled                                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_corner_samp_dyadic(WORD32 i4_x, WORD32 i4_y, WORD32 i4_xD, WORD32 i4_yD,
                              UWORD8 *pu1_refarray_y, UWORD8 *pu1_refarray_cb,
                              UWORD8 *pu1_refarray_cr)
{
    WORD32 i4_ref_xD, i4_ref_yD;
    WORD32 i4_c_ref_xD, i4_c_ref_yD;
    WORD32 i4_xc, i4_yc;
    WORD32 i4_c_xc, i4_c_yc;
    WORD32 i4_samp1, i4_samp2;
    UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

    i4_ref_xD = i4_x - i4_xD;
    i4_ref_yD = i4_y - i4_yD;
    i4_xc = i4_ref_xD + SIGN(i4_xD);
    i4_yc = i4_ref_yD + SIGN(i4_yD);

    /* Luma */
    pu1_tmp_src = pu1_refarray_y + (i4_yc * DYADIC_REF_W_Y);
    i4_samp1 = pu1_tmp_src[i4_ref_xD];
    pu1_tmp_src = pu1_refarray_y + (i4_ref_yD * DYADIC_REF_W_Y);
    i4_samp2 = pu1_tmp_src[i4_xc];
    pu1_tmp_dst = pu1_tmp_src;
    pu1_tmp_dst[i4_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;

    /* Chroma */
    i4_c_ref_xD = i4_ref_xD >> 1;
    i4_c_ref_yD = i4_ref_yD >> 1;
    i4_c_xc = i4_c_ref_xD + SIGN(i4_xD);
    i4_c_yc = i4_c_ref_yD + SIGN(i4_yD);

    /* Cb */
    pu1_tmp_src = pu1_refarray_cb + (i4_c_yc * DYADIC_REF_W_C);
    i4_samp1 = pu1_tmp_src[i4_c_ref_xD];
    pu1_tmp_src = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
    i4_samp2 = pu1_tmp_src[i4_c_xc];
    pu1_tmp_dst = pu1_tmp_src;
    pu1_tmp_dst[i4_c_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;

    /* Cr */
    pu1_tmp_src = pu1_refarray_cr + (i4_c_yc * DYADIC_REF_W_C);
    i4_samp1 = pu1_tmp_src[i4_c_ref_xD];
    pu1_tmp_src = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
    i4_samp2 = pu1_tmp_src[i4_c_xc];
    pu1_tmp_dst = pu1_tmp_src;
    pu1_tmp_dst[i4_c_ref_xD] = (i4_samp1 + i4_samp2 + 1) >> 1;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_fill_non_ava_pixel                                  */
/*                                                                           */
/*  Description   :  This function does the core pixel level processing      */
/*                    while filling the non available pixel                  */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  i4_refarray_wd : width of the reference array            */
/*                  i4_refarray_ht : height of the reference array           */
/*                  ps_mb_coord  : current mb coord structure                */
/*                  i4_chroma_flag : chroam processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : based on the map buffer values the non available pixels  */
/*                   are filled using border extension algorithm             */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_fill_non_avail_pixel(void *pv_map_ctxt, UWORD8 *pu1_refarray_1, UWORD8 *pu1_refarray_2,
                                WORD32 i4_refarray_stride, mb_coord_t *ps_mb_coord,
                                WORD32 i4_chroma_flag, UWORD8 u1_avail_map[4][4])
{
    /* --------------------------------------------------------------------- */
    /* Index Variables                                                         */
    /* --------------------------------------------------------------------- */
    intra_samp_map_ctxt_t *ps_map_ctxt;
    ref_mb_map_t *ps_x_off_len;
    ref_mb_map_t *ps_y_off_len;
    WORD32 i4_x, i4_y;
    WORD32 i4_corner_pixel_available;

    /* --------------------------------------------------------------------- */
    /* Local Pointer Declaration for Segment lookup                             */
    /* --------------------------------------------------------------------- */
    seg_lookup_desc_t *ps_segments_x;
    seg_lookup_desc_t *ps_segments_y;
    seg_description_t *ps_seg_desc_x, *ps_seg_desc_y;
    seg_description_t *ps_seg_x_tmp, *ps_seg_y_tmp;
    UWORD8 u1_num_sgmts_x, u1_num_sgmts_y;

    /* --------------------------------------------------------------------- */
    /* Temp Variables for Mapping context                                     */
    /* --------------------------------------------------------------------- */
    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refmb_wd;
    WORD32 i4_refmb_ht;
    WORD32 i4_mbaddr_x;
    WORD32 i4_mbaddr_y;
    WORD32 i4_xr_index, i4_yr_index;
    WORD32 i4_j, i4_i;
    WORD32 i4_cur_x;
    UWORD32 u4_lookup_4bit, u4_lookup_5bit, u4_4thbit;
    WORD32 i4_pad_size;
    WORD32 i4_x_min;
    WORD32 i4_y_min;
    WORD32 i4_x_start_pos, i4_y_start_pos;
    ref_min_max_map_t *ps_x_min_max;
    ref_min_max_map_t *ps_y_min_max;
    UWORD8 *pu1_ref_idx_x, *pu1_ref_idx_y;
    ftype_intra_samp_padding *pf_intra_samp_padding;
    ftype_intra_samp_padding **pf_intra_samp_lookup;

    ps_map_ctxt = (intra_samp_map_ctxt_t *) pv_map_ctxt;
    ps_x_min_max = ps_map_ctxt->ps_x_min_max;
    ps_y_min_max = ps_map_ctxt->ps_y_min_max;
    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;
    i4_mbaddr_y = ps_mb_coord->u2_mb_y;
    i4_mbaddr_x = ps_mb_coord->u2_mb_x;
    i4_x_offset = ps_x_off_len[i4_mbaddr_x].i2_offset;
    i4_y_offset = ps_y_off_len[i4_mbaddr_y].i2_offset;
    i4_refmb_wd = (MB_WIDTH >> i4_chroma_flag) - 1;
    i4_refmb_ht = (MB_HEIGHT >> i4_chroma_flag) - 1;

    if(0 == i4_chroma_flag)
    {
        pf_intra_samp_lookup = gpf_lookup_fxns_luma;
    }
    else
    {
        pf_intra_samp_lookup = gpf_lookup_fxns_chroma;
    }

    /* get the min and max positions */
    i4_x_min = ps_x_min_max[i4_mbaddr_x].i2_min_pos;
    i4_y_min = ps_y_min_max[i4_mbaddr_y].i2_min_pos;

    /* get the start position of the MB in reference layer */
    i4_pad_size = 2 >> i4_chroma_flag;
    i4_x_start_pos = (i4_x_min - i4_pad_size);
    i4_y_start_pos = (i4_y_min - i4_pad_size);
    i4_xr_index = (i4_x_start_pos + i4_x_offset) & i4_refmb_wd;
    i4_yr_index = (i4_y_start_pos + i4_y_offset) & i4_refmb_ht;

    /* Find the number of segments in x and y direction */
    ps_segments_x = (ps_map_ctxt->ps_seg_lookup_horz + i4_xr_index);
    ps_segments_y = (ps_map_ctxt->ps_seg_lookup_vert + i4_yr_index);
    u1_num_sgmts_x = ps_segments_x->u1_num_segments;
    u1_num_sgmts_y = ps_segments_y->u1_num_segments;
    ps_seg_desc_x = ps_segments_x->s_segments;
    ps_seg_desc_y = ps_segments_y->s_segments;
    pu1_ref_idx_x = ps_map_ctxt->pu1_refarray_x_idx;
    pu1_ref_idx_y = ps_map_ctxt->pu1_refarray_y_idx;
    i4_cur_x = pu1_ref_idx_x[i4_x_start_pos];
    u4_4thbit = ps_segments_x->u4_start_pos;

    for(i4_j = 0; i4_j < u1_num_sgmts_y; i4_j++)
    {
        UWORD8 i4_idx_a, i4_idx_b;
        UWORD8 u1_seg_ht, u1_seg_wd;
        UWORD8 u1_mb_adjoin_x, u1_mb_adjoin_y;
        WORD8 i1_nearst_mb_bdry_x, i1_nearst_mb_bdry_y;
        UWORD32 u4_num_valid_segs;
        WORD32 i4_idx_a_plus_ny, i4_idx_b_plus_nx, i4_index;
        WORD8 i1_yd_index, i1_xd_index;

        ps_seg_y_tmp = &ps_seg_desc_y[i4_j];
        i4_y = i4_y_start_pos + ps_seg_y_tmp->u1_seg_off;
        u1_seg_ht = ps_seg_y_tmp->u1_seg_dim;
        i1_yd_index = ps_seg_y_tmp->i1_dist_idx;
        i1_nearst_mb_bdry_y = ps_seg_y_tmp->i1_nearst_mb_bdry;
        u1_mb_adjoin_y = ps_seg_y_tmp->u1_mb_adjoin;
        i4_idx_a = pu1_ref_idx_y[i4_y];
        i4_idx_a_plus_ny = (i4_idx_a + i1_nearst_mb_bdry_y);

        /* Pack the availabilities of the next three horizontal MBs in 3bit
           format and 4th bit indicating if the start position is greater than the mb_width/2 */
        u4_lookup_4bit = u4_4thbit | u1_avail_map[i4_idx_a][i4_cur_x + 2] << 2 |
                         u1_avail_map[i4_idx_a][i4_cur_x + 1] << 1 |
                         u1_avail_map[i4_idx_a][i4_cur_x];

        u4_num_valid_segs = gu4_valid_segs_lookup[u4_lookup_4bit];
        i4_i = isvcd_left_most_bit_detect(u4_num_valid_segs);
        u4_num_valid_segs <<= (i4_i + 1);

        for(; i4_i < u1_num_sgmts_x; i4_i++)
        {
            ps_seg_x_tmp = &ps_seg_desc_x[i4_i];
            i4_x = i4_x_start_pos + ps_seg_x_tmp->u1_seg_off;
            i4_idx_b = pu1_ref_idx_x[i4_x];
            u1_seg_wd = ps_seg_x_tmp->u1_seg_dim;
            i1_xd_index = ps_seg_x_tmp->i1_dist_idx;
            i1_nearst_mb_bdry_x = ps_seg_x_tmp->i1_nearst_mb_bdry;
            u1_mb_adjoin_x = ps_seg_x_tmp->u1_mb_adjoin;
            i4_idx_b_plus_nx = (i4_idx_b + i1_nearst_mb_bdry_x);

            /* Find the avalability of (x,y-Yd),(x-Xd,y),(x-Xd,y-Yd) and pack it to 3 bits */
            u4_lookup_5bit = u1_avail_map[i4_idx_a_plus_ny][i4_idx_b_plus_nx] << 2 |
                             u1_avail_map[i4_idx_a_plus_ny][i4_idx_b] << 1 |
                             u1_avail_map[i4_idx_a][i4_idx_b_plus_nx] | u1_mb_adjoin_x |
                             u1_mb_adjoin_y;

            i4_corner_pixel_available = u1_avail_map[i4_idx_a_plus_ny][i4_idx_b_plus_nx];

            /* Function pointer table from lookup to get Left,Top,Bottom,Right,Diagonal padding */
            if(u4_lookup_5bit > 31)
            {
                u4_lookup_5bit = 0;
            }
            pf_intra_samp_padding = pf_intra_samp_lookup[u4_lookup_5bit];

            if(pf_intra_samp_padding != NULL)
            {
                pf_intra_samp_padding(i4_x, i4_y, i1_xd_index, i1_yd_index, u1_seg_wd, u1_seg_ht,
                                      pu1_refarray_1, pu1_refarray_2, i4_refarray_stride,
                                      u1_mb_adjoin_x, u1_mb_adjoin_y, i4_corner_pixel_available);
            }

            /* increment to the next unavailable segment */
            i4_index = isvcd_left_most_bit_detect(u4_num_valid_segs);
            u4_num_valid_segs <<= (i4_index + 1);
            i4_i += i4_index;

        } /* end of loop over ref array width */

    }     /* end of loop over ref array height */
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_ref_layer_mbtype                                */
/*                                                                           */
/*  Description   : This function is used to find the mb type of the         */
/*                    corresponding MB in the reference layer                */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra samp context                  */
/*                  pi1_ref_mb_modes : ref mb modes buffer pointer           */
/*                  i4_ref_mode_stride : mb mode buffer stride               */
/*                  i4_x_ref : reference location X                          */
/*                  i4_y_ref : reference location Y                          */
/*                  pi4_mb_type : pointer to store the mb type               */
/*                  i4_chroma_flag : chroma flag                             */
/*  Globals       : none                                                     */
/*  Processing    : it derives the bit corresponding to reference MB and     */
/*                  stores the mbtype as INTRA if the bit is set             */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD8 isvcd_get_ref_layer_mbtype(WORD8 *pi1_ref_mb_modes, WORD32 *pi4_mb_type,
                                 WORD8 i1_curr_slice_id, WORD8 i1_cons_intr_samp_flag)
{
    WORD8 i1_intra_slice_id;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
    WORD8 i1_mb_mode;

    ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes;
    i1_mb_mode = ps_inter_lyr_mb_prms->i1_mb_mode;

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        /* INTER */
        *pi4_mb_type = SVC_INTER_MB;
        i1_intra_slice_id = -1;
    }
    else
    {
        /* INTRA */
        *pi4_mb_type = SVC_INTRA_MB;
        i1_intra_slice_id = ps_inter_lyr_mb_prms->i1_slice_id;

        if(1 == i1_cons_intr_samp_flag)
        {
            /* check for different slice idc */
            if(ps_inter_lyr_mb_prms->i1_slice_id != i1_curr_slice_id)
            {
                /* store the mode as INTER (not available for upsampling) */
                *pi4_mb_type = SVC_INTER_MB;
            }
        }
    }

    /* if contarained intra flag is 1 then check for same mb mode */
    return (i1_intra_slice_id);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_reflayer_construction                               */
/*                                                                           */
/*  Description   :  This function constructs the reference array buffer     */
/*                    used for intra resampling of a component in an MB      */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt: intra sampling context               */
/*                  pu1_inp : input (reference layer data)                   */
/*                  i4_inp_stride : input buffer stride                      */
/*                  ps_ref_mb_mode_map : ref layer mb mode buffer desc       */
/*                  pi4_refarr_wd : pointer to store the reference array WD  */
/*                  pi4_refarr_ht : pointer to store the reference array HT  */
/*                  pi4_x_offset : pointer to store the reference X offset   */
/*                  pi4_y_offset : pointer to store the reference Y offset   */
/*                  ps_coord     : mb co-ordinate structure                  */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it fills the reference layer data if they are falling in */
/*                  INTRA MB region. If all the pixels are not filled  it    */
/*                  calls the border extension algorithm to fill them        */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*****************************************************************************/
WORD32 isvcd_reflayer_construction(void *pv_intra_samp_ctxt, UWORD8 *pu1_inp_1, UWORD8 *pu1_inp_2,
                                   WORD32 i4_inp_stride, WORD32 i4_refarray_stride,
                                   mem_element_t *ps_ref_mb_mode_map, mb_coord_t *ps_coord,
                                   WORD32 i4_chroma_flag)
{
    WORD32 i4_x, i4_y;

    /* --------------------------------------------------------------------- */
    /* Context and reference layer realted varaibles                         */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_map_ctxt_t *ps_map_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD8 *pi1_ref_mb_modes, *pi1_ref_mb_modes_bkp_1;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    ref_mb_map_t *ps_x_off_len;
    ref_mb_map_t *ps_y_off_len;
    WORD32 i4_mbaddr_y;
    WORD32 i4_mbaddr_x;
    WORD32 i4_mb_ht, i4_mb_wd;
    UWORD8 u1_map_buf[4][4] = {0}; /*!< 4x4 mb grid buffer to store the mb availablity */
    /* --------------------------------------------------------------------- */
    /* Temp Variables for Mapping context                                     */
    /* --------------------------------------------------------------------- */
    WORD32 i4_ref_wd;
    WORD32 i4_ref_ht;
    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refarray_wd;
    WORD32 i4_refarray_ht;
    WORD32 i4_mb_type;
    WORD8 i1_cons_intr_samp_flag;
    WORD8 i1_slice_id = 0;
    WORD32 i4_mb_wd_sft, i4_mb_ht_sft;

    /* --------------------------------------------------------------------- */
    /* Local Pointer Declaration for arrays in Mapping context                 */
    /* --------------------------------------------------------------------- */

    WORD32 i4_unfill_check;
    UWORD8 *pu1_refarray_1, *pu1_refarray_2;

    UNUSED(pu1_inp_2);

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode_map->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode_map->i4_element_size;

    /* get the condtrained intra sampling flag */
    i1_cons_intr_samp_flag = ps_lyr_ctxt->i1_constrained_intra_rsmpl_flag;

    if(NULL == pi1_ref_mb_modes)
    {
        return NOT_OK;
    }

    /* --------------------------------------------------------------------- */
    /* Based on Chroma and Luma, extracting the context information struct     */
    /* --------------------------------------------------------------------- */
    if(1 == i4_chroma_flag)
        ps_map_ctxt = &ps_lyr_ctxt->s_chroma_map_ctxt;
    else
        ps_map_ctxt = &ps_lyr_ctxt->s_luma_map_ctxt;

    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;

    /* --------------------------------------------------------------------- */
    /* Deriving the parameters required for further processing                 */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_base_width = ps_lyr_ctxt->i4_ref_width;
        WORD32 i4_base_height = ps_lyr_ctxt->i4_ref_height;

        i4_ref_wd = i4_base_width >> i4_chroma_flag;
        i4_ref_ht = i4_base_height >> i4_chroma_flag;
        i4_mb_wd_sft = (MB_WIDTH_SHIFT - i4_chroma_flag);
        i4_mb_ht_sft = (MB_HEIGHT_SHIFT - i4_chroma_flag);
    }

    /* --------------------------------------------------------------------- */
    /* Local variables based on the MB address                                 */
    /* --------------------------------------------------------------------- */
    i4_mbaddr_y = ps_coord->u2_mb_y;
    i4_mbaddr_x = ps_coord->u2_mb_x;
    i4_x_offset = ps_x_off_len[i4_mbaddr_x].i2_offset;
    i4_y_offset = ps_y_off_len[i4_mbaddr_y].i2_offset;
    i4_refarray_wd = ps_x_off_len[i4_mbaddr_x].i2_length;
    i4_refarray_ht = ps_y_off_len[i4_mbaddr_y].i2_length;
    i4_mb_wd = (MB_WIDTH >> i4_chroma_flag);
    i4_mb_ht = (MB_HEIGHT >> i4_chroma_flag);

    /* --------------------------------------------------------------------- */
    /* Derivation of ref slice MB idc                                         */
    /* --------------------------------------------------------------------- */
    if(1 == i1_cons_intr_samp_flag)
    {
        WORD32 i4_x_min, i4_x_max;
        WORD32 i4_y_min, i4_y_max;
        ref_min_max_map_t *ps_x_min_max;
        ref_min_max_map_t *ps_y_min_max;

        ps_x_min_max = ps_map_ctxt->ps_x_min_max;
        ps_y_min_max = ps_map_ctxt->ps_y_min_max;

        /* get the min and max positions */
        i4_x_min = ps_x_min_max[i4_mbaddr_x].i2_min_pos;
        i4_x_max = ps_x_min_max[i4_mbaddr_x].i2_max_pos;
        i4_y_min = ps_y_min_max[i4_mbaddr_y].i2_min_pos;
        i4_y_max = ps_y_min_max[i4_mbaddr_y].i2_max_pos;

        /* default initialization */
        i4_mb_type = SVC_INTER_MB;

        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_mb_x, i4_mb_y;

            i4_y_ref = (i4_y_min + 1) + i4_y_offset;
            i4_x_ref = (i4_x_min + 1) + i4_x_offset;
            i4_mb_x = (i4_x_ref >> i4_mb_wd_sft);
            i4_mb_y = (i4_y_ref >> i4_mb_ht_sft);
            pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;

            /* get the location of the byte which has the current mb mode */
            pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_element_size);
            pi1_ref_mb_modes += (i4_mb_x * i4_element_size);
        }

        for(i4_y = (i4_y_min + 1); i4_y <= (i4_y_max - 1);)
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_distleftX, i4_rangeX;
            WORD32 i4_disttopY, i4_rangeY;

            i4_y_ref = (i4_y + i4_y_offset);
            i4_disttopY = (i4_y_ref) & (i4_mb_ht - 1);
            i4_rangeY = (i4_mb_ht - i4_disttopY);

            pi1_ref_mb_modes_bkp_1 = pi1_ref_mb_modes;

            for(i4_x = (i4_x_min + 1); i4_x <= (i4_x_max - 1);)
            {
                i4_x_ref = (i4_x + i4_x_offset);
                i4_distleftX = (i4_x_ref) & (i4_mb_wd - 1);
                i4_rangeX = (i4_mb_wd - i4_distleftX);

                /* get the referecne layer mb type */
                i1_slice_id =
                    isvcd_get_ref_layer_mbtype(pi1_ref_mb_modes_bkp_1, &i4_mb_type, i1_slice_id, 0);
                if(SVC_INTRA_MB == i4_mb_type)
                {
                    /* if an Intra MB is returned then break the loop */
                    break;
                }

                i4_x += i4_rangeX;
                pi1_ref_mb_modes_bkp_1 += i4_element_size;
            } /* end of loop in horizontal direction */

            if(SVC_INTRA_MB == i4_mb_type)
            {
                /* if an Intra MB is returned then break the loop */
                break;
            }

            i4_y += i4_rangeY;
            pi1_ref_mb_modes += (i4_ref_mode_stride * i4_element_size);

        } /* end of loop in vertical direction */
    }
    else
    {
        /* set to non valid value */
        i1_slice_id = -1;
    }

    i4_unfill_check = 0;

    /* --------------------------------------------------------------------- */
    /* Copying the data from recon buffer to refSample Array.                */
    /* NOTE: The copying of the data from recon buffer to refSample Array    */
    /*       can be optimized by bring in data at N-MB level,thus taking     */
    /*       advantage of the overlapping data which now gets copied every MB*/
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_x_ref_start, i4_x_ref_end;
        WORD32 i4_y_ref_start, i4_y_ref_end;
        WORD32 i4_rangeW, i4_rangeH;
        WORD32 i4_offset;
        UWORD8 *pu1_src, *pu1_dst;
        UWORD8 *pu1_dst1, *pu1_dst2;

        /* Copy (refW x refH) dimension into reference sample array */
        i4_x_ref_start = MAX(0, MIN((i4_ref_wd - 1), i4_x_offset));
        i4_x_ref_end = MAX(0, MIN((i4_ref_wd - 1), (i4_refarray_wd - 1) + i4_x_offset));
        i4_y_ref_start = MAX(0, MIN((i4_ref_ht - 1), i4_y_offset));
        i4_y_ref_end = MAX(0, MIN((i4_ref_ht - 1), (i4_refarray_ht - 1) + i4_y_offset));

        /* find the actual data to be copied */
        i4_rangeW = (i4_x_ref_end - i4_x_ref_start + 1);
        i4_rangeH = (i4_y_ref_end - i4_y_ref_start + 1);

        /* get the reconbuffer pointer and ref sample array pointer */
        i4_offset =
            (i4_x_ref_start - i4_x_offset) + ((i4_y_ref_start - i4_y_offset) * i4_refarray_stride);

        if(0 == i4_chroma_flag)
        {
            pu1_refarray_1 = ps_ctxt->pu1_refarray_buffer;
            pu1_refarray_2 = NULL;
            pu1_src = pu1_inp_1;
            pu1_dst = pu1_refarray_1 + i4_offset;

            /* Copy luma data into refsample array */
            isvcd_copy_data(pu1_src, i4_inp_stride, pu1_dst, i4_refarray_stride, i4_rangeW,
                            i4_rangeH);
        }
        else
        {
            pu1_refarray_1 = ps_ctxt->pu1_refarray_buffer;
            pu1_refarray_2 = ps_ctxt->pu1_refarray_cb;
            pu1_src = pu1_inp_1;
            pu1_dst1 = pu1_refarray_1 + i4_offset;
            pu1_dst2 = pu1_refarray_2 + i4_offset;
            isvcd_copy_data_semiplanr(pu1_src, i4_inp_stride, pu1_dst1, pu1_dst2,
                                      i4_refarray_stride, i4_rangeW, i4_rangeH);
        }
    }

    /* --------------------------------------------------------------------- */
    /* Loop to fill ref sample array and corresponding map for interpolation */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_i, i4_j;
        UWORD8 *pu1_ref_idx_x, *pu1_ref_idx_y;
        WORD32 i4_x_ref;
        WORD32 i4_y_ref;
        WORD32 i4_mb_x, i4_mb_y;

        i4_y_ref = i4_y_offset;
        i4_x_ref = i4_x_offset;
        i4_mb_x = (i4_x_ref >> i4_mb_wd_sft);
        i4_mb_y = (i4_y_ref >> i4_mb_ht_sft);
        pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;

        /* get the location of the byte which has the current mb mode */
        pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_element_size);
        pi1_ref_mb_modes += (i4_mb_x * i4_element_size);
        pu1_ref_idx_x = ps_map_ctxt->pu1_refarray_x_idx;
        pu1_ref_idx_y = ps_map_ctxt->pu1_refarray_y_idx;

        i4_j = 0;
        for(i4_y = 0; i4_y < i4_refarray_ht;)
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_distleftX, i4_rangeX;
            WORD32 i4_disttopY, i4_rangeY;

            i4_y_ref = i4_y + i4_y_offset;
            i4_disttopY = (i4_y_ref) & (i4_mb_ht - 1);
            i4_rangeY = (i4_mb_ht - i4_disttopY);

            /* find the y-index lookup */
            memset(pu1_ref_idx_y, i4_j, i4_rangeY);
            pu1_ref_idx_y += i4_rangeY;

            i4_i = 0;
            pi1_ref_mb_modes_bkp_1 = pi1_ref_mb_modes;
            for(i4_x = 0; i4_x < i4_refarray_wd;)
            {
                i4_x_ref = i4_x + i4_x_offset;
                i4_distleftX = (i4_x_ref) & (i4_mb_wd - 1);
                i4_rangeX = (i4_mb_wd - i4_distleftX);

                if(0 == i4_j)
                {
                    /* find the x-index lookup */
                    memset(pu1_ref_idx_x, i4_i, i4_rangeX);
                    pu1_ref_idx_x += i4_rangeX;
                }

                /* get the referecne layer mb type */
                isvcd_get_ref_layer_mbtype(pi1_ref_mb_modes_bkp_1, &i4_mb_type, i1_slice_id,
                                           i1_cons_intr_samp_flag);

                if(SVC_INTRA_MB == i4_mb_type)
                {
                    u1_map_buf[i4_j][i4_i] = 1;
                }
                else
                {
                    i4_unfill_check = 1;
                }

                i4_x = i4_x + i4_rangeX;
                i4_i++;
                pi1_ref_mb_modes_bkp_1 += i4_element_size;

            } /* end of loop over ref array width */

            i4_j++;
            i4_y = i4_y + i4_rangeY;
            pi1_ref_mb_modes += (i4_ref_mode_stride * i4_element_size);
        } /* end of loop over ref array height */
    }

    /* --------------------------------------------------------------------- */
    /* Calling boundary extension algorithm to fill unfilled pixels             */
    /* --------------------------------------------------------------------- */
    if(i4_unfill_check == 1)
    {
        isvcd_fill_non_avail_pixel(ps_map_ctxt, pu1_refarray_1, pu1_refarray_2, i4_refarray_stride,
                                   ps_coord, i4_chroma_flag, u1_map_buf);
    }
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_reflayer_construction_dyadic                        */
/*                                                                           */
/*  Description   :  This function constructs the reference array buffer     */
/*                   for dyadic cases used for intra resampling of a         */
/*                     component in an MB                                    */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_mb_mode_map : ref layer mb mode buffer desc       */
/*                  pu1_inp_luma : luma input (reference layer data)         */
/*                    pu1_inp_chroma : chroma input (reference layer data)   */
/*                  i4_inp_luma_stride : luma input buffer stride            */
/*                    i4_inp_chroma_stride : chroma input buffer stride      */
/*                  i4_top : indicates whether the core 8x8 reference block  */
/*                             is one of 0 and 1 or one of 2 and 3           */
/*                    i4_left : indicates whether the core 8x8 ref block     */
/*                             is one of 0 and 2 or one of 1 and 3           */
/*                    ps_ref_mb_coord : coordinates of the reference MB      */
/*  Globals       : none                                                     */
/*  Processing    : it fills the reference layer data if they are falling in */
/*                  INTRA MB region. If all the pixels are not filled  it    */
/*                  calls the border extension algorithm to fill them        */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_reflayer_construction_dyadic(void *pv_intra_samp_ctxt,
                                          mem_element_t *ps_ref_mb_mode_map, UWORD8 *pu1_inp_luma,
                                          UWORD8 *pu1_inp_chroma, WORD32 i4_inp_luma_stride,
                                          WORD32 i4_inp_chroma_stride, WORD32 i4_top,
                                          WORD32 i4_left, UWORD16 u2_mb_x, UWORD16 u2_mb_y)
{
    /* Index variables */
    WORD32 i4_x, i4_y;
    WORD32 i4_x0, i4_y0;
    WORD32 i4_xc0, i4_yc0;
    WORD32 i4_ref_xD, i4_ref_yD;
    WORD32 i4_c_ref_xD, i4_c_ref_yD;

    /* --------------------------------------------------------------------- */
    /* Context and reference layer related variables                         */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD8 *pi1_ref_mb_modes;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    WORD32 i4_mbaddr_y;
    WORD32 i4_mbaddr_x;

    /* --------------------------------------------------------------------- */
    /* Temp Variables for Mapping context                                     */
    /* --------------------------------------------------------------------- */
    WORD32 i4_refarray_wd_luma, i4_refarray_wd_chroma;
    WORD32 i4_refarray_ht_luma, i4_refarray_ht_chroma;
    WORD32 i4_avlblty;
    WORD8 i1_cons_intr_samp_flag;
    WORD8 i1_slice_id;
    WORD8 i1_corner_samp_avlbl_flag;
    UWORD8 u1_ny_avlblty;

    /* --------------------------------------------------------------------- */
    /* Local Pointer Declaration for arrays in Mapping context                 */
    /* --------------------------------------------------------------------- */
    UWORD8 *pu1_refarray_luma;
    UWORD8 *pu1_refarray_cb, *pu1_refarray_cr;

    /* --------------------------------------------------------------------- */
    /* Derivation of local variables                                         */
    /* --------------------------------------------------------------------- */
    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode_map->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode_map->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode_map->i4_element_size;

    /* --------------------------------------------------------------------- */
    /* get the constrained intra resampling flag                             */
    /* --------------------------------------------------------------------- */
    i1_cons_intr_samp_flag = ps_lyr_ctxt->i1_constrained_intra_rsmpl_flag;
    if(NULL == pi1_ref_mb_modes)
    {
        return NOT_OK;
    }

    pu1_refarray_luma = ps_ctxt->pu1_refarray_buffer;
    pu1_refarray_cb = ps_ctxt->pu1_refarray_cb;
    pu1_refarray_cr = ps_ctxt->pu1_refarray_cr;

    /* --------------------------------------------------------------------- */
    /* Get the coordinates of the reference layer MB                         */
    /* --------------------------------------------------------------------- */
    i4_mbaddr_x = u2_mb_x;
    i4_mbaddr_y = u2_mb_y;

    /* --------------------------------------------------------------------- */
    /* Getting the size of the valid area of ref array to be brought in         */
    /* --------------------------------------------------------------------- */
    i4_refarray_wd_luma = 20;
    i4_refarray_ht_luma = 20;
    i4_refarray_wd_chroma = i4_refarray_wd_luma >> 1;
    i4_refarray_ht_chroma = i4_refarray_ht_luma >> 1;

    /* --------------------------------------------------------------------- */
    /* Derivation of ref slice MB idc                                         */
    /* --------------------------------------------------------------------- */
    if(1 == i1_cons_intr_samp_flag)
    {
        inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
        WORD8 *pi1_ref_mb_mode_tmp;
        WORD8 i1_mb_mode;

        /* get the location of the byte which has the current mb mode */
        pi1_ref_mb_mode_tmp = pi1_ref_mb_modes;
        pi1_ref_mb_mode_tmp += (i4_mbaddr_y * i4_ref_mode_stride * i4_element_size);
        pi1_ref_mb_mode_tmp += (i4_mbaddr_x * i4_element_size);
        ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_mode_tmp;
        i1_mb_mode = ps_inter_lyr_mb_prms->i1_mb_mode;

        /* The reference layer MB should be intra */
        UNUSED(i1_mb_mode);

        i1_slice_id = ps_inter_lyr_mb_prms->i1_slice_id;
    }
    else
    {
        /* set to non valid value */
        i1_slice_id = -1;
    }

    /* --------------------------------------------------------------------- */
    /* Bring in the reference array                                          */
    /* --------------------------------------------------------------------- */
    {
        UWORD8 *pu1_src, *pu1_dst;
        WORD32 i4_src_stride, i4_dst_stride;

        /* Copy luma */
        i4_src_stride = i4_inp_luma_stride;
        i4_dst_stride = DYADIC_REF_W_Y;
        pu1_src = pu1_inp_luma;
        pu1_dst = pu1_refarray_luma;
        isvcd_copy_data(pu1_src, i4_src_stride, pu1_dst, i4_dst_stride, i4_refarray_wd_luma,
                        i4_refarray_ht_luma);
        // Semi planar
        i4_src_stride = i4_inp_chroma_stride;
        i4_dst_stride = DYADIC_REF_W_C;
        pu1_src = pu1_inp_chroma;
        isvcd_copy_data_semiplanr(pu1_src, i4_src_stride, pu1_refarray_cb, pu1_refarray_cr,
                                  i4_dst_stride, i4_refarray_wd_chroma, i4_refarray_ht_chroma);
    }

    /* --------------------------------------------------------------------- */
    /* Get the availability of 5 neighboring MBs                             */
    /* --------------------------------------------------------------------- */
    {
        /* mb_x + left, mb_y + top */
        isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                           i4_mbaddr_x + i4_left, i4_mbaddr_y + i4_top, &i4_avlblty,
                                           i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty = i4_avlblty;

        /* mb_x + left, mb_y */
        isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                           i4_mbaddr_x + i4_left, i4_mbaddr_y, &i4_avlblty,
                                           i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 1);

        /* mb_x, mb_y + top */
        isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                           i4_mbaddr_x, i4_mbaddr_y + i4_top, &i4_avlblty,
                                           i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 2);

        /* mb_x - left, mb_y + top */
        isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                           i4_mbaddr_x - i4_left, i4_mbaddr_y + i4_top, &i4_avlblty,
                                           i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 3);

        /* mb_x + left, mb_y - top */
        isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                           i4_mbaddr_x + i4_left, i4_mbaddr_y - i4_top, &i4_avlblty,
                                           i1_slice_id, i1_cons_intr_samp_flag);
        u1_ny_avlblty += (i4_avlblty << 4);
    }

    /* --------------------------------------------------------------------- */
    /* Filling the unavailable samples, if any                                 */
    /* --------------------------------------------------------------------- */
    if(0x7 == u1_ny_avlblty)
    {
        /* All are available, exit */
        return OK;
    }

    if(!(u1_ny_avlblty & 0x7))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Set the 4 corner samples to (x-xD,y-yD) */
        i4_x0 = 9 + (i4_left << 3) + i4_left;
        i4_y0 = 9 + (i4_top << 3) + i4_top;

        i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
        i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;
        pu1_tmp_dst1[i4_x0] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst2[i4_x0] = pu1_tmp_src[i4_ref_xD];
        pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];

        /* Set the corner sample of Cb and Cr to (x-xD,y-yD) */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;
        i4_c_ref_yD = i4_ref_yD >> 1;
        i4_c_ref_xD = i4_ref_xD >> 1;
        pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
        pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];
    }

    if(!(u1_ny_avlblty & 0x5))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Copy (x0,ref_yD), (x0+1,ref_yD), ..., (x0+7,ref_yD) to */
        /* (x0,y0), (x0+1,y0), ..., (x0+7,y0) and   */
        /* (x0,y0+1), (x0+1,y0+1), ..., (x0+7,y0+1) */
        i4_x0 = 2;
        i4_y0 = 9 + (i4_top << 3) + i4_top;
        if(i4_left > 0)
        {
            i4_x0 += 8;
        }
        i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;

        for(i4_x = i4_x0; i4_x < i4_x0 + 8; i4_x++)
        {
            pu1_tmp_dst1[i4_x] = pu1_tmp_src[i4_x];
            pu1_tmp_dst2[i4_x] = pu1_tmp_src[i4_x];
        }

        /* Cb and Cr copy */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;
        i4_c_ref_yD = i4_ref_yD >> 1;
        pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);

        for(i4_x = i4_xc0; i4_x < i4_xc0 + 4; i4_x++)
        {
            pu1_tmp_dst1[i4_x] = pu1_tmp_src1[i4_x];
            pu1_tmp_dst2[i4_x] = pu1_tmp_src2[i4_x];
        }
    }

    if(!(u1_ny_avlblty & 0x3))
    {
        UWORD8 *pu1_tmp_src, *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        /* Copy (ref_xD,y0) to (x0,y0) and (x0+1,y0); */
        /* copy (ref_xD,y0+1) to (x0,y0+1) and (x0+1,y0+1); ... ;*/
        /* copy (ref_xD,y0+7) to (x0,y0+7) and (x0+1,y0+7) */
        i4_x0 = 9 + (i4_left << 3) + i4_left;
        i4_y0 = 2;
        if(i4_top > 0)
        {
            i4_y0 += 8;
        }
        i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

        pu1_tmp_src = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
        pu1_tmp_dst1 = pu1_tmp_src;

        for(i4_y = i4_y0; i4_y < i4_y0 + 8; i4_y++)
        {
            pu1_tmp_dst1[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_src += DYADIC_REF_W_Y;
            pu1_tmp_dst1 += DYADIC_REF_W_Y;
        }

        /* Cb and Cr copy */
        i4_xc0 = i4_x0 >> 1;
        i4_yc0 = i4_y0 >> 1;
        i4_c_ref_xD = i4_ref_xD >> 1;
        pu1_tmp_src1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst1 = pu1_tmp_src1;
        pu1_tmp_src2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
        pu1_tmp_dst2 = pu1_tmp_src2;

        for(i4_y = i4_yc0; i4_y < i4_yc0 + 4; i4_y++)
        {
            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];
            pu1_tmp_src1 += DYADIC_REF_W_C;
            pu1_tmp_src2 += DYADIC_REF_W_C;
            pu1_tmp_dst1 += DYADIC_REF_W_C;
            pu1_tmp_dst2 += DYADIC_REF_W_C;
        }
    }

    if(!(u1_ny_avlblty & 0x4))
    {
        if(!(u1_ny_avlblty & 0x8))
        {
            /* (mb_x-left,mb_y+top) not available */
            UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

            i4_x0 = 9 - i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;

            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            /* Copy (x0,ref_yD) and (x0+1,ref_yD) to (x0,y0) and (x0+1,y0), and */
            /* to (x0,y0+1) and (x0+1,y0+1) */
            pu1_tmp_src = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
            pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_x0];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_x0 + 1];
            pu1_tmp_dst += DYADIC_REF_W_Y;
            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_x0];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_x0 + 1];

            /* Cb copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            pu1_tmp_src = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_xc0];

            /* Cr copy */
            pu1_tmp_src = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_xc0];

        } /* if (mb_x-left,mb_y+top) not available */
        else
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride,
                                               i4_element_size, i4_mbaddr_x - i4_left, i4_mbaddr_y,
                                               &i4_avlblty, i1_slice_id, i1_cons_intr_samp_flag);
            i1_corner_samp_avlbl_flag = i4_avlblty;

            i4_x0 = 9 - i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);
            i4_ref_xD = i4_x0 - (i4_left * 7) - (i4_left >> 1);
            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;
            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Fill corner sample if not available */
            if(!i1_corner_samp_avlbl_flag)
            {
                isvcd_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                         pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call diagonal construction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvcd_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                    DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                            DYADIC_REF_W_C);

            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                            DYADIC_REF_W_C);
        }
    }

    if(!(u1_ny_avlblty & 0x2))
    {
        if(!(u1_ny_avlblty & 0x10))
        {
            UWORD8 *pu1_tmp_src, *pu1_tmp_dst;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 - i4_top;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

            /* Copy (ref_xD,y0) to (x0,y0), (x0+1,y0), and  */
            /* copy (ref_xD,y0+1) to (x0,y0+1), (x0+1,y0+1) */
            pu1_tmp_src = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst = pu1_tmp_src;
            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_src += DYADIC_REF_W_Y;
            pu1_tmp_dst += DYADIC_REF_W_Y;
            pu1_tmp_dst[i4_x0] = pu1_tmp_src[i4_ref_xD];
            pu1_tmp_dst[i4_x0 + 1] = pu1_tmp_src[i4_ref_xD];

            /* Cb copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;
            pu1_tmp_src = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_tmp_src;
            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_c_ref_xD];

            /* Cr copy */
            pu1_tmp_src = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst = pu1_tmp_src;
            pu1_tmp_dst[i4_xc0] = pu1_tmp_src[i4_c_ref_xD];

        } /* if (mb_x+left,mb_y-top) not available */
        else
        {
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            isvcd_get_ref_layer_avlblty_dyadic(pi1_ref_mb_modes, i4_ref_mode_stride,
                                               i4_element_size, i4_mbaddr_x, i4_mbaddr_y - i4_top,
                                               &i4_avlblty, i1_slice_id, i1_cons_intr_samp_flag);
            i1_corner_samp_avlbl_flag = i4_avlblty;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 - i4_top;
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
            i4_ref_yD = i4_y0 - (i4_top * 7) - (i4_top >> 1);
            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;
            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            if(!i1_corner_samp_avlbl_flag)
            {
                isvcd_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                         pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call diagonal consrtuction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvcd_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                    DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                            DYADIC_REF_W_C);

            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                            DYADIC_REF_W_C);
        }
    }

    if(u1_ny_avlblty & 1)
    {
        if(!(u1_ny_avlblty & 2))
        {
            /* (mb_x+left,mb_y) is unavailable */
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;
            UWORD8 *pu1_tmp_dst;
            UWORD8 u1_filled_samp;

            i1_corner_samp_avlbl_flag = (u1_ny_avlblty & 4) >> 2;

            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 2;
            i4_ref_yD = 1;
            if(i4_top > 0)
            {
                i4_y0 += 8;
                i4_ref_yD = 18;
            }

            i4_ref_xD = i4_x0 - (i4_left) - (i4_left >> 1);
            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Fill corner sample if unavailable */
            if(!i1_corner_samp_avlbl_flag)
            {
                isvcd_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                         pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call the diagonal construction for the 8 rows */
            if(i4_top == i4_left)
            {
                /* if top * left = 1 (x0,y0) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0, i4_y0, i4_xD, i4_yD,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);

                /* (x0,y0+1), ..., (x0,y0+7) and */
                /* (x0+1,y0), ..., (x0+1,y0+6)   */
                for(i4_y = i4_y0 + 1; i4_y < i4_y0 + 8; i4_y++)
                {
                    i4_yD++;
                    u1_filled_samp = isvcd_diagonal_construct_dyadic(
                        i4_x0, i4_y, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
                    pu1_tmp_dst[i4_x0 + 1] = u1_filled_samp;
                    pu1_tmp_dst += DYADIC_REF_W_Y;
                }

                /* (x0+1,y0+7) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(
                    i4_x0 + 1, i4_y0 + 7, i4_xD + 1, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
            }
            else
            {
                /* top * left = -1 (x0+1,y0) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0 + 1, i4_y0, i4_xD + 1, i4_yD,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);

                /* (x0,y0), ..., (x0,y0+6) and   */
                /* (x0+1,y0+1), ..., (x0+1,y0+7) */
                for(i4_y = i4_y0; i4_y < i4_y0 + 7; i4_y++)
                {
                    u1_filled_samp = isvcd_diagonal_construct_dyadic(
                        i4_x0, i4_y, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);

                    pu1_tmp_dst += DYADIC_REF_W_Y;
                    pu1_tmp_dst[i4_x0 + 1] = u1_filled_samp;
                    i4_yD++;
                }

                /* (x0,y0+7) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0, i4_y0 + 7, i4_xD, i4_yD,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);
            }

            /* For Cb and Cr */
            for(i4_y = i4_yc0; i4_y < i4_yc0 + 4; i4_y++)
            {
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_xc0, i4_y, i4_c_xD, i4_c_yD,
                                                                 pu1_refarray_cb, DYADIC_REF_W_C);
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_xc0, i4_y, i4_c_xD, i4_c_yD,
                                                                 pu1_refarray_cr, DYADIC_REF_W_C);
                i4_c_yD++;
            }

        } /* (mb_x+left,mb_y) is unavailable */

        if(!(u1_ny_avlblty & 4))
        {
            /* (mb_x,mb_y+top) is unavailable */
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;
            UWORD8 *pu1_tmp_dst;
            UWORD8 u1_filled_samp;

            i1_corner_samp_avlbl_flag = (u1_ny_avlblty & 2) >> 1;
            i4_y0 = 9 + (i4_top << 3) + (i4_top);
            i4_x0 = 2;
            i4_ref_xD = 1;
            if(i4_left > 0)
            {
                i4_x0 += 8;
                i4_ref_xD = 18;
            }

            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);
            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            if(!i1_corner_samp_avlbl_flag)
            {
                isvcd_corner_samp_dyadic(i4_x0, i4_y0, i4_xD, i4_yD, pu1_refarray_luma,
                                         pu1_refarray_cb, pu1_refarray_cr);
            }

            /* Call the diagonal construction for the 2 rows */
            if(i4_top == i4_left)
            {
                /* if top * left = 1 (x0,y0) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0, i4_y0, i4_xD, i4_yD,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + ((i4_y0 + 1) * DYADIC_REF_W_Y);

                /* (x0+1,y0), ..., (x0+7,y0) and */
                /* (x0,y0+1), ..., (x0+6,y0+1)   */
                for(i4_x = i4_x0 + 1; i4_x < i4_x0 + 8; i4_x++)
                {
                    i4_xD++;
                    u1_filled_samp = isvcd_diagonal_construct_dyadic(
                        i4_x, i4_y0, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);
                    pu1_tmp_dst[i4_x - 1] = u1_filled_samp;
                }

                /* (x0+7,y0+1) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(
                    i4_x0 + 7, i4_y0 + 1, i4_xD, i4_yD + 1, pu1_refarray_luma, DYADIC_REF_W_Y);
            }
            else
            {
                /* top * left = -1 */
                /* (x0,y0+1) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0, i4_y0 + 1, i4_xD, i4_yD + 1,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);

                pu1_tmp_dst = pu1_refarray_luma + ((i4_y0 + 1) * DYADIC_REF_W_Y);

                /* (x0,y0), ..., (x0,y0+6) and   */
                /* (x0+1,y0+1), ..., (x0+1,y0+7) */
                for(i4_x = i4_x0; i4_x < i4_x0 + 7; i4_x++)
                {
                    u1_filled_samp = isvcd_diagonal_construct_dyadic(
                        i4_x, i4_y0, i4_xD, i4_yD, pu1_refarray_luma, DYADIC_REF_W_Y);

                    pu1_tmp_dst[i4_x + 1] = u1_filled_samp;
                    i4_xD++;
                }

                /* (x0+7,y0) */
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x0 + 7, i4_y0, i4_xD, i4_yD,
                                                                 pu1_refarray_luma, DYADIC_REF_W_Y);
            }

            /* For Cb and Cr */
            for(i4_x = i4_xc0; i4_x < i4_xc0 + 4; i4_x++)
            {
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x, i4_yc0, i4_c_xD, i4_c_yD,
                                                                 pu1_refarray_cb, DYADIC_REF_W_C);
                u1_filled_samp = isvcd_diagonal_construct_dyadic(i4_x, i4_yc0, i4_c_xD, i4_c_yD,
                                                                 pu1_refarray_cr, DYADIC_REF_W_C);
                i4_c_xD++;
            }

        } /* (mb_x,mb_y+top) is unavailable */
    }     /* if (mb_x+left,mb_y+top) not available */
    else
    {
        UWORD8 *pu1_tmp_dst1, *pu1_tmp_dst2;
        UWORD8 *pu1_tmp_src1, *pu1_tmp_src2;

        if(0x02 == (u1_ny_avlblty & 0x6))
        {
            /* (mb_x+left,mb_y) available, (mb_x,mb_y+top) unavailable */
            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);

            /* Copy (x0,ref_yD), (x0+1,ref_yD) to  */
            /* (x0,y0), (x0+1,y0), and (x0,y0+1), (x0+1,y0+1) */
            pu1_tmp_src1 = pu1_refarray_luma + (i4_ref_yD * DYADIC_REF_W_Y);
            pu1_tmp_dst1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst2 = pu1_tmp_dst1 + DYADIC_REF_W_Y;
            pu1_tmp_dst1[i4_x0] = pu1_tmp_src1[i4_x0];
            pu1_tmp_dst2[i4_x0] = pu1_tmp_src1[i4_x0];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src1[i4_x0 + 1];
            pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src1[i4_x0 + 1];

            /* Cb and Cr copy */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            pu1_tmp_src1 = pu1_refarray_cb + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_src2 = pu1_refarray_cr + (i4_c_ref_yD * DYADIC_REF_W_C);
            pu1_tmp_dst2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_xc0];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_xc0];

        } /* if (mb_x+left,mb_y) available, (mb_x,mb_y+top) unavailable */
        else if(0x04 == (u1_ny_avlblty & 0x6))
        {
            /* (mb_x+left,mb_y) unavailable, (mb_x,mb_y+top) available */
            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);

            /* Copy (ref_xD,y0) to (x0,y0) and (x0+1,y0) */
            /* copy (ref_xD,y0+1) to (x0,y0+1) and (x0+1,y0+1) */
            pu1_tmp_src1 = pu1_refarray_luma + (i4_y0 * DYADIC_REF_W_Y);
            pu1_tmp_dst1 = pu1_tmp_src1;
            pu1_tmp_src2 = pu1_tmp_src1 + DYADIC_REF_W_Y;
            pu1_tmp_dst2 = pu1_tmp_src2;

            pu1_tmp_dst1[i4_x0] = pu1_tmp_src1[i4_ref_xD];
            pu1_tmp_dst1[i4_x0 + 1] = pu1_tmp_src1[i4_ref_xD];
            pu1_tmp_dst2[i4_x0] = pu1_tmp_src2[i4_ref_xD];
            pu1_tmp_dst2[i4_x0 + 1] = pu1_tmp_src2[i4_ref_xD];

            /* Copy Cb and Cr */
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;

            pu1_tmp_src1 = pu1_refarray_cb + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst1 = pu1_tmp_src1;
            pu1_tmp_src2 = pu1_refarray_cr + (i4_yc0 * DYADIC_REF_W_C);
            pu1_tmp_dst2 = pu1_tmp_src2;

            pu1_tmp_dst1[i4_xc0] = pu1_tmp_src1[i4_c_ref_xD];
            pu1_tmp_dst2[i4_xc0] = pu1_tmp_src2[i4_c_ref_xD];

        } /* if (mb_x+left,mb_y) unavailable, (mb_x,mb_y+top) available */
        else if(0x6 == (u1_ny_avlblty & 0x6))
        {
            /* (mb_x+left,mb_y) available, (mb_x,mb_y+top) available */
            WORD32 i4_xD, i4_yD;
            WORD32 i4_c_xD, i4_c_yD;

            i4_y0 = 9 + (i4_top << 3) + i4_top;
            i4_x0 = 9 + (i4_left << 3) + i4_left;
            i4_ref_xD = i4_x0 - i4_left - (i4_left >> 1);
            i4_ref_yD = i4_y0 - i4_top - (i4_top >> 1);
            i4_xD = i4_x0 - i4_ref_xD;
            i4_yD = i4_y0 - i4_ref_yD;
            i4_xc0 = i4_x0 >> 1;
            i4_yc0 = i4_y0 >> 1;
            i4_c_ref_xD = i4_ref_xD >> 1;
            i4_c_ref_yD = i4_ref_yD >> 1;
            i4_c_xD = i4_xc0 - i4_c_ref_xD;
            i4_c_yD = i4_yc0 - i4_c_ref_yD;

            /* Call diagonal construction for luma */
            for(i4_y = i4_y0; i4_y < i4_y0 + 2; i4_y++)
            {
                for(i4_x = i4_x0; i4_x < i4_x0 + 2; i4_x++)
                {
                    isvcd_diagonal_construct_dyadic(i4_x, i4_y, i4_xD, i4_yD, pu1_refarray_luma,
                                                    DYADIC_REF_W_Y);
                    i4_xD++;
                }
                i4_yD++;
                i4_xD -= 2;
            }

            /* Call diagonal construction for chroma */
            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cb,
                                            DYADIC_REF_W_C);

            isvcd_diagonal_construct_dyadic(i4_xc0, i4_yc0, i4_c_xD, i4_c_yD, pu1_refarray_cr,
                                            DYADIC_REF_W_C);

        } /* if (mb_x+left,mb_y) available, (mb_x,mb_y+top) available */
    }     /* (mb_x+left,mb_y+top) available */

    return OK;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_base_luma_dyadic                        */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  intra resampling for dyadic scaling ratios               */
/*  Inputs        : pu1_inp_buf : ptr to the 12x12 reference sample buffer   */
/*                    pi2_tmp_filt_buf : ptr to the 12x16 buffer to hold the */
/*                        vertically interpolated data                       */
/*                  pu1_out_buf : output buffer pointer                      */
/*                  i4_out_stride : output buffer stride                     */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_base_luma_dyadic(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                        UWORD8 *pu1_out_buf, WORD32 i4_out_stride)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1, i4_samp_2, i4_samp_3;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp, *pu1_out;
    WORD16 *pi2_tmp;

    /* Filter coefficient values for phase 4 */
    i4_coeff_0 = -3;
    i4_coeff_1 = 28;
    i4_coeff_2 = 8;
    i4_coeff_3 = -1;
    i4_filt_stride = 12;
    i4_src_stride = DYADIC_REF_W_Y;
    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    pu1_out = pu1_out_buf;

    /* Vertical interpolation */
    for(i4_x = 0; i4_x < 12; i4_x++)
    {
        /* y = 0, y_phase = 12 */
        i4_samp_0 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_2 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_3 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        /* since y_phase 12 for y = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_3;
        i4_rslt_1 += i4_samp_1 * i4_coeff_2;
        i4_rslt_1 += i4_samp_2 * i4_coeff_1;
        i4_rslt_1 += i4_samp_3 * i4_coeff_0;

        /* Store the output */
        pi2_tmp[i4_x] = i4_rslt_1;
        /* Increment the output ptr */
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 15; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = i4_samp_2;
            i4_samp_2 = i4_samp_3;
            i4_samp_3 = pu1_inp[i4_x];

            /* y_phase is 4 for odd values of y */
            /* and 12 for even values of y    */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;
            i4_rslt_1 += i4_samp_2 * i4_coeff_2;
            i4_rslt_1 += i4_samp_3 * i4_coeff_3;
            i4_rslt_2 = i4_samp_0 * i4_coeff_3;
            i4_rslt_2 += i4_samp_1 * i4_coeff_2;
            i4_rslt_2 += i4_samp_2 * i4_coeff_1;
            i4_rslt_2 += i4_samp_3 * i4_coeff_0;

            /* Storing the results */
            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;

            /* Incrementing the pointers */
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* y = 15, y_phase = 4 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = i4_samp_2;
        i4_samp_2 = i4_samp_3;
        i4_samp_3 = pu1_inp[i4_x];

        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;
        i4_rslt_1 += i4_samp_2 * i4_coeff_2;
        i4_rslt_1 += i4_samp_3 * i4_coeff_3;

        /* Store the output */
        pi2_tmp[i4_x] = i4_rslt_1;

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;
    }

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 16; i4_y++)
    {
        /* x = 0, x_phase = 12 */
        i4_samp_0 = *pi2_tmp++;
        i4_samp_1 = *pi2_tmp++;
        i4_samp_2 = *pi2_tmp++;
        i4_samp_3 = *pi2_tmp++;

        /* since x_phase 12 for x = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_3;
        i4_rslt_1 += i4_samp_1 * i4_coeff_2;
        i4_rslt_1 += i4_samp_2 * i4_coeff_1;
        i4_rslt_1 += i4_samp_3 * i4_coeff_0;
        i4_rslt_1 += 512;

        i4_rslt_1 >>= 10;

        /* Store the output */
        pu1_out[0] = CLIPUCHAR(i4_rslt_1);

        for(i4_x = 1; i4_x < 15; i4_x += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = i4_samp_2;
            i4_samp_2 = i4_samp_3;
            i4_samp_3 = *pi2_tmp++;

            /* x_phase is 4 for odd values of x */
            /* and 12 for even values of x    */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;
            i4_rslt_1 += i4_samp_2 * i4_coeff_2;
            i4_rslt_1 += i4_samp_3 * i4_coeff_3;
            i4_rslt_1 += 512;

            i4_rslt_2 = i4_samp_0 * i4_coeff_3;
            i4_rslt_2 += i4_samp_1 * i4_coeff_2;
            i4_rslt_2 += i4_samp_2 * i4_coeff_1;
            i4_rslt_2 += i4_samp_3 * i4_coeff_0;
            i4_rslt_2 += 512;

            i4_rslt_1 >>= 10;
            i4_rslt_2 >>= 10;

            /* Store the output */
            pu1_out[i4_x] = CLIPUCHAR(i4_rslt_1);
            pu1_out[i4_x + 1] = CLIPUCHAR(i4_rslt_2);
        }

        /* x = 15 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = i4_samp_2;
        i4_samp_2 = i4_samp_3;
        i4_samp_3 = *pi2_tmp++;

        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;
        i4_rslt_1 += i4_samp_2 * i4_coeff_2;
        i4_rslt_1 += i4_samp_3 * i4_coeff_3;
        i4_rslt_1 += 512;

        i4_rslt_1 >>= 10;

        /* Store the output */
        pu1_out[i4_x] = CLIPUCHAR(i4_rslt_1);

        /* Increment the output ptr */
        pu1_out += i4_out_stride;

    } /* End of loop over y */
} /* isvcd_interpolate_base_luma_dyadic */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_1                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                    chroma for the following ref_lyr_chroma_phase_y_plus1  */
/*                    chroma_phase_y_plus1:                                  */
/*                        ref_lyr        cur_lyr                             */
/*                            0            0                                 */
/*                            1            0                                 */
/*                            1            1                                 */
/*                            1            2                                 */
/*                            2            1                                 */
/*                            2            2                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_1(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                         WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_inp = pu1_inp_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;

    /* Vertical interpolation */
    for(i4_x = 0; i4_x < 6; i4_x++)
    {
        /* y = 0, y_phase = phase_0 */
        i4_samp_0 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        /* since y_phase = phase_0 for y = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;

        /* Store the output */
        pi2_tmp[i4_x] = i4_rslt_1;

        /* Increment the output ptr */
        pi2_tmp += i4_filt_stride;

        for(i4_y = 1; i4_y < 7; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = pu1_inp[i4_x];

            /* y_phase is phase_1 for odd values of y */
            /* and phase_0 for even values of y          */
            i4_rslt_1 = i4_samp_0 * i4_coeff_2;
            i4_rslt_1 += i4_samp_1 * i4_coeff_3;
            i4_rslt_2 = i4_samp_0 * i4_coeff_0;
            i4_rslt_2 += i4_samp_1 * i4_coeff_1;

            /* Storing the results */
            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;

            /* Incrementing the pointers */
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* y = 7, y_phase = phase_1 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = pu1_inp[i4_x];

        i4_rslt_1 = i4_samp_0 * i4_coeff_2;
        i4_rslt_1 += i4_samp_1 * i4_coeff_3;

        /* Store the output */
        pi2_tmp[i4_x] = i4_rslt_1;

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;

    } /* End of loop over x */
} /* isvcd_vert_interpol_chroma_dyadic_1 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_2                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                    chroma for the following ref_lyr_chroma_phase_y_plus1  */
/*                    chroma_phase_y_plus1:                                  */
/*                        ref_lyr        cur_lyr                             */
/*                            0            1                                 */
/*                            0            2                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_2(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                         WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf + i4_src_stride;

    /* Vertical interpolation */
    for(i4_x = 0; i4_x < 6; i4_x++)
    {
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        for(i4_y = 0; i4_y < 8; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = pu1_inp[i4_x];

            /* y_phase is phase_1 for odd values of y and phase_0 for even values of y */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;

            i4_rslt_2 = i4_samp_0 * i4_coeff_2;
            i4_rslt_2 += i4_samp_1 * i4_coeff_3;

            /* Storing the results */
            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;

            /* Incrementing the pointers */
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf + i4_src_stride;
        pi2_tmp = pi2_tmp_filt_buf;

    } /* End of loop over x */
} /* isvcd_vert_interpol_chroma_dyadic_2 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_vert_interpol_chroma_dyadic_3                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  vertical intra resampling for dyadic scaling ratios for  */
/*                    chroma for the following ref_lyr_chroma_phase_y_plus1  */
/*                    chroma_phase_y_plus1:                                  */
/*                        ref_lyr        cur_lyr                             */
/*                            2            0                                 */
/*  Inputs        : pu1_inp_buf : ptr to the 6x6 reference sample buffer     */
/*                    pi2_tmp_filt_buf : ptr to the 6x8 buffer to hold the   */
/*                        vertically interpolated data                       */
/*                    i4_phase_0 : y phase for even values of y              */
/*                    i4_phase_1 : y phase for odd values of y               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : vertically resampled samples                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_vert_interpol_chroma_dyadic_3(UWORD8 *pu1_inp_buf, WORD16 *pi2_tmp_filt_buf,
                                         WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_filt_stride, i4_src_stride;
    UWORD8 *pu1_inp;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pi2_tmp = pi2_tmp_filt_buf;
    i4_filt_stride = 6;
    i4_src_stride = DYADIC_REF_W_C;
    pu1_inp = pu1_inp_buf;

    /* Vertical interpolation */
    for(i4_x = 0; i4_x < 6; i4_x++)
    {
        i4_samp_1 = pu1_inp[i4_x];
        pu1_inp += i4_src_stride;

        for(i4_y = 0; i4_y < 8; i4_y += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = pu1_inp[i4_x];

            /* y_phase is phase_1 for odd values of y */
            /* and phase_0 for even values of y          */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;

            i4_rslt_2 = i4_samp_0 * i4_coeff_2;
            i4_rslt_2 += i4_samp_1 * i4_coeff_3;

            /* Storing the results */
            pi2_tmp[i4_x] = i4_rslt_1;
            pi2_tmp += i4_filt_stride;
            pi2_tmp[i4_x] = i4_rslt_2;

            /* Incrementing the pointers */
            pi2_tmp += i4_filt_stride;
            pu1_inp += i4_src_stride;

        } /* End of loop over y */

        /* Reinitializing the ptrs */
        pu1_inp = pu1_inp_buf;
        pi2_tmp = pi2_tmp_filt_buf;

    } /* End of loop over x */
} /* isvcd_vert_interpol_chroma_dyadic_3 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_1                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  horizontal intra resampling for dyadic scaling ratios for*/
/*                    chroma with following ref_lyr_chroma_phase_x_plus1_flag*/
/*                    and chroma_phase_x_plus1_flag:                         */
/*                        ref_lyr        cur_lyr                             */
/*                            0            0                                 */
/*                            1            0                                 */
/*                            1            1                                 */
/*  Inputs        : pi2_tmp_filt_buf : ptr to the 6x8 buffer containing the  */
/*                        vertically interpolated data                       */
/*                    pu1_out_buf : pointer to the output buffer             */
/*                    i4_out_stride : output buffer stride                   */
/*                    i4_phase_0 : x phase for even values of x              */
/*                    i4_phase_1 : x phase for odd values of x               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : resampled samples                                        */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_horz_interpol_chroma_dyadic_1(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                         WORD32 i4_out_stride, WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_dst_stride;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf;
    i4_dst_stride = i4_out_stride;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y++)
    {
        /* x = 0, x_phase = phase_0 */
        i4_samp_0 = *pi2_tmp++;
        i4_samp_1 = *pi2_tmp++;

        /* since x_phase = phase_0 for x = 0 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_0;
        i4_rslt_1 += i4_samp_1 * i4_coeff_1;

        /* Round to 8-bit value */
        i4_rslt_1 += 32;
        i4_rslt_1 >>= 6;

        /* Store the output */
        pu1_out[0] = i4_rslt_1;

        for(i4_x = 1; i4_x < 7; i4_x += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = *pi2_tmp++;

            /* x_phase is phase_1 for odd values of x and phase_0 for even values of x */
            i4_rslt_1 = i4_samp_0 * i4_coeff_2;
            i4_rslt_1 += i4_samp_1 * i4_coeff_3;
            i4_rslt_2 = i4_samp_0 * i4_coeff_0;
            i4_rslt_2 += i4_samp_1 * i4_coeff_1;

            /* Rounding to 8-bit values */
            i4_rslt_1 += 32;
            i4_rslt_1 >>= 6;
            i4_rslt_2 += 32;
            i4_rslt_2 >>= 6;

            /* Storing the results */
            pu1_out[2 * i4_x] = i4_rslt_1;
            pu1_out[2 * (i4_x + 1)] = i4_rslt_2;

        } /* End of loop over y */

        /* y = 7, y_phase = phase_1 */
        i4_samp_0 = i4_samp_1;
        i4_samp_1 = *pi2_tmp++;

        /* since x_phase = phase_1 for x = 7 */
        i4_rslt_1 = i4_samp_0 * i4_coeff_2;
        i4_rslt_1 += i4_samp_1 * i4_coeff_3;

        /* Round to 8-bit value */
        i4_rslt_1 += 32;
        i4_rslt_1 >>= 6;

        /* Store the output */
        pu1_out[2 * 7] = i4_rslt_1;

        /* Incrementing the output ptr */
        pu1_out += i4_dst_stride;

    } /* End of loop over x */
} /* isvcd_horz_interpol_chroma_dyadic_1 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_horz_interpol_chroma_dyadic_2                       */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                  horizontal intra resampling for dyadic scaling ratios for*/
/*                    chroma with following ref_lyr_chroma_phase_x_plus1_flag*/
/*                    and chroma_phase_x_plus1_flag:                         */
/*                        ref_lyr        cur_lyr                             */
/*                            0            1                                 */
/*  Inputs        : pi2_tmp_filt_buf : ptr to the 6x8 buffer containing the  */
/*                        vertically interpolated data                       */
/*                    pu1_out_buf : pointer to the output buffer             */
/*                    i4_out_stride : output buffer stride                   */
/*                    i4_phase_0 : x phase for even values of x              */
/*                    i4_phase_1 : x phase for odd values of x               */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction          */
/*  Outputs       : resampled samples                                        */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_horz_interpol_chroma_dyadic_2(WORD16 *pi2_tmp_filt_buf, UWORD8 *pu1_out_buf,
                                         WORD32 i4_out_stride, WORD32 i4_phase_0, WORD32 i4_phase_1)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_coeff_0, i4_coeff_1, i4_coeff_2, i4_coeff_3;
    WORD32 i4_samp_0, i4_samp_1;
    WORD32 i4_rslt_1, i4_rslt_2;
    WORD32 i4_dst_stride;
    UWORD8 *pu1_out;
    WORD16 *pi2_tmp;

    i4_coeff_0 = 8 - i4_phase_0;
    i4_coeff_1 = i4_phase_0;
    i4_coeff_2 = 8 - i4_phase_1;
    i4_coeff_3 = i4_phase_1;

    pu1_out = pu1_out_buf;
    pi2_tmp = pi2_tmp_filt_buf + 1;
    i4_dst_stride = i4_out_stride;

    /* Horizontal interpolation */
    for(i4_y = 0; i4_y < 8; i4_y++)
    {
        /* x = 0, x_phase = phase_0 */
        i4_samp_1 = *pi2_tmp++;

        for(i4_x = 0; i4_x < 8; i4_x += 2)
        {
            i4_samp_0 = i4_samp_1;
            i4_samp_1 = *pi2_tmp++;

            /* x_phase is phase_1 for odd values of x */
            /* and phase_0 for even values of x          */
            i4_rslt_1 = i4_samp_0 * i4_coeff_0;
            i4_rslt_1 += i4_samp_1 * i4_coeff_1;

            i4_rslt_2 = i4_samp_0 * i4_coeff_2;
            i4_rslt_2 += i4_samp_1 * i4_coeff_3;

            /* Rounding to 8-bit values */
            i4_rslt_1 += 32;
            i4_rslt_1 >>= 6;
            i4_rslt_2 += 32;
            i4_rslt_2 >>= 6;

            /* Storing the results */
            pu1_out[2 * i4_x] = i4_rslt_1;
            pu1_out[2 * (i4_x + 1)] = i4_rslt_2;

        } /* End of loop over x */

        /* Incrementing the ptrs */
        pi2_tmp += 1;
        pu1_out += i4_dst_stride;

    } /* End of loop over y */
} /* isvcd_horz_interpol_chroma_dyadic_2 */

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_intra_resamp_mb_dyadic                               */
/*                                                                           */
/*  Description   : MB level function which performs the intra resampling    */
/*                  of data of an MB (luma and chroma inclusive) for dyadic  */
/*                    scaling ratios                                         */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_mb_mode_map : ref layer mb mode map buff desc     */
/*                  ps_curr_luma : current layer out luma buffer desc        */
/*                  ps_curr_chroma : current layer out chroma buffer desc    */
/*                  x,y : current mb coorinate                               */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interpolation function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           creation                             */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_intra_resamp_mb_dyadic(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                                    mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode_map,
                                    mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                                    mb_coord_t *ps_mb_coord, void *pv_svc_dec)
{
    /* --------------------------------------------------------------------- */
    /* I/O buffer params                                                     */
    /* --------------------------------------------------------------------- */
    UWORD8 *pu1_inp_luma, *pu1_inp_chroma;
    UWORD8 *pu1_out_luma, *pu1_out_chroma;
    UWORD8 *pu1_out_cb, *pu1_out_cr;
    UWORD8 *pu1_refarray_luma, *pu1_refarray_cb, *pu1_refarray_cr;
    WORD16 *pi2_tmp_filt_buf;
    WORD32 i4_inp_luma_stride, i4_inp_chroma_stride;
    WORD32 i4_out_luma_stride, i4_out_chroma_stride;
    UWORD16 u2_mb_x_ref, u2_mb_y_ref;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    /* --------------------------------------------------------------------- */
    /* Intra resampling ctxt pointers                                         */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;

    /* --------------------------------------------------------------------- */
    /* reference and current layer MB coordinates                             */
    /* --------------------------------------------------------------------- */
    WORD32 i4_scaled_mb_x, i4_scaled_mb_y;
    WORD32 i4_top, i4_left;
    WORD32 ret;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    /* --------------------------------------------------------------------- */
    /* Pointer derivation                                                     */
    /* --------------------------------------------------------------------- */
    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    /* --------------------------------------------------------------------- */
    /* MB coordinate derivation                                                 */
    /* --------------------------------------------------------------------- */
    i4_scaled_mb_x =
        ps_mb_coord->u2_mb_x - (ps_svc_slice_params->i4_scaled_ref_layer_left_offset >> 4);
    i4_scaled_mb_y =
        ps_mb_coord->u2_mb_y - (ps_svc_slice_params->i4_scaled_ref_layer_top_offset >> 4);

    if(i4_scaled_mb_x & 0x1)
    {
        i4_left = 1;
    }
    else
    {
        i4_left = -1;
    }
    if(i4_scaled_mb_y & 0x1)
    {
        i4_top = 1;
    }
    else
    {
        i4_top = -1;
    }

    u2_mb_x_ref = (i4_scaled_mb_x >> 1);
    u2_mb_y_ref = (i4_scaled_mb_y >> 1);

    /* --------------------------------------------------------------------- */
    /* Reference Array Consrtuction - luma and chroma                         */
    /* --------------------------------------------------------------------- */
    pu1_inp_luma = (UWORD8 *) ps_ref_luma->pv_buffer;
    pu1_inp_chroma = (UWORD8 *) ps_ref_chroma->pv_buffer;
    i4_inp_luma_stride = ps_ref_luma->i4_num_element_stride;
    i4_inp_chroma_stride = ps_ref_chroma->i4_num_element_stride;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_reflayer_construction_dyadic(
        pv_intra_samp_ctxt, ps_ref_mb_mode_map, pu1_inp_luma, pu1_inp_chroma, i4_inp_luma_stride,
        i4_inp_chroma_stride, i4_top, i4_left, u2_mb_x_ref, u2_mb_y_ref);

    if(ret != OK)
    {
        return NOT_OK;
    }
    /* --------------------------------------------------------------------- */
    /* LUMA INTERPOLATION                                                       */
    /* --------------------------------------------------------------------- */
    pu1_refarray_luma = ps_ctxt->pu1_refarray_buffer;
    if(1 == i4_top)
    {
        pu1_refarray_luma += (DYADIC_REF_W_Y << 3);
    }
    if(1 == i4_left)
    {
        pu1_refarray_luma += 8;
    }
    pu1_out_luma = (UWORD8 *) ps_curr_luma->pv_buffer;
    i4_out_luma_stride = ps_curr_luma->i4_num_element_stride;
    pi2_tmp_filt_buf = (WORD16 *) ps_ctxt->pi4_temp_interpolation_buffer;

    ps_ctxt->pf_interpolate_base_luma_dyadic(pu1_refarray_luma, pi2_tmp_filt_buf, pu1_out_luma,
                                             i4_out_luma_stride);

    /* --------------------------------------------------------------------- */
    /* CHROMA INTERPOLATION                                                     */
    /* --------------------------------------------------------------------- */
    pu1_out_chroma = (UWORD8 *) ps_curr_chroma->pv_buffer;
    i4_out_chroma_stride = ps_curr_chroma->i4_num_element_stride;

    /* CB */
    pu1_out_cb = pu1_out_chroma;
    pu1_refarray_cb = ps_ctxt->pu1_refarray_cb;

    if(1 == i4_top)
    {
        pu1_refarray_cb += (DYADIC_REF_W_C << 2);
    }
    if(1 == i4_left)
    {
        pu1_refarray_cb += 4;
    }

    /* Vertical interpolation */
    ps_lyr_ctxt->pf_vert_chroma_interpol(pu1_refarray_cb, pi2_tmp_filt_buf,
                                         ps_lyr_ctxt->i4_y_phase_0, ps_lyr_ctxt->i4_y_phase_1);

    /* Horizontal interpolation */
    ps_lyr_ctxt->pf_horz_chroma_interpol(pi2_tmp_filt_buf, pu1_out_cb, i4_out_chroma_stride,
                                         ps_lyr_ctxt->i4_x_phase_0, ps_lyr_ctxt->i4_x_phase_1);

    /* CR */
    pu1_out_cr = pu1_out_chroma + 1;
    pu1_refarray_cr = ps_ctxt->pu1_refarray_cr;

    if(1 == i4_top)
    {
        pu1_refarray_cr += (DYADIC_REF_W_C << 2);
    }
    if(1 == i4_left)
    {
        pu1_refarray_cr += 4;
    }

    /* Vertical interpolation */
    ps_lyr_ctxt->pf_vert_chroma_interpol(pu1_refarray_cr, pi2_tmp_filt_buf,
                                         ps_lyr_ctxt->i4_y_phase_0, ps_lyr_ctxt->i4_y_phase_1);

    /* Horizontal interpolation */
    ps_lyr_ctxt->pf_horz_chroma_interpol(pi2_tmp_filt_buf, pu1_out_cr, i4_out_chroma_stride,
                                         ps_lyr_ctxt->i4_x_phase_0, ps_lyr_ctxt->i4_x_phase_1);
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_intra_base                              */
/*                                                                           */
/*  Description   : This function takes the reference array buffer & performs*/
/*                    interpolation of a component to find the intra         */
/*                     resampled value                                       */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  pu1_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_intra_base(void *pv_intra_samp_ctxt, UWORD8 *pu1_out, WORD32 i4_out_stride,
                                  WORD32 i4_refarray_wd, WORD32 i4_mb_x, WORD32 i4_mb_y,
                                  WORD32 i4_chroma_flag, WORD32 i4_refarray_flag)
{
    /* --------------------------------------------------------------------- */
    /* Index Parameters                                                         */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_map_ctxt_t *ps_map_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;

    WORD32 i4_temp_array_ht;
    WORD32 *pi4_interp_buff;
    WORD32 *pi4_interp_buff_temp;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD32 i4_x_min, i4_x_max;
    ref_min_max_map_t *ps_x_min_max;
    UWORD8 *pu1_refarray = NULL;

    /* --------------------------------------------------------------------- */
    /* Extracting pointers from the  context                                  */
    /* --------------------------------------------------------------------- */
    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    if(0 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_buffer;
    }
    else if(1 == i4_refarray_flag)
    {
        pu1_refarray = ps_ctxt->pu1_refarray_cb;
    }

    /* --------------------------------------------------------------------- */
    /* LUMA    or CHROMA */
    /* --------------------------------------------------------------------- */
    if(1 == i4_chroma_flag)
        ps_map_ctxt = &(ps_lyr_ctxt->s_chroma_map_ctxt);
    else
        ps_map_ctxt = &(ps_lyr_ctxt->s_luma_map_ctxt);

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;
    ps_x_min_max = ps_map_ctxt->ps_x_min_max;
    i4_frm_mb_y = i4_mb_y * i4_mb_ht;
    i4_frm_mb_x = i4_mb_x * i4_mb_wd;

    /* get the min and max positions */
    i4_x_min = ps_x_min_max[i4_mb_x].i2_min_pos;
    i4_x_max = ps_x_min_max[i4_mb_x].i2_max_pos;

    /* --------------------------------------------------------------------- */
    /* Projected frame level pointers                                        */
    /* --------------------------------------------------------------------- */
    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

    /* --------------------------------------------------------------------- */
    /* Pointers and Dimenstion of the temporary buffer                         */
    /* --------------------------------------------------------------------- */
    i4_temp_array_ht = i4_mb_ht;
    pi4_interp_buff = ps_ctxt->pi4_temp_interpolation_buffer;
    pi4_interp_buff_temp = pi4_interp_buff;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation in vertical direction                             */
    /* --------------------------------------------------------------------- */
    for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
    {
        for(i4_x = (i4_x_min - 1); i4_x <= (i4_x_max + 2); i4_x++)
        {
            UWORD8 *pu1_refarray_temp;
            WORD32 i4_y_ref;
            WORD32 i4_y_phase;
            /* ------------------------------------------------------------ */
            /* Finding the offset                                            */
            /* ------------------------------------------------------------ */
            i4_y_ref = ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos;
            i4_y_phase = ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            pu1_refarray_temp = pu1_refarray + i4_x + (i4_y_ref * i4_refarray_wd);

            /* ------------------------------------------------------------ */
            /* Check for Luma/Chroma Processing                                */
            /* ------------------------------------------------------------ */
            if(0 == i4_chroma_flag)
            {
                *(pi4_interp_buff + i4_x) =
                    (g_ai1_interp_filter_luma[i4_y_phase]) *
                        (*(pu1_refarray_temp - i4_refarray_wd)) +

                    (g_ai1_interp_filter_luma[16 + i4_y_phase]) * (*(pu1_refarray_temp)) +

                    (g_ai1_interp_filter_luma[32 + i4_y_phase]) *
                        (*(pu1_refarray_temp + i4_refarray_wd)) +

                    (g_ai1_interp_filter_luma[48 + i4_y_phase]) *
                        (*(pu1_refarray_temp + (2 * i4_refarray_wd)));
            }
            else
            {
                *(pi4_interp_buff + i4_x) =
                    (g_au1_interp_filter_chroma[i4_y_phase]) * (*(pu1_refarray_temp)) +

                    (g_au1_interp_filter_chroma[16 + i4_y_phase]) *
                        (*(pu1_refarray_temp + i4_refarray_wd));
            }

        } /* end of loop over array  width */
        pi4_interp_buff = pi4_interp_buff + i4_refarray_wd;
    }     /* end of loop over temp array height*/

    pi4_interp_buff = pi4_interp_buff_temp;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation in horizontal direction                         */
    /* --------------------------------------------------------------------- */
    for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            WORD32 i4_x_ref;
            WORD32 i4_x_phase;
            /* ------------------------------------------------------------- */
            /* Finding the offset                                             */
            /* ------------------------------------------------------------- */
            i4_x_ref = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            i4_x_phase = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;
            pi4_interp_buff_temp = pi4_interp_buff + i4_x_ref;

            /* ------------------------------------------------------------- */
            /* Check for Luma/Chroma Processing                                 */
            /* ------------------------------------------------------------- */
            if(0 == i4_chroma_flag)
            {
                *(pu1_out + i4_x + (i4_y * i4_out_stride)) =

                    CLIPUCHAR(
                        ((g_ai1_interp_filter_luma[i4_x_phase]) * (*(pi4_interp_buff_temp - 1)) +
                         (g_ai1_interp_filter_luma[16 + i4_x_phase]) * (*(pi4_interp_buff_temp)) +
                         (g_ai1_interp_filter_luma[32 + i4_x_phase]) *
                             (*(pi4_interp_buff_temp + 1)) +
                         (g_ai1_interp_filter_luma[48 + i4_x_phase]) *
                             (*(pi4_interp_buff_temp + 2)) +
                         512) >>
                        10);
            }
            else
            {
                *(pu1_out + (2 * i4_x) + (i4_y * i4_out_stride)) = CLIPUCHAR(
                    ((g_au1_interp_filter_chroma[i4_x_phase]) * (*(pi4_interp_buff_temp)) +
                     (g_au1_interp_filter_chroma[16 + i4_x_phase]) * (*(pi4_interp_buff_temp + 1)) +
                     512) >>
                    10);
            }

        } /* end of loop over array width */
        pi4_interp_buff = pi4_interp_buff + i4_refarray_wd;
    }     /* end of loop over MB height */

    return;
} /* End of Interpolation Function */
/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_intra_resamp_mb                                      */
/*                                                                           */
/*  Description   : MB level function whcih perform the intra resampling     */
/*                  of data of an MB (luma and chroma insclusive)            */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt : intra sampling context              */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_mb_mode_map : ref layer mb mode map buff desc     */
/*                  ps_curr_luma : current layer out luma buffer desc        */
/*                  ps_curr_chroma : current layer out chroma buffer desc    */
/*                  ps_mb_coord : current mb coorinate                       */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interpolation function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_intra_resamp_mb(void *pv_intra_samp_ctxt, mem_element_t *ps_ref_luma,
                             mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode_map,
                             mem_element_t *ps_curr_luma, mem_element_t *ps_curr_chroma,
                             mb_coord_t *ps_mb_coord)
{
    /* --------------------------------------------------------------------- */
    /* I/O buffer params                                                     */
    /* --------------------------------------------------------------------- */
    intra_sampling_ctxt_t *ps_ctxt;
    UWORD8 *pu1_inp_luma, *pu1_inp_cb, *pu1_inp_cr;
    UWORD8 *pu1_out_luma, *pu1_out_cb, *pu1_out_cr;
    WORD32 i4_inp_stride;
    WORD32 i4_out_stride;
    WORD32 i4_refarray_stride;

    /* --------------------------------------------------------------------- */
    /* ref sample array and corresponding parametrs                             */
    /* --------------------------------------------------------------------- */
    WORD32 i4_chroma_flag, ret;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;
    /* needs to be populated at proper place */
    i4_refarray_stride = ps_ctxt->i4_refarray_stride;

    /* --------------------------------------------------------------------- */
    /* LUMA PROCESSING                                                        */
    /* --------------------------------------------------------------------- */
    pu1_inp_luma = (UWORD8 *) ps_ref_luma->pv_buffer;
    pu1_out_luma = (UWORD8 *) ps_curr_luma->pv_buffer;
    i4_inp_stride = ps_ref_luma->i4_num_element_stride;
    i4_out_stride = ps_curr_luma->i4_num_element_stride;
    i4_chroma_flag = 0;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_reflayer_construction(pv_intra_samp_ctxt, pu1_inp_luma, NULL, i4_inp_stride,
                                      i4_refarray_stride, ps_ref_mb_mode_map, ps_mb_coord,
                                      i4_chroma_flag);

    if(ret != OK)
    {
        return NOT_OK;
    }
    /* ---- Interpolation process for Intra_Base prediction     ------ */
    ps_ctxt->pf_interpolate_intra_base(pv_intra_samp_ctxt, pu1_out_luma, i4_out_stride,
                                       i4_refarray_stride, ps_mb_coord->u2_mb_x,
                                       ps_mb_coord->u2_mb_y, i4_chroma_flag, 0);

    /* --------------------------------------------------------------------- */
    /* CHROMA PROCESSING                                                       */
    /* --------------------------------------------------------------------- */
    /* CB */
    i4_inp_stride = ps_ref_chroma->i4_num_element_stride;
    pu1_inp_cb = (UWORD8 *) ps_ref_chroma->pv_buffer;
    pu1_inp_cr = pu1_inp_cb + 1;
    i4_chroma_flag = 1;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_reflayer_construction(pv_intra_samp_ctxt, pu1_inp_cb, pu1_inp_cr, i4_inp_stride,
                                      i4_refarray_stride, ps_ref_mb_mode_map, ps_mb_coord,
                                      i4_chroma_flag);

    if(ret != OK)
    {
        return NOT_OK;
    }
    i4_out_stride = ps_curr_chroma->i4_num_element_stride;
    pu1_out_cb = (UWORD8 *) ps_curr_chroma->pv_buffer;
    pu1_out_cr = pu1_out_cb + 1;

    /* ---- Cb Interpolation process for Intra_Base prediction     ------ */
    ps_ctxt->pf_interpolate_intra_base(pv_intra_samp_ctxt, pu1_out_cb, i4_out_stride,
                                       i4_refarray_stride, ps_mb_coord->u2_mb_x,
                                       ps_mb_coord->u2_mb_y, i4_chroma_flag, 0);

    /* ---- Cr Interpolation process for Intra_Base prediction     ------ */
    ps_ctxt->pf_interpolate_intra_base(pv_intra_samp_ctxt, pu1_out_cr, i4_out_stride,
                                       i4_refarray_stride, ps_mb_coord->u2_mb_x,
                                       ps_mb_coord->u2_mb_y, i4_chroma_flag, 1);
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_intra_resamp_generate_segment_lookup                */
/*                                                                           */
/*  Description   : This function generates segment lookup used to derive    */
/*                    segments which have to be be intra resampled           */
/*                                                                           */
/*  Inputs        : pv_lookup_table : look up table                          */
/*                  i4_dimension    : dimension of the block which is used in*/
/*                                      resampling process.                  */
/*                  i4_mb_size        : size of the mb                       */
/*  Globals       : None                                                     */
/*  Processing    : This function generates segment lookup used to derive    */
/*                    segments which have to be be intra resampled           */
/*  Outputs          : none                                                  */
/*  Returns          : none                                                  */
/*                                                                           */
/*  Issues       : None                                                      */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Creation                             */
/*                                                                           */
/*****************************************************************************/
void isvcd_intra_resamp_generate_segment_lookup(seg_lookup_desc_t *ps_seg_lookup_table,
                                                WORD32 i4_dimension, WORD32 i4_mb_size,
                                                WORD32 i4_shift_val)
{
    WORD32 i4_x;
    WORD32 i4_position, i4_dist_prev_mb, i4_dist_next_mb;
    UWORD8 u1_seg_dim;
    UWORD8 u1_num_sgmts;
    WORD32 i4_block_size = i4_mb_size >> 1;
    UWORD8 u1_offset = 0;
    seg_lookup_desc_t *ps_segments;
    seg_description_t *ps_seg_desc;

    memset(ps_seg_lookup_table, 0, i4_mb_size * sizeof(seg_lookup_desc_t));

    for(i4_x = 0; i4_x < i4_mb_size; i4_x++)
    {
        ps_segments = &ps_seg_lookup_table[i4_x];
        ps_seg_desc = ps_segments->s_segments;
        i4_position = i4_x;

        if(i4_x >= i4_block_size)
        {
            /* set the fourth bit so that later it can be directly OR ed */
            ps_segments->u4_start_pos = 8;
        }
        else
        {
            ps_segments->u4_start_pos = 0;
        }

        u1_num_sgmts = 0;
        u1_offset = 0;

        while(i4_position < (i4_x + i4_dimension))
        {
            /* check and fill the nearest mb boundry flag */
            if((i4_position & (i4_mb_size - 1)) < i4_block_size)
            {
                ps_seg_desc->i1_nearst_mb_bdry = -1;
            }
            else
            {
                ps_seg_desc->i1_nearst_mb_bdry = 1;
            }

            /* find the distance from the previous MB for start of segment*/
            i4_dist_prev_mb = (i4_position & (i4_mb_size - 1));

            ps_seg_desc->i1_dist_idx =
                ((i4_dist_prev_mb >= i4_mb_size >> 1) ? (i4_mb_size - i4_dist_prev_mb)
                                                      : -(i4_dist_prev_mb + 1));

            /* find the size of the segment */
            u1_seg_dim = (i4_block_size - (i4_position & (i4_block_size - 1)));
            i4_position += u1_seg_dim;
            if(i4_position > (i4_x + i4_dimension))
            {
                i4_position = (i4_x + i4_dimension);
                u1_seg_dim = (i4_position & (i4_block_size - 1));
            }

            /* find the distance from the next MB for end of segment */
            i4_dist_next_mb = (i4_position & (i4_mb_size - 1));
            ps_seg_desc->u1_seg_dim = u1_seg_dim;
            ps_seg_desc->u1_seg_off = u1_offset;

            /* check if the segment has a adjoining MB edge */
            if(i4_dist_prev_mb == 0)
            {
                if(0 == u1_num_sgmts)
                {
                    ps_seg_desc->u1_mb_adjoin = 0;
                }
                else
                {
                    ps_seg_desc->u1_mb_adjoin = 1 << i4_shift_val;
                }
            }
            else if(i4_dist_next_mb == 0)
            {
                if(i4_position == (i4_x + i4_dimension))
                {
                    ps_seg_desc->u1_mb_adjoin = 0;
                }
                else
                {
                    ps_seg_desc->u1_mb_adjoin = 1 << i4_shift_val;
                }
            }
            else
            {
                ps_seg_desc->u1_mb_adjoin = 0;
            }

            /* Updations */
            u1_offset += u1_seg_dim;
            u1_num_sgmts++;
            ps_seg_desc++;
        }

        /* fill the number of segments for this position */
        ps_segments->u1_num_segments = u1_num_sgmts;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_intra_resamp_populate_list                            */
/*                                                                           */
/*  Description   : This is a seq or frame level init function which fills   */
/*                  all offsets, projected locations arrays based on         */
/*                  the two resolutions  and cropping parameters             */
/*  Inputs        : refer ot doxygen comments below                          */
/*  Globals       : none                                                     */
/*  Processing    : it projects the locations and computes the values        */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2009   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_intra_resamp_populate_list(intra_samp_map_ctxt_t *ps_map_ctxt,
                                      res_prms_t *ps_curr_res_prms, res_prms_t *ps_ref_res_prms,
                                      WORD32 i4_chroma_flag, svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    /* --------------------------------------------------------------------- */
    /* Local variables required for finding the mapping between the layers     */
    /* --------------------------------------------------------------------- */
    UWORD32 i4_shift_x, i4_shift_y, i4_scale_x, i4_scale_y;
    WORD32 i4_offset_x, i4_offset_y;
    WORD32 i4_add_x, i4_add_y, i4_delta_x, i4_delta_y, i4_refphase_x, i4_refphase_y;
    WORD32 i4_phase_x, i4_phase_y, i4_sub_wd, i4_sub_ht, i4_mb_wd, i4_mb_ht;
    WORD32 i4_horz_dim, i4_vert_dim, i4_tmp;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    /* --------------------------------------------------------------------- */
    /* Local Pointer Declaration for arrays in Mapping context                 */
    /* --------------------------------------------------------------------- */
    ref_mb_map_t *ps_x_off_len, *ps_y_off_len;
    UWORD32 i4_ref_wd, i4_ref_ht, i4_scaled_wd, i4_scaled_ht, i4_curr_lyr_width, i4_curr_lyr_height;

    /* --------------------------------------------------------------------- */
    /* Local Flag Declaration                                                 */
    /* --------------------------------------------------------------------- */
    WORD32 i4_ref_layer_field_pic_flag, i4_field_pic_flag, i4_frame_mbs_only_flag;
    WORD32 i4_ref_layer_frame_Mbs_only_flag, i4_field_Mb_flag, i4_bot_field_flag;

    /* --------------------------------------------------------------------- */
    /* Cropping Parameters Declaration                                         */
    /* --------------------------------------------------------------------- */
    WORD32 i4_scaled_ref_layer_left_offset, i4_scaled_ref_layer_top_offset;
    WORD32 i4_scaled_ref_layer_right_offset, i4_scaled_ref_layer_bottom_offset;
    dec_seq_params_t *ps_sps;
    dec_svc_seq_params_t *ps_subset_sps;
    ps_sps = ps_dec->ps_cur_sps;
    ps_subset_sps = ps_svc_lyr_dec->ps_cur_subset_sps;

    /* --------------------------------------------------------------------- */
    /* Hardcoding flag information    (assuming no field support) */
    /* --------------------------------------------------------------------- */
    i4_ref_layer_field_pic_flag = SVCD_FALSE;
    i4_field_pic_flag = SVCD_FALSE;
    i4_frame_mbs_only_flag = SVCD_TRUE;
    i4_field_Mb_flag = SVCD_FALSE;
    i4_bot_field_flag = SVCD_FALSE;
    i4_ref_layer_frame_Mbs_only_flag = SVCD_TRUE;
    i4_horz_dim = 0;
    i4_vert_dim = 0;

    /* --------------------------------------------------------------------- */
    /* Pointer and Paramater are intialized    - Chroma and Luma */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_base_width;
        WORD32 i4_base_height;
        WORD32 i4_ref_layer_chroma_phase_x_plus1_flag;
        WORD32 i4_ref_layer_chroma_phase_y_plus1;
        WORD32 i4_chroma_phase_x_plus1_flag;
        WORD32 i4_chroma_phase_y_plus1;

        /* ------------------------------------------------------------- */
        /* HARD CODED FOR 420                                             */
        /* ------------------------------------------------------------- */
        WORD32 i4_sub_wd_chroma = 2;
        WORD32 i4_sub_ht_chroma = 2;

        i4_base_width = ps_ref_res_prms->i4_res_width;
        i4_base_height = ps_ref_res_prms->i4_res_height;

        i4_ref_layer_chroma_phase_x_plus1_flag =
            ps_curr_res_prms->i1_ref_lyr_chroma_phase_x_plus1_flag;

        i4_ref_layer_chroma_phase_y_plus1 = ps_curr_res_prms->i1_ref_lyr_chroma_phase_y_plus1;
        i4_chroma_phase_x_plus1_flag = ps_subset_sps->s_sps_svc_ext.u1_chroma_phase_x_plus1_flag;
        i4_chroma_phase_y_plus1 = ps_subset_sps->s_sps_svc_ext.u1_chroma_phase_y_plus1;
        i4_scaled_ref_layer_bottom_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_bot;
        i4_scaled_ref_layer_left_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_left;
        i4_scaled_ref_layer_top_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_top;
        i4_scaled_ref_layer_right_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_rt;

        /* ----------------------------------------------------------------- */
        /* Computing Effective Frame Dimensions                                 */
        /* ------------------------------------------------------------------*/
        i4_ref_wd = (i4_base_width >> i4_chroma_flag);
        i4_ref_ht = (i4_base_height >> i4_chroma_flag) * (1 + i4_ref_layer_field_pic_flag);

        i4_scaled_wd = ps_curr_res_prms->u2_scaled_ref_width;
        i4_scaled_ht = ps_curr_res_prms->u2_scaled_ref_height;
        i4_scaled_wd = (i4_scaled_wd >> i4_chroma_flag);
        i4_scaled_ht = (i4_scaled_ht >> i4_chroma_flag) * (1 + i4_field_pic_flag);

        if(1 == i4_chroma_flag)
        {
            i4_refphase_x = i4_ref_layer_chroma_phase_x_plus1_flag - 1;
            i4_refphase_y = i4_ref_layer_chroma_phase_y_plus1 - 1;
            i4_phase_x = i4_chroma_phase_x_plus1_flag - 1;
            i4_phase_y = i4_chroma_phase_y_plus1 - 1;
            i4_sub_wd = i4_sub_wd_chroma;
            i4_sub_ht = i4_sub_ht_chroma;
            i4_mb_wd = MB_WIDTH >> 1;
            i4_mb_ht = MB_HEIGHT >> 1;
        }
        else
        {
            i4_refphase_x = 0;
            i4_refphase_y = 0;
            i4_phase_x = 0;
            i4_phase_y = 0;
            i4_sub_wd = 1;
            i4_sub_ht = 1;
            i4_mb_wd = MB_WIDTH;
            i4_mb_ht = MB_HEIGHT;
        }
    }

    /* --------------------------------------------------------------------- */
    /* Derive shift x and y based on level idd                               */
    /* --------------------------------------------------------------------- */
    if(ps_sps->u1_level_idc <= 30)
    {
        i4_shift_x = 16;
        i4_shift_y = 16;
    }
    else
    {
        i4_shift_x = 31 - isvcd_get_ceil_log2(i4_ref_wd);
        i4_shift_y = 31 - isvcd_get_ceil_log2(i4_ref_ht);
    }

    /* --------------------------------------------------------------------- */
    /* The following condition is not true in our case for time being         */
    /* --------------------------------------------------------------------- */
    if((SVCD_FALSE == i4_frame_mbs_only_flag) || (SVCD_FALSE == i4_ref_layer_frame_Mbs_only_flag))
    {
        i4_phase_y = i4_phase_y + 4 * i4_bot_field_flag;

        if(1 == i4_ref_layer_frame_Mbs_only_flag)
            i4_refphase_y = (2 * i4_refphase_y) + 2;
        else
            i4_refphase_y = i4_refphase_y + (4 * i4_bot_field_flag);
    }

    /* --------------------------------------------------------------------- */
    /* Dx and Dy Computation - Ratio of the base and enhance layer width     */
    /* --------------------------------------------------------------------- */
    i4_scale_x = ((i4_ref_wd << i4_shift_x) + (i4_scaled_wd >> 1)) / (i4_scaled_wd);
    i4_scale_y = ((i4_ref_ht << i4_shift_y) + (i4_scaled_ht >> 1)) / (i4_scaled_ht);

    i4_offset_x = i4_scaled_ref_layer_left_offset / i4_sub_wd;
    i4_add_x = (((i4_ref_wd * (2 + i4_phase_x)) << (i4_shift_x - 2)) + (i4_scaled_wd >> 1)) /
                   i4_scaled_wd +
               (1 << (i4_shift_x - 5));
    i4_delta_x = 4 * (2 + i4_refphase_x);

    if((SVCD_TRUE == i4_frame_mbs_only_flag) && (SVCD_TRUE == i4_ref_layer_frame_Mbs_only_flag))
    {
        i4_offset_y = i4_scaled_ref_layer_top_offset / i4_sub_ht;
        i4_add_y = (((i4_ref_ht * (2 + i4_phase_y)) << (i4_shift_y - 2)) + (i4_scaled_ht >> 1)) /
                       i4_scaled_ht +
                   (1 << (i4_shift_y - 5));
        i4_delta_y = 4 * (2 + i4_refphase_y);
    }
    else
    {
        i4_offset_y = i4_scaled_ref_layer_top_offset / (2 * i4_sub_ht);
        i4_add_y = (((i4_ref_ht * (2 + i4_phase_y)) << (i4_shift_y - 3)) + (i4_scaled_ht >> 1)) /
                       i4_scaled_ht +
                   (1 << (i4_shift_y - 5));
        i4_delta_y = 2 * (2 + i4_refphase_y);
    }

    /* --------------------------------------------------------------------- */
    /* Intializing Local Pointers    - Chroma and Luma                         */
    /* --------------------------------------------------------------------- */
    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;
    i4_curr_lyr_width = ps_curr_res_prms->i4_res_width >> i4_chroma_flag;
    i4_curr_lyr_height = ps_curr_res_prms->i4_res_height >> i4_chroma_flag;

    /* --------------------------------------------------------------------- */
    /* Dyadic Scaling Ratios Handling                                         */
    /* --------------------------------------------------------------------- */
    if(1 == ps_curr_res_prms->u1_dyadic_flag)
    {
        WORD32 i4_refarray_wd, i4_x_offset;
        WORD32 i4_refarray_ht, i4_y_offset;
        WORD32 i4_crp_wd_lt, i4_crp_ht_top;
        WORD32 i4_crp_wd_rt, i4_crp_ht_bot;
        WORD32 i4_ref_lyr_wd, i4_ref_lyr_ht;
        WORD32 i4_ref_x, i4_ref_y;
        WORD32 i4_ofst;
        WORD32 i4_i, i4_j;

        /* Hard coded for dyadic case */
        i4_refarray_wd = 20 >> i4_chroma_flag;
        i4_ofst = -2 >> i4_chroma_flag;
        i4_crp_wd_lt = i4_scaled_ref_layer_left_offset >> i4_chroma_flag;
        i4_crp_wd_rt = i4_scaled_ref_layer_right_offset >> i4_chroma_flag;
        i4_ref_lyr_wd = (i4_curr_lyr_width >> 1);

        i4_ref_x = 0;
        for(i4_i = 0; i4_i < (WORD32) i4_curr_lyr_width; i4_i += (i4_mb_wd << 1))
        {
            i4_x_offset = MAX(i4_ofst, (i4_ref_x + i4_ofst));
            i4_x_offset = MIN(i4_x_offset, (i4_ref_lyr_wd - i4_ofst));
            ps_x_off_len->i2_offset = i4_x_offset;
            ps_x_off_len->i2_length = i4_refarray_wd;
            ps_x_off_len++;
            ps_x_off_len->i2_offset = i4_x_offset;
            ps_x_off_len->i2_length = i4_refarray_wd;
            ps_x_off_len++;
            if(i4_i >= i4_crp_wd_lt)
            {
                if(i4_i <= (WORD32) (i4_curr_lyr_width - i4_crp_wd_rt))
                {
                    i4_ref_x += i4_mb_wd;
                }
            }
        }

        i4_refarray_ht = 20 >> i4_chroma_flag;
        i4_crp_ht_top = i4_scaled_ref_layer_top_offset >> i4_chroma_flag;
        i4_crp_ht_bot = i4_scaled_ref_layer_bottom_offset >> i4_chroma_flag;
        i4_ref_lyr_ht = (i4_curr_lyr_height >> 1);

        i4_ref_y = 0;
        for(i4_j = 0; i4_j < (WORD32) i4_curr_lyr_height; i4_j += (i4_mb_ht << 1))
        {
            i4_y_offset = MAX(i4_ofst, (i4_ref_y + i4_ofst));
            i4_y_offset = MIN(i4_y_offset, (i4_ref_lyr_ht - i4_ofst));
            ps_y_off_len->i2_offset = i4_y_offset;
            ps_y_off_len->i2_length = i4_refarray_ht;
            ps_y_off_len++;
            ps_y_off_len->i2_offset = i4_y_offset;
            ps_y_off_len->i2_length = i4_refarray_ht;
            ps_y_off_len++;
            if(i4_j >= i4_crp_ht_top)
            {
                if(i4_j <= (WORD32) (i4_curr_lyr_height - i4_crp_ht_bot))
                {
                    i4_ref_y += i4_mb_ht;
                }
            }
        }
        /* No need to process further, return */
        return;
    } /* If dyadic path */

    /* Proposed Algo for Optimization */
    {
        WORD32 i4_max, i4_min;
        ref_min_max_map_t *ps_x_min_max;
        ref_min_max_map_t *ps_y_min_max;
        WORD32 i4_i, i4_j;

        ps_x_min_max = ps_map_ctxt->ps_x_min_max;
        ps_y_min_max = ps_map_ctxt->ps_y_min_max;
        /* ----------------------------------------------------------------- */
        /* Computation of offsetX refArrayW Xmin and Xmax Lists                 */
        /* ----------------------------------------------------------------- */
        for(i4_i = 0; i4_i < (WORD32) i4_curr_lyr_width; i4_i = i4_i + i4_mb_wd)
        {
            WORD32 i4_refarray_wd, i4_xr_index;
            WORD32 i4_x_refmin16;
            WORD32 i4_x_refmax16;
            WORD32 i4_x_offset;

            i4_x_refmin16 = (WORD64) (((WORD64) (i4_i - i4_offset_x) * i4_scale_x + i4_add_x) >>
                                      ((WORD32) (i4_shift_x - 4))) -
                            i4_delta_x;

            i4_x_refmax16 =
                (WORD64) (((WORD64) (i4_i + i4_mb_wd - 1 - i4_offset_x) * i4_scale_x + i4_add_x) >>
                          ((WORD32) (i4_shift_x - 4))) -
                i4_delta_x;

            /* ------------------------------------------------------------- */
            /* Modified AC205                                                 */
            /* Minimum width required - So adding 2 pixels on each side         */
            /* ------------------------------------------------------------- */
            i4_refarray_wd = ((i4_x_refmax16 + 15) >> 4) - (i4_x_refmin16 >> 4) + 1 + 4;

            /* ------------------------------------------------------------- */
            /* Setting the offset 2 pixels before                             */
            /* ------------------------------------------------------------- */
            i4_x_offset = (i4_x_refmin16 >> 4) - 2;

            /* ------------------------------------------------------------- */
            /* Modifying the values based on the location                     */
            /* ------------------------------------------------------------- */
            i4_min = i4_x_offset;
            i4_xr_index = i4_min - ((i4_min / i4_mb_wd) * i4_mb_wd);

            if(i4_xr_index < (i4_mb_wd >> 1))
            {
                i4_refarray_wd = i4_refarray_wd + (i4_mb_wd >> 1);
                i4_x_offset = i4_x_offset - (i4_mb_wd >> 1);
            }

            i4_max = ((i4_x_refmax16 + 15) >> 4) + 2;
            i4_xr_index = i4_max - ((i4_max / i4_mb_wd) * i4_mb_wd);

            if(i4_xr_index >= (i4_mb_wd >> 1)) i4_refarray_wd = i4_refarray_wd + (i4_mb_wd >> 1);

            /* ------------------------------------------------------------- */
            /* Filling the arrays with offset, min, max and refArray dim     */
            /* ------------------------------------------------------------- */
            ps_x_off_len->i2_offset = i4_x_offset;
            ps_x_off_len->i2_length = i4_refarray_wd;

            ps_x_min_max->i2_min_pos = (i4_x_refmin16 >> 4) - i4_x_offset;
            ps_x_min_max->i2_max_pos = ((i4_x_refmax16 + 15) >> 4) - i4_x_offset;

            i4_tmp = (WORD32) (ps_x_min_max->i2_max_pos - ps_x_min_max->i2_min_pos) +
                     (4 >> i4_chroma_flag);
            if(i4_tmp > i4_horz_dim)
            {
                i4_horz_dim = i4_tmp;
            }

            /* increment the pointer */
            ps_x_off_len++;
            ps_x_min_max++;
        } /* end of loop over scaled width */

        /* ----------------------------------------------------------------- */
        /* Computation of offsetY refArrayH Ymin and Ymax Lists                 */
        /* ----------------------------------------------------------------- */
        for(i4_j = 0; i4_j < (WORD32) i4_curr_lyr_height; i4_j = i4_j + i4_mb_ht)
        {
            WORD32 i4_refarray_ht, i4_yr_index;
            WORD32 i4_y_refmin16;
            WORD32 i4_y_refmax16;
            WORD32 i4_y_offset;

            i4_y_refmin16 = (WORD64) (((WORD64) (i4_j - i4_offset_y) * i4_scale_y + i4_add_y) >>
                                      ((WORD32) (i4_shift_y - 4))) -
                            i4_delta_y;

            i4_y_refmax16 =
                (WORD64) (((WORD64) (i4_j + i4_mb_ht - 1 - i4_offset_y) * i4_scale_y + i4_add_y) >>
                          ((WORD32) (i4_shift_y - 4))) -
                i4_delta_y;

            /* ------------------------------------------------------------- */
            /* Modified AC205                                                 */
            /* Minimum width required - So adding 2 pixels on each side            */
            /* ------------------------------------------------------------- */
            i4_refarray_ht = ((i4_y_refmax16 + 15) >> 4) - (i4_y_refmin16 >> 4) + 1 + 4;

            /* ------------------------------------------------------------- */
            /* Setting the offset 2 pixels before                             */
            /* ------------------------------------------------------------- */
            i4_y_offset = (i4_y_refmin16 >> 4) - 2;

            /* ------------------------------------------------------------- */
            /* Modifying the values based on the location                     */
            /* ------------------------------------------------------------- */
            i4_min = i4_y_offset;
            i4_yr_index = i4_min - ((i4_min / i4_mb_ht) * i4_mb_ht);
            if(i4_yr_index < (i4_mb_ht >> 1))
            {
                i4_refarray_ht = i4_refarray_ht + (i4_mb_ht >> 1);
                i4_y_offset = i4_y_offset - (i4_mb_ht >> 1);
            }

            i4_max = ((i4_y_refmax16 + 15) >> 4) + 2;
            i4_yr_index = i4_max - ((i4_max / i4_mb_ht) * i4_mb_ht);
            if(i4_yr_index >= (i4_mb_ht >> 1)) i4_refarray_ht = i4_refarray_ht + (i4_mb_ht >> 1);

            /* ------------------------------------------------------------- */
            /* Filling the arrays with offset, min, max and refArray dim     */
            /* ------------------------------------------------------------- */
            ps_y_off_len->i2_offset = i4_y_offset;
            ps_y_off_len->i2_length = i4_refarray_ht;
            ps_y_min_max->i2_min_pos = (i4_y_refmin16 >> 4) - i4_y_offset;
            ps_y_min_max->i2_max_pos = ((i4_y_refmax16 + 15) >> 4) - i4_y_offset;

            i4_tmp = (WORD32) (ps_y_min_max->i2_max_pos - ps_y_min_max->i2_min_pos) +
                     (4 >> i4_chroma_flag);
            if(i4_tmp > i4_vert_dim)
            {
                i4_vert_dim = i4_tmp;
            }

            /* increment the pointer */
            ps_y_off_len++;
            ps_y_min_max++;
        } /* end of loop over scaled height */
    }

    /* --------------------------------------------------------------------- */
    /* Computation of Xref and Xphase List as per standard                     */
    /* --------------------------------------------------------------------- */
    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;

    {
        ref_pixel_map_t *ps_x_pos_phase;
        WORD32 i4_xc;
        WORD32 i4_offset_x_index;

        ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;

        for(i4_xc = 0; i4_xc < (WORD32) i4_curr_lyr_width; i4_xc++)
        {
            WORD32 i4_x_offset;
            WORD32 i4_x_ref16;

            i4_offset_x_index = i4_xc / i4_mb_wd;

            i4_x_offset = ps_x_off_len[i4_offset_x_index].i2_offset;

            i4_x_ref16 = (WORD64) (((WORD64) (i4_xc - i4_offset_x) * i4_scale_x + i4_add_x) >>
                                   ((WORD32) (i4_shift_x - 4))) -
                         i4_delta_x;

            /* store the values */
            ps_x_pos_phase->i2_ref_pos = (i4_x_ref16 >> 4) - i4_x_offset;
            ps_x_pos_phase->i2_phase = i4_x_ref16 & 15;

            /* increment the pointer */
            ps_x_pos_phase++;

        } /* end of loop over scaled width */
    }

    /* --------------------------------------------------------------------- */
    /* Computation of Yref and Yphase List as per standard                     */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_yc;
        ref_pixel_map_t *ps_y_pos_phase;

        ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

        for(i4_yc = 0; i4_yc < (WORD32) i4_curr_lyr_height; i4_yc++)
        {
            WORD32 i4_y_offset;
            WORD32 i4_y_ref16;
            WORD32 i4_offset_y_index;

            i4_offset_y_index = i4_yc / i4_mb_ht;

            i4_y_offset = ps_y_off_len[i4_offset_y_index].i2_offset;

            if((SVCD_FALSE == i4_frame_mbs_only_flag) ||
               (SVCD_FALSE == i4_ref_layer_frame_Mbs_only_flag))
            {
                i4_yc = i4_yc >> (1 - i4_field_Mb_flag);
            }

            i4_y_ref16 = (WORD64) (((WORD64) (i4_yc - i4_offset_y) * i4_scale_y + i4_add_y) >>
                                   ((WORD32) (i4_shift_y - 4))) -
                         i4_delta_y;
            ps_y_pos_phase->i2_ref_pos = (i4_y_ref16 >> 4) - i4_y_offset;
            ps_y_pos_phase->i2_phase = i4_y_ref16 & 15;

            /* increment the pointer */
            ps_y_pos_phase++;

        } /* end of loop over scaled height */
    }

    /* --------------------------------------------------------------------- */
    /* Computation of Corresponding Diagonal Location                         */
    /* --------------------------------------------------------------------- */
    {
        WORD16 *pi2_xd_index;
        WORD16 *pi2_yd_index;
        WORD16 *pi2_ya_index;
        WORD32 i4_i, i4_j;

        pi2_xd_index = ps_map_ctxt->pi2_xd_index;
        pi2_yd_index = ps_map_ctxt->pi2_yd_index;
        pi2_ya_index = ps_map_ctxt->pi2_ya_index;

        for(i4_i = 0; i4_i < i4_mb_wd; i4_i++)
        {
            *(pi2_xd_index + i4_i) = ((i4_i >= i4_mb_wd >> 1) ? (i4_i - i4_mb_wd) : (i4_i + 1));

        } /* end of loop over MB width */

        for(i4_j = 0; i4_j < i4_mb_ht; i4_j++)
        {
            *(pi2_yd_index + i4_j) = ((i4_j >= i4_mb_ht >> 1) ? (i4_j - i4_mb_ht) : (i4_j + 1));

            *(pi2_ya_index + i4_j) =
                *(pi2_yd_index + i4_j) - (((i4_mb_ht >> 1) + 1) * (SIGN(*(pi2_yd_index + i4_j))));

        } /* end of loop over MB height */
    }

    /* generate the lookup to generate horizontal segments */
    isvcd_intra_resamp_generate_segment_lookup(ps_map_ctxt->ps_seg_lookup_horz, i4_horz_dim,
                                               i4_mb_wd, 3);

    /* generate the lookup to generate vertical segments */
    isvcd_intra_resamp_generate_segment_lookup(ps_map_ctxt->ps_seg_lookup_vert, i4_vert_dim,
                                               i4_mb_ht, 4);

    return;
} /* end of function "isvcd_intra_resamp_populate_list"*/

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_populate_res_prms                                   */
/*                                                                           */
/*  Description   :this function populates the current layer params          */
/*                 from the base layer and decoder context                   */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_populate_res_prms(void *pv_svc_dec)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    res_prms_t *ps_curr_lyr_res_prms;
    svc_dec_lyr_struct_t *ps_svc_ref_lyr_dec;
    ps_svc_ref_lyr_dec = ps_svc_lyr_dec->ps_dec_svc_ref_layer;
    ps_curr_lyr_res_prms = &ps_svc_lyr_dec->s_res_prms;

    ps_curr_lyr_res_prms->i4_res_width = ps_dec->u2_pic_wd;
    ps_curr_lyr_res_prms->i4_res_height = ps_dec->u2_pic_ht;
    ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_left_offset << 1;
    ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_top_offset << 1;
    ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_rt =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_right_offset << 1;
    ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_bot =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_bottom_offset << 1;
    ps_curr_lyr_res_prms->u2_scaled_ref_width =
        (ps_dec->u2_frm_wd_in_mbs << 4) - (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left +
                                           ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_rt);

    ps_curr_lyr_res_prms->u2_scaled_ref_height =
        (ps_dec->u2_frm_ht_in_mbs << 4) - (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top +
                                           ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_bot);

    ps_curr_lyr_res_prms->u1_cropping_change_flag = 0;
    if(2 == ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc)
    {
        ps_curr_lyr_res_prms->u1_cropping_change_flag = 1;

        ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left =
            ps_svc_lyr_dec->s_svc_slice_params.i4_scaled_ref_layer_left_offset << 1;

        ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top =
            ps_svc_lyr_dec->s_svc_slice_params.i4_scaled_ref_layer_top_offset << 1;
        ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_rt =
            ps_svc_lyr_dec->s_svc_slice_params.i4_scaled_ref_layer_right_offset << 1;
        ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_bot =
            ps_svc_lyr_dec->s_svc_slice_params.i4_scaled_ref_layer_bottom_offset << 1;
        ps_curr_lyr_res_prms->u2_scaled_ref_width =
            (ps_dec->u2_frm_wd_in_mbs << 4) -
            (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left +
             ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_rt);

        ps_curr_lyr_res_prms->u2_scaled_ref_height =
            (ps_dec->u2_frm_ht_in_mbs << 4) -
            (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top +
             ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_bot);

        return NOT_OK;
    }

    ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_TRUE;

    ps_curr_lyr_res_prms->u1_disable_inter_lyr_dblk_filter_idc =
        ps_svc_lyr_dec->s_svc_slice_params.u4_disable_inter_layer_deblk_filter_idc;
    ps_curr_lyr_res_prms->i1_inter_lyr_alpha_c0_offset =
        ps_svc_lyr_dec->s_svc_slice_params.i4_inter_layer_slice_alpha_c0_offset_div2;
    ps_curr_lyr_res_prms->i1_inter_lyr_beta_offset =
        ps_svc_lyr_dec->s_svc_slice_params.i4_inter_layer_slice_beta_offset_div2;
    ps_curr_lyr_res_prms->i1_constrained_intra_rsmpl_flag =
        ps_svc_lyr_dec->s_svc_slice_params.u1_constrained_intra_resampling_flag;
    ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_x_plus1_flag =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_x_plus1_flag;
    ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_y_plus1 =
        ps_svc_lyr_dec->ps_cur_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_y_plus1;
    ps_curr_lyr_res_prms->u1_direct_8x8_inference_flag =
        ps_dec->ps_cur_sps->u1_direct_8x8_inference_flag;

    ps_curr_lyr_res_prms->u1_remap_req_flag = 1;
    ps_curr_lyr_res_prms->u1_dyadic_flag = ps_svc_lyr_dec->u1_dyadic_flag;

    /* Derive the reference layer width and height */

    if(SVCD_TRUE != ps_svc_lyr_dec->u1_base_res_flag)
    {
        WORD32 i4_ref_lyr_width;
        WORD32 i4_ref_lyr_ht;
        WORD32 i4_dyadic_flag = SVCD_FALSE;
        i4_ref_lyr_width = ps_svc_ref_lyr_dec->s_res_prms.i4_res_width;
        i4_ref_lyr_ht = ps_svc_ref_lyr_dec->s_res_prms.i4_res_height;

        /* set the Restricted Spatial Resolution change flag */
        ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_TRUE;

        if(0 == ((ps_curr_lyr_res_prms->u2_scaled_ref_width == i4_ref_lyr_width) ||
                 (ps_curr_lyr_res_prms->u2_scaled_ref_width == (i4_ref_lyr_width << 1))))
        {
            ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_FALSE;
        }

        if(0 == ((ps_curr_lyr_res_prms->u2_scaled_ref_height == i4_ref_lyr_ht) ||
                 (ps_curr_lyr_res_prms->u2_scaled_ref_height == (i4_ref_lyr_ht << 1))))
        {
            ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_FALSE;
        }

        if(0 != (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left & 15))
        {
            ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_FALSE;
        }

        if(0 != (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top & 15))
        {
            ps_curr_lyr_res_prms->u1_rstrct_res_change_flag = SVCD_FALSE;
        }

        /* populate the dyadic status */

        if((ps_curr_lyr_res_prms->u2_scaled_ref_width == (i4_ref_lyr_width << 1)) &&
           (ps_curr_lyr_res_prms->u2_scaled_ref_height == (i4_ref_lyr_ht << 1)))
        {
            i4_dyadic_flag = SVCD_TRUE;
        }
        else if((ps_curr_lyr_res_prms->u2_scaled_ref_width != ((i4_ref_lyr_width * 3) >> 1)) ||
                (ps_curr_lyr_res_prms->u2_scaled_ref_height != ((i4_ref_lyr_ht * 3) >> 1)))
        {
            ps_curr_lyr_res_prms->u1_dyadic_flag = i4_dyadic_flag;
            ps_svc_lyr_dec->u1_dyadic_flag = ps_curr_lyr_res_prms->u1_dyadic_flag;
            return NOT_OK;
        }

        /* check if cropping is MB aligned */
        if(SVCD_TRUE == i4_dyadic_flag)
        {
            if((0 != (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_left & 15)) ||
               (0 != (ps_curr_lyr_res_prms->s_ref_lyr_scaled_offset.i2_top & 15)))
            {
                i4_dyadic_flag = SVCD_FALSE;
            }
        }

        ps_curr_lyr_res_prms->u1_dyadic_flag = i4_dyadic_flag;
        ps_svc_lyr_dec->u1_dyadic_flag = ps_curr_lyr_res_prms->u1_dyadic_flag;
    }

    {
        inter_lyr_mb_prms_t *ps_tmp_prms, *ps_tmp_prms_2;
        inter_lyr_mb_prms_t *ps_ref_mb_prms;
        WORD32 i4_stride;
        WORD32 i4_ht_in_mbs, i4_wd_in_mbs;
        WORD32 i4_i;

        /* Derive the reference mb mode map */

        ps_ref_mb_prms = ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start;
        i4_stride = ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride;

        i4_ht_in_mbs = ps_dec->u2_frm_ht_in_mbs;
        i4_wd_in_mbs = ps_dec->u2_frm_wd_in_mbs;

        /* Set the first border row to 0xFF */
        ps_tmp_prms = (ps_ref_mb_prms - 1 - i4_stride);

        memset(ps_svc_lyr_dec->ps_inter_lyr_mb_prms_base, -1,
               ps_svc_lyr_dec->u4_inter_lyr_mb_prms_size);
        memset(ps_svc_lyr_dec->pu1_svc_base_mode_flag, 0,
               ps_svc_lyr_dec->i4_frm_svc_base_mode_cabac_size);

        for(i4_i = 0; i4_i < (i4_wd_in_mbs + 2); i4_i++)
        {
            ps_tmp_prms->i1_mb_mode = (WORD8) 0xFF;
            ps_tmp_prms += 1;
        }

        /* Set the left and right border pixels of each row to 0 */
        ps_tmp_prms = ps_ref_mb_prms - 1;

        for(i4_i = 0; i4_i < i4_ht_in_mbs; i4_i++)
        {
            ps_tmp_prms->i1_mb_mode = (WORD8) 0xFF;
            ps_tmp_prms_2 = ps_tmp_prms + (i4_wd_in_mbs + 1);
            ps_tmp_prms_2->i1_mb_mode = (WORD8) 0xFF;
            ps_tmp_prms += i4_stride;
        }

        /* Set the last border row to 0xFF */
        for(i4_i = 0; i4_i < (i4_wd_in_mbs + 2); i4_i++)
        {
            ps_tmp_prms->i1_mb_mode = (WORD8) 0xFF;
            ps_tmp_prms += 1;
        }
    }

    /* reset residual luma, chroma buffer*/
    memset(ps_svc_lyr_dec->pi2_il_residual_resample_luma_base, 0,
           ps_svc_lyr_dec->u4_residual_resample_luma_size);
    memset(ps_svc_lyr_dec->pi2_il_residual_resample_chroma_base, 0,
           ps_svc_lyr_dec->u4_residual_resample_chroma_size);

    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_crop_wnd_flag_res_int                               */
/*                                                                           */
/*  Description   : This routine computes the crop window flag for entire    */
/*                  dependency layer and places it in the crop window flag   */
/*                  buffer                                                   */
/*  Inputs        : 1. ECD context structure                                 */
/*                  2. Crop offset structure                                 */
/*  Globals       : None                                                     */
/*  Processing    : For Mbs within the crop window, flag set to 1 and for    */
/*                  others it is set to 0                                    */
/*                                                                           */
/*  Outputs       : Updates crop window flag buffer                          */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
void isvcd_crop_wnd_flag_res_int(void *pv_svc_dec)
{
    UWORD8 *pu1_crop_wnd_flag;
    WORD32 i4_num_mbs;
    WORD32 i4_crop_mbs_x;
    WORD32 i4_crop_mbs_y;
    WORD32 i4_cnt;
    WORD32 i4_left_offset, i4_rt_offset;
    WORD32 i4_top_offset, i4_bot_offset;
    WORD32 i4_frm_wd_in_mbs;
    WORD32 i4_frm_ht_in_mbs;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    res_prms_t *ps_res_prms;

    ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    ps_dec = &ps_svc_lyr_dec->s_dec;
    ps_res_prms = &ps_svc_lyr_dec->s_res_prms;
    i4_frm_wd_in_mbs = ps_dec->u2_frm_wd_in_mbs;
    i4_frm_ht_in_mbs = ps_dec->u2_frm_ht_in_mbs;

    /* Initializations */
    pu1_crop_wnd_flag = ps_svc_lyr_dec->pu1_crop_wnd_flag;
    i4_num_mbs = i4_frm_wd_in_mbs * i4_frm_ht_in_mbs;

    /* bottom most layer in a resolution */
    if(ps_res_prms->s_ref_lyr_scaled_offset.i2_left >= 0)
    {
        /* check for offset greater than 0 */
        i4_left_offset = (ps_res_prms->s_ref_lyr_scaled_offset.i2_left + 15) >> 4;
    }
    else
    {
        /* if negative set it to 0*/
        i4_left_offset = 0;
    }

    if(ps_res_prms->s_ref_lyr_scaled_offset.i2_rt >= 0)
    {
        /* check for offset greater than 0 */
        i4_rt_offset =
            (ps_res_prms->i4_res_width - ps_res_prms->s_ref_lyr_scaled_offset.i2_rt) >> 4;
    }
    else
    {
        /* if negative set it to framewidth in MBs */
        i4_rt_offset = (ps_res_prms->i4_res_width >> 4);
    }

    if(ps_res_prms->s_ref_lyr_scaled_offset.i2_top >= 0)
    {
        /* check for offset greater than 0 */
        i4_top_offset = (ps_res_prms->s_ref_lyr_scaled_offset.i2_top + 15) >> 4;
    }
    else
    {
        /* if negative set it to 0 */
        i4_top_offset = 0;
    }

    if(ps_res_prms->s_ref_lyr_scaled_offset.i2_bot >= 0)
    {
        /* check for offset greater than 0 */
        i4_bot_offset =
            (ps_res_prms->i4_res_height - ps_res_prms->s_ref_lyr_scaled_offset.i2_bot) >> 4;
    }
    else
    {
        /* if negative set it to frameheight in MBs */
        i4_bot_offset = (ps_res_prms->i4_res_height >> 4);
    }

    i4_crop_mbs_x = i4_rt_offset - i4_left_offset;
    i4_crop_mbs_y = i4_bot_offset - i4_top_offset;

    /* Set crop window flag to 0 for all mbs */
    memset(pu1_crop_wnd_flag, 0, i4_num_mbs);

    pu1_crop_wnd_flag += (i4_frm_wd_in_mbs * i4_top_offset);
    pu1_crop_wnd_flag += i4_left_offset;
    /* Loop over MBs in crop window */
    for(i4_cnt = 0; i4_cnt < i4_crop_mbs_y; i4_cnt++)
    {
        memset(pu1_crop_wnd_flag, 1, i4_crop_mbs_x);
        pu1_crop_wnd_flag += i4_frm_wd_in_mbs;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_intra_resamp_res_init                                 */
/*                                                                           */
/*  Description   : this function calculates the scale factors and initialise*/
/*                  the context structure                                    */
/*                                                                           */
/*  Inputs        : pv_intra_samp_ctxt: handle to private structure          */
/*                  ps_curr_lyr_res_prms: pointer to current resolution      */
/*                                               params                      */
/*                  ps_ref_lyr_res_prms : pointer to ref resolution params   */
/*  Globals       : none                                                     */
/*  Processing    : it stores the layer dimensions                           */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 08 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_intra_resamp_res_init(void *pv_svc_dec)
{
    intra_sampling_ctxt_t *ps_ctxt;
    intra_samp_lyr_ctxt *ps_lyr_ctxt;
    dec_svc_seq_params_t *ps_cur_subset_sps;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    void *pv_intra_samp_ctxt = ps_svc_lyr_dec->pv_intra_sample_ctxt;
    res_prms_t *ps_curr_lyr_res_prms = &ps_svc_lyr_dec->s_res_prms;
    ref_mb_map_t **pps_luma_map_horz = &ps_svc_lyr_dec->ps_intsam_luma_map_horz;
    ref_mb_map_t **pps_chroma_map_horz = &ps_svc_lyr_dec->ps_intsam_chroma_map_horz;
    ref_mb_map_t **pps_luma_map_vert = &ps_svc_lyr_dec->ps_intsam_luma_map_vert;
    ref_mb_map_t **pps_chroma_map_vert = &ps_svc_lyr_dec->ps_intsam_chroma_map_vert;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    ps_cur_subset_sps = ps_svc_lyr_dec->ps_cur_subset_sps;

    ps_ctxt = (intra_sampling_ctxt_t *) pv_intra_samp_ctxt;

    /* if called for base resolution store default values */
    if(SVCD_TRUE == ps_svc_lyr_dec->u1_base_res_flag)
    {
        *pps_luma_map_horz = NULL;
        *pps_chroma_map_horz = NULL;
        *pps_luma_map_vert = NULL;
        *pps_chroma_map_vert = NULL;
        ps_ctxt->i4_res_lyr_id = -1;
        ps_ctxt->i4_ref_width = ps_dec->u2_pic_wd;
        ps_ctxt->i4_ref_height = ps_dec->u2_pic_ht;

        /* Note: The stride option is provided for bringing in data at NMB */
        /* level. Hence to set a NMB level stride refSample array buffer   */
        /* have to be increased                                            */
        ps_ctxt->i4_refarray_stride = REF_ARRAY_WIDTH;
        return OK;
    }

    /* derive the current sps */
    /* store the res id appropriately */
    ps_ctxt->i4_res_lyr_id = ps_svc_lyr_dec->u1_layer_id - 1;

    /* store the resolution params */
    ps_ctxt->ps_res_prms = ps_curr_lyr_res_prms;

    /* get the current layer ctxt */
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_svc_lyr_dec->u1_layer_id - 1];

    ps_ctxt->i4_res_lyr_id = ps_svc_lyr_dec->u1_layer_id - 1;
    /* get the width and heights */
    ps_lyr_ctxt->i4_curr_width = ps_dec->u2_pic_wd;
    ps_lyr_ctxt->i4_curr_height = ps_dec->u2_pic_ht;
    ps_lyr_ctxt->i4_ref_width = ps_ctxt->i4_ref_width;
    ps_lyr_ctxt->i4_ref_height = ps_ctxt->i4_ref_height;
    ps_lyr_ctxt->i1_constrained_intra_rsmpl_flag =
        ps_svc_slice_params->u1_constrained_intra_resampling_flag;

    /* store the structure pointer containing projected locations */
    *pps_luma_map_horz = ps_lyr_ctxt->s_luma_map_ctxt.ps_x_offset_length;
    *pps_chroma_map_horz = ps_lyr_ctxt->s_chroma_map_ctxt.ps_x_offset_length;
    *pps_luma_map_vert = ps_lyr_ctxt->s_luma_map_ctxt.ps_y_offset_length;
    *pps_chroma_map_vert = ps_lyr_ctxt->s_chroma_map_ctxt.ps_y_offset_length;

    /* check for recomputation of mapping required */
    if(SVCD_TRUE == ps_curr_lyr_res_prms->u1_remap_req_flag)
    {
        res_prms_t s_ref_res_prms = {0};
        WORD32 i4_chroma_x_phase, i4_chroma_y_phase;
        WORD32 i4_ref_chroma_x_phase, i4_ref_chroma_y_phase;
        WORD32 i4_x_phase_0, i4_x_phase_1;
        WORD32 i4_y_phase_0, i4_y_phase_1;
        WORD32 i4_vert_flag;

        /* store the reference layer resolution width and height */
        s_ref_res_prms.i4_res_width = ps_ctxt->i4_ref_width;
        s_ref_res_prms.i4_res_height = ps_ctxt->i4_ref_height;

        /* call the frame level projections calculation function */
        isvcd_intra_resamp_populate_list(&ps_lyr_ctxt->s_luma_map_ctxt, ps_curr_lyr_res_prms,
                                         &s_ref_res_prms, 0, ps_svc_lyr_dec);

        isvcd_intra_resamp_populate_list(&ps_lyr_ctxt->s_chroma_map_ctxt, ps_curr_lyr_res_prms,
                                         &s_ref_res_prms, 1, ps_svc_lyr_dec);

        /* Compute the chroma xPhase and yPhase values */
        if(1 == ps_curr_lyr_res_prms->u1_dyadic_flag)
        {
            i4_ref_chroma_x_phase = ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_x_plus1_flag;
            i4_ref_chroma_y_phase = ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_y_plus1;
            i4_chroma_x_phase = ps_cur_subset_sps->s_sps_svc_ext.u1_chroma_phase_x_plus1_flag;
            i4_chroma_y_phase = ps_cur_subset_sps->s_sps_svc_ext.u1_chroma_phase_y_plus1;

            i4_x_phase_0 = i4_chroma_x_phase - (i4_ref_chroma_x_phase << 1);
            i4_x_phase_1 = (3 + i4_x_phase_0) & 0x7;
            i4_x_phase_0 += 7;
            i4_x_phase_0 &= 0x7;
            i4_y_phase_0 = i4_chroma_y_phase - (i4_ref_chroma_y_phase << 1);
            i4_y_phase_1 = (3 + i4_y_phase_0) & 0x7;
            i4_y_phase_0 += 7;
            i4_y_phase_0 &= 0x7;

            ps_lyr_ctxt->i4_x_phase_0 = i4_x_phase_0;
            ps_lyr_ctxt->i4_x_phase_1 = i4_x_phase_1;
            ps_lyr_ctxt->i4_y_phase_0 = i4_y_phase_0;
            ps_lyr_ctxt->i4_y_phase_1 = i4_y_phase_1;

            /* Choose the appropriate chroma interpolation functions */
            if((0 == i4_ref_chroma_x_phase) && (1 == i4_chroma_x_phase))
            {
                ps_lyr_ctxt->pf_horz_chroma_interpol = ps_ctxt->pf_horz_chroma_interpol[1];
            }
            else
            {
                ps_lyr_ctxt->pf_horz_chroma_interpol = ps_ctxt->pf_horz_chroma_interpol[0];
            }

            i4_vert_flag = 0;
            if(0 == i4_ref_chroma_y_phase)
            {
                if((1 == i4_chroma_y_phase) || (2 == i4_chroma_y_phase))
                {
                    i4_vert_flag = 1;
                }
            }
            else if((2 == i4_ref_chroma_y_phase) && (0 == i4_chroma_y_phase))
            {
                i4_vert_flag = 2;
            }

            if(1 == i4_vert_flag)
            {
                ps_lyr_ctxt->pf_vert_chroma_interpol = ps_ctxt->pf_vert_chroma_interpol[1];
            }
            else if(2 == i4_vert_flag)
            {
                ps_lyr_ctxt->pf_vert_chroma_interpol = ps_ctxt->pf_vert_chroma_interpol[2];
            }
            else
            {
                ps_lyr_ctxt->pf_vert_chroma_interpol = ps_ctxt->pf_vert_chroma_interpol[0];
            }
        }
    }
    else
    {
        /* should take false value */
        if(SVCD_FALSE != ps_curr_lyr_res_prms->u1_remap_req_flag)
        {
            return NOT_OK;
        }
    }

    /* store the current layer width and height to context */
    ps_ctxt->i4_ref_width = ps_curr_lyr_res_prms->i4_res_width;
    ps_ctxt->i4_ref_height = ps_curr_lyr_res_prms->i4_res_height;

    return OK;
}