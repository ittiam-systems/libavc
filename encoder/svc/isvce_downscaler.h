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
*  isvce_downscaler.h
*
* @brief
*  Contains downscaler functions required by the SVC encoder
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_get_downscaler_data_size()
*  - isvce_get_downscaler_padding_dims()
*  - isvce_isvce_process_ctxt_t_downscaler()
*  - isvce_get_downscaler_normalized_filtered_pixel()
*  - isvce_horizontal_downscale_and_transpose()
*  - isvce_process_downscaler()
*  - isvce_initialize_downscaler()
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_DOWNSCALER_H_
#define _ISVCE_DOWNSCALER_H_

#include "ih264_typedefs.h"
#include "iv2.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_defs.h"

typedef struct
{
    /**
     * pointer to the state of downscaler
     */
    void *pv_scaler_state;

    /**
     * scaling factor between the dimensions of two consecutive SVC layers
     */
    DOUBLE d_scaling_factor;

    /**
     * Num spatial layers
     */
    UWORD8 u1_num_spatial_layers;

} downscaler_ctxt_t;

typedef struct
{
    UWORD8 u1_left_pad_size;

    UWORD8 u1_right_pad_size;

    UWORD8 u1_top_pad_size;

    UWORD8 u1_bottom_pad_size;

} padding_dims_t;

/**
*******************************************************************************
*
* @brief
*   initializes the downscaler context
*
* @par Description:
*   initializes the downscaler context for the given scaling factor
*   with padding size, filter size, etc.
*
* @param[in] ps_scaler
*   pointer downscaler context
*
* @param[in] ps_mem_rec
*   pointer to memory allocated to downscaler process
*
* @param[in] d_scaling_factor
*   scaling reatio of width/ height between two consecutive SVC layers
*
* @param[in] u1_num_spatial_layers
*   scaling reatio of width/ height between two consecutive SVC layers
*
* @param[in] u4_wd
*   width of the input
*
* @param[in] u4_ht
*   height of the input
*
* @param[in] e_arch
*   architecure type
*
* @returns
*
* @remarks
*  when ARM intrinsics are added, update should be done here
*
*******************************************************************************
*/

extern void isvce_initialize_downscaler(downscaler_ctxt_t *ps_scaler, iv_mem_rec_t *ps_mem_rec,
                                        DOUBLE d_scaling_factor, UWORD8 u1_num_spatial_layers,
                                        UWORD32 u4_in_width, UWORD32 u4_in_height,
                                        IV_ARCH_T e_arch);

/**
*******************************************************************************
*
* @brief
*   gets the memory size required for downscaler
*
* @par Description:
*   returns the memory required by the downscaler context and state structs
*   for allocation.
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

extern UWORD32 isvce_get_downscaler_data_size(UWORD8 u1_num_spatial_layers, DOUBLE d_scaling_factor,
                                              UWORD32 u4_width, UWORD32 u4_height);

/**
*******************************************************************************
*
* @brief
*   processes downscaler
*
* @par Description:
*   calls the function for padding and scaling
*
* @param[in] ps_scaler
*  pointer to downdownscaler context
*
* @param[in] ps_src_buf_props
*  pointer to source buffer props struct
*
* @param[in] u4_blk_wd
*  width of the block to be processed
*
* @param[in] u4_blk_ht
*  height of the block to be processed
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

extern void isvce_process_downscaler(downscaler_ctxt_t *ps_scaler,
                                     yuv_buf_props_t *ps_src_buf_props,
                                     yuv_buf_props_t *ps_dst_buf_props, UWORD32 u4_blk_wd,
                                     UWORD32 u4_blk_ht);

/**
*******************************************************************************
*
* @brief
*   gets the padding size required for filtering
*
* @par Description:
*   gets the padding size required for filtering
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

extern void isvce_get_downscaler_padding_dims(padding_dims_t *ps_pad_dims);

#endif
