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
*  isvce_interface_structs.h
*
* @brief
*  Contains struct definition used for interface objects such as input,
*  output, and rec
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_INTERFACE_STRUCTS_H_
#define _ISVCE_INTERFACE_STRUCTS_H_

#include "isvc_structs.h"

typedef struct isvce_raw_inp_buf_t
{
    /** Descriptor of raw buffer                                     */
    iv_raw_buf_t s_raw_buf;

    /** Lower 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_low;

    /** Upper 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_high;

    /** Flag to indicate if the current buffer is last buffer */
    UWORD32 u4_is_last;

    /** Flag to indicate if mb info is sent along with input buffer     */
    UWORD32 u4_mb_info_type;

    /** Flag to indicate the size of mb info structure                  */
    UWORD32 u4_mb_info_size;

    /** Buffer containing mb info if isvce_mb_info_type is non-zero           */
    void *pv_mb_info;

    /** Flag to indicate if pic info is sent along with input buffer     */
    UWORD32 u4_pic_info_type;

    /** Buffer containing pic info if isvce_mb_info_type is non-zero           */
    void *pv_pic_info;

    /** SEI CCV params flag                                              */
    UWORD8 u1_sei_ccv_params_present_flag;

    /** SEI CCV params info                                              */
    sei_ccv_params_t s_sei_ccv;

} isvce_raw_inp_buf_t;

typedef struct
{
    /** Descriptor of bitstream buffer                                     */
    iv_bits_buf_t as_bits_buf[MAX_NUM_SPATIAL_LAYERS];

    /** Lower 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_low;

    /** Upper 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_high;

    /** Flag to indicate if the current buffer is last buffer */
    UWORD32 u4_is_last;

} isvce_out_buf_t;

typedef struct
{
    /** Descriptor of picture buffer                                     */
    svc_au_buf_t s_pic_buf;

    /** Lower 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_low;

    /** Upper 32bits of time stamp corresponding to the above buffer */
    UWORD32 u4_timestamp_high;

    /** Flag to indicate if the current buffer is last buffer */
    UWORD32 u4_is_last;

    /** Picture count corresponding to current picture */
    WORD32 i4_pic_cnt;

} isvce_rec_buf_t;

#endif
