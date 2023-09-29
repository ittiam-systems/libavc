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
*  ih264e_time_stamp.h
*
* @brief
*  This file contains function declarations used for managing input and output
*  frame time stamps
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_TIME_STAMP_H_
#define _IH264E_TIME_STAMP_H_

/*****************************************************************************/
/* Structures                                                                */
/*****************************************************************************/

/**
 * Parameters for Src/Tgt frames that are encoded
 */
typedef struct frame_time_t
{
    /* common time base(=LCM) between source and target frame rate (in ticks)*/
    WORD32 common_time_base;

    /* number of ticks between two source frames */
    UWORD32 u4_src_frm_time_incr;

    /* number of ticks between two target frames */
    UWORD32 u4_tgt_frm_time_incr;

    /* Source frame time - measured as modulo of common time base
     and incremented by src_frm_time_incr */
    UWORD32 u4_src_frm_time;

    /* Target frame time - measured as modulo of common time base
     and incremented by tgt_frm_time_incr */
    UWORD32 u4_tgt_frm_time;

    /* Number of frames not to be skipped while maintaining
     tgt_frm_rate due to delta_time_stamp  */
    UWORD32 u4_num_frms_dont_skip;

}frame_time_t;

typedef struct frame_time_t *frame_time_handle;

/**
 *  Parameters that go in the bitstream based on tgt_frm_rate
 *   1) Initialize the vop_time_incr_res with the max_frame_rate (in frames per 1000 bits)
 *      - To represent all kinds of frame rates
 *   2) Decide the vop_time_incr based on the source frame rate
 *      - The decoder would like to know which source frame is encoded i.e. the source time
 *    id of the target frame encoded and there by adjusting its time of delay
 *   3) vop_time increments every source frame and whenever a frame is encoded (target frame),
 *      the encoder queries the vop time of the source frame and sends it in the bit stream.
 *   4) Since the Source frame skip logic is taken care by the frame_time module, whenever the
 *      encoder queries the time stamp module (which gets updated outside the encoder) the
 *      time stamp module would have the source time
 */
typedef struct time_stamp_t
{
    /*vop_time_incr_res is a integer that indicates
     the number of evenly spaced subintervals, called ticks,
     within one modulo time. */
    UWORD32 u4_vop_time_incr_res;

    /* number of bits to represent vop_time_incr_res */
    UWORD32 u4_vop_time_incr_range;

    /* The number of ticks elapsed between two source vops */
    UWORD32 u4_vop_time_incr;

    /* incremented by vop_time_incr for every source frame.
     Represents the time offset after a modulo_time_base = 1 is sent
     in bit stream*/
    UWORD32 u4_vop_time;

    /* A temporary buffer to copy of vop time and modulo time base
     is stored since update is called before query (get time stamp) and
     so these extra variables cur_tgt_vop_time,  */
    UWORD32 u4_cur_tgt_vop_time;

    UWORD32 u4_prev_tgt_vop_time;

    /* This variable is set to 1 if we scale max frame rate by a factor of 2.
     For mpeg4 standard, we just have 16bits and we can't accommodate more than 60000 as frame rate.
     So we scale it and work with it */
    WORD32 is_max_frame_rate_scaled;

} time_stamp_t;

typedef struct time_stamp_t *time_stamp_handle;

/*****************************************************************************/
/* Function declarations                                                     */
/*****************************************************************************/

void ih264e_init_frame_time(frame_time_t *ps_frame_time,
                            UWORD32 u4_src_frm_rate,
                            UWORD32 u4_tgt_frm_rate);

UWORD8 ih264e_should_src_be_skipped(frame_time_t *ps_frame_time,
                                    UWORD32 u4_delta_time_stamp,
                                    UWORD32 *pu4_frm_not_skipped_for_dts);

void ih264e_init_time_stamp(time_stamp_handle time_stamp,
                            UWORD32 max_frm_rate,
                            UWORD32 src_frm_rate);

void ih264e_update_time_stamp(time_stamp_handle time_stamp);

WORD32 ih264e_frame_time_get_init_free_memtab(frame_time_handle *pps_frame_time,
                                              itt_memtab_t *ps_memtab,
                                              ITT_FUNC_TYPE_E e_func_type);

WORD32 ih264e_time_stamp_get_init_free_memtab(time_stamp_handle *pps_time_stamp,
                                              itt_memtab_t *ps_memtab,
                                              ITT_FUNC_TYPE_E e_func_type);

WORD32 ih264e_frame_time_get_src_frame_rate(frame_time_t *ps_frame_time);

WORD32 ih264e_frame_time_get_tgt_frame_rate(frame_time_t *ps_frame_time);

WORD32 ih264e_frame_time_get_src_ticks(frame_time_t *ps_frame_time);

WORD32 ih264e_frame_time_get_tgt_ticks(frame_time_t *ps_frame_time);

WORD32 ih264e_frame_time_get_src_time(frame_time_t *frame_time);

WORD32 ih264e_frame_time_get_tgt_time(frame_time_t *frame_time);

void ih264e_frame_time_update_src_frame_rate(frame_time_t *ps_frame_time,
                                             WORD32 src_frm_rate);

void ih264e_frame_time_update_tgt_frame_rate(frame_time_t *ps_frame_time,
                                             WORD32 tgt_frm_rate);

void ih264_time_stamp_update_frame_rate(time_stamp_t *ps_time_stamp,
                                        UWORD32 src_frm_rate);

#endif /*_IH264E_TIME_STAMP_H_ */

