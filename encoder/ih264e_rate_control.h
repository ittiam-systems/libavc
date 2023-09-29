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
*  ih264e_rate_control.h
*
* @brief
*  This file contains declarations of api functions for h264 rate control
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_RATE_CONTROL_H_
#define _IH264E_RATE_CONTROL_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_rc_init(void *pv_rc_api,
                    void *pv_frame_time,
                    void *pv_time_stamp,
                    void *pv_pd_frm_rate,
                    UWORD32 u4_max_frm_rate,
                    UWORD32 u4_src_frm_rate,
                    UWORD32 u4_tgt_frm_rate,
                    rc_type_e e_rate_control_type,
                    UWORD32 u4_avg_bit_rate,
                    UWORD32 u4_peak_bit_rate,
                    UWORD32 u4_max_delay,
                    UWORD32 u4_intra_frame_interval,
                    WORD32  i4_inter_frm_int,
                    UWORD8 *pu1_init_qp,
                    WORD32 i4_max_inter_frm_int,
                    UWORD8 *pu1_min_max_qp,
                    UWORD8 u1_profile_level);

picture_type_e ih264e_rc_get_picture_details(void *pv_rc_api,
                                             WORD32 *pi4_pic_id,
                                             WORD32 *pi4_pic_disp_order_no);

WORD32 ih264e_update_rc_framerates(void *ps_rate_control_api,
                                   void *ps_pd_frm_rate,
                                   void *ps_time_stamp,
                                   void *ps_frame_time);

void ih264e_update_rc_mb_info(frame_info_t *ps_frame_info, void *pv_proc);

void ih264e_rc_get_buffer_status(void *pv_rc_api,
                                 WORD32 i4_total_frame_bits,
                                 picture_type_e e_pic_type,
                                 WORD32 *pi4_num_bits_to_prevent_vbv_underflow,
                                 UWORD8 *pu1_is_enc_buf_overflow,
                                 UWORD8 *pu1_is_enc_buf_underflow);

WORD32 ih264e_rc_post_enc(void *ps_rate_control_api,
                         frame_info_t *ps_frame_info,
                         void *ps_pd_frm_rate,
                         void *ps_time_stamp,
                         void *ps_frame_time,
                         WORD32 i4_total_mb_in_frame,
                         picture_type_e *pe_vop_coding_type,
                         WORD32 i4_is_first_frame,
                         WORD32 *pi4_is_post_encode_skip,
                         UWORD8 u1_frame_qp,
                         WORD32 *pi4_num_intra_in_prev_frame,
                         WORD32 *pi4_avg_activity);

void ih264e_update_rc_bits_info(frame_info_t *ps_frame_info, void *pv_entropy);

#endif /* _IH264E_RATE_CONTROL_H_ */

