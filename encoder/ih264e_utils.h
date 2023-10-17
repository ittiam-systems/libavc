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
*  ih264e_utils.h
*
* @brief
*  Contains declarations of miscellaneous utility functions used by the encoder
*
* @author
*  Harish
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_UTILS_H_
#define _IH264E_UTILS_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

WORD32 ih264e_input_queue_update(codec_t *ps_codec,
                                 ive_video_encode_ip_t *ps_ive_ip,
                                 inp_buf_t *ps_enc_buff);

WORD32 ih264e_get_min_level(WORD32 wd, WORD32 ht);

WORD32 ih264e_get_lvl_idx(WORD32 level);

WORD32 ih264e_get_dpb_size(WORD32 level, WORD32 pic_size);

WORD32 ih264e_get_total_pic_buf_size(WORD32 pic_size, WORD32 level,
                                     WORD32 horz_pad, WORD32 vert_pad,
                                     WORD32 num_ref_frames,
                                     WORD32 num_reorder_frames);

WORD32 ih264e_get_pic_mv_bank_size(WORD32 num_luma_samples);

IH264E_ERROR_T ih264e_pic_buf_mgr_add_bufs(codec_t *ps_codec);

IH264E_ERROR_T ih264e_mv_buf_mgr_add_bufs(codec_t *ps_codec);

void ih264e_init_quant_params(process_ctxt_t *ps_proc, int qp);

IH264E_ERROR_T ih264e_init_air_map(codec_t *ps_codec);

void ih264e_speed_preset_side_effects(codec_t *ps_codec);

IH264E_ERROR_T ih264e_codec_init(codec_t *ps_codec,
                                 UWORD32 u4_timestamp_low,
                                 UWORD32 u4_timestamp_high);

IH264E_ERROR_T ih264e_pic_init(codec_t *ps_codec, inp_buf_t *ps_inp_buf);

#endif /* _IH264E_UTILS_H_ */
