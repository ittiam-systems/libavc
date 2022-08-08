/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

/*****************************************************************************/
/*                                                                           */
/*  File Name         : ih264d_api_utils.h                                   */
/*                                                                           */
/*  Description       : This file contains function prototypes for non-API   */
/*                      functions defined in ih264d_api.c                    */
/*****************************************************************************/

#ifndef _IH264D_API_UTILS_H_
#define _IH264D_API_UTILS_H_

#include "ih264_typedefs.h"
#include "iv.h"

#define MIN_IN_BUFS 1
#define MIN_OUT_BUFS_420 3
#define MIN_OUT_BUFS_422ILE 1
#define MIN_OUT_BUFS_RGB565 1
#define MIN_OUT_BUFS_420SP 2

extern WORD32 ih264d_create(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_delete(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_video_decode(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_set_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_set_num_cores(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_set_processor(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_set_degrade(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_set_flush_mode(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern WORD32 ih264d_get_buf_info(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

extern UWORD32 ih264d_get_outbuf_size(WORD32 u4_pic_wd, UWORD32 u4_pic_ht, UWORD8 u1_chroma_format,
                                      UWORD32 *pu4_buf_size);

extern void ih264d_export_sei_params(ivd_sei_decode_op_t *ps_sei_decode_op, dec_struct_t *ps_dec);

extern UWORD32 ih264d_map_error(UWORD32 i4_err_status);

extern void ih264d_signal_decode_thread(dec_struct_t *ps_dec);

extern void ih264d_signal_bs_deblk_thread(dec_struct_t *ps_dec);

extern WORD32 ih264d_deblock_display(dec_struct_t *ps_dec);

#endif
