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
/*  File Name         : imvcd_api_utils.h                                    */
/*                                                                           */
/*  Description       : Utility functions used by 'imvcd_api.c'              */
/*                                                                           */
/*****************************************************************************/

#ifndef _IMVCD_API_UTILS_H_
#define _IMVCD_API_UTILS_H_

#include "iv.h"
#include "ih264d_structs.h"
#include "imvcd_structs.h"

extern IV_API_CALL_STATUS_T imvcd_check_dec_handle(iv_obj_t *ps_handle);

extern IV_API_CALL_STATUS_T imvcd_check_create_structs(imvcd_create_ip_t *ps_ip,
                                                       imvcd_create_op_t *ps_op);

extern IV_API_CALL_STATUS_T imvcd_check_ctl_structs(void *pv_ip, void *pv_op);

extern IV_API_CALL_STATUS_T imvcd_check_decode_structs(iv_obj_t *ps_dec_hdl,
                                                       imvcd_video_decode_ip_t *ps_ip,
                                                       imvcd_video_decode_op_t *ps_op);

extern void imvcd_au_init(iv_obj_t *ps_dec_hdl, imvcd_video_decode_ip_t *ps_ip,
                          imvcd_video_decode_op_t *ps_op);

extern void imvcd_view_init(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern IV_API_CALL_STATUS_T imvcd_bitstream_buf_alloc(dec_struct_t *ps_view_ctxt,
                                                      UWORD16 u2_num_views);

extern void imvcd_bitsteam_buf_free(dec_struct_t *ps_view_ctxt);

extern void imvcd_convert_to_app_disp_buf(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                          iv_yuv_buf_t *ps_view_disp_bufs);
#endif
