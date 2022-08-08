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
/*  File Name         : imvcd_utils.c                                        */
/*                                                                           */
/*  Description       : MVCD Utility functions used by 'imvcd_api.c'         */
/*                                                                           */
/*****************************************************************************/

#ifndef _IMVCD_UTILS_H_
#define _IMVCD_UTILS_H_

#include <stdbool.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "imvc_defs.h"
#include "ih264d_mvpred.h"
#include "ih264d_structs.h"
#include "imvcd_structs.h"

#define SWAP(x, y, data_type)                 \
    {                                         \
        data_type temp;                       \
        memcpy(&temp, &y, sizeof(data_type)); \
        memcpy(&y, &x, sizeof(data_type));    \
        memcpy(&x, &temp, sizeof(data_type)); \
    }

extern IV_API_CALL_STATUS_T imvcd_get_next_display_au_buf(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern UWORD32 imvcd_get_num_mbs_in_level(UWORD8 u1_level_idc);

extern WORD32 imvcd_allocate_dynamic_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern WORD16 imvcd_free_dynamic_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern WORD32 imvcd_init_au_buffers(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern WORD32 imvcd_init_au_mv_pred_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern void imvcd_convert_au_buf_to_view_buf(mvc_au_buffer_t *ps_au_buf, pic_buffer_t *ps_view_buf,
                                             UWORD16 u2_view_order_id, UWORD16 u2_view_id);

extern void imvcd_init_ref_idx_to_ref_buf_map(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern void imvcd_ivp_buf_copier(mvc_au_buffer_t *ps_au_buf_src, mvc_au_buffer_t *ps_au_buf_dst,
                                 mvc_au_mv_pred_t *ps_au_mv_data_src,
                                 mvc_au_mv_pred_t *ps_au_mv_data_dst, UWORD16 u2_src_view_id,
                                 UWORD16 u2_dst_view_id);

/* Function defined in 'ih264d_utils.c' and declared nowhere else */
extern WORD32 ih264d_init_dec_mb_grp(dec_struct_t *ps_dec);

extern void ih264d_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec);

extern void ih264d_get_implicit_weights(dec_struct_t *ps_dec);

extern void imvcd_free_ref_bufs(mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                                mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr,
                                WORD32 i4_pic_buf_id);

extern void imvcd_release_all_ref_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt, WORD32 i4_num_bufs);

extern void imvcd_free_ref_and_io_bufs(mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                                       mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr,
                                       WORD32 i4_pic_buf_id);

extern void imvcd_release_all_ref_and_io_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt, WORD32 i4_num_bufs);

extern bool is_header_decoded(WORD32 i4_header_decoded, AVC_EXT_NALU_ID_T e_nalu_id);

extern bool is_mvc_nalu(AVC_EXT_NALU_ID_T e_nalu_id);

extern bool is_slice_nalu_type(AVC_EXT_NALU_ID_T e_nalu_id);

extern nalu_mvc_ext_t *imvcd_get_cur_nalu_mvc_ext(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern nalu_mvc_ext_t *imvcd_get_nalu_mvc_ext(nalu_mvc_ext_t *ps_nalu_mvc_exts,
                                              UWORD16 u2_num_views_decoded, UWORD16 u2_view_id);

extern ref_pic_list_mod_data_t *imvcd_get_cur_ref_pic_list_mod_data(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern subset_sps_t *imvcd_get_valid_subset_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern void imvcd_modulate_max_disp_seq(dec_struct_t *ps_view_ctxt);

extern mv_pred_t imvcd_get_default_mv_pred(void);

extern UWORD32 imvcd_get_max_num_ivp_refs(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern bool imvcd_is_idr_au(mvc_dec_ctxt_t *ps_mvcd_ctxt);

extern coordinates_t imvcd_get_buf_pad_dims(bool b_is_chroma);

extern WORD32 imvcd_get_ref_pic_pad_offset(WORD32 i4_stride, bool b_is_chroma);

extern UWORD32 imvcd_get_next_bits(dec_bit_stream_t *ps_bitstream);

extern void imvcd_set_view_buf_id_to_buf_map(dec_struct_t *ps_view_ctxt);

#endif
