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
*  ih264e_process.h
*
* @brief
*  Contains functions for codec thread
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_PROCESS_H_
#define _IH264E_PROCESS_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

IH264E_ERROR_T ih264e_generate_sps_pps(codec_t *ps_codec);

IH264E_ERROR_T ih264e_init_entropy_ctxt(process_ctxt_t *ps_proc);

IH264E_ERROR_T ih264e_entropy(process_ctxt_t *ps_proc);

IH264E_ERROR_T ih264e_pack_header_data(process_ctxt_t *ps_proc);

WORD32 ih264e_update_proc_ctxt(process_ctxt_t *ps_proc);

IH264E_ERROR_T ih264e_init_proc_ctxt(process_ctxt_t *ps_proc);

IH264E_ERROR_T ih264e_pad_recon_buffer(process_ctxt_t *ps_proc,
                                       UWORD8 *pu1_curr_pic_luma,
                                       UWORD8 *pu1_curr_pic_chroma,
                                       WORD32 i4_mb_x,
                                       WORD32 i4_mb_y,
                                       WORD32 i4_pad_ht);

IH264E_ERROR_T ih264e_halfpel_generation(process_ctxt_t *ps_proc,
                                         UWORD8 *pu1_curr_pic_luma,
                                         WORD32 i4_mb_x,
                                         WORD32 i4_mb_y);

WORD32 ih264e_process(process_ctxt_t *ps_proc);

WORD32 ih264e_update_rc_post_enc(codec_t *ps_codec, WORD32 ctxt_sel, WORD32 pic_cnt);

WORD32 ih264e_process_thread(void *pv_proc);

#endif /* _IH264E_PROCESS_H_ */
