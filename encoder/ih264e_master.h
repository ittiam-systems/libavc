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
*  ih264e_master.h
*
* @brief
*  Contains declarations of functions used by master thread
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_MASTER_H_
#define _IH264E_MASTER_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_join_threads(codec_t *ps_codec);

void ih264e_compute_quality_stats(process_ctxt_t *ps_proc);

IH264E_ERROR_T ih264e_wait_for_thread(UWORD32 sleep_us);

WORD32 ih264e_encode(iv_obj_t *ps_codec_obj, void *pv_api_ip, void *pv_api_op);

IH264E_ERROR_T ih264e_codec_update_config(codec_t *ps_codec, cfg_params_t *ps_cfg);

#endif /* _IH264E_MASTER_H_ */
