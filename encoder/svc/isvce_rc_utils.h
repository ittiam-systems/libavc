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
*  isvce_rc_utils.h
*
* @brief
*  Contains get gpp function required by the SVC encoder
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_rc_utils_init()
*  - isvce_get_rc_utils_data_size()
*  - isvce_compute_gpp()
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_RC_UTILS_H_
#define _ISVCE_RC_UTILS_H_

#include "ih264_typedefs.h"
#include "isvc_structs.h"

typedef struct
{
    /**
     * pointer to the state of rc utils
     */
    void *pv_rc_utils_state;

} svc_rc_utils_ctxt_t;

/**
*******************************************************************************
*
* @brief
*   initializes the rc utils context
*
* @par Description:
*   initializes the rc utils context
*
* @param[in] ps_svc_rc_utils_ctxt
*   pointer to svc rc utils context
*
* @param[in] ps_mem_rec
*   pointer to memory allocated to compute gpp process
*
* @param[in] e_arch
*   architecure type
*
* @returns
*
* @remarks
*  none
*
*******************************************************************************
*/

extern void isvce_rc_utils_init(svc_rc_utils_ctxt_t *ps_svc_rc_utils_ctxt, iv_mem_rec_t *ps_mem_rec,
                                IV_ARCH_T e_arch);

/**
*******************************************************************************
*
* @brief
*   gets the memory size required for compute gpp
*
* @par Description:
*   returns the memory required by the rc utils context and state structs
*   for allocation.
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

extern UWORD32 isvce_get_rc_utils_data_size();

/**
*******************************************************************************
*
* @brief
*   compute gpp process
*
* @par Description:
*   calls the function to compute gpp
*
* @param[in] ps_svc_rc_utils_ctxt
*  pointer to svc rc utils context
*
* @param[in] ps_input_buf
*  pointer to yuv buffer properties
*
* @returns
*  calculated gpp value
*
* @remarks
*  none
*
*******************************************************************************
*/

extern DOUBLE isvce_compute_gpp(svc_rc_utils_ctxt_t *ps_svc_rc_utils_ctxt,
                                yuv_buf_props_t *ps_input_buf);

#endif
