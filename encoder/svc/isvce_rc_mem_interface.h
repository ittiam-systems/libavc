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
******************************************************************************
* @file
*  isvce_rc_mem_interface.h
*
* @brief
*  This file contains function declaration and structures for rate control
*  memtabs
*
* @author
*  ittiam
*
* @remarks
*  The rate control library is a global library across various codecs. It
*  anticipates certain structures definitions. Those definitions are to be
*  imported from global workspace. Instead of that, the structures needed for
*  rc library are copied in to this file and exported to rc library. If the
*  structures / enums / ... in the global workspace change, this file also needs
*  to be modified accordingly.
*
******************************************************************************
*/
#ifndef _ISVCE_RC_MEM_INTERFACE_H_
#define _ISVCE_RC_MEM_INTERFACE_H_

#include "ih264e_rc_mem_interface.h"

/**
 ***************************************************************************
 * Enum to hold mem records in RC
 ****************************************************************************
 */
typedef enum RC_MEM_TYPES_T
{
    RC_MEM_FRAME_TIME,

    RC_MEM_TIME_STAMP,

    RC_MEM_FRAME_RATE,

    RC_MEM_API_L0,

    RC_MEM_API_L1,

    RC_MEM_API_L2,

    RC_MEM_CNT

    /*
     * Do not add anything below
     */
} RC_MEM_TYPES_T;

extern WORD32 isvce_get_rate_control_mem_tab(void *pv_rate_control, iv_mem_rec_t *ps_mem,
                                             ITT_FUNC_TYPE_E e_func_type);

#endif
