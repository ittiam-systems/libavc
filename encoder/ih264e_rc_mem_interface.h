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
******************************************************************************
* @file
*  ih264e_rc_mem_interface.h
*
* @brief
*  Get memory requirements of rate control library
*
* @author
*  ittiam
*
******************************************************************************
*/

#ifndef _IH264E_RC_MEM_INTERFACE_H_
#define _IH264E_RC_MEM_INTERFACE_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/
WORD32 ih264e_get_rate_control_mem_tab(void *pv_rate_control,
                                       iv_mem_rec_t *ps_mem,
                                       ITT_FUNC_TYPE_E e_func_type);


#endif /* _IH264E_RC_MEM_INTERFACE_H_ */

