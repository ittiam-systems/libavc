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
*  isvcd_function_selector.c
*
* @brief
*  Contains functions to initialize function pointers used in svc
*
* @author
*  Kishore
*
* @par List of Functions:
* @remarks
*  None
*
*******************************************************************************
*/
/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "isvcd_structs.h"
#include "ih264d_function_selector.h"
#include "isvcd_function_selector.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the intra/inter/transform/deblk function pointers of
 * codec context
 *
 * @par Description: the current routine initializes the function pointers of
 * codec context basing on the architecture in use for svc
 *
 * @param[in] ps_svc_lyr_dec
 *  svc dec layer context pointer
 *
 * @returns  none
 *
 * @remarks none
 *
 *******************************************************************************
 */
void isvcd_init_function_ptr(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    isvcd_init_function_ptr_generic(ps_svc_lyr_dec);
    switch(ps_svc_lyr_dec->s_dec.e_processor_arch)
    {
        case ARCH_X86_GENERIC:
            isvcd_init_function_ptr_generic(ps_svc_lyr_dec);
            break;
        case ARCH_X86_SSSE3:
            ih264d_init_function_ptr_ssse3(&ps_svc_lyr_dec->s_dec);
            break;
        case ARCH_X86_SSE42:
        default:
            ih264d_init_function_ptr_ssse3(&ps_svc_lyr_dec->s_dec);
            isvcd_init_function_ptr_sse42(ps_svc_lyr_dec);
            break;
    }
}
