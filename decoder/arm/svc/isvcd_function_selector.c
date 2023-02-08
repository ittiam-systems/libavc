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
*  Contains functions to initialize function pointers used in hevc
*
* @author
*  Kishore
*
* @par List of Functions:
*  - isvcd_init_function_ptr()
*
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

#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"
#include "ih264_defs.h"
#include "ih264_size_defs.h"
#include "ih264_error.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "ih264d_structs.h"
#include "ih264d_function_selector.h"
#include "isvcd_structs.h"
#include "isvcd_function_selector.h"

/**
*******************************************************************************
*
* @brief Initialize the intra/inter/transform/deblk function pointers of
* codec context
*
* @par Description: the current routine initializes the function pointers of
* codec context basing arm v8/x86/a9q architecture
*
* @param[in] ps_codec
*  ps_codec context pointer
*
* @returns  none
*
* @remarks none
*
*******************************************************************************
*/
void isvcd_init_function_ptr(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_codec = &ps_svc_lyr_dec->s_dec;
    IVD_ARCH_T e_proc_arch = ps_codec->e_processor_arch;
    isvcd_init_function_ptr_generic(ps_svc_lyr_dec);
    switch(e_proc_arch)
    {
#if defined(ARMV8)
        case ARCH_ARMV8_GENERIC:
        default:
            ih264d_init_function_ptr_av8(ps_codec);
            isvcd_init_function_ptr_neonintr(ps_svc_lyr_dec);
            break;
#elif !defined(DISABLE_NEON)
        case ARCH_ARM_A5:
        case ARCH_ARM_A7:
        case ARCH_ARM_A9:
        case ARCH_ARM_A15:
        case ARCH_ARM_A9Q:
        default:
            ih264d_init_function_ptr_a9q(ps_codec);
            isvcd_init_function_ptr_neonintr(ps_svc_lyr_dec);
            break;
#else
        default:
#endif
        case ARCH_ARM_NONEON:
            break;
    }
}
