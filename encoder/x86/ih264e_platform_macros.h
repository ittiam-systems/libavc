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
*  ih264e_platform_macros.h
*
* @brief
*  Contains platform specific routines used for codec context intialization
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_PLATFORM_MACROS_H_
#define _IH264E_PLATFORM_MACROS_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_init_function_ptr_generic(codec_t *ps_codec);
void ih264e_init_function_ptr_ssse3(codec_t *ps_codec);
void ih264e_init_function_ptr_sse42(codec_t *ps_codec);
void ih264e_init_function_ptr(void *pv_codec);
IV_ARCH_T ih264e_default_arch(void);

#endif /* _IH264E_PLATFORM_MACROS_H_ */
