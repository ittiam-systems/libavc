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
 *  isvcd_function_selector.h
 *
 * @brief
 *  Structure definitions used in the decoder
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_FUNCTION_SELECTOR_H_
#define _ISVCD_FUNCTION_SELECTOR_H_

void isvcd_init_function_ptr(svc_dec_lyr_struct_t *ps_codec);

void isvcd_init_function_ptr_generic(svc_dec_lyr_struct_t *ps_codec);

void isvcd_init_function_ptr_neonintr(svc_dec_lyr_struct_t *ps_codec);

void isvcd_init_function_ptr_sse42(svc_dec_lyr_struct_t *ps_codec);
#endif /* _ISVCD_FUNCTION_SELECTOR_H_ */
