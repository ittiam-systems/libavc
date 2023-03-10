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
*  isvce_encode.h
*
* @brief
*  Contains functions for encode API
*
*******************************************************************************
*/

#ifndef _ISVCE_ENCODE_H_
#define _ISVCE_ENCODE_H_

#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"

extern WORD32 isvce_encode(iv_obj_t *ps_codec_obj, void *pv_api_ip, void *pv_api_op);

#endif
