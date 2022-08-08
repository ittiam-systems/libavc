/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

/*****************************************************************************/
/*                                                                           */
/*  File Name         : imvcd_slice_functions.h                              */
/*                                                                           */
/*  Description       : Functions for MVC Slice parsing, etc.                */
/*                                                                           */
/*****************************************************************************/

#ifndef _IMVCD_SLICE_FUNCTIONS_H_
#define _IMVCD_SLICE_FUNCTIONS_H_

#include "ih264_typedefs.h"
#include "imvcd_structs.h"

extern WORD32 imvcd_parse_decode_slice(mvc_dec_ctxt_t *ps_mvcd_ctxt);
#endif
