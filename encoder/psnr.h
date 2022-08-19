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
 *****************************************************************************/

#ifndef __PSNR_H__
#define __PSNR_H__

#define MAX_PSNR 100.0

void get_sse(UWORD8 *pu1_src_luma, UWORD8 *pu1_rec_luma, UWORD8 *pu1_src_chroma,
             UWORD8 *pu1_rec_chroma, WORD32 src_strd, WORD32 rec_strd, WORD32 width, WORD32 height,
             DOUBLE *pd_sse);

DOUBLE sse_to_psnr(DOUBLE samples, DOUBLE sse);

#endif
