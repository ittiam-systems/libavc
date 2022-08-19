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

/* System Include Files */
#include <math.h>

#include "ih264_macros.h"
#include "ih264_typedefs.h"
#include "psnr.h"

void get_sse(UWORD8 *pu1_src_luma, UWORD8 *pu1_rec_luma, UWORD8 *pu1_src_chroma,
             UWORD8 *pu1_rec_chroma, WORD32 src_strd, WORD32 rec_strd, WORD32 width, WORD32 height,
             DOUBLE pd_sse[3])
{
    WORD32 i, j;

    for(j = 0; j < height; j++)
    {
        for(i = 0; i < width; i++)
        {
            WORD32 diff = pu1_src_luma[i] - pu1_rec_luma[i];
            pd_sse[0] += diff * diff;
        }
        pu1_src_luma += src_strd;
        pu1_rec_luma += rec_strd;
    }

    for(j = 0; j < height / 2; j++)
    {
        for(i = 0; i < width / 2; i++)
        {
            WORD32 diff = pu1_src_chroma[i * 2] - pu1_rec_chroma[i * 2];
            pd_sse[1] += diff * diff;
            diff = pu1_src_chroma[i * 2 + 1] - pu1_rec_chroma[i * 2 + 1];
            pd_sse[2] += diff * diff;
        }
        pu1_src_chroma += src_strd;
        pu1_rec_chroma += rec_strd;
    }
}

DOUBLE sse_to_psnr(DOUBLE samples, DOUBLE sse)
{
    DOUBLE psnr;
    if(samples <= 0) return -1;
    if (sse<=0) return MAX_PSNR;
    psnr = 10.0 * (log10(samples) + 2*log10(255) - log10(sse));
    psnr = MIN(MAX_PSNR, psnr);
    return psnr;
}
