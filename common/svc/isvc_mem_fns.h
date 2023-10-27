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
*  isvc_mem_fns.h
*
* @brief
*  Function declarations used for memory functions
*
* @author
*  Ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/
#ifndef _ISVC_MEM_FNS_H_
#define _ISVC_MEM_FNS_H_

#include "ih264_typedefs.h"

typedef void *FT_MEM_ALLOC(UWORD32 u4_size);

typedef void FT_MEM_FREE(void *pv_mem);

typedef void FT_MEMCPY(UWORD8 *pu1_dst, UWORD8 *pu1_src, UWORD32 num_bytes);

typedef void FT_COPY_2D(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 *pu1_src,
                        WORD32 i4_src_stride, WORD32 i4_blk_wd, WORD32 i4_blk_ht);

typedef void FT_MEMSET_2D(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 u1_val, WORD32 i4_blk_wd,
                          WORD32 i4_blk_ht);

typedef void FT_MEMSET(UWORD8 *pu1_dst, UWORD8 value, UWORD32 num_bytes);

typedef void FT_MEMSET_16BIT(UWORD16 *pu2_dst, UWORD16 value, UWORD32 num_words);

typedef void FT_16BIT_INTERLEAVED_COPY(WORD16 *pi2_src, WORD16 *pi2_dst, WORD32 src_strd,
                                       WORD32 dst_strd, WORD32 ht, WORD32 wd);

typedef void FT_16BIT_INTERLEAVED_MEMSET(WORD16 *pi2_src, WORD32 i4_src_strd, WORD16 i2_value,
                                         WORD32 i4_wd, WORD32 i4_ht);

typedef UWORD8 FT_NONZERO_CHECKER(UWORD8 *pu1_data, WORD32 i4_data_strd, UWORD32 u4_wd,
                                  UWORD32 u4_ht);

/* C function declarations */
extern FT_MEMCPY ih264_memcpy_mul_8;
extern FT_MEMSET ih264_memset_mul_8;
extern FT_MEMSET_16BIT ih264_memset_16bit;
extern FT_MEMSET_16BIT ih264_memset_16bit_mul_8;
extern FT_COPY_2D isvc_copy_2d;
extern FT_MEMSET_2D isvc_memset_2d;
extern FT_NONZERO_CHECKER isvc_is_nonzero_blk;
extern FT_MEM_ALLOC isvc_memory_alloc;
extern FT_MEM_FREE isvc_memory_free;

/* A9 Q function declarations */
extern FT_MEMCPY ih264_memcpy_mul_8_a9q;
extern FT_MEMSET ih264_memset_mul_8_a9q;
extern FT_MEMSET_16BIT ih264_memset_16bit_a9q;
extern FT_MEMSET_16BIT ih264_memset_16bit_mul_8_a9q;

/* AV8 function declarations */
extern FT_MEMCPY ih264_memcpy_mul_8_av8;
extern FT_MEMSET ih264_memset_mul_8_av8;
extern FT_MEMSET_16BIT ih264_memset_16bit_av8;
extern FT_MEMSET_16BIT ih264_memset_16bit_mul_8_av8;

/* NEON function declarations */
extern FT_MEMSET_2D isvc_memset_2d_neon;

/* SSSE3 variants */
extern FT_MEMCPY ih264_memcpy_mul_8_ssse3;
extern FT_MEMSET ih264_memset_mul_8_ssse3;
extern FT_MEMSET_16BIT ih264_memset_16bit_mul_8_ssse3;
extern FT_COPY_2D isvc_copy_2d_ssse3;

/* SSE4.2 variants */
extern FT_MEMSET_2D isvc_memset_2d_sse42;

#endif
