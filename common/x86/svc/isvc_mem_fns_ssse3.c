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
 *  isvc_mem_fns_atom_intr.c
 *
 * @brief
 *  Functions used for memory operations
 *
 * @author
 *  Ittiam
 *
 * @par List of Functions:
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
#include <assert.h>

#include "ih264_typedefs.h"
#include "isvc_mem_fns.h"

#include <immintrin.h>

/**
********************************************************************************
*  @brief  copies a 2d blk from one location to another
*
*  @param[out] pu1_dst : dst pointer
*
*  @param[in] i4_dst_stride: stride of destination
*
*  @param[in] pu1_src : src ptr
*
*  @param[in] i4_src_stride: stride of src
*
*  @param[in] i4_blk_wd : blk width
*
*  @param[in] i4_blk_ht : blk height
*
*  @return void
********************************************************************************
*/
void isvc_copy_2d_ssse3(UWORD8 *pu1_dst, WORD32 i4_dst_stride, UWORD8 *pu1_src,
                        WORD32 i4_src_stride, WORD32 i4_blk_wd, WORD32 i4_blk_ht)
{
    WORD32 i, j;
    /* all 128 bit registers are named with a suffix mxnb, where m is the */
    /* number of n bits packed in the register                            */

    if(((i4_blk_wd % 4) != 0) || ((i4_blk_ht % 4) != 0))
    {
        isvc_copy_2d(pu1_dst, i4_dst_stride, pu1_src, i4_src_stride, i4_blk_wd, i4_blk_ht);

        return;
    }

    if(0 == (i4_blk_wd & 31)) /* wd multiple of 32 case */
    {
        __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;
        __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;

        if(0 == (i4_blk_ht & 7)) /* ht multiple of 8 case */
        {
            __m128i src8_16x8b, src9_16x8b, src10_16x8b, src11_16x8b;
            __m128i src12_16x8b, src13_16x8b, src14_16x8b, src15_16x8b;

            for(i = 0; i < i4_blk_ht; i += 8)
            {
                for(j = 0; j < i4_blk_wd; j += 32)
                {
                    src0_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src));  // i = 0
                    src1_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3
                    src4_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 4 * i4_src_stride));  // i = 4
                    src5_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 5 * i4_src_stride));  // i = 5
                    src6_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 6 * i4_src_stride));  // i = 6
                    src7_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 7 * i4_src_stride));  // i = 7
                    /* Add 16 as offset */
                    src8_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 16));  // i = 0
                    src9_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + i4_src_stride + 16));  // i = 1
                    src10_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride + 16));  // i = 2
                    src11_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride + 16));  // i = 3
                    src12_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 4 * i4_src_stride + 16));  // i = 4
                    src13_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 5 * i4_src_stride + 16));  // i = 5
                    src14_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 6 * i4_src_stride + 16));  // i = 6
                    src15_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 7 * i4_src_stride + 16));  // i = 7

                    _mm_storeu_si128((__m128i *) (pu1_dst), src0_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + i4_dst_stride), src1_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 4 * i4_dst_stride), src4_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 5 * i4_dst_stride), src5_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 6 * i4_dst_stride), src6_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 7 * i4_dst_stride), src7_16x8b);

                    _mm_storeu_si128((__m128i *) (pu1_dst + 16), src8_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + i4_dst_stride + 16), src9_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride + 16), src10_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride + 16), src11_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 4 * i4_dst_stride + 16), src12_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 5 * i4_dst_stride + 16), src13_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 6 * i4_dst_stride + 16), src14_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 7 * i4_dst_stride + 16), src15_16x8b);

                    pu1_src += 32;
                    pu1_dst += 32;
                }

                pu1_src = pu1_src - i4_blk_wd + 8 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 8 * i4_dst_stride;
            }
        }
        else /* ht multiple of 4 case */
        {
            for(i = 0; i < i4_blk_ht; i += 4)
            {
                for(j = 0; j < i4_blk_wd; j += 32)
                {
                    src0_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src));  // i = 0
                    src1_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3
                    /* Add 16 as offset */
                    src4_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 16));  // i = 0
                    src5_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + i4_src_stride + 16));  // i = 1
                    src6_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride + 16));  // i = 2
                    src7_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride + 16));  // i = 3

                    _mm_storeu_si128((__m128i *) (pu1_dst), src0_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + i4_dst_stride), src1_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 16), src4_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + i4_dst_stride + 16), src5_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride + 16), src6_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride + 16), src7_16x8b);

                    pu1_src += 32;
                    pu1_dst += 32;
                }

                pu1_src = pu1_src - i4_blk_wd + 4 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 4 * i4_dst_stride;
            }
        }
    }
    else if(0 == (i4_blk_wd & 15)) /* wd multiple of 16 case */
    {
        __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;

        if(0 == (i4_blk_ht & 7)) /* ht multiple of 8 case */
        {
            __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;

            for(i = 0; i < i4_blk_ht; i += 8)
            {
                for(j = 0; j < i4_blk_wd; j += 16)
                {
                    src0_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3
                    src4_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 4 * i4_src_stride));  // i = 4
                    src5_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 5 * i4_src_stride));  // i = 5
                    src6_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 6 * i4_src_stride));  // i = 6
                    src7_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 7 * i4_src_stride));  // i = 7

                    _mm_storeu_si128((__m128i *) (pu1_dst + 0 * i4_dst_stride), src0_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 1 * i4_dst_stride), src1_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 4 * i4_dst_stride), src4_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 5 * i4_dst_stride), src5_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 6 * i4_dst_stride), src6_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 7 * i4_dst_stride), src7_16x8b);

                    pu1_src += 16;
                    pu1_dst += 16;
                }

                pu1_src = pu1_src - i4_blk_wd + 8 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 8 * i4_dst_stride;
            }
        }
        else /* ht multiple of 4 case */
        {
            for(i = 0; i < i4_blk_ht; i += 4)
            {
                for(j = 0; j < i4_blk_wd; j += 16)
                {
                    src0_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadu_si128((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3

                    _mm_storeu_si128((__m128i *) (pu1_dst + 0 * i4_dst_stride), src0_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 1 * i4_dst_stride), src1_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storeu_si128((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);

                    pu1_src += 16;
                    pu1_dst += 16;
                }

                pu1_src = pu1_src - i4_blk_wd + 4 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 4 * i4_dst_stride;
            }
        }
    }
    else if(0 == (i4_blk_wd & 7)) /* wd multiple of 8 case */
    {
        __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;

        if(0 == (i4_blk_ht & 7)) /* ht multiple of 8 case */
        {
            __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;

            for(i = 0; i < i4_blk_ht; i += 8)
            {
                for(j = 0; j < i4_blk_wd; j += 8)
                {
                    src0_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3
                    src4_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 4 * i4_src_stride));  // i = 4
                    src5_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 5 * i4_src_stride));  // i = 5
                    src6_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 6 * i4_src_stride));  // i = 6
                    src7_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 7 * i4_src_stride));  // i = 7

                    _mm_storel_epi64((__m128i *) (pu1_dst + 0 * i4_dst_stride), src0_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 1 * i4_dst_stride), src1_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 4 * i4_dst_stride), src4_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 5 * i4_dst_stride), src5_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 6 * i4_dst_stride), src6_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 7 * i4_dst_stride), src7_16x8b);

                    pu1_src += 8;
                    pu1_dst += 8;
                }

                pu1_src = pu1_src - i4_blk_wd + 8 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 8 * i4_dst_stride;
            }
        }
        else /* ht multiple of 4 case */
        {
            for(i = 0; i < i4_blk_ht; i += 4)
            {
                for(j = 0; j < i4_blk_wd; j += 8)
                {
                    src0_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3

                    _mm_storel_epi64((__m128i *) (pu1_dst + 0 * i4_dst_stride), src0_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 1 * i4_dst_stride), src1_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 2 * i4_dst_stride), src2_16x8b);
                    _mm_storel_epi64((__m128i *) (pu1_dst + 3 * i4_dst_stride), src3_16x8b);

                    pu1_src += 8;
                    pu1_dst += 8;
                }

                pu1_src = pu1_src - i4_blk_wd + 4 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 4 * i4_dst_stride;
            }
        }
    }
    else /* wd multiple of 4 case */
    {
        __m128i src0_16x8b, src1_16x8b, src2_16x8b, src3_16x8b;
        WORD32 src0, src1, src2, src3;
        if(0 == (i4_blk_ht & 7)) /* ht multiple of 8 case */
        {
            __m128i src4_16x8b, src5_16x8b, src6_16x8b, src7_16x8b;
            WORD32 src4, src5, src6, src7;

            for(i = 0; i < i4_blk_ht; i += 8)
            {
                for(j = 0; j < i4_blk_wd; j += 4)
                {
                    src0_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3
                    src4_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 4 * i4_src_stride));  // i = 4
                    src5_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 5 * i4_src_stride));  // i = 5
                    src6_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 6 * i4_src_stride));  // i = 6
                    src7_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 7 * i4_src_stride));  // i = 7

                    src0 = _mm_cvtsi128_si32(src0_16x8b);
                    src1 = _mm_cvtsi128_si32(src1_16x8b);
                    src2 = _mm_cvtsi128_si32(src2_16x8b);
                    src3 = _mm_cvtsi128_si32(src3_16x8b);
                    src4 = _mm_cvtsi128_si32(src4_16x8b);
                    src5 = _mm_cvtsi128_si32(src5_16x8b);
                    src6 = _mm_cvtsi128_si32(src6_16x8b);
                    src7 = _mm_cvtsi128_si32(src7_16x8b);

                    *(WORD32 *) (&pu1_dst[0 * i4_dst_stride]) = src0;
                    *(WORD32 *) (&pu1_dst[1 * i4_dst_stride]) = src1;
                    *(WORD32 *) (&pu1_dst[2 * i4_dst_stride]) = src2;
                    *(WORD32 *) (&pu1_dst[3 * i4_dst_stride]) = src3;
                    *(WORD32 *) (&pu1_dst[4 * i4_dst_stride]) = src4;
                    *(WORD32 *) (&pu1_dst[5 * i4_dst_stride]) = src5;
                    *(WORD32 *) (&pu1_dst[6 * i4_dst_stride]) = src6;
                    *(WORD32 *) (&pu1_dst[7 * i4_dst_stride]) = src7;

                    pu1_src += 4;
                    pu1_dst += 4;
                }

                pu1_src = pu1_src - i4_blk_wd + 8 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 8 * i4_dst_stride;
            }
        }
        else /* ht multiple of 4 case */
        {
            for(i = 0; i < i4_blk_ht; i += 4)
            {
                for(j = 0; j < i4_blk_wd; j += 4)
                {
                    src0_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 0 * i4_src_stride));  // i = 0
                    src1_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 1 * i4_src_stride));  // i = 1
                    src2_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 2 * i4_src_stride));  // i = 2
                    src3_16x8b =
                        _mm_loadl_epi64((__m128i *) (pu1_src + 3 * i4_src_stride));  // i = 3

                    src0 = _mm_cvtsi128_si32(src0_16x8b);
                    src1 = _mm_cvtsi128_si32(src1_16x8b);
                    src2 = _mm_cvtsi128_si32(src2_16x8b);
                    src3 = _mm_cvtsi128_si32(src3_16x8b);

                    *(WORD32 *) (&pu1_dst[0 * i4_dst_stride]) = src0;
                    *(WORD32 *) (&pu1_dst[1 * i4_dst_stride]) = src1;
                    *(WORD32 *) (&pu1_dst[2 * i4_dst_stride]) = src2;
                    *(WORD32 *) (&pu1_dst[3 * i4_dst_stride]) = src3;

                    pu1_src += 4;
                    pu1_dst += 4;
                }

                pu1_src = pu1_src - i4_blk_wd + 4 * i4_src_stride;
                pu1_dst = pu1_dst - i4_blk_wd + 4 * i4_dst_stride;
            }
        }
    }
}
