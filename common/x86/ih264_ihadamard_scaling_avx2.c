/******************************************************************************
+ *
+ * Copyright (C) 2015 The Android Open Source Project
+ *
+ * Licensed under the Apache License, Version 2.0 (the "License");
+ * you may not use this file except in compliance with the License.
+ * You may obtain a copy of the License at:
+ *
+ * http://www.apache.org/licenses/LICENSE-2.0
+ *
+ * Unless required by applicable law or agreed to in writing, software
+ * distributed under the License is distributed on an "AS IS" BASIS,
+ * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
+ * See the License for the specific language governing permissions and
+ * limitations under the License.
+ *
+ *****************************************************************************
+ * Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
+*/
/**
+ *******************************************************************************
+ * @file
+ *  ih264_ihadamard_scaling_avx2.c
+ *
+ * @brief
+ *  Contains definition of functions for h264 inverse hadamard 4x4 transform and scaling
+ *
+ * @author
+ *  Priyanka
+ *
+ *  @par List of Functions:
+ *  - ih264_ihadamard_scaling_4x4_avx2()
+ *
+ * @remarks
+ *
+ *******************************************************************************
+ */
/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "ih264_structs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include <immintrin.h>

/*
+ ********************************************************************************
+ *
+ * @brief This function performs a 4x4 inverse hadamard transform on the 4x4 DC coefficients
+ * of a 16x16 intra prediction macroblock, and then performs scaling.
+ * prediction buffer
+ *
+ * @par Description:
+ *  The DC coefficients pass through a 2-stage inverse hadamard transform.
+ *  This inverse transformed content is scaled to based on Qp value.
+ *
+ * @param[in] pi2_src
+ *  input 4x4 block of DC coefficients
+ *
+ * @param[out] pi2_out
+ *  output 4x4 block
+ *
+ * @param[in] pu2_iscal_mat
+ *  pointer to scaling list
+ *
+ * @param[in] pu2_weigh_mat
+ *  pointer to weight matrix
+ *
+ * @param[in] u4_qp_div_6
+ *  Floor (qp/6)
+ *
+ * @param[in] pi4_tmp
+ * temporary buffer of size 1*16
+ *
+ * @returns none
+ *
+ * @remarks none
+ *
+ *******************************************************************************
*/

#include <stdint.h>
#include <string.h>

#include <stdio.h>

#ifdef __ANDROID__
#include "log/log.h"
#include <cutils/log.h>
#endif


void ih264_ihadamard_scaling_4x4_avx2(WORD16* pi2_src,
                                      WORD16* pi2_out,
                                      const UWORD16 *pu2_iscal_mat,
                                      const UWORD16 *pu2_weigh_mat,
                                      UWORD32 u4_qp_div_6,
                                      WORD32* pi4_tmp)
{ 
    __m256i src,r0_r1,r2_r3,r3_r2,r1_r3,r0_r2;
    __m256i src_r0_r1, src_r2_r3;
    __m256i temp0, temp1,tmp0, tmp1, tmp2, tmp3;
    __m256i add_rshift = _mm256_set1_epi32((u4_qp_div_6 < 6) ? (1 << (5 - u4_qp_div_6)) : 0);
    __m256i mult_val = _mm256_set1_epi32(pu2_iscal_mat[0] * pu2_weigh_mat[0]);
    __m256i zero =  _mm256_setzero_si256();

    __m128i t0 ,t1;
    UNUSED (pi4_tmp);

    src_r0_r1 = _mm256_loadu_si256((__m256i *) (pi2_src)); //a00 a01 a02 a03 a10 a11 a12 a13 -- the source matrix 0th,1st row

    temp0 = _mm256_unpacklo_epi64(src_r0_r1, zero);  
    temp1 = _mm256_unpackhi_epi64(src_r0_r1, zero);        // b0 b1 b2..         d0 d1...
    temp0 = _mm256_unpacklo_epi16(temp0, temp1); 
    tmp0 =  _mm256_permute2x128_si256(temp0,zero,0x20);    //tmp0 tmp3
    tmp1 =  _mm256_permute2x128_si256(temp0,zero,0x31);    //tmp1  tmp2

    temp0 = _mm256_unpacklo_epi32(tmp0, tmp1);             //a0 c0 a1 c1 a2 c2 a3 c3 a0 c0 a1 c1 b0 d0 b1 c1
    temp1 = _mm256_unpackhi_epi32(tmp0, tmp1);
    
    temp1 = _mm256_shuffle_epi32(temp1,0b01001110);
    tmp0 = _mm256_add_epi16(temp0, temp1);
    tmp1 = _mm256_sub_epi16(temp0, temp1);

    temp0 =   _mm256_unpacklo_epi64(tmp0, tmp1);
    temp1 =   _mm256_unpackhi_epi64(tmp0, tmp1);
    tmp0 = _mm256_add_epi16(temp0, temp1);
    tmp1 = _mm256_sub_epi16(temp0, temp1);

    temp0 = _mm256_unpacklo_epi32(tmp0, tmp1);             //a0 c0 a1 c1 a2 c2 a3 c3 a0 c0 a1 c1 b0 d0 b1 c1
    temp1 = _mm256_unpackhi_epi32(tmp0, tmp1); 
    
    temp0 = _mm256_unpacklo_epi64(tmp0, tmp1);             //a0 c0 a1 c1 a2 c2 a3 c3 a0 c0 a1 c1 b0 d0 b1 c1
    temp1 = _mm256_unpackhi_epi64(tmp0, tmp1);
   
    tmp0 = _mm256_unpacklo_epi16(temp0, temp1);             //a0 c0 a1 c1 a2 c2 a3 c3 a0 c0 a1 c1 b0 d0 b1 c1
    tmp1 = _mm256_unpackhi_epi16(temp0, temp1);

    temp0 = _mm256_unpacklo_epi32(tmp0, tmp1);             //a0 c0 a1 c1 a2 c2 a3 c3 a0 c0 a1 c1 b0 d0 b1 c1
    temp1 = _mm256_unpackhi_epi32(tmp0, tmp1);
    
    temp1 = _mm256_shuffle_epi32(temp1, _MM_SHUFFLE(1, 0, 3, 2));
    tmp0 = _mm256_add_epi16(temp0, temp1);
    tmp1 = _mm256_sub_epi16(temp0, temp1);
    temp0 =   _mm256_unpacklo_epi64(tmp0, tmp1);
    temp1 =   _mm256_unpackhi_epi64(tmp0, tmp1);
    tmp0 = _mm256_add_epi16(temp0, temp1);
    tmp1 = _mm256_sub_epi16(temp0, temp1);

    temp0 =   _mm256_unpacklo_epi64(tmp0, tmp1);
    temp1 =   _mm256_unpackhi_epi64(tmp0, tmp1);


    r0_r1 =_mm256_cvtepi16_epi32(_mm256_castsi256_si128(temp0));
    r2_r3 = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(temp1));

    src_r0_r1 = _mm256_mullo_epi32(r0_r1, mult_val);
    src_r2_r3 = _mm256_mullo_epi32(r2_r3, mult_val);

    //Scaling
    if(u4_qp_div_6 >= 6)
    {
        src_r0_r1 = _mm256_slli_epi32(src_r0_r1, u4_qp_div_6 - 6);
        src_r2_r3 = _mm256_slli_epi32(src_r2_r3, u4_qp_div_6 - 6);
    }
    else
    {
        temp0  = _mm256_add_epi32(src_r0_r1, add_rshift);
        temp1  = _mm256_add_epi32(src_r2_r3, add_rshift);
        src_r0_r1 = _mm256_srai_epi32(temp0, 6 - u4_qp_div_6);
        src_r2_r3 = _mm256_srai_epi32(temp1, 6 - u4_qp_div_6);
    }

    src = _mm256_packs_epi32(src_r0_r1, src_r2_r3);
    _mm256_storeu_si256((__m256i *) (&pi2_out[0]), src);
}

