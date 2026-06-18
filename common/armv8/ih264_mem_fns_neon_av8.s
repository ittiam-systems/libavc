//******************************************************************************
//*
//* Copyright (C) 2015 The Android Open Source Project
//*
//* Licensed under the Apache License, Version 2.0 (the "License");
//* you may not use this file except in compliance with the License.
//* You may obtain a copy of the License at:
//*
//* http://www.apache.org/licenses/LICENSE-2.0
//*
//* Unless required by applicable law or agreed to in writing, software
//* distributed under the License is distributed on an "AS IS" BASIS,
//* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//* See the License for the specific language governing permissions and
//* limitations under the License.
//*
//*****************************************************************************
//* Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
//*/
///**
// *******************************************************************************
// * @file
// *  ih264_mem_fns_neon.s
// *
// * @brief
// *  Contains function definitions for memory manipulation
// *
// * @author
// *     Naveen SR
// *
// * @par List of Functions:
// *  - ih264_memset_16bit_mul_8_av8()
// *  - ih264_memset_16bit_av8()
// *
// * @remarks
// *  None
// *
// *******************************************************************************
//*/

.include "ih264_neon_macros.s"
.text
//void ih264_memset_16bit_mul_8(UWORD16 *pu2_dst,
//                                      UWORD16 value,
//                                      UWORD32 num_words)
//**************Variables Vs Registers*************************
//    x0 => *pu2_dst
//    w1 => value
//    w2 => num_words


    .global ih264_memset_16bit_mul_8_av8

ENTRY ih264_memset_16bit_mul_8_av8

// Assumptions: num_words is either 8, 16 or 32

    // Memset 8 words
    dup       v0.4h, w1
loop_memset_16bit_mul_8:
    st1       {v0.4h}, [x0], #8
    st1       {v0.4h}, [x0], #8

    subs      w2, w2, #8
    bne       loop_memset_16bit_mul_8

    EXIT_FUNC
    ret



//void ih264_memset_16bit(UWORD16 *pu2_dst,
//                       UWORD16 value,
//                       UWORD32 num_words)
//**************Variables Vs Registers*************************
//    x0 => *pu2_dst
//    w1 => value
//    w2 => num_words



    .global ih264_memset_16bit_av8

ENTRY ih264_memset_16bit_av8
    subs      w2, w2, #8
    blt       arm_memset_16bit
    dup       v0.4h, w1
loop_neon_memset_16bit:
    // Memset 8 words
    st1       {v0.4h}, [x0], #8
    st1       {v0.4h}, [x0], #8

    subs      w2, w2, #8
    bge       loop_neon_memset_16bit
    cmn       w2, #8
    beq       end_func3

arm_memset_16bit:
    add       w2, w2, #8

loop_arm_memset_16bit:
    strh      w1, [x0], #2
    subs      w2, w2, #1
    bne       loop_arm_memset_16bit
    EXIT_FUNC
    ret

end_func3:
    EXIT_FUNC
    ret



