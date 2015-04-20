@/******************************************************************************
@ *
@ * Copyright (C) 2015 The Android Open Source Project
@ *
@ * Licensed under the Apache License, Version 2.0 (the "License");
@ * you may not use this file except in compliance with the License.
@ * You may obtain a copy of the License at:
@ *
@ * http://www.apache.org/licenses/LICENSE-2.0
@ *
@ * Unless required by applicable law or agreed to in writing, software
@ * distributed under the License is distributed on an "AS IS" BASIS,
@ * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
@ * See the License for the specific language governing permissions and
@ * limitations under the License.
@ *
@ *****************************************************************************
@ * Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
@*/
@**
@ *******************************************************************************
@ * @file
@ *  ih264_itrans_recon_neon_a9.s
@ *
@ * @brief
@ *  Contains function definitions for single stage  inverse transform
@ *
@ *
@ * @par List of Functions:
@ *  - ih264_itrans_recon_4x4_a9()
@ *
@ * @remarks
@ *  None
@ *
@ *******************************************************************************
@*
@**
@ *******************************************************************************
@ *
@ * @brief
@ *  This function performs Inverse transform type Ci4 for 4*4 block
@ *
@ * @par Description:
@ *  Performs inverse transform Ci4 and adds the residue to get the
@ *  reconstructed block
@ *
@ * @param[in] pi16_levelBlock
@ *  Input 4x4 coefficients
@ *
@ * @param[in] puc_predBuffer
@ *  Prediction 4x4 block
@ *
@ * @param[out] puc_reconPic
@ *  Output 4x4 block
@ *
@ * @param[in] ui16_picWidth
@ *  Input stride
@ *
@ * @param[in] pred_strd
@ *  Prediction stride
@ *
@ * @param[in] dst_strd
@ *  Output Stride
@ *
@ * @param[in] zero_cols
@ *  Zero columns in pi2_src
@ *
@ * @returns  Void
@ *
@ * @remarks
@ *  None
@ *
@ *
@ *******************************************************************************
@ *
@void ih264_itrans_recon_4x4(
@       WORD16 *pi2_src,
@       UWORD8 *pu1_pred,
@       UWORD8 *pu1_recon,
@       WORD32 src_strd,
@       WORD32 pred_strd,
@       WORD32 dst_strd,
@       UWORD32 q_lev,          //quantizer level
@       WORD32 *pi4_tmp)
@**************Variables Vs Registers*****************************************
@r0 => *pi2_src
@r1 => *pu1_pred
@r2 => *pu1_recon
@r3 =>  src_strd
@r4 =>  pred_strd
@r5 =>  dst_strd
@r6 =>  q_lev
@r7 =>  *pi4_tmp

.text
.p2align 2


    .global ih264_itrans_recon_4x4_a9

ih264_itrans_recon_4x4_a9:
    stmfd         sp!, {r4-r12, r14}    @stack stores the values of the arguments
    lsl           r3, r3, #1

    vld1.16       d0, [r0], r3          @0th row pi2_src_tmp[0]
    ldr           r4, [sp, #40]         @Loads pred_strd

    vld1.16       d1, [r0], r3          @I row pi2_src_tmp[0]
    ldr           r5, [sp, #44]         @Loads *dst_strd

    vld1.16       d2, [r0], r3          @II row pi2_src_tmp[0]

    vld1.16       d3, [r0]              @III row pi2_src_tmp[0]
    ldr           r7, [sp, #52]         @Loads *pi4_tmp

    vpush         {d8-d15}

    vtrn.16       d0, d1                @Transpose to get all the 0th element in the single D register
    vtrn.16       d2, d3
    vtrn.32       d0, d2
    vtrn.32       d1, d3                @D0 --> pi2_src_tmp[0], D1 --> pi2_src_tmp[1]
                                        @D2 --> pi2_src_tmp[2], D3 --> pi2_src_tmp[3]

    vaddl.s16     q3, d0, d2            @x0 = (pi2_src_tmp[0] +  pi2_src_tmp[2])
    vsubl.s16     q4, d0, d2            @x1 = (pi2_src_tmp[0] -  pi2_src_tmp[2])
    vshr.s16      d4, d1, #1            @pi2_src_tmp[1] >> 1
    vshr.s16      d5, d3, #1            @pi2_src_tmp[3] >> 1

    vsubl.s16     q5, d4, d3            @x2 = D_SHIFT(pi2_src_tmp[1],1,shft) -  pi2_src_tmp[3]

    vaddl.s16     q6, d1, d5            @x3 = pi2_src_tmp[1] + D_SHIFT(pi2_src_tmp[3],1,shft)

    vadd.s32      q8, q4, q5            @x1 + x2
    vsub.s32      q9, q4, q5            @x1 - x2

    vadd.s32      q7, q3, q6            @x0 + x3
    vsub.s32      q10, q3, q6           @x0 - x3

    vtrn.32       q7, q8                @Transpose the register to have the adjacent values

    vtrn.32       q9, q10
    vadd.s32      d6, d14, d15          @x0(0,1) = (pi4_tblk[0,1]     +  pi4_tblk[8,9])

    vsub.s32      d7, d14, d15          @x1(0,1) = (pi4_tblk[0,1]     -  pi4_tblk[8,9])

    vshr.s32      d4, d16, #1           @pi4_tblk[4,5] >> 1
    vshr.s32      d5, d17, #1           @pi4_tblk[12,13] >> 1

    vsub.s32      d8, d4, d17           @x2(0,1) = D_SHIFT(pi4_tblk[4,5],1,shft) -  pi4_tblk[12,13]
    vadd.s32      d9, d16, d5           @x3(0,1) =  pi4_tblk[4,5] + D_SHIFT(pi4_tblk[12,13],1,shft)

    vadd.s32      d10, d18, d19         @x0(2,3) = (pi4_tblk[2,3]     +  pi4_tblk[10,11])
    vsub.s32      d11, d18, d19         @x1(2,3) = (pi4_tblk[2,3]     -  pi4_tblk[10,11])
    vshr.s32      d4, d20, #1           @pi4_tblk[6,7] >> 1
    vshr.s32      d5, d21, #1           @pi4_tblk[14,15] >> 1

    vld1.32       d30[0], [r1], r4      @I row Load pu1_pred buffer
    vsub.s32      d12, d4, d21          @x2(2,3) = D_SHIFT(pi4_tblk[6,7],1,shft) -  pi4_tblk[14,15]

    vmovl.u8      q15, d30              @I row Convert 8 bit pred buffer to 16 bit
    vadd.s32      d13, d20, d5          @x3(2,3) =  pi4_tblk[6,7] + D_SHIFT(pi4_tblk[14,15],1,shft)

    vadd.s32      d16, d6, d9           @I row i_macro(0,1) = x0(0,1) + x3(0,1)

    vld1.32       d28[0], [r1], r4      @II row Load pu1_pred buffer
    vadd.s32      d17, d10, d13         @I row i_macro(2,3) = x0(2,3) + x3(2,3)

    vqrshrn.s32   d16, q8, #6           @I row i_macro = D_SHIFT(i_macro,6,shft)

    vmovl.u8      q14, d28              @II row Convert 8 bit pred buffer to 16 bit
    vadd.u16      d16, d16, d30         @I row i_macro += *pu1_pred_tmp

    vqmovun.s16   d16, q8               @I row CLIP_U8(i_macro)
    vadd.s32      d18, d7, d8           @II row i_macro(0,1) = x1(0,1) + x2(0,1)

    vld1.32       d26[0], [r1], r4      @III row Load pu1_pred buffer
    vadd.s32      d19, d11, d12         @II row i_macro(2,3) = x1(2,3) + x2(2,3)

    vqrshrn.s32   d18, q9, #6           @II row i_macro = D_SHIFT(i_macro,6,shft)

    vmovl.u8      q13, d26              @III row Convert 8 bit pred buffer to 16 bit
    vadd.u16      d18, d18, d28         @II row i_macro += *pu1_pred_tmp

    vst1.32       d16[0], [r2], r5      @I row store the value
    vsub.s32      d20, d7, d8           @III row i_macro(0,1) = x1(0,1) - x2(0,1)

    vqmovun.s16   d18, q9               @II row CLIP_U8(i_macro)
    vsub.s32      d21, d11, d12         @III row i_macro(2,3) = x1(2,3) - x2(2,3)

    vld1.32       d24[0], [r1], r4      @IV row Load pu1_pred buffer
    vqrshrn.s32   d20, q10, #6          @III row i_macro = D_SHIFT(i_macro,6,shft)

    vmovl.u8      q12, d24              @IV row Convert 8 bit pred buffer to 16 bit
    vadd.u16      d20, d20, d26         @III row i_macro += *pu1_pred_tmp

    vqmovun.s16   d20, q10              @III row CLIP_U8(i_macro)
    vsub.s32      d22, d6, d9           @IV row i_macro(0,1) = x0(0,1) - x3(0,1)

    vst1.32       d18[0], [r2], r5      @II row store the value
    vsub.s32      d23, d10, d13         @IV row i_macro(2,3) = x0(2,3) - x3(2,3)

    vqrshrn.s32   d22, q11, #6          @IV row i_macro = D_SHIFT(i_macro,6,shft)

    vst1.32       d20[0], [r2], r5      @III row store the value
    vadd.u16      d22, d22, d24         @IV row i_macro += *pu1_pred_tmp

    vqmovun.s16   d22, q11              @IV row CLIP_U8(i_macro)
    vst1.32       d22[0], [r2], r5      @IV row store the value


    vpop          {d8-d15}
    ldmfd         sp!, {r4-r12, r15}    @Reload the registers from SP




