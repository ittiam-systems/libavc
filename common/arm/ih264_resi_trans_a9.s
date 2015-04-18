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
@/**
@*******************************************************************************
@* @file
@*  ih264_resi_trans_a9.s
@*
@* @brief
@*  Contains function definitions for residual and forward trans
@*
@* @author
@*  Ittiam
@*
@* @par List of Functions:
@*  ih264_resi_trans_4x4_a9
@*  ih264_resi_trans_8x8_a9
@* @remarks
@*  None
@*
@*******************************************************************************


.text
.p2align 2
@*****************************************************************************
@*
@* Function Name     : ih264_resi_trans_4x4_a9
@* Description       : This function does cf4 of H264 followed by and approximate scaling
@*
@* Arguments         :
@                       R0 :pointer to src buffer
@                       R1 :pointer to pred buffer
@                       R2 :pointer to dst buffer
@                       R3 :src_stride
@                       STACk :pred_stride,dst_stride

@* Values Returned   : NONE
@*
@* Register Usage    :
@* Stack Usage       :
@* Cycles            : Around
@* Interruptiaility  : Interruptable
@*
@* Known Limitations
@*   \Assumptions    :
@*
@* Revision History  :
@*         DD MM YYYY    Author(s)   Changes
@*         30 12 2009    100633      First version
@*
@*****************************************************************************


    .global ih264_resi_trans_4x4_a9
    .extern g_scal_coff_h264_4x4
g_scal_coff_h264_4x4_addr:
    .long g_scal_coff_h264_4x4 - 4x4lbl - 8

ih264_resi_trans_4x4_a9:

    @R0 :pointer to src buffer
    @R1 :pointer to pred buffer
    @R2 :pointer to dst buffer
    @R3 :src_stride
    @STACk :pred_stride,dst_stride

    push          {r4-r12, lr}          @push all the variables first

    mov           r6, sp
    add           r6, r6, #40           @decrement stack pointer,to accomodate two variables
    ldmfd         r6, {r4-r5}           @load the strides into registers
                                        @R4 pred_stride
                                        @R5 dst_stride


    @we have to give the stride as post inrement in VLDR1
    @but since thr stride is from end of row 1 to start of row 2,
    @we need to add the size of the curent row to strides ie we need to add 4 to it (4 bytes)
    @ADD R3,#4
    @ADD R4,#4
    @ADD R5,#4
    @in case of dst the stride represnts 16 bit ie 2*8bits
    @hence we need to add #4 to it and thenm multiply by 2
    @--------------------function loading done------------------------

    @lets find residual
    @data is like 1a -> d0[1:31]  d0[32:64]
    @                    a b c d   # # # #
    vld1.u8       d30, [r0], r3         @load 4 pixels of row1 current buffer
    vld1.u8       d31, [r1], r4         @load 4 pixels of row1 pred buffer
    @ data is like 1a -> q4[1:63]  q4[64:148]
    @                    d8[1:63]  d9[1:63]
    @                    a b c d   # # # #

    vld1.u8       d28, [r0], r3         @load row 2 of src to d28[0]
    vld1.u8       d29, [r1], r4         @load row2 of pred to d29[0]

    vld1.u8       d26, [r0], r3         @load row 3 of src to d26[0]
    vsubl.u8      q0, d30, d31          @curr - pred for row one

    vld1.u8       d27, [r1], r4         @load row 3of pred t0 d27[0]
    vsubl.u8      q1, d28, d29          @find row 2 of src -pred to d0

    vld1.u8       d24, [r0], r3         @load row 4 of src to d24[0]

    vld1.u8       d25, [r1], r4         @load row 4 of src tp d25[0]
    vsubl.u8      q2, d26, d27          @load src-pred row 3 to d[2]

    lsl           r5, r5, #2            @ multiply dst stride by since we are storing 32 bit values
    ldr           r6, g_scal_coff_h264_4x4_addr
4x4lbl:
    add           r6, r6, pc            @  load the address of global array

    vsubl.u8      q3, d24, d25          @load row 4 of src - pred to q6

    @after this
    @D0  -> 1a
    @D2 -> 2a
    @D4 -> 3a
    @D6 -> 4a

    @transpose the matrix so that we can do the horizontal transform first
    @#1 #2  #3  #4
    @a  b   c   d       ---- D0
    @e  f   g   h       -----D2
    @i  j   k   l       -----D4
    @m  n   o   p       -----D6
    @transpose the inner 2x2 blocks
    vtrn.16       d0, d2
    vld1.s16      {q10}, [r6]!          @   load the scaling values 0-7;
    vtrn.16       d4, d6
    @a  e   c   g
    @b  f   d   h
    @i  m   k   o
    @j  n   l   p
    vtrn.32       d0, d4
    vtrn.32       d2, d6
    @a  e   i   m  #1  -- D0 --- x4
    @b  f   j   n  #2  -- D2 --- x5
    @c  g   k   o  #3  -- D4 ----x6
    @d  h   l   p  #4  -- D6 ----x7

    @we have loaded the residuals into the registers , now we need to add and subtract them
    @let us do the horiz transform first

    vsub.s16      d5, d2, d4            @x2 = x5-x6
    vsub.s16      d7, d0, d6            @x3 = x4-x7;

    vadd.s16      d3, d2, d4            @x1 = x5+x6
    vadd.s16      d1, d0, d6            @x0 = x4+x7


    vshl.s16      d31, d7, #1           @
    vshl.s16      d30, d5, #1           @

    vadd.s16      d0, d1, d3            @x0 + x1;
    vsub.s16      d4, d1, d3            @x0 - x1;

    vadd.s16      d2, d31, d5           @U_SHIFT(x3,1,shft) + x2;
    vsub.s16      d6, d7, d30           @x3 - U_SHIFT(x2,1,shft);

    @taking transform again so as to make do vert transform
    vtrn.16       d0, d2
    vtrn.16       d4, d6

    vtrn.32       d0, d4
    vtrn.32       d2, d6

    @let us do vertical transform
    @same code as horiz

    vadd.s16      d1, d0, d6            @x0 = x4+x7
    vadd.s16      d3, d2, d4            @x1 = x5+x6
    vsub.s16      d7, d0, d6            @x3 = x4-x7;
    vsub.s16      d5, d2, d4            @x2 = x5-x6


@Since we are going to do scal / quant or whatever, we are going to divide by
@a 32 bit number. So we have to expand the values

    @VADDL.S16 Q12,D1,D3;x0 + x1
    @VSUBL.S16 Q14,D1,D3;x0 - x1

    @VSHL.S16  D8,D5,#1;
    @VSHL.S16  D9,D7,#1;

    @VADDL.S16 Q13,D9,D5 ; + x2
    @VSUBL.S16 Q15,D7,D8 ;x3 - U_SHIFT(x2,1,shft)

@scaling follows

@now we need to do the scaling,so load the scaling matrix
@mutliplying by the scaling coeffient; store the results from q5-q8 ;

    vadd.s16      d24, d3, d1           @x4 = x0 + x1
    vsub.s16      d28, d1, d3           @x6 = x0 - x1

    vshl.s16      d0, d7, #1            @ U_SHIFT(x3,1,shft)
    vmull.s16     q4, d24, d20          @x4*s0

    vshl.s16      d2, d5, #1            @ U_SHIFT(x2,1,shft)

    vadd.s16      d26, d0, d5           @x5 = U_SHIFT(x3,1,shft) + x2
    vmull.s16     q5, d26, d21          @x5*s1

    vst1.s32      {q4}, [r2], r5        @save 4 pixels of row1 current buffer and increment pointer by stride

    vld1.s16      {q10}, [r6]           @load 8-16 scaling coeffcients

    vsub.s16      d30, d7, d2           @x7 = x3 - U_SHIFT(x2,1,shft)

    vmull.s16     q6, d28, d20          @x6*s2
    vst1.s32      {q5}, [r2], r5

    vmull.s16     q7, d30, d21          @x7*s3


    vst1.s32      {q6}, [r2], r5
    vst1.s32      {q7}, [r2]

    pop           {r4-r12, pc}          @pop back all variables




@*****************************************************************************
@* Function Name     : ih264_resi_trans_8x8_a9
@* Description       : This function does cf8 followd by an approximate normalization of H264
@*
@* Arguments         :
@*                      R0 :pointer to src buffer
@                       R1 :pointer to pred buffer
@                       R2 :pointer to dst buffer
@                       R3 :src_stride
@                       STACk :pred_stride,dst_st
@*
@*
@* Values Returned   : NONE
@*
@* Register Usage    :
@* Stack Usage       :
@* Cycles            : Around
@* Interruptiaility  : Interruptable
@*
@* Known Limitations
@*   \Assumptions    :
@*
@* Revision History  :
@*         DD MM YYYY    Author(s)   Changes
@*         30 12 2009    100633      First version
@*
@*****************************************************************************


    .global ih264_resi_trans_8x8_a9
    .extern g_scal_coff_h264_8x8
g_scal_coff_h264_8x8_addr:
    .long g_scal_coff_h264_8x8 - 8x8lbl - 8


ih264_resi_trans_8x8_a9:

    @R0 :pointer to src buffer
    @R1 :pointer to pred buffer
    @R2 :pointer to dst buffer
    @R3 :src_stride
    @STACk :pred_stride,dst_stride

    push          {r4-r12, lr}          @push all the variables first

    mov           r6, sp
    add           r6, r6, #40           @decrement stack pointer,to accomodate two variables
    ldmfd         r6, {r4-r5}           @load the strides into registers
                                        @R4 pred_stride
                                        @R5 dst_stride

    @we have to give the stride as post inrement in vst1
    @in case of dst the stride represnts 16 bit ie 2*8bits
    @hence we need to add #4 to it and thenm multiply by 2
    @--------------------function loading done------------------------

    @lets find residual
    @data is like 1a -> d0[1:31]  d0[32:64]
    @                    a b c d   # # # #
    vld1.u8       d30, [r0], r3         @load 4 pixels of row1 current buffer
    vld1.u8       d31, [r1], r4         @load 4 pixels of row1 pred buffer

    vld1.u8       d28, [r0], r3         @src  rw2
    vld1.u8       d29, [r1], r4         @pred rw2
    vsubl.u8      q0, d30, d31          @src-pred rw1

    vld1.u8       d26, [r0], r3
    vld1.u8       d27, [r1], r4
    vsubl.u8      q1, d28, d29

    vld1.u8       d24, [r0], r3
    vld1.u8       d25, [r1], r4
    vsubl.u8      q2, d26, d27

    vld1.u8       d22, [r0], r3
    vld1.u8       d23, [r1], r4
    vsubl.u8      q3, d24, d25

    vld1.u8       d20, [r0], r3
    vld1.u8       d21, [r1], r4
    vsubl.u8      q4, d22, d23

    vld1.u8       d18, [r0], r3
    vld1.u8       d19, [r1], r4
    vsubl.u8      q5, d20, d21

    vld1.u8       d16, [r0], r3
    vld1.u8       d17, [r1], r4
    vsubl.u8      q6, d18, d19

    lsl           r5, r5, #2


    vsubl.u8      q7, d16, d17

    @after this
    @Q0 -> 1a
    @Q1 -> 2a
    @Q2 -> 3a
    @Q3 -> 4a
    @Q4 -> 5a
    @Q5 -> 6a
    @Q6 -> 7a
    @Q7 -> 8a

    @transpose the matrix so that we can do the horizontal transform first

    @transpose the inner 2x2 blocks
    vtrn.16       q0, q1
    vtrn.16       q2, q3
    vtrn.16       q4, q5
    vtrn.16       q6, q7

    @transpose the inner 4x4 blocks
    vtrn.32       q0, q2
    vtrn.32       q1, q3

    vtrn.32       q4, q6
    vtrn.32       q5, q7

    @transpose the outer 8x8 blocks
    vswp          d1, d8
    vswp          d7, d14
    vswp          d3, d10
    vswp          d5, d12
    @transpose done

@@this point we will have data in Q0-Q7
@Q7 will be populated within 2 clock cycle
@all others are availabe @ this clock cycle

    @we have loaded the residuals into the registers , now we need to add and subtract them
    @let us do the horiz transform first

    vadd.s16      q8, q0, q7            @      a0 = r0 + r7;
    vadd.s16      q9, q1, q6            @      a1 = r1 + r6;
    vadd.s16      q10, q2, q5           @     a2 = r2 + r5;
    vadd.s16      q11, q3, q4           @     a3 = r3 + r4;

    vsub.s16      q12, q0, q7           @     b0 = r0 - r7;
    vsub.s16      q13, q1, q6           @     b1 = r1 - r6;
    vsub.s16      q15, q3, q4           @     b3 = r3 - r4;
    vsub.s16      q14, q2, q5           @     b2 = r2 - r5;

    vadd.s16      q1, q8, q11           @     a4 = a0 + a3;
    vadd.s16      q3, q9, q10           @     a5 = a1 + a2;
    vsub.s16      q7, q9, q10           @     a7 = a1 - a2;
    vsub.s16      q5, q8, q11           @     a6 = a0 - a3;

    ldr           r6, g_scal_coff_h264_8x8_addr
8x8lbl:
    add           r6, r6, pc            @  load the address of global array

    vadd.s16      q0, q1, q3            @      pi2_res[0] = a4 + a5;
    vshr.s16      q8, q7, #1            @      pi2_res[2] = a6 + D_SHIFT(a7,1,shft);

    vsub.s16      q4, q1, q3            @      pi2_res[4] = a4 - a5;

    vadd.s16      q2, q5, q8            @


    vshr.s16      q9, q5, #1            @      pi2_res[6] = D_SHIFT(a6,1,shft) - a7;
    vsub.s16      q6, q9, q7            @

@do not change Q0,Q2.Q4,Q6 they contain results
@Q1,Q3,Q5,Q7 TO STORE RESULTS
@Q8 Q9 Q10 Q11 USE @WILL

    vshr.s16      q1, q12, #1           @     D_SHIFT(b0,1,shft)
    vshr.s16      q3, q13, #1           @     D_SHIFT(b1,1,shft)
    vshr.s16      q5, q14, #1           @     D_SHIFT(b2,1,shft)
    vshr.s16      q7, q15, #1           @     D_SHIFT(b3,1,shft)

    vadd.s16      q8, q1, q12           @     (D_SHIFT(b0,1,shft) + b0);
    vadd.s16      q9, q3, q13           @     (D_SHIFT(b1,1,shft) + b1);
    vadd.s16      q10, q5, q14          @    (D_SHIFT(b2,1,shft) + b2);
    vadd.s16      q11, q7, q15          @    (D_SHIFT(b3,1,shft) + b3);

    vadd.s16      q1, q14, q8           @     b2 + (D_SHIFT(b0,1,shft) + b0);
    vsub.s16      q5, q15, q9           @     b3 - (D_SHIFT(b1,1,shft) + b1);
    vadd.s16      q3, q15, q10          @    b3 + (D_SHIFT(b2,1,shft) + b2);
    vsub.s16      q7, q11, q14          @    -b2 + (D_SHIFT(b3,1,shft) + b3);

    vadd.s16      q8, q13, q1           @     b4 = b1 + b2 + (D_SHIFT(b0,1,shft) + b0);
    vsub.s16      q9, q12, q3           @     b5 = b0 - b3 - (D_SHIFT(b2,1,shft) + b2);
    vadd.s16      q10, q12, q5          @    b6 = b0 + b3 - (D_SHIFT(b1,1,shft) + b1);
    vadd.s16      q11, q13, q7          @    b7 = b1 - b2 + (D_SHIFT(b3,1,shft) + b3);

    vshr.s16      q15, q8, #2           @     D_SHIFT(b4,2,shft)
    vshr.s16      q14, q9, #2           @     D_SHIFT(b5,2,shft);
    vshr.s16      q13, q10, #2          @    D_SHIFT(b6,2,shft);
    vshr.s16      q12, q11, #2          @    D_SHIFT(b7,2,shft);


    vadd.s16      q3, q9, q13           @     pi2_res[3] = b5 + D_SHIFT(b6,2,shft);
    vsub.s16      q5, q10, q14          @    pi2_res[5] = b6 - D_SHIFT(b5,2,shft);
    vadd.s16      q1, q8, q12           @     pi2_res[1] = b4 + D_SHIFT(b7,2,shft);
    vsub.s16      q7, q15, q11          @    pi2_res[7] = D_SHIFT(b4,2,shft) - b7;

    @------------horiz transform done-------------------------
    @results are in Q0-Q7
    @all other neon registes can be used at will

@doing vertical transform
@code exact copy of horiz transform above

    @transpose the inner 2x2 blocks
    vtrn.16       q0, q1
    vtrn.16       q2, q3
    vtrn.16       q4, q5
    vtrn.16       q6, q7

    @transpose the inner 4x4 blocks
    vtrn.32       q0, q2
    vtrn.32       q1, q3

    vtrn.32       q4, q6
    vtrn.32       q5, q7

    @transpose the outer 8x8 blocks
    vswp          d1, d8
    vswp          d3, d10
    vswp          d5, d12
    vswp          d7, d14

    @transpose done

    vadd.s16      q8, q0, q7            @      a0 = r0 + r7;
    vadd.s16      q9, q1, q6            @      a1 = r1 + r6;
    vadd.s16      q10, q2, q5           @     a2 = r2 + r5;
    vadd.s16      q11, q3, q4           @     a3 = r3 + r4;

    vsub.s16      q12, q0, q7           @     b0 = r0 - r7;
    vsub.s16      q13, q1, q6           @     b1 = r1 - r6;
    vsub.s16      q14, q2, q5           @     b2 = r2 - r5;
    vsub.s16      q15, q3, q4           @     b3 = r3 - r4;

    vadd.s16      q1, q8, q11           @     a4 = a0 + a3;
    vadd.s16      q3, q9, q10           @     a5 = a1 + a2;
    vsub.s16      q5, q8, q11           @     a6 = a0 - a3;
    vsub.s16      q7, q9, q10           @     a7 = a1 - a2;


    vadd.s16      q0, q1, q3            @      pi2_res[0] = a4 + a5;

    vshr.s16      q8, q7, #1            @      pi2_res[2] = a6 + D_SHIFT(a7,1,shft);
    @DSHIFT_TO_0 Q8,Q7,#1,#0
    vadd.s16      q2, q5, q8            @

    vsub.s16      q4, q1, q3            @      pi2_res[4] = a4 - a5;

    vshr.s16      q9, q5, #1            @      pi2_res[6] = D_SHIFT(a6,1,shft) - a7;
    vsub.s16      q6, q9, q7            @

@do not change Q0,Q2.Q4,Q6 they contain results
@Q1,Q3,Q5,Q7 TO STORE RESULTS
@Q8 Q9 Q10 Q11 USE @WILL

    vshr.s16      q1, q12, #1           @     D_SHIFT(b0,1,shft)
    vshr.s16      q3, q13, #1           @     D_SHIFT(b1,1,shft)
    vshr.s16      q5, q14, #1           @     D_SHIFT(b2,1,shft)
    vshr.s16      q7, q15, #1           @     D_SHIFT(b3,1,shft)


    vadd.s16      q8, q1, q12           @     (D_SHIFT(b0,1,shft) + b0);
    vadd.s16      q9, q3, q13           @     (D_SHIFT(b1,1,shft) + b1);
    vadd.s16      q10, q5, q14          @    (D_SHIFT(b2,1,shft) + b2);
    vadd.s16      q11, q7, q15          @    (D_SHIFT(b3,1,shft) + b3);

    vadd.s16      q1, q14, q8           @     b2 + (D_SHIFT(b0,1,shft) + b0);
    vadd.s16      q3, q15, q10          @    b3 + (D_SHIFT(b2,1,shft) + b2);
    vsub.s16      q5, q15, q9           @     b3 - (D_SHIFT(b1,1,shft) + b1);
    vsub.s16      q7, q11, q14          @    -b2 + (D_SHIFT(b3,1,shft) + b3);

    vadd.s16      q8, q13, q1           @     b4 = b1 + b2 + (D_SHIFT(b0,1,shft) + b0);
    vsub.s16      q9, q12, q3           @     b5 = b0 - b3 - (D_SHIFT(b2,1,shft) + b2);
    vadd.s16      q10, q12, q5          @    b6 = b0 + b3 - (D_SHIFT(b1,1,shft) + b1);
    vadd.s16      q11, q13, q7          @    b7 = b1 - b2 + (D_SHIFT(b3,1,shft) + b3);

    vshr.s16      q15, q8, #2           @     D_SHIFT(b4,2,shft)
    vshr.s16      q14, q9, #2           @     D_SHIFT(b5,2,shft);
    vshr.s16      q13, q10, #2          @    D_SHIFT(b6,2,shft);
    vshr.s16      q12, q11, #2          @    D_SHIFT(b7,2,shft);


@since we are going to scal by small values, we need not expand the guys to 32 bit bit values
    vsub.s16      q5, q10, q14          @    pi2_res[5] = b6 - D_SHIFT(b5,2,shft);
    vsub.s16      q7, q15, q11          @    pi2_res[7] = D_SHIFT(b4,2,shft) - b7;
    vadd.s16      q3, q9, q13           @     pi2_res[3] = b5 + D_SHIFT(b6,2,shft);
    vadd.s16      q1, q8, q12           @     pi2_res[1] = b4 + D_SHIFT(b7,2,shft);

    @------------vert transform done-------------------------
    @results are in Q0-Q7
    @all other neon registes can be used at will

    @scaling
    @since the 8x8 scaling matrix repeats in 1x4,1x4 block ,
    @we need only load 4 values for each row and in total 4 rows
    vld1.s16      {q14-q15}, [r6]       @

    @since we need to get a 32 bit o/p for two 16 bit multiplications
    @we need a VMULL instruction
@-----------------------------first and second row

    vmull.s16     q8, d0, d28           @scale the first row first 4 elem
    vmull.s16     q9, d28, d1           @scale the second row last 4 elemts

    vmull.s16     q10, d2, d29          @ scale second row first 4 elem
    vmull.s16     q11, d29, d3          @scale the second row last 4 elem
    vmull.s16     q12, d4, d30          @scale third row first  4 elem

    vst1.s32      {q8, q9}, [r2], r5    @ write the first row complete

    vmull.s16     q13, d30, d5          @scale the third row last 4 elem
    vmull.s16     q8, d6, d31           @scale the fourth row first 4 elem


    vst1.s32      {q10, q11}, [r2], r5  @store the second row complete

@------------------------------- 3rd and 4th row

    vmull.s16     q9, d31, d7           @scale the fourth row second column

    vst1.s32      {q12, q13}, [r2], r5  @store the third row complete

    vmull.s16     q10, d8, d28          @scale the 5th row fisrst 4 elms
    vmull.s16     q11, d28, d9          @scale the 5th row second 4 elems

    vmull.s16     q12, d10, d29         @scale the 6th row first4 elements


    vst1.s32      {q8, q9}, [r2], r5    @store fifth row

@--------------------------------5th and 6th row

    vmull.s16     q13, d29, d11         @scale 6th row sendond 4 elems

    vmull.s16     q8, d12, d30          @scale 7th rw first 4 elms

    vst1.s32      {q10, q11}, [r2], r5  @store 6th row second 4 elements

    vmull.s16     q9, d30, d13          @scale 7th rw second 4 elms
    vmull.s16     q10, d14, d31         @scale 8th rw forst 4 elms


    vst1.s32      {q12, q13}, [r2], r5  @store 6th row

@----------------------------------7th and 8th row
    vmull.s16     q11, d31, d15         @scale 8th row second 4 elms

    vst1.s32      {q8, q9}, [r2], r5    @store 7th row
    vst1.s32      {q10, q11}, [r2], r5  @store 8th row

@----------------------------------done writing

    pop           {r4-r12, pc}          @pop back all variables






