/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
*  ih264e_globals.c
*
* @brief
*  Contains definitions of global variables used across the encoder
*
* @author
*  ittiam
*
* @par List of functions
*
*
* @remarks
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264e_defs.h"
#include "ih264e_globals.h"

/*****************************************************************************/
/* Extern global definitions                                                 */
/*****************************************************************************/

/**
******************************************************************************
* @brief  lamda for varying quantizer scales that would be used to
* compute the RD cost while deciding on the MB modes.
* input  : qp
* output : lambda
* @remarks lambda = 0.85 * pow(2, (qp - 12)/3), when SSD is used as metric
* for computing distortion (Bit rate estimation for cost function of H.264/
* AVC by Mohd Golam Sarwer et. al.)  If the use of distortion metric is SAD
* rather than SSD in the stage of encoding, consider sqrt(lambda) simply to
* adjust lambda for the lack of squaring operation in the error computation
* (from rate distortion optimization for video compression by sullivan).
******************************************************************************
*/
const UWORD16 gu2_qp_lambda[52]=
{
       0,      0,      0,      0,      0,      0,      0,      1,
       1,      1,      1,      1,      1,      1,      1,      1,
       1,      2,      2,      2,      2,      3,      3,      3,
       4,      4,      5,      5,      6,      7,      7,      8,
       9,     10,     12,     13,     15,     17,     19,     21,
      23,     26,     30,     33,     37,     42,     47,     53,
      59,     66,     74,     83,
};

/**
******************************************************************************
* @brief  Lamda for varying quantizer scales that would be used to
* compute the RD cost while deciding on the MB modes.
* input  : qp
* output : lambda
* @remarks lambda = pow(2, (qp - 12)/6)
******************************************************************************
*/
const UWORD8 gu1_qp0[52]=
{
       0,      0,      0,      0,      0,      0,      0,      0,
       0,      0,      0,      0,      1,      1,      1,      1,
       2,      2,      2,      2,      3,      3,      3,      4,
       4,      4,      5,      6,      6,      7,      8,      9,
      10,     11,     13,     14,     16,     18,     20,     23,
      25,     29,     32,     36,     40,     45,     51,     57,
      64,     72,     81,     91,
};

/**
******************************************************************************
* @brief  unsigned exp. goulumb codelengths to assign cost to a coefficient of
* mb types.
* input  : Integer
* output : codelength
* @remarks Refer sec. 9-1 in h264 specification
******************************************************************************
*/
const UWORD8 u1_uev_codelength[32] =
{
     1,      3,      3,      5,      5,      5,      5,      7,
     7,      7,      7,      7,      7,      7,      7,      9,
     9,      9,      9,      9,      9,      9,      9,      9,
     9,      9,      9,      9,      9,      9,      9,      11,
};


/**
******************************************************************************
* @brief  Look up table to assign cost to a coefficient of a residual block
* basing on its surrounding coefficients
* input  : Numbers of T1's
* output : coeff_cost
* @remarks Refer Section 2.3 Elimination of single coefficients in inter
* macroblocks in document JVT-O079
******************************************************************************
*/
const UWORD8 gu1_coeff_cost[6] =
{
     3, 2, 2, 1, 1, 1
};

/**
******************************************************************************
* @brief  Indices map to raster scan for luma 4x4 block
* input  : scan index
* output : scan location
* @remarks None
******************************************************************************
*/
const UWORD8 gu1_luma_scan_order[16] =
{
     0,  1,  4,  8,  5,  2,  3,  6,  9,  12, 13, 10, 7,  11, 14, 15
};

/**
******************************************************************************
* @brief  Indices map to raster scan for chroma AC block
* input  : scan index
* output : scan location
* @remarks None
******************************************************************************
*/
const UWORD8 gu1_chroma_scan_order[15] =
{
     1,  4,  8,  5,  2,  3,  6,  9,  12, 13, 10, 7,  11, 14, 15
};

/**
******************************************************************************
* @brief  Indices map to raster scan for luma 4x4 dc block
* input  : scan index
* output : scan location
* @remarks : None
******************************************************************************
*/
const UWORD8 gu1_luma_scan_order_dc[16] =
{
     0, 1,  4,  8,  5,  2,  3,  6,  9,  12, 13, 10, 7,  11, 14, 15
};

/**
******************************************************************************
* @brief  Indices map to raster scan for chroma 2x2 dc block
* input  : scan index
* output : scan location
* @remarks None
******************************************************************************
*/
const UWORD8 gu1_chroma_scan_order_dc[4] =
{
     0, 1,  2,  3
};

/**
******************************************************************************
* @brief  choice of motion vectors to be used during mv prediction
* input  : formatted reference idx comparison metric
* output : mv prediction has to be median or a simple straight forward selec
* tion from neighbors.
* @remarks If only one of the candidate blocks has a reference frame equal to
    the current block then use the same block as the final predictor. A simple
    look up table to assist this mv prediction condition
******************************************************************************
*/
const WORD8 gi1_mv_pred_condition[8] =
{
     -1,    0,    1,    -1,    2,    -1,    -1,    -1
};

/**
******************************************************************************
* @brief  maps the h264 quantizer to the mpeg2 quantizer scale
* input  : h264 qp
* output : equivalent mpeg 2 qp
* @remarks mpeg2qscale = 2 ^ [((h264qp - 12) / 6) + 1]
******************************************************************************
*/
const UWORD8 gau1_h264_to_mpeg2_qmap[H264_QP_ELEM] =
{
       1,      1,      1,      1,      1,      1,      1,      1,
       2,      2,      2,      2,      3,      3,      3,      4,
       4,      4,      5,      6,      6,      7,      8,      9,
      10,     11,     13,     14,     16,     18,     20,     23,
      25,     29,     32,     36,     40,     45,     51,     57,
      64,     72,     81,     91,    102,    114,    128,    144,
     161,    181,    203,    228,
};

/**
******************************************************************************
* @brief  maps the mpeg2 quantizer to the h264 quantizer scale
* input  : mpeg2 qp
* output : equivalent h264qp
* @remarks  MPEG-2 dequantization: (2*QFij + k)*Wij*qscale/32
*      k = 0 (for intra)  k = sign(QFij)
*   H.264 dequantization: (QFij*R(QP%6,i,j))>>(6 - QP/6)
*
*   Excluding the portion of R(QP%6,i,j) that is due to
*   the DCT scale factors, the 6 entries after dividing by 64 (2^6)
*   correspond to dequant values of
*   2.5, 2.8125, 3.125, 3.5625, 3.9375, 4.4375.
*   (a=0.5 b=sqrt(2/5) - refer to JVT-B038.doc)
*
*   Assuming that h264Qp=12 corresponds to MPEG2 qscale of 2
*   (the actual mapping seems to be to MPEG2 qscale of 2.5),
*   and the fact that the effective h264 quantizer changes by
*   a factor of 2 for every 6 steps, the following mapping is
*   obtained:
*    h264qp = 6*(log2(mpeg2qscale/2)) + 12.
*
*   Note that the quant matrix entry assumed for the above
*   equality is 16. Hence when the mpeg2 quant matrix entries
*   are all 16, this lookup can be used as is (which is the
*   default inter quant matrix in mpeg-2).
******************************************************************************
*/
const UWORD8 gau1_mpeg2_to_h264_qmap[MPEG2_QP_ELEM] =
{
       0,      4,     10,     14,     16,     18,     20,     21,     22,     23,     24,     25,     26,     26,     27,     27,
      28,     29,     29,     29,     30,     30,     31,     31,     32,     32,     32,     33,     33,     33,     33,     34,
      34,     34,     35,     35,     35,     35,     35,     36,     36,     36,     36,     37,     37,     37,     37,     37,
      38,     38,     38,     38,     38,     38,     39,     39,     39,     39,     39,     39,     39,     40,     40,     40,
      40,     40,     40,     40,     41,     41,     41,     41,     41,     41,     41,     41,     41,     42,     42,     42,
      42,     42,     42,     42,     42,     42,     43,     43,     43,     43,     43,     43,     43,     43,     43,     43,
      44,     44,     44,     44,     44,     44,     44,     44,     44,     44,     44,     44,     45,     45,     45,     45,
      45,     45,     45,     45,     45,     45,     45,     45,     45,     46,     46,     46,     46,     46,     46,     46,
      46,     46,     46,     46,     46,     46,     46,     46,     47,     47,     47,     47,     47,     47,     47,     47,
      47,     47,     47,     47,     47,     47,     47,     47,     47,     48,     48,     48,     48,     48,     48,     48,
      48,     48,     48,     48,     48,     48,     48,     48,     48,     48,     48,     49,     49,     49,     49,     49,
      49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,     49,
};

