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
* @remarks
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* User Include Files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264e_defs.h"
#include "ih264e_globals.h"

/*****************************************************************************/
/* Global definitions                                                        */
/*****************************************************************************/

/**
******************************************************************************
* @brief  lambda for varying quantizer scales that would be used to
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
const UWORD8 gu1_qp_lambdaIP[52]=
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
* @brief  lambda for varying quantizer scales that would be used to
* compute the RD cost while deciding on the MB modes.
* input  : qp
* output : lambda
* @remarks lambda = max(2, min(4, pow(2, (qp - 12)/6))) * gu1_qp_lambdaIP[]
******************************************************************************
*/
const UWORD8 gu1_qp_lambdaB[52]=
{
        0,       0,       0,       0,       1,       1,       1,       1,
        1,       1,       1,       1,       1,       1,       2,       2,
        2,       2,       3,       3,       3,       4,       4,       5,
        5,       6,       7,       8,      10,      11,      13,      15,
       17,      20,      22,      26,      30,      33,      37,      42,
       47,      53,      59,      66,      74,      83,      94,     105,
      118,     132,     149,     167,
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
* @brief  Translation of Qstep <-> QP
* Qstep(QP) = Qstep(QP%6) * (2 ^ floor(QP/6))
* Qstep(QP, n = 0) {0.625, 0.6875, 0.8125, 0.875, 1.0, 1.125} for QP[0-5]
* @remarks QPRange[0-51] & QstepRange[1 - 224].
******************************************************************************
*/
const UWORD8 gau1_h264_to_mpeg2_qmap[H264_QP_ELEM] =
{
     1,    1,    1,    1,    1,    1,    1,    1,
     2,    2,    2,    2,    3,    3,    3,    4,
     4,    5,    5,    6,    7,    7,    8,    9,
    10,   11,   13,   14,   16,   18,   20,   22,
    26,   28,   32,   36,   40,   44,   52,   56,
    64,   72,   80,   88,  104,  112,  128,  144,
    160,  176,  208,  224,
};
const UWORD8 gau1_mpeg2_to_h264_qmap[MPEG2_QP_ELEM] =
{
     0,    4,   10,   13,   16,   18,   19,   21,
    22,   23,   24,   25,   25,   26,   27,   27,
    28,   28,   29,   30,   30,   30,   31,   31,
    31,   32,   32,   32,   33,   33,   33,   34,
    34,   34,   34,   35,   35,   36,   36,   36,
    36,   36,   36,   37,   37,   37,   37,   37,
    37,   38,   38,   38,   38,   38,   38,   39,
    39,   39,   39,   39,   39,   40,   40,   40,
    40,   40,   40,   40,   40,   41,   41,   41,
    41,   42,   42,   42,   42,   42,   42,   42,
    42,   42,   42,   42,   42,   43,   43,   43,
    43,   43,   43,   43,   43,   43,   43,   43,
    43,   44,   44,   44,   44,   44,   44,   44,
    44,   44,   44,   44,   44,   45,   45,   45,
    45,   45,   45,   45,   45,   45,   45,   45,
    45,   46,   46,   46,   46,   46,   46,   46,
    46,   46,   46,   46,   46,   46,   46,   46,
    46,   47,   47,   47,   47,   47,   47,   47,
    47,   48,   48,   48,   48,   48,   48,   48,
    48,   48,   48,   48,   48,   48,   48,   48,
    48,   48,   48,   48,   48,   48,   48,   48,
    48,   49,   49,   49,   49,   49,   49,   49,
    49,   49,   49,   49,   49,   49,   49,   49,
    49,   49,   49,   49,   49,   49,   49,   49,
    49,   50,   50,   50,   50,   50,   50,   50,
    50,   50,   50,   50,   50,   50,   50,   50,
    50,   50,   50,   50,   50,   50,   50,   50,
    50,   51,   51,   51,   51,   51,   51,   51,
    51,   51,   51,   51,   51,   51,   51,   51,
    51,   51,   51,   51,   51,   51,   51,   51,
    51,   52,   52,   52,   52,   52,   52,   52,
    52,   52,   52,   52,   52,   52,   52,   52,
};

