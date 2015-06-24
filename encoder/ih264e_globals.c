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


/*******************************************************************************
 * Translation of MPEG QP to H264 QP
 ******************************************************************************/
/*
 * Note : RC library models QP and bits assuming the QP to be MPEG2.
 *        Since MPEG qp varies linearly, when the relationship is computed,
 *        it learns that delta(qp) => delta(bits). Now what we are doing by the
 *        transation of qp is that
 *              QPrc = a + b*2^(QPen)
 *        By not considering the weight matrix in both MPEG and H264 we in effect
 *        only changing the relation to
 *              QPrc = c + d*2^(QPen)
 *        This will only entatil changin the RC model parameters, and this will
 *        not affect rc relation at all
 *
 *
 * We have MPEG qp which varies from 0-228. The quantization factor has a linear
 * relation ship with the size of quantized values
 *
 * We also have H264 Qp, which varies such that for a change in QP of 6 , we
 * double the corresponding scaling factor. Hence the scaling is linear in terms
 * of 2^(QPh/6)
 *
 * Now we want to have translation between QPm and QPh. Hence we can write
 *
 * QPm = a + b*2^(QPh/6)
 *
 * Appling boundary condition that
 *      1) QPm = 1 if QPh = 0
 *      2) QPm = 228 if QPh = 51,
 *
 * we will have
 *  a = -0.372, b = 0.628
 *
 * Hence the relatiohship is
 *  QPm = a + b*2^(Qph/6)
 *  QPh = 6*log((Qpm - a)/b)
 *
 *
 * Unrounded values for gau1_h264_to_mpeg2_qmap[H264_QP_ELEM] =
 *
 *   0.33291     0.41923     0.51613     0.62489     0.74697     0.88400
 *   1.03781     1.21046     1.40425     1.62178     1.86594     2.14000
 *   2.44762     2.79292     3.18050     3.61555     4.10388     4.65200
 *   5.26725     5.95784     6.73301     7.60310     8.57975     9.67600
 *   10.90650    12.28769    13.83802    15.57821    17.53150    19.72400
 *   22.18500    24.94737    28.04804    31.52841    35.43500    39.82000
 *   44.74199    50.26675    56.46807    63.42882    71.24200    80.01200
 *   89.85599    100.90549   113.30814   127.22965   142.85601   160.39600
 *   180.08398   202.18299   226.98829
 *
 * Unrounded values for gau1_mpeg2_to_h264_qmap[MPEG2_QP_ELEM]
 *
 *  -4.5328    6.7647   11.5036   14.5486   16.7967   18.5797   20.0575
 *  21.3193   22.4204   23.3971   24.2747   25.0715   25.8010   26.4738
 *  27.0981   27.6804   28.2259   28.7391   29.2236   29.6824   30.1181
 *  30.5329   30.9287   31.3072   31.6699   32.0180   32.3526   32.6748
 *  32.9854   33.2852   33.5750   33.8554   34.1270   34.3904   34.6460
 *  34.8942   35.1355   35.3703   35.5989   35.8216   36.0387   36.2505
 *  36.4572   36.6591   36.8564   37.0494   37.2381   37.4228   37.6036
 *  37.7807   37.9543   38.1244   38.2913   38.4550   38.6157   38.7735
 *  38.9284   39.0806   39.2302   39.3772   39.5218   39.6640   39.8039
 *  39.9416   40.0771   40.2106   40.3420   40.4714   40.5990   40.7247
 *  40.8486   40.9707   41.0911   41.2099   41.3271   41.4427   41.5568
 *  41.6694   41.7806   41.8903   41.9987   42.1057   42.2115   42.3159
 *  42.4191   42.5211   42.6219   42.7216   42.8201   42.9175   43.0138
 *  43.1091   43.2033   43.2965   43.3887   43.4799   43.5702   43.6596
 *  43.7480   43.8356   43.9223   44.0081   44.0930   44.1772   44.2605
 *  44.3431   44.4248   44.5058   44.5861   44.6656   44.7444   44.8224
 *  44.8998   44.9765   45.0525   45.1279   45.2026   45.2766   45.3501
 *  45.4229   45.4951   45.5667   45.6378   45.7082   45.7781   45.8474
 *  45.9162   45.9844   46.0521   46.1193   46.1859   46.2521   46.3177
 *  46.3829   46.4475   46.5117   46.5754   46.6386   46.7014   46.7638
 *  46.8256   46.8871   46.9481   47.0087   47.0689   47.1286   47.1880
 *  47.2469   47.3054   47.3636   47.4213   47.4787   47.5357   47.5923
 *  47.6486   47.7045   47.7600   47.8152   47.8700   47.9245   47.9787
 *  48.0325   48.0859   48.1391   48.1919   48.2444   48.2966   48.3485
 *  48.4000   48.4513   48.5022   48.5529   48.6033   48.6533   48.7031
 *  48.7526   48.8018   48.8508   48.8995   48.9478   48.9960   49.0438
 *  49.0914   49.1388   49.1858   49.2327   49.2792   49.3256   49.3716
 *  49.4175   49.4630   49.5084   49.5535   49.5984   49.6430   49.6875
 *  49.7317   49.7756   49.8194   49.8629   49.9062   49.9493   49.9922
 *  50.0348   50.0773   50.1196   50.1616   50.2034   50.2451   50.2865
 *  50.3278   50.3688   50.4097   50.4503   50.4908   50.5311   50.5712
 *  50.6111   50.6508   50.6904   50.7298   50.7690   50.8080   50.8468
 *  50.8855   50.9240   50.9623   51.0004   51.0384
 *
 */

const UWORD8 gau1_h264_to_mpeg2_qmap[H264_QP_ELEM] =
{
                 1,     1,     1,     1,     1,     1,
                 1,     1,     1,     2,     2,     2,
                 2,     3,     3,     4,     4,     5,
                 5,     6,     7,     8,     9,     10,
                 11,    12,    14,    16,    18,    20,
                 22,    25,    28,    32,    35,    40,
                 45,    50,    56,    63,    71,    80,
                 90,    101,   113,   127,   143,   160,
                 180,   202,   227
};

const UWORD8 gau1_mpeg2_to_h264_qmap[MPEG2_QP_ELEM] =
{
                 0,    7,    12,   15,   17,   19,   20,
                 21,   22,   23,   24,   25,   26,   26,
                 27,   28,   28,   29,   29,   30,   30,
                 31,   31,   31,   32,   32,   32,   33,
                 33,   33,   34,   34,   34,   34,   35,
                 35,   35,   35,   36,   36,   36,   36,
                 36,   37,   37,   37,   37,   37,   38,
                 38,   38,   38,   38,   38,   39,   39,
                 39,   39,   39,   39,   40,   40,   40,
                 40,   40,   40,   40,   40,   41,   41,
                 41,   41,   41,   41,   41,   41,   42,
                 42,   42,   42,   42,   42,   42,   42,
                 42,   43,   43,   43,   43,   43,   43,
                 43,   43,   43,   43,   43,   44,   44,
                 44,   44,   44,   44,   44,   44,   44,
                 44,   44,   45,   45,   45,   45,   45,
                 45,   45,   45,   45,   45,   45,   45,
                 45,   45,   46,   46,   46,   46,   46,
                 46,   46,   46,   46,   46,   46,   46,
                 46,   46,   47,   47,   47,   47,   47,
                 47,   47,   47,   47,   47,   47,   47,
                 47,   47,   47,   47,   47,   48,   48,
                 48,   48,   48,   48,   48,   48,   48,
                 48,   48,   48,   48,   48,   48,   48,
                 48,   48,   49,   49,   49,   49,   49,
                 49,   49,   49,   49,   49,   49,   49,
                 49,   49,   49,   49,   49,   49,   49,
                 49,   49,   50,   50,   50,   50,   50,
                 50,   50,   50,   50,   50,   50,   50,
                 50,   50,   50,   50,   50,   50,   50,
                 50,   50,   50,   50,   50,   51,   51,
                 51,   51,   51,   51,   51,   51,   51,
                 51,   51,   51,   51,   51
};

