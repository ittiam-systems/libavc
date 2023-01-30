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
*  isvce_error.h
*
* @brief
*  SVC specific error codes
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_ERROR_H_
#define _ISVCE_ERROR_H_

#include "ih264e_error.h"

typedef enum ISVCE_ERRORS_T
{
    /**Invalid SVC params */
    IH264E_INVALID_SVC_PARAMS = IH264E_CODEC_ERROR_START + 0x100,

    /**Invalid num_temporal_layers */
    IH264E_INVALID_NUM_TEMPORAL_LAYERS = IH264E_CODEC_ERROR_START + 0x101,

    /**Invalid num_spatial_layers */
    IH264E_INVALID_NUM_SPATIAL_LAYERS = IH264E_CODEC_ERROR_START + 0x102,

    /**Invalid spatial_res_ratio */
    IH264E_INVALID_SPATIAL_RES_RATIO = IH264E_CODEC_ERROR_START + 0x103,

    /** Weighted prediction not supported */
    IH264E_WEIGHTED_PRED_NOT_SUPPORTED = IH264E_CODEC_ERROR_START + 0x104,

    /** CABAC entropy mode not supported for SVC */
    IH264E_CABAC_NOT_SUPPORTED = IH264E_CODEC_ERROR_START + 0x105,

    /**Invalid input dimensions */
    IH264E_INVALID_SVC_INPUT_DIMENSIONS = IH264E_CODEC_ERROR_START + 0x106,

    /** Invalid init QP */
    IH264E_INVALID_DYN_INIT_QP = IH264E_CODEC_ERROR_START + 0x107,

} ISVCE_ERRORS_T;

#endif
