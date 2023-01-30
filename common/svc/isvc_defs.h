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
*  isvc_defs.h
*
* @brief
*  Contains macro defintions, and other typedefs used for SVC encoding
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVC_DEFS_H_
#define _ISVC_DEFS_H_

#define MAX_NUM_TEMPORAL_LAYERS 3

#define MAX_NUM_SPATIAL_LAYERS 3

#define MAX_VUI_EXT_NUM_ENTRIES (MAX_NUM_TEMPORAL_LAYERS * MAX_NUM_SPATIAL_LAYERS)

#define SVC_INTER_MB (1 << 0) /*!< Intra MBs other than IPCM and I_BL */

#define SVC_INTRA_MB (1 << 1) /*!< P or B MBs decoded or inferred*/

#define SVC_IPCM_MB (1 << 2) /*!< IPCM_MB  decoder or inferred*/

#define SVC_IBL_MB (1 << 3) /*!< I_BL MB always inferred */

#define SVC_INTRA_INTER_MB                                         \
    (1 << 4) /*!< Intra Inter MB will have an alternate prediction \
                process*/

#define MB_WIDTH_SHIFT 4

#define MB_HEIGHT_SHIFT 4

#define UV 1

#define NUM_SP_COMPONENTS 2

#define NUM_COMPONENTS 3

#define SVC_EXTRACT_MB_MODE(x) ((x) &0x1F)

#define GET_BIT_TX_SIZE(x, y) ((x) & (1 << (7 - (y))))

typedef enum SVC_PROFILES_T
{
    IH264_SCALABLE_BASELINE = 83,
    IH264_SCALABLE_HIGH_PROFILE = 86
} SVC_PROFILES_T;

typedef enum PRED_MODE_T
{
    L0 = 0,
    L1 = 1,
    BI = 2,
    NUM_PRED_DIRS = 2,
    INVALID_PRED_MODE = 4,
} PRED_MODE_T;

#endif
