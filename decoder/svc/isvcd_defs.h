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
 *  isvcd_defs.h
 *
 * @brief
 *  Type definitions used in the code
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_DEFS_H_
#define _ISVCD_DEFS_H_

#include <stdint.h>

typedef enum
{
    ERROR_SVC_FIELD_PIC_UNSUPPORTED = 0xC1,
    ERROR_SVC_INV_SCAN_IDX = 0xC2,
    ERROR_SVC_INV_NAL_UNIT = 0xC3,
    ERROR_SVC_INV_SUBSET_SPS = 0xC4
} isvcd_decoder_error_code_t;

#define FLUSH 2

#define SCALABLE_BASELINE_PROFILE_IDC 83
#define SCALABLE_HIGH_PROFILE_IDC 86
#define SCALABLE_HIGH_INTRA_IDC 118

#define SPS_EXTENSION_NAL 13
#define PREFIX_UNIT_NAL 14
#define SUBSET_SPS_NAL 15
#define CODED_SLICE_EXTENSION_NAL 20

#define EP_SLICE 5
#define EB_SLICE 6
#define EI_SLICE 7

#define D_INTRA_IBL 16

#define CAB_INFERRED 0xFF

#define MAX_TOTAL_LYRS (MAX_QUALITY_ID + 1) * (MAX_DEPENDENCY_ID + 1) * 16

#define MAX_QUALITY_LYRS 5 /* ReqRename */
/*!< Maximum number of layers with same
     dependency id
 */
#define MAX_DEPENDENCY_LYRS 6 /* ReqRename */
/*!< Maximum number of layers with
     different dependency id
 */

/** Maximum number of layers without spatial resolution change */
#define MAX_NUM_LYRS_IN_RES 5

/** Maximum number of dependency layers in a resolution */
#define MAX_DEP_LYRS_IN_RES 3

/* Maximum number of spatial resolution layers */
#define MAX_NUM_RES_LYRS 3

/* Maximum number of layers in an access unit */
#define MAX_NUM_LAYERS MAX_NUM_LYRS_IN_RES *MAX_NUM_RES_LYRS

#define MAX_NUM_PIC_BUFS (32 + 1)

/*SVC Standard Specific Macros*/
#define MAX_QUALITY_ID 0
#define MAX_DEPENDENCY_ID 4
#define MAX_TEMPORAL_ID 7
#define MAX_PRIORITY_ID 63
#define MAX_REF_DEP_ID ((MAX_DEPENDENCY_ID << 4) | MAX_QUALITY_ID)

#define BASE_LAYER 0
#define MEDIAL_ENHANCEMENT_LAYER 1
#define TARGET_LAYER 2

#define MB_INFER 250

#define SVC_INTER_MB (1 << 0)       /*!< Intra MBs other than IPCM and I_BL */
#define SVC_INTRA_MB (1 << 1)       /*!< P or B MBs decoded or inferred*/
#define SVC_IPCM_MB (1 << 2)        /*!< IPCM_MB  decoder or inferred*/
#define SVC_IBL_MB (1 << 3)         /*!< I_BL MB always inferred */
#define SVC_INTRA_INTER_MB (1 << 4) /*!< Intra Inter MB */

#endif                              /*_ISVCD_DEFS_H_*/
