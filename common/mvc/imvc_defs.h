/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

#ifndef _IMVC_DEFS_H_
#define _IMVC_DEFS_H_

#define MAX_NUM_VIEWS 6

#define LOG2_MAX_NUM_VIEWS 3

#define MAX_NUM_IVP_REFS MAX_NUM_VIEWS

#define MAX_NUM_LEVEL_VALUES_SIGNALLED 1

#define MAX_NUM_OPERATING_POINTS 1

#define MULTIVIEW_HIGH_PROFILE_IDC 118

#define NUM_OF_ZERO_BYTES_BEFORE_START_CODE 2

#define EMULATION_PREVENTION_BYTE 0x03

#define MIN_H264_QP 0

#define MAX_H264_QP 51

#define NUM_SP_COMPONENTS 2

#define NUM_COMPONENTS 3

#define FORCEINLINE __attribute__((always_inline)) inline

typedef void *FT_ALIGNED_ALLOC(void *pv_mem_ctxt, WORD32 i4_alignment, WORD32 i4_size);

typedef void FT_ALIGNED_FREE(void *pv_mem_ctxt, void *pv_buf);

typedef enum COMPONENT_TYPES_T
{
    Y = 0,
    UV = 1,
    U = 1,
    V = 2,
} COMPONENT_TYPES_T;

typedef enum AVC_EXT_NALU_ID_T
{
    UNSPEC_0 = 0,

    SLICE_NON_IDR = 1,

    SLICE_DPA = 2,

    SLICE_DPB = 3,

    SLICE_DPC = 4,

    SLICE_IDR = 5,

    SEI = 6,

    SPS = 7,

    PPS = 8,

    AUD = 9,

    EOSEQ = 10,

    EOSTR = 11,

    FILLER = 12,

    SPSE = 13,

    PREFIX_NAL = 14,

    SUBSET_SPS = 15,

    AUX_PIC = 19,

    CODED_SLICE_EXTENSION = 20,

    UNSPEC_31 = 24

} AVC_EXT_NALU_ID_T;

typedef enum SLICE_TYPES_T
{
    PSLICE = 0,
    BSLICE = 1,
    ISLICE = 2,
    SPSLICE = 3,
    SISLICE = 4,
    MAXSLICE_TYPE,
} SLICE_TYPES_T;

#endif
