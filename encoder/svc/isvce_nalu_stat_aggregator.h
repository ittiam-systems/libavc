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
*  isvce_nalu_stat_aggregator.h
*
* @brief
*  Contains objects used for aggregating nalu statistics
*
*******************************************************************************
*/

#ifndef _ISVCE_NALU_STAT_AGGREGATOR_H_
#define _ISVCE_NALU_STAT_AGGREGATOR_H_

#include <stdbool.h>

#include "ih264_typedefs.h"
#include "isvce.h"
#include "isvc_defs.h"
#include "isvce_defs.h"

/* Macros */
/* +1 for '\0' */
#define MAX_BYTES_PER_NALU_INFO (45 + 1)

/* SPS + (MAX_NUM_SPATIAL_LAYERS - 1) * SUBSET_SPS +
 * MAX_NUM_SPATIAL_LAYERS * PPS + */
/* 1 PREFIX_NALU + 1 SLICE_[NON|]IDR + (MAX_NUM_SPATIAL_LAYERS - 1) *
 * CODED_SLICE_EXTENSION */
#define MAX_NALU_PER_LAYER 10

/* Structs */
typedef struct nalu_info_t
{
    NAL_UNIT_TYPE_T e_nalu_type;

    WORD64 i8_num_bits;

    bool b_is_vcl_nal;

    bool b_is_idr;

    UWORD8 u1_spatial_layer_id;

    UWORD8 u1_temporal_layer_id;

    UWORD8 u1_num_slices;
} nalu_info_t;

typedef struct nalu_descriptors_t
{
    nalu_info_t as_nalu_info[MAX_NALU_PER_LAYER];

    UWORD8 u1_num_nalus;

} nalu_descriptors_t;

/* Function declarations */
static FORCEINLINE UWORD32 isvce_get_nalu_info_buf_size(UWORD8 u1_num_spatial_layers)
{
    return MAX_NALU_PER_LAYER * u1_num_spatial_layers * MAX_BYTES_PER_NALU_INFO;
}

extern void isvce_nalu_info_au_init(nalu_descriptors_t *ps_nalu_descriptor,
                                    UWORD8 u1_num_spatial_layers);

extern void isvce_nalu_info_csv_translator(nalu_descriptors_t *ps_nalu_descriptor,
                                           isvce_nalu_info_buf_t *ps_csv_buf);

extern nalu_info_t *isvce_get_next_nalu_info_buf(nalu_descriptors_t *ps_nalu_descriptor);

extern void isvce_nalu_info_buf_init(nalu_info_t *ps_nalu_info, WORD64 i8_init_bytes,
                                     NAL_UNIT_TYPE_T e_nalu_type, UWORD8 u1_spatial_layer_id,
                                     UWORD8 u1_temporal_layer_id, UWORD8 u1_num_slices,
                                     bool b_is_idr);

extern void isvce_update_nalu_count(nalu_descriptors_t *ps_nalu_descriptor);

#endif
