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
*  isvc_svc_ilp_mv_private_defs.h
*
* @brief
*  Contains datatype and macro definitions used exclusively in
*  ILP MV derivations
*
*******************************************************************************
*/

#ifndef _ISVCE_ILP_MV_PRIVATE_DEFS_H_
#define _ISVCE_ILP_MV_PRIVATE_DEFS_H_

#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"

/* Structs */
/* Offsets, etc used for resLayer MV upsampling */
/* Derived as per 'G.8.6.1.1' for all MB's once during init */
typedef struct ilp_mv_mb_state_t
{
    coordinates_t as_pu_positions[MAX_PU_IN_MB_COL][MAX_PU_IN_MB_ROW];

    coordinates_t as_mb_positions[MAX_PU_IN_MB_COL][MAX_PU_IN_MB_ROW];
} ilp_mv_mb_state_t;

typedef struct ilp_mv_layer_state_t
{
    layer_resampler_props_t *ps_props;

    ilp_mv_mb_state_t *ps_mb_states;

    coordinates_t s_mv_scale;

} ilp_mv_layer_state_t;

typedef struct ilp_mv_state_t
{
    /* Array of size numSpatialLayers */
    ilp_mv_layer_state_t *ps_layer_state;

} ilp_mv_state_t;

#endif
