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
*  isvce_nalu_stat_aggregator.c
*
* @brief
*  Contains objects used for aggregating nalu statistics
*
*******************************************************************************
*/
#include <stdio.h>
#include <string.h>

#include "ih264_typedefs.h"
#include "iv2.h"
#include "isvce_structs.h"
#include "isvce_nalu_stat_aggregator.h"

void isvce_nalu_info_au_init(nalu_descriptors_t *ps_nalu_descriptor, UWORD8 u1_num_spatial_layers)
{
    WORD32 i;

    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        ps_nalu_descriptor[i].u1_num_nalus = 0;
    }
}

void isvce_nalu_info_csv_translator(nalu_descriptors_t *ps_nalu_descriptor,
                                    isvce_nalu_info_buf_t *ps_csv_buf)
{
    char ac_csv_string[MAX_BYTES_PER_NALU_INFO];
    WORD32 i;

    WORD64 i8_num_bytes_available = ps_csv_buf->u4_buf_size - ps_csv_buf->u4_num_bytes;

    for(i = 0; i < ps_nalu_descriptor->u1_num_nalus; i++)
    {
        if(ps_nalu_descriptor->as_nalu_info[i].b_is_vcl_nal)
        {
            snprintf(ac_csv_string, MAX_BYTES_PER_NALU_INFO, "%d,%u,%d,%d,%d,%d,%d\n",
                     ps_nalu_descriptor->as_nalu_info[i].e_nalu_type,
                     (UWORD32) (ps_nalu_descriptor->as_nalu_info[i].i8_num_bits / 8),
                     ps_nalu_descriptor->as_nalu_info[i].u1_spatial_layer_id,
                     ps_nalu_descriptor->as_nalu_info[i].u1_temporal_layer_id,
                     ps_nalu_descriptor->as_nalu_info[i].b_is_idr, 1, 1);
        }
        else
        {
            snprintf(ac_csv_string, MAX_BYTES_PER_NALU_INFO, "%d,%u,%d,%d,%d,%d,%d\n",
                     ps_nalu_descriptor->as_nalu_info[i].e_nalu_type,
                     (UWORD32) (ps_nalu_descriptor->as_nalu_info[i].i8_num_bits / 8), -1, -1, -1,
                     -1, -1);
        }

        snprintf((char *) (ps_csv_buf->pu1_buf + ps_csv_buf->u4_num_bytes), i8_num_bytes_available,
                 "%s", ac_csv_string);

        ps_csv_buf->u4_num_bytes = (UWORD32) strlen((char *) ps_csv_buf->pu1_buf);
        i8_num_bytes_available = ps_csv_buf->u4_buf_size - ps_csv_buf->u4_num_bytes;

        ASSERT(i8_num_bytes_available >= 0);
    }
}

nalu_info_t *isvce_get_next_nalu_info_buf(nalu_descriptors_t *ps_nalu_descriptor)
{
    return &ps_nalu_descriptor->as_nalu_info[ps_nalu_descriptor->u1_num_nalus];
}

void isvce_nalu_info_buf_init(nalu_info_t *ps_nalu_info, WORD64 i8_init_bits,
                              NAL_UNIT_TYPE_T e_nalu_type, UWORD8 u1_spatial_layer_id,
                              UWORD8 u1_temporal_layer_id, UWORD8 u1_num_slices, bool b_is_idr)
{
    ps_nalu_info->e_nalu_type = e_nalu_type;
    ps_nalu_info->i8_num_bits = i8_init_bits;
    ps_nalu_info->b_is_idr = b_is_idr;

    switch(e_nalu_type)
    {
        case NAL_SLICE_NON_IDR:
        case NAL_SLICE_IDR:
        case NAL_CODED_SLICE_EXTENSION:
        {
            ps_nalu_info->b_is_vcl_nal = true;
            ps_nalu_info->u1_spatial_layer_id = u1_spatial_layer_id;
            ps_nalu_info->u1_temporal_layer_id = u1_temporal_layer_id;
            ps_nalu_info->u1_num_slices = u1_num_slices;

            break;
        }
        default:
        {
            ps_nalu_info->b_is_vcl_nal = false;

            break;
        }
    }
}

void isvce_update_nalu_count(nalu_descriptors_t *ps_nalu_descriptor)
{
    ps_nalu_descriptor->u1_num_nalus++;
}
