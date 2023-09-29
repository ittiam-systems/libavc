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
******************************************************************************
* @file
*  ih264e_rc_mem_interface.c
*
* @brief
*  This file contains api function definitions for rate control memtabs
*
* @author
*  ittiam
*
* List of Functions
*  - ih264e_map_rc_mem_recs_to_itt_api
*  - ih264e_map_itt_mem_rec_to_rc_mem_rec
*  - ih264e_get_rate_control_mem_tab
*
* @remarks
*  none
*
*******************************************************************************
*/


/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System Include Files */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* User Include Files */
#include "ih264e_config.h"
#include "ih264_typedefs.h"

#include "ih264_debug.h"
#include "ih264_defs.h"
#include "ih264_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_structs.h"
#include "ih264_size_defs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "ih264_cabac_tables.h"

#include "ime_defs.h"
#include "ime_distortion_metrics.h"
#include "ime_structs.h"

#include "irc_mem_req_and_acq.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "irc_rate_control_api.h"
#include "irc_common.h"
#include "irc_rd_model.h"
#include "irc_est_sad.h"
#include "irc_fixed_point_error_bits.h"
#include "irc_vbr_storage_vbv.h"
#include "irc_picture_type.h"
#include "irc_bit_allocation.h"
#include "irc_mb_model_based.h"
#include "irc_cbr_buffer_control.h"
#include "irc_vbr_str_prms.h"
#include "irc_rate_control_api_structs.h"

#include "ih264e.h"
#include "ih264e_error.h"
#include "ih264e_defs.h"
#include "ih264e_time_stamp.h"
#include "ih264e_modify_frm_rate.h"
#include "ih264e_rate_control.h"
#include "ih264e_bitstream.h"
#include "ih264e_cabac_structs.h"
#include "ih264e_structs.h"


/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
******************************************************************************
*
* @brief This function maps rc mem records structure to encoder lib mem records
*  structure
*
* @par   Description
*  This function maps rc mem records structure to encoder lib mem records
*  structure
*
* @param[in]   ps_mem
*  pointer to encoder lib mem records
*
* @param[in]   rc_memtab
*  pointer to rc mem records
*
* @param[in]   num_mem_recs
*  number of memory records
*
* @return      void
*
******************************************************************************
*/
void ih264e_map_rc_mem_recs_to_itt_api(iv_mem_rec_t *ps_mem,
                                       itt_memtab_t *rc_memtab,
                                       UWORD32 num_mem_recs)
{
    UWORD32 j;
    UWORD32 Size, align;

    for (j = 0; j < num_mem_recs; j++)
    {
        Size = rc_memtab->u4_size;
        align = rc_memtab->i4_alignment;

        /* we always ask for external persistent cacheable memory */
        FILL_MEMTAB(ps_mem, j, Size, align, IV_EXTERNAL_CACHEABLE_PERSISTENT_MEM);

        rc_memtab++;
    }
}

/**
*******************************************************************************
*
* @brief This function maps encoder lib mem records structure to RC memory
* records structure
*
* @par   Description
*  This function maps encoder lib mem records structure to RC memory
*  records structure
*
* @param[in] ps_mem
*  pointer to encoder lib mem records
*
* @param[in] rc_memtab
*  pointer to rc mem records
*
* @param[in] num_mem_recs
*  Number of memory records

* @returns none
*
* @remarks
*
*******************************************************************************
*/
void ih264e_map_itt_mem_rec_to_rc_mem_rec(iv_mem_rec_t *ps_mem,
                                          itt_memtab_t *rc_memtab,
                                          UWORD32 num_mem_recs)
{
    UWORD32 i;

    for (i = 0; i < num_mem_recs; i++)
    {
        rc_memtab->i4_alignment = ps_mem->u4_mem_alignment;
        rc_memtab->u4_size = ps_mem->u4_mem_size;
        rc_memtab->pv_base = ps_mem->pv_base;

        /* only DDR memory is available */
        rc_memtab->e_mem_region = DDR;
        rc_memtab->e_usage = PERSISTENT;

        rc_memtab++;
        ps_mem++;
    }
}

/**
******************************************************************************
*
* @brief Get/Init memtabs for rate control
*
* @par   Description
*  This routine is used to Get/init memtabs for rate control
*
* @param[in] pv_rate_control
*  pointer to rate control context (handle)
*
* @param[in] ps_mem
*  pointer to encoder lib mem records
*
* @param[in] e_func_type
*  enum that dictates fill memory records or Init memory records
*
* @return total number of mem records
*
******************************************************************************
*/
WORD32 ih264e_get_rate_control_mem_tab(void *pv_rate_control,
                                       iv_mem_rec_t  *ps_mem,
                                       ITT_FUNC_TYPE_E e_func_type)
{
    itt_memtab_t as_itt_memtab[NUM_RC_MEMTABS];
    WORD32 i4_num_memtab = 0, j = 0;
    void *refptr2[4];
    void **refptr1[4];
    rate_control_ctxt_t *ps_rate_control = pv_rate_control;

    for (j = 0; j < 4; j++)
        refptr1[j] = &(refptr2[j]);

    j = 0;

    if (e_func_type == USE_BASE || e_func_type == FILL_BASE)
    {
        refptr1[1] = &ps_rate_control->pps_frame_time;
        refptr1[2] = &ps_rate_control->pps_time_stamp;
        refptr1[3] = &ps_rate_control->pps_pd_frm_rate;
        refptr1[0] = &ps_rate_control->pps_rate_control_api;
    }

    /* Get the total number of memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab((rate_control_api_t **)refptr1[0], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    ih264e_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    /* Fill the memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab((rate_control_api_t **)refptr1[0],as_itt_memtab+j,e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    ih264e_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Frame time Module */
    i4_num_memtab = ih264e_frame_time_get_init_free_memtab((frame_time_t **)refptr1[1], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    ih264e_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    /* Fill the memtabs used by Frame time Module */
    i4_num_memtab = ih264e_frame_time_get_init_free_memtab((frame_time_t **)refptr1[1], as_itt_memtab+j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    ih264e_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Time stamp Module */
    i4_num_memtab = ih264e_time_stamp_get_init_free_memtab((time_stamp_t **)refptr1[2], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    ih264e_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    /* Fill the memtabs used by Time Stamp Module */
    i4_num_memtab = ih264e_time_stamp_get_init_free_memtab((time_stamp_t **)refptr1[2], as_itt_memtab+j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    ih264e_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Frame rate Module */
    i4_num_memtab = ih264e_pd_frm_rate_get_init_free_memtab((pd_frm_rate_t **)refptr1[3], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    ih264e_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    /* Fill the memtabs used by Frame Rate Module */
    i4_num_memtab = ih264e_pd_frm_rate_get_init_free_memtab((pd_frm_rate_t **)refptr1[3], as_itt_memtab+j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    ih264e_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab+j, i4_num_memtab);
    j += i4_num_memtab;

    return j; /* Total MemTabs Needed by Rate Control Module */
}
