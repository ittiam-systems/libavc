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
******************************************************************************
* @file
*  isvce_rc_mem_interface.c
*
* @brief
*  This file contains api function definitions for rate control memtabs
*
* @author
*  ittiam
*
* List of Functions
*  - fill_memtab()
*  - use_or_fill_base()
*  - isvce_map_rc_mem_recs_to_itt_api()
*  - isvce_map_itt_mem_rec_to_rc_mem_rec()
*  - isvce_get_rate_control_mem_tab()
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>

/* User Include Files */
#include "ih264e_config.h"
#include "ih264_typedefs.h"
#include "ih264_size_defs.h"
#include "iv2.h"
#include "ive2.h"
#include "ime_distortion_metrics.h"
#include "ime_defs.h"
#include "ime_structs.h"
#include "isvce.h"
#include "ithread.h"
#include "isvc_defs.h"
#include "ih264_debug.h"
#include "isvc_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_error.h"
#include "isvc_structs.h"
#include "isvc_trans_quant_itrans_iquant.h"
#include "isvc_inter_pred_filters.h"
#include "isvc_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "isvc_common_tables.h"
#include "ih264_list.h"
#include "isvc_cabac_tables.h"
#include "ih264e_error.h"
#include "isvce_defs.h"
#include "ih264e_bitstream.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "isvce_rate_control.h"
#include "isvce_cabac_structs.h"
#include "isvce_structs.h"
#include "ih264e_master.h"
#include "ih264_buf_mgr.h"
#include "ih264_dpb_mgr.h"
#include "isvce_utils.h"
#include "ih264e_platform_macros.h"
#include "ih264_cavlc_tables.h"
#include "ih264e_statistics.h"
#include "ih264e_trace.h"
#include "ih264e_fmt_conv.h"
#include "isvce_cavlc.h"
#include "ih264e_rc_mem_interface.h"
#include "isvce_rc_mem_interface.h"
#include "ih264e_time_stamp.h"
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
#include "irc_rate_control_api.h"
#include "irc_rate_control_api_structs.h"
#include "ih264e_modify_frm_rate.h"

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
static void isvce_map_rc_mem_recs_to_itt_api(iv_mem_rec_t *ps_mem, itt_memtab_t *rc_memtab,
                                             UWORD32 num_mem_recs)
{
    UWORD32 j;
    UWORD32 Size, align;

    for(j = 0; j < num_mem_recs; j++)
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
static void isvce_map_itt_mem_rec_to_rc_mem_rec(iv_mem_rec_t *ps_mem, itt_memtab_t *rc_memtab,
                                                UWORD32 num_mem_recs)
{
    UWORD32 i;

    for(i = 0; i < num_mem_recs; i++)
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
* @brief Get memtabs for rate control
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
WORD32 isvce_get_rate_control_mem_tab(void *pv_rate_control, iv_mem_rec_t *ps_mem,
                                      ITT_FUNC_TYPE_E e_func_type)
{
    itt_memtab_t as_itt_memtab[NUM_SVCE_RC_MEMTABS];
    WORD32 i4_num_memtab = 0, j = 0;
    void *refptr2[RC_MEM_CNT];
    void **refptr1[RC_MEM_CNT];
    isvce_rate_control_ctxt_t *ps_rate_control = pv_rate_control;

    for(j = 0; j < RC_MEM_CNT; j++) refptr1[j] = &(refptr2[j]);

    j = 0;

    if(e_func_type == USE_BASE || e_func_type == FILL_BASE)
    {
        refptr1[RC_MEM_FRAME_TIME] = &ps_rate_control->pps_frame_time;
        refptr1[RC_MEM_TIME_STAMP] = &ps_rate_control->pps_time_stamp;
        refptr1[RC_MEM_FRAME_RATE] = &ps_rate_control->pps_pd_frm_rate;
        refptr1[RC_MEM_API_L0] = &ps_rate_control->apps_rate_control_api[0];
        refptr1[RC_MEM_API_L1] = &ps_rate_control->apps_rate_control_api[1];
        refptr1[RC_MEM_API_L2] = &ps_rate_control->apps_rate_control_api[2];
    }

    /* Get the total number of memtabs used by Frame time Module */
    i4_num_memtab = ih264e_frame_time_get_init_free_memtab(
        (frame_time_t **) refptr1[RC_MEM_FRAME_TIME], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Frame time Module */
    i4_num_memtab = ih264e_frame_time_get_init_free_memtab(
        (frame_time_t **) refptr1[RC_MEM_FRAME_TIME], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Time stamp Module */
    i4_num_memtab = ih264e_time_stamp_get_init_free_memtab(
        (time_stamp_t **) refptr1[RC_MEM_TIME_STAMP], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Time Stamp Module */
    i4_num_memtab = ih264e_time_stamp_get_init_free_memtab(
        (time_stamp_t **) refptr1[RC_MEM_TIME_STAMP], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Frame rate Module */
    i4_num_memtab = ih264e_pd_frm_rate_get_init_free_memtab(
        (pd_frm_rate_t **) refptr1[RC_MEM_FRAME_RATE], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Frame Rate Module */
    i4_num_memtab = ih264e_pd_frm_rate_get_init_free_memtab(
        (pd_frm_rate_t **) refptr1[RC_MEM_FRAME_RATE], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L0], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L0], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L1], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L1], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    /* Get the total number of memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L2], NULL, GET_NUM_MEMTAB);
    /* Few extra steps during init */
    isvce_map_itt_mem_rec_to_rc_mem_rec((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    /* Fill the memtabs used by Rate Controller */
    i4_num_memtab = irc_rate_control_num_fill_use_free_memtab(
        (rate_control_api_t **) refptr1[RC_MEM_API_L2], as_itt_memtab + j, e_func_type);
    /* Mapping ittiam memtabs to App. memtabs */
    isvce_map_rc_mem_recs_to_itt_api((&ps_mem[j]), as_itt_memtab + j, i4_num_memtab);
    j += i4_num_memtab;

    return j; /* Total MemTabs Needed by Rate Control Module */
}
