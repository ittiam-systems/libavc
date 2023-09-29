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
*  ih264e_modify_frm_rate.h
*
* @brief
*  Handle source frame rate pulldown
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_MODIFY_FRM_RATE_H_
#define _IH264E_MODIFY_FRM_RATE_H_

/*****************************************************************************/
/* Constant Definitions                                                      */
/*****************************************************************************/

#define MAX_NUM_FRAME   120


/*****************************************************************************/
/* Structures                                                                */
/*****************************************************************************/
typedef struct pd_frm_rate_t
{
    /*
     * The input frame rate set in the encoder (per 1000 sec)
     */
    UWORD32 u4_input_frm_rate;

    /*
     * Frame rate of current frame due to pull down
     */
    UWORD32 u4_cur_frm_rate[MAX_NUM_FRAME];

    /*
     * current frame num in the above buffer
     */
    UWORD32 u4_frm_num;

    /*
     * Total number of frames encoded.
     * if greater than input frame rate stays at input frame rate
     */
    UWORD32 u4_tot_frm_encoded;

}pd_frm_rate_t;

typedef struct pd_frm_rate_t *pd_frm_rate_handle;


/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

WORD32 ih264e_pd_frm_rate_get_init_free_memtab(pd_frm_rate_handle *pps_pd_frm_rate,
                                               itt_memtab_t *ps_memtab,
                                               ITT_FUNC_TYPE_E e_func_type);

void ih264e_init_pd_frm_rate(pd_frm_rate_handle ps_pd_frm_rate,
                             UWORD32 u4_input_frm_rate);

void ih264e_update_pd_frm_rate(pd_frm_rate_handle ps_pd_frm_rate,
                               UWORD32 u4_cur_frm_rate);

UWORD32 ih264e_get_pd_avg_frm_rate(pd_frm_rate_handle ps_pd_frm_rate);

#endif /* _IH264E_MODIFY_FRM_RATE_H_ */
