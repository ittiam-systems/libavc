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
*  irc_mem_req_and_acq.h
*
* @brief
*  This file contains interface definitions for allocating rate control
*  memtabs
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _RC_MEM_REQ_AND_ACQ_H_
#define _RC_MEM_REQ_AND_ACQ_H_

/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/

#define FILL_MEMTAB(m_pv_mem_rec, m_j, m_mem_size, m_align, m_type)      \
{                                                                        \
    m_pv_mem_rec[m_j].u4_size = sizeof(iv_mem_rec_t);                    \
    m_pv_mem_rec[m_j].u4_mem_size = m_mem_size;                          \
    m_pv_mem_rec[m_j].u4_mem_alignment = m_align;                        \
    m_pv_mem_rec[m_j].e_mem_type = m_type;                               \
}

/*****************************************************************************/
/* Enums                                                                     */
/*****************************************************************************/
typedef enum
{
    ALIGN_BYTE = 1,
    ALIGN_WORD16 = 2,
    ALIGN_WORD32 = 4,
    ALIGN_WORD64 = 8,
    ALIGN_128_BYTE = 128
}ITT_MEM_ALIGNMENT_TYPE_E;

typedef enum
{
    SCRATCH = 0,
    PERSISTENT = 1,
    WRITEONCE  = 2
}ITT_MEM_USAGE_TYPE_E;

typedef enum
{
    L1D = 0,
    SL2 = 1,
    DDR = 3
}ITT_MEM_REGION_E;

typedef enum
{
    GET_NUM_MEMTAB = 0,
    FILL_MEMTAB = 1,
    USE_BASE = 2,
    FILL_BASE =3
}ITT_FUNC_TYPE_E;

/*****************************************************************************/
/* Structures                                                                */
/*****************************************************************************/

typedef struct
{
    /* Size in bytes */
    UWORD32 u4_size;

    /* Alignment in bytes */
    WORD32 i4_alignment;

    /* decides which memory region to be placed */
    ITT_MEM_REGION_E e_mem_region;

    /* memory is scratch or persistent */
    ITT_MEM_USAGE_TYPE_E e_usage;

    /* Base pointer for allocated memory */
    void *pv_base;
} itt_memtab_t;

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void fill_memtab(itt_memtab_t *ps_mem_tab,
                 WORD32 u4_size,
                 WORD32 i4_alignment,
                 ITT_MEM_USAGE_TYPE_E e_usage,
                 ITT_MEM_REGION_E e_mem_region);

WORD32 use_or_fill_base(itt_memtab_t *ps_mem_tab,
                        void **ptr_to_be_filled,
                        ITT_FUNC_TYPE_E e_func_type);


#endif // _RC_MEM_REQ_AND_ACQ_H_

