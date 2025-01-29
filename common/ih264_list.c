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
*  ih264_list.c
*
* @brief
*  Contains functions for buf queue
*
* @author
*  ittiam
*
* @par List of Functions:
*  - ih264_list_size
*  - ih264_list_lock
*  - ih264_list_unlock
*  - ih264_list_yield
*  - ih264_list_free
*  - ih264_list_init
*  - ih264_list_reset
*  - ih264_list_deinit
*  - ih264_list_terminate
*  - ih264_list_queue
*  - ih264_list_dequeue
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
#include "ih264_typedefs.h"
#include "ithread.h"
#include "ih264_debug.h"
#include "ih264_macros.h"
#include "ih264_error.h"
#include "ih264_list.h"
#include "ih264_platform_macros.h"

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
*******************************************************************************
*
* @brief Returns size for job queue context.
*
* @par Description
*  Returns size for job queue context.
*
* @param[in] num_entries
*  max number of jobs that can be queued
*
* @param[in] entry_size
*  memory needed for a single job
*
* @returns Size of the job queue context
*
* @remarks
*
*******************************************************************************
*/
WORD32 ih264_list_size(WORD32 num_entries, WORD32 entry_size)
{
    WORD32 size;
    WORD32 clz;

    size = sizeof(list_t);
    size += ithread_get_mutex_lock_size();

    /* Use next power of two number of entries*/
    clz = CLZ(num_entries);
    num_entries = 1 << (32 - clz);

    size += num_entries * entry_size;
    return size;
}

/**
*******************************************************************************
*
* @brief Locks the list context
*
* @par Description
*  Locks the list context by calling ithread_mutex_lock()
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @returns IH264_FAIL if mutex lock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_lock(list_t *ps_list)
{
    WORD32 retval;

    retval = ithread_mutex_lock(ps_list->pv_mutex);
    if(retval)
        return IH264_FAIL;
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief Unlocks the list context
*
* @par Description
*  Unlocks the list context by calling ithread_mutex_unlock()
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @returns IH264_FAIL if mutex unlock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_unlock(list_t *ps_list)
{
    WORD32 retval;

    retval = ithread_mutex_unlock(ps_list->pv_mutex);
    if(retval)
        return IH264_FAIL;
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief Yields the thread
*
* @par Description
*  Unlocks the list context by calling  ih264_list_unlock(), ithread_yield()
*  and then ih264_list_lock(). List is unlocked before to ensure its
*  access by other threads. If unlock is not done before calling yield then
*  no other thread can access the list functions and update list.
*
* @param[in] ps_list
*  pointer to Job Queue context
*
* @returns IH264_FAIL if mutex lock unlock or yield fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_yield(list_t *ps_list)
{
    IH264_ERROR_T ret;

    ret = ih264_list_unlock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    ithread_yield();

    if(ps_list->i4_yield_interval_us > 0)
        ithread_usleep(ps_list->i4_yield_interval_us);

    ret = ih264_list_lock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief free the list context
*
* @par Description
*  Frees the list context
*
* @param[in] ps_list
*  pointer to Job Queue context
*
* @returns IH264_FAIL if mutex desttroy fails else IH264_SUCCESS
*
* @remarks
*  Since it will be called only once by master thread this is not thread safe.
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_free(list_t *ps_list)
{
    WORD32 ret;

    ret = ithread_mutex_destroy(ps_list->pv_mutex);
    if(0 == ret)
        return IH264_SUCCESS;
    return IH264_FAIL;
}

/**
*******************************************************************************
*
* @brief Initialize the buf queue
*
* @par Description
*  Initializes the list context and sets write and read pointers to start of
*  buf queue buffer
*
* @param[in] pv_buf
*  Memory for job queue context
*
* @param[in] buf_size
*  Size of the total memory allocated
*
* @param[in] num_entries
*  max number of jobs that can be queued
*
* @param[in] entry_size
*  memory needed for a single job
*
* @param[in] yield_interval_us
*  Thread sleep duration
*
* @returns Pointer to job queue context
*
* @remarks
*  Since it will be called only once by master thread this is not thread safe.
*
*******************************************************************************
*/
void* ih264_list_init(void *pv_buf,
                      WORD32 buf_size,
                      WORD32 num_entries,
                      WORD32 entry_size,
                      WORD32 yield_interval_us)
{
    list_t *ps_list = (list_t *)pv_buf;
    UWORD8 *pu1_buf = (UWORD8 *)pv_buf;

    pu1_buf += sizeof(list_t);
    buf_size -= sizeof(list_t);

    ps_list->pv_mutex = pu1_buf;
    pu1_buf += ithread_get_mutex_lock_size();
    buf_size -= ithread_get_mutex_lock_size();

    if (buf_size <= 0)
      return NULL;

    ithread_mutex_init(ps_list->pv_mutex);

    /* Ensure num_entries is power of two */
    ASSERT(0 == (num_entries & (num_entries - 1)));

    /* Ensure remaining buffer is large enough to hold given number of entries */
    ASSERT((num_entries * entry_size) <= buf_size);

    ps_list->pv_buf_base = pu1_buf;
    ps_list->i4_terminate = 0;
    ps_list->i4_entry_size = entry_size;
    ps_list->i4_buf_rd_idx = 0;
    ps_list->i4_buf_wr_idx = 0;
    ps_list->i4_log2_buf_max_idx = 32 - CLZ(num_entries);
    ps_list->i4_buf_max_idx = num_entries;
    ps_list->i4_yield_interval_us = yield_interval_us;

    return ps_list;
}

/**
*******************************************************************************
*
* @brief Resets the list context
*
* @par Description
*  Resets the list context by initializing buf queue context elements
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @returns IH264_FAIL if lock unlock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_reset(list_t *ps_list)
{
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_list_lock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    ps_list->i4_terminate  = 0;
    ps_list->i4_buf_rd_idx = 0;
    ps_list->i4_buf_wr_idx = 0;

    ret = ih264_list_unlock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return ret;
}

/**
*******************************************************************************
*
* @brief De-initializes the list context
*
* @par Description
*  De-initializes the list context by calling ih264_list_reset() and then
*  destroying the mutex created
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @returns IH264_FAIL if lock unlock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_deinit(list_t *ps_list)
{
    WORD32 retval;
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_list_reset(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    retval = ithread_mutex_destroy(ps_list->pv_mutex);
    if(retval)
        return IH264_FAIL;
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief Terminates the list
*
* @par Description
*  Terminates the list by setting a flag in context.
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @returns IH264_FAIL if lock unlock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_terminate(list_t *ps_list)
{
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_list_lock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    ps_list->i4_terminate = 1;

    ret = ih264_list_unlock(ps_list);
    RETURN_IF((ret != IH264_SUCCESS), ret);
    return ret;
}

/**
*******************************************************************************
*
* @brief Adds a job to the queue
*
* @par Description
*  Adds a buffer to the queue and updates write address to next location.
*
* @param[in] ps_list
*  Pointer to job queue context
*
* @param[in] pv_buf
*  Pointer to the location that contains details of the job to be added
*
* @param[in] blocking
*  To signal if the write is blocking or non-blocking.
*
* @returns IH264_SUCCESS on success and IH264_FAIL on fail
*
* @remarks
*  Job Queue buffer is assumed to be allocated to handle worst case number of
*  buffers. Wrap around is not supported
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_queue(list_t *ps_list, void *pv_buf, WORD32 blocking)
{
    IH264_ERROR_T ret = IH264_SUCCESS;
    IH264_ERROR_T rettmp;
    WORD32 diff;
    void *pv_buf_wr;
    volatile WORD32 *pi4_wr_idx, *pi4_rd_idx;
    WORD32 buf_size = ps_list->i4_entry_size;


    rettmp = ih264_list_lock(ps_list);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    while(1)
    {
        /* Ensure wr idx does not go beyond rd idx by more than number of entries
         */
        pi4_wr_idx = &ps_list->i4_buf_wr_idx;
        pi4_rd_idx = &ps_list->i4_buf_rd_idx;
        diff = *pi4_wr_idx - *pi4_rd_idx;

        if(diff < ps_list->i4_buf_max_idx)
        {
            WORD32 wr_idx;
            wr_idx = ps_list->i4_buf_wr_idx & (ps_list->i4_buf_max_idx - 1);
            pv_buf_wr = (UWORD8 *)ps_list->pv_buf_base + wr_idx * buf_size;

            memcpy(pv_buf_wr, pv_buf, buf_size);
            ps_list->i4_buf_wr_idx++;
            break;
        }
        else
        {
            /* wr is ahead, so wait for rd to consume */
            if(blocking)
            {
                ih264_list_yield(ps_list);
            }
            else
            {
                ret = IH264_FAIL;
                break;
            }
        }
    }
    ps_list->i4_terminate = 0;

    rettmp = ih264_list_unlock(ps_list);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    return ret;
}

/**
*******************************************************************************
*
* @brief Gets next job from the job queue
*
* @par   Description
*  Gets next job from the job queue and updates rd address to next location.
*  If it is a blocking call and if there is no new buf then this functions
*  unlocks the mutex and calls yield and then locks it back and continues
*  till a buf is available or terminate is set
*
* @param[in] ps_list
*  Pointer to Job Queue context
*
* @param[out] pv_buf
*  Pointer to the location that contains details of the buf to be written
*
* @param[in] blocking
*  To signal if the read is blocking or non-blocking.
*
* @returns
*
* @remarks
*  Job Queue buffer is assumed to be allocated to handle worst case number of
*  buffers. Wrap around is not supported
*
*******************************************************************************
*/
IH264_ERROR_T ih264_list_dequeue(list_t *ps_list, void *pv_buf, WORD32 blocking)
{
    IH264_ERROR_T ret = IH264_SUCCESS;
    IH264_ERROR_T rettmp;
    WORD32 buf_size = ps_list->i4_entry_size;
    WORD32 diff;
    void *pv_buf_rd;
    volatile WORD32 *pi4_wr_idx, *pi4_rd_idx;

    rettmp = ih264_list_lock(ps_list);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    while(1)
    {
        /* Ensure wr idx is ahead of rd idx and
         * wr idx does not go beyond rd idx by more than number of entries
         */
        pi4_wr_idx = &ps_list->i4_buf_wr_idx;
        pi4_rd_idx = &ps_list->i4_buf_rd_idx;
        diff = *pi4_wr_idx - *pi4_rd_idx;

        if(diff > 0)
        {
            WORD32 rd_idx;
            rd_idx = ps_list->i4_buf_rd_idx & (ps_list->i4_buf_max_idx - 1);
            pv_buf_rd = (UWORD8 *)ps_list->pv_buf_base + rd_idx * buf_size;

            memcpy(pv_buf, pv_buf_rd, buf_size);
            ps_list->i4_buf_rd_idx++;
            break;
        }
        else
        {
            /* If terminate is signaled then break */
            if(ps_list->i4_terminate)
            {
                ret = IH264_FAIL;
                break;
            }
            /* wr is ahead, so wait for rd to consume */
            if(blocking)
            {
                ih264_list_yield(ps_list);
            }
            else
            {
                ret = IH264_FAIL;
                break;
            }
        }
    }

    rettmp = ih264_list_unlock(ps_list);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    return ret;
}

#ifdef KEEP_THREADS_ACTIVE
/**
*******************************************************************************
*
* @brief
*   Gets the number of jobs
*
* @par   Description
*   Gets the number of jobs to be processed in job queue context.
*
* @param[in] ps_list
*   Job Queue context
*
* @returns 0 if lock unlock fails else number of jobs to be processed
*
* @remarks
*
*******************************************************************************
*/
WORD32 ih264_get_job_count_in_list(list_t *ps_list)
{
    WORD32 jobs = 0;
    RETURN_IF((ih264_list_lock(ps_list) != IH264_SUCCESS), 0);
    jobs = ps_list->i4_buf_wr_idx - ps_list->i4_buf_rd_idx;
    RETURN_IF((ih264_list_unlock(ps_list) != IH264_SUCCESS), 0);
    return jobs;
}
#endif /* KEEP_THREADS_ACTIVE */
