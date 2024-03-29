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
*  ih264_buf_mgr.c
*
* @brief
*  Contains function definitions for buffer management
*
* @author
*  ittiam
*
* @par List of Functions:
*  - ih264_buf_mgr_size
*  - ih264_buf_mgr_lock
*  - ih264_buf_mgr_unlock
*  - ih264_buf_mgr_yield
*  - ih264_buf_mgr_free
*  - ih264_buf_mgr_init
*  - ih264_buf_mgr_add
*  - ih264_buf_mgr_get_next_free
*  - ih264_buf_mgr_check_free
*  - ih264_buf_mgr_release
*  - ih264_buf_mgr_set_status
*  - ih264_buf_mgr_get_status
*  - ih264_buf_mgr_get_buf
*  - ih264_buf_mgr_get_bufid
*  - ih264_buf_mgr_get_num_active_buf
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
#include <stdlib.h>

/* User Include Files */
#include "ih264_typedefs.h"
#include "ithread.h"
#include "ih264_macros.h"
#include "ih264_error.h"
#include "ih264_defs.h"
#include "ih264_buf_mgr.h"


/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
*******************************************************************************
*
* @brief Returns size for buffer queue context. Does not include the size of the
*  actual buffers that are stored here.
*
* @par Description
*  Returns size for buf queue context. Does not include the size of the
*  actual buffers that are stored here.
*
* @returns Size of the buf queue context
*
* @remarks none
*
*******************************************************************************
*/
WORD32 ih264_buf_mgr_size(void)
{
    WORD32 size;

    size = sizeof(buf_mgr_t);
    size += ithread_get_mutex_lock_size();

    return size;
}

/**
*******************************************************************************
*
* @brief Locks the buffer manager context
*
* @par Description
*  Locks the buffer manager context by calling ithread_mutex_lock()
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns IH264_FAIL if mutex lock fails else IH264_SUCCESS
*
* @remarks none
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_lock(buf_mgr_t *ps_buf_mgr)
{
    WORD32 retval;

    retval = ithread_mutex_lock(ps_buf_mgr->pv_mutex);
    if(retval)
        return IH264_FAIL;
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief Unlocks the buffer manager context
*
* @par Description
*  Unlocks the buffer manager context by calling ithread_mutex_unlock()
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns IH264_FAIL if mutex unlock fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_unlock(buf_mgr_t *ps_buf_mgr)
{
    WORD32 retval;

    retval = ithread_mutex_unlock(ps_buf_mgr->pv_mutex);
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
*  Unlocks the buf_mgr context by calling ih264_buf_mgr_unlock(),
*  ithread_yield() and then ih264_buf_mgr_lock(). buf_mgr is unlocked before to
*  ensure the buf_mgr can be accessed by other threads. If unlock is not done
*  before calling yield then no other thread can access the buf_mgr functions
*  and update buf_mgr.
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns IH264_FAIL if mutex lock unlock or yield fails else IH264_SUCCESS
*
* @remarks
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_yield(buf_mgr_t *ps_buf_mgr)
{
    IH264_ERROR_T ret;

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    //ithread_usleep(10);
    ithread_yield();

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);
    return IH264_SUCCESS;
}

/**
*******************************************************************************
*
* @brief frees the buffer manager context
*
* @par Description
*  Frees the buffer manager context
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns IH264_FAIL if mutex destory fails else IH264_SUCCESS
*
* @remarks
*  Since it will be called only once by master thread this is not thread safe.
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_free(buf_mgr_t *ps_buf_mgr)
{
    WORD32 ret;

    ret = ithread_mutex_destroy(ps_buf_mgr->pv_mutex);
    if(0 == ret)
        return IH264_SUCCESS;
    return IH264_FAIL;
}

/**
*******************************************************************************
*
* @brief Buffer manager context initialization.
*
* @par Description:
*  Initializes the buffer manager structure
*
* @param[in] pv_buf
*  Pointer to the buffer manager
*
* @returns none
*
* @remarks none
*
*******************************************************************************
*/
void* ih264_buf_mgr_init(void *pv_buf)
{
    WORD32 id;
    UWORD8 *pu1_buf = (UWORD8 *)pv_buf;
    buf_mgr_t *ps_buf_mgr = (buf_mgr_t *)pu1_buf;

    pu1_buf += sizeof(buf_mgr_t);
    ps_buf_mgr->pv_mutex = pu1_buf;

    pu1_buf += ithread_get_mutex_lock_size();
    ithread_mutex_init(ps_buf_mgr->pv_mutex);

    ps_buf_mgr->i4_max_buf_cnt = BUF_MGR_MAX_CNT;
    ps_buf_mgr->i4_active_buf_cnt = 0;
    for(id = 0; id < BUF_MGR_MAX_CNT; id++)
    {
        ps_buf_mgr->au4_status[id] = 0;
        ps_buf_mgr->apv_ptr[id] = NULL;
    }

    return ps_buf_mgr;
}

/**
*******************************************************************************
*
* @brief Buffer manager reset function.
*
* @par Description:
*  resets the buffer manager structure
*
* @param[in] pv_buf
*  Pointer to the buffer manager
*
* @returns none
*
* @remarks none
*
*******************************************************************************
*/
void ih264_buf_mgr_reset(void *pv_buf_mgr)
{
    WORD32 id;
    buf_mgr_t *ps_buf_mgr;

    ps_buf_mgr = (buf_mgr_t *) pv_buf_mgr;

    for(id = 0; id < BUF_MGR_MAX_CNT; id++)
    {
        ps_buf_mgr->au4_status[id] = 0;
    }

    return;
}

/**
*******************************************************************************
*
* @brief Adds and increments the buffer and buffer count.
*
* @par Description:
*  Adds a buffer to the buffer manager if it is not already present and
*  increments the active buffer count
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] pv_ptr
*  Pointer to the buffer to be added
*
* @param[in] buf_id
*  buffer id
*
* @returns  Returns 0 on success, -1 otherwise
*
* @remarks
*  none
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_add(buf_mgr_t *ps_buf_mgr,
                                void *pv_ptr,
                                WORD32 buf_id)
{
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    /* Check if buffer ID is within allowed range */
    if(buf_id >= ps_buf_mgr->i4_max_buf_cnt)
    {
        ret = ih264_buf_mgr_unlock(ps_buf_mgr);
        RETURN_IF((ret != IH264_SUCCESS), ret);
        return IH264_FAIL;
    }

    /* Check if the current ID is being used to hold some other buffer */
    if((ps_buf_mgr->apv_ptr[buf_id] != NULL) &&
       (ps_buf_mgr->apv_ptr[buf_id] != pv_ptr))
    {
        ret = ih264_buf_mgr_unlock(ps_buf_mgr);
        RETURN_IF((ret != IH264_SUCCESS), ret);
        return IH264_FAIL;
    }
    ps_buf_mgr->apv_ptr[buf_id] = pv_ptr;
    ps_buf_mgr->i4_active_buf_cnt++;

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return ret;
}

/**
*******************************************************************************
*
* @brief Gets the next free buffer.
*
* @par Description:
*  Returns the next free buffer available and sets the corresponding status
*  to Busy
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] pi4_buf_id
*  Pointer to the id of the free buffer
*
* @returns  Pointer to the free buffer
*
* @remarks
*  none
*
*******************************************************************************
*/
void* ih264_buf_mgr_get_next_free(buf_mgr_t *ps_buf_mgr, WORD32 *pi4_buf_id)
{
    WORD32 id;
    void *pv_ret_ptr = NULL;
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), NULL);

    for(id = 0; id < ps_buf_mgr->i4_active_buf_cnt; id++)
    {
        /* Check if the buffer is non-null and status is zero */
        if((ps_buf_mgr->au4_status[id] == 0) && (ps_buf_mgr->apv_ptr[id]))
        {
            *pi4_buf_id = id;
            /* mark buffer as busy before returning */
            ps_buf_mgr->au4_status[id] = 1;
            pv_ret_ptr = ps_buf_mgr->apv_ptr[id];
            break;
        }
    }
    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), NULL);

    return pv_ret_ptr;
}

/**
*******************************************************************************
*
* @brief Checks the buffer manager for free buffers available.
*
* @par Description:
*  Checks if there are any free buffers available
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns  Returns IH264_SUCCESS if available, IH264_FAIL otherwise
*
* @remarks
*  none
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_check_free(buf_mgr_t *ps_buf_mgr)
{
    WORD32 id;
    IH264_ERROR_T ret = IH264_FAIL;
    IH264_ERROR_T rettmp = IH264_SUCCESS;

    rettmp = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    for(id = 0; id < ps_buf_mgr->i4_active_buf_cnt; id++)
    {
        if((ps_buf_mgr->au4_status[id] == 0) &&
           (ps_buf_mgr->apv_ptr[id]))
        {
            ret = IH264_SUCCESS;
            break;
        }
    }
    rettmp = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((rettmp != IH264_SUCCESS), rettmp);

    return ret;
}

/**
*******************************************************************************
*
* @brief Release buffer to buffer manager
*
* @par Description:
*  Clear the status bit basing on the mask. If the buffer is no longer used by
*  the codec then release back to the pool
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] buf_id
*  ID of the buffer status to be released
*
* @param[in] mask
*  Contains the bits that are to be reset
*
* @returns  IH264_SUCCESS if success, IH264_FAIL otherwise
*
* @remarks
*  none
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_release(buf_mgr_t *ps_buf_mgr,
                                    WORD32 buf_id,
                                    UWORD32 mask)
{
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    /* If the given id is pointing to an id which is not yet added */
    if(buf_id >= ps_buf_mgr->i4_active_buf_cnt)
    {
        ret = ih264_buf_mgr_unlock(ps_buf_mgr);
        RETURN_IF((ret != IH264_SUCCESS), ret);
        return IH264_FAIL;
    }

    ps_buf_mgr->au4_status[buf_id] &= ~mask;

    /* If both the REF and DISP are zero, DEC is set to zero */
    if(ps_buf_mgr->au4_status[buf_id] == 1)
    {
        ps_buf_mgr->au4_status[buf_id] = 0;
    }

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return ret;
}

/**
*******************************************************************************
*
* @brief Sets the status bit.
*
* @par Description:
*  sets the status bits that the mask contains (status corresponding to the id)
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] buf_id
*  ID of the buffer whose status needs to be modified
*
* @param[in] mask
*  Contains the bits that are to be set
*
* @returns  IH264_SUCCESS if success, IH264_FAIL otherwise
*
* @remarks
*  none
*
*******************************************************************************
*/
IH264_ERROR_T ih264_buf_mgr_set_status(buf_mgr_t *ps_buf_mgr,
                                       WORD32 buf_id,
                                       UWORD32 mask)
{
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    if(buf_id >= ps_buf_mgr->i4_active_buf_cnt)
    {
        ret = ih264_buf_mgr_unlock(ps_buf_mgr);
        RETURN_IF((ret != IH264_SUCCESS), ret);
        return IH264_FAIL;
    }

    if((ps_buf_mgr->au4_status[buf_id] & mask) != 0)
    {
        ret = ih264_buf_mgr_unlock(ps_buf_mgr);
        RETURN_IF((ret != IH264_SUCCESS), ret);
        return IH264_FAIL;
    }

    ps_buf_mgr->au4_status[buf_id] |= mask;
    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return ret;
}

/**
*******************************************************************************
*
* @brief Returns the status of the buffer.
*
* @par Description:
*  Returns the status of the buffer corresponding to the id
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] buf_id
*  ID of the buffer status required
*
* @returns  Status of the buffer corresponding to the id
*
* @remarks
*  none
*
*******************************************************************************
*/
WORD32 ih264_buf_mgr_get_status(buf_mgr_t *ps_buf_mgr, WORD32 buf_id)
{
    IH264_ERROR_T ret = IH264_SUCCESS;
    UWORD32 status;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    status = ps_buf_mgr->au4_status[buf_id];

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return status;
}

/**
*******************************************************************************
*
* @brief Gets the buffer from the buffer manager
*
* @par Description:
*  Returns the pointer to the buffer corresponding to the id
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] buf_id
*  ID of the buffer required
*
* @returns  Pointer to the buffer required
*
* @remarks
*  none
*
*******************************************************************************
*/
void* ih264_buf_mgr_get_buf(buf_mgr_t *ps_buf_mgr, WORD32 buf_id)
{
    IH264_ERROR_T ret = IH264_SUCCESS;
    void *pv_buf;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), NULL);

    pv_buf = ps_buf_mgr->apv_ptr[buf_id];

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), NULL);

    return pv_buf;
}

/**
*******************************************************************************
*
* @brief Gets the buffer id from the buffer manager if the buffer is added to
*  the buffer manager
*
* @par Description:
*  Returns the buffer id corresponding to the given buffer if it exists
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @param[in] pv_buf
*  Pointer to the buffer
*
* @returns  Buffer id if exists, else -1
*
* @remarks
*  none
*
*******************************************************************************
*/
WORD32 ih264_buf_mgr_get_bufid(buf_mgr_t *ps_buf_mgr, void *pv_buf)
{
    WORD32 id;
    WORD32 buf_id = -1;
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    for(id = 0; id < ps_buf_mgr->i4_active_buf_cnt; id++)
    {
        if(ps_buf_mgr->apv_ptr[id] == pv_buf)
        {
            buf_id = id;
            break;
        }
    }
    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return buf_id;
}

/**
*******************************************************************************
*
* @brief Gets the no.of active buffer
*
* @par Description:
*  Return the number of active buffers in the buffer manager
*
* @param[in] ps_buf_mgr
*  Pointer to the buffer manager
*
* @returns  number of active buffers
*
* @remarks
*  none
*
*******************************************************************************
*/
UWORD32 ih264_buf_mgr_get_num_active_buf(buf_mgr_t *ps_buf_mgr)
{
    UWORD32 u4_buf_cnt = 0;
    IH264_ERROR_T ret = IH264_SUCCESS;

    ret = ih264_buf_mgr_lock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    u4_buf_cnt = ps_buf_mgr->i4_active_buf_cnt;

    ret = ih264_buf_mgr_unlock(ps_buf_mgr);
    RETURN_IF((ret != IH264_SUCCESS), ret);

    return u4_buf_cnt;
}
