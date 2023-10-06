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
*  ithread.h
*
* @brief
*  Wrapper functions for thread operations. Library uses these functions during
*  multi-thread execution. These calls can be mapped to POSIX, win threads, ...
*  basing on platform
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _ITHREAD_H_
#define _ITHREAD_H_

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/
UWORD32 ithread_get_handle_size(void);

UWORD32 ithread_get_mutex_lock_size(void);

WORD32  ithread_create(void *thread_handle, void *attribute, void *strt, void *argument);

void    ithread_exit(void *val_ptr);

WORD32  ithread_join(void *thread_id, void ** val_ptr);

WORD32  ithread_get_mutex_struct_size(void);

WORD32  ithread_mutex_init(void *mutex);

WORD32  ithread_mutex_destroy(void *mutex);

WORD32  ithread_mutex_lock(void *mutex);

WORD32  ithread_mutex_unlock(void *mutex);

void    ithread_yield(void);

void    ithread_sleep(UWORD32 u4_time);

void    ithread_msleep(UWORD32 u4_time_ms);

void    ithread_usleep(UWORD32 u4_time_us);

UWORD32 ithread_get_sem_struct_size(void);

WORD32  ithread_sem_init(void *sem,WORD32 pshared,UWORD32 value);

WORD32  ithread_sem_post(void *sem);

WORD32  ithread_sem_wait(void *sem);

WORD32  ithread_sem_destroy(void *sem);

WORD32  ithread_set_affinity(WORD32 core_id);

void    ithread_set_name(CHAR *pc_thread_name);

WORD32  ithread_get_cond_struct_size(void);

WORD32  ithread_cond_init(void *cond);

WORD32  ithread_cond_destroy(void *cond);

WORD32  ithread_cond_wait(void *cond, void *mutex);

WORD32  ithread_cond_signal(void *cond);

#endif /* _ITHREAD_H_ */
