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
/*****************************************************************************/
/*                                                                           */
/*  File Name         : ithread.c                                            */
/*                                                                           */
/*  Description       : Contains abstraction for threads, mutex and semaphores*/
/*                                                                           */
/*  List of Functions :                                                      */
/*                                                                           */
/*  Issues / Problems : None                                                 */
/*                                                                           */
/*  Revision History  :                                                      */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         07 09 2012   Harish          Initial Version                      */
/*****************************************************************************/
/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/
#include <string.h>
#include "ih264_typedefs.h"



#include "ithread.h"
#include <sys/types.h>


#define UNUSED(x) ((void)(x))

//#define PTHREAD_AFFINITY
//#define SYSCALL_AFFINITY

#ifdef PTHREAD_AFFINITY
#define _GNU_SOURCE
#define __USE_GNU
#endif

#include <pthread.h>
#include <sched.h>
#include <semaphore.h>
#include <unistd.h>
#ifdef PTHREAD_AFFINITY
#include <sys/prctl.h>
#endif


UWORD32 ithread_get_handle_size(void)
{
    return sizeof(pthread_t);
}

UWORD32 ithread_get_mutex_lock_size(void)
{
    return sizeof(pthread_mutex_t);
}


WORD32 ithread_create(void *thread_handle, void *attribute, void *strt, void *argument)
{
    UNUSED(attribute);
    return pthread_create((pthread_t *)thread_handle, NULL,(void *(*)(void *)) strt, argument);
}

WORD32 ithread_join(void *thread_handle, void ** val_ptr)
{
    pthread_t *pthread_handle   = (pthread_t *)thread_handle;
    UNUSED(val_ptr);
    return pthread_join(*pthread_handle, NULL);
}

WORD32 ithread_get_mutex_struct_size(void)
{
    return(sizeof(pthread_mutex_t));
}
WORD32 ithread_mutex_init(void *mutex)
{
    return pthread_mutex_init((pthread_mutex_t *) mutex, NULL);
}

WORD32 ithread_mutex_destroy(void *mutex)
{
    return pthread_mutex_destroy((pthread_mutex_t *) mutex);
}

WORD32 ithread_mutex_lock(void *mutex)
{
    return pthread_mutex_lock((pthread_mutex_t *)mutex);
}

WORD32 ithread_mutex_unlock(void *mutex)
{
    return pthread_mutex_unlock((pthread_mutex_t *)mutex);
}

void ithread_yield(void)
{
    sched_yield();
}

void ithread_sleep(UWORD32 u4_time)
{
    usleep(u4_time * 1000 * 1000);
}

void ithread_msleep(UWORD32 u4_time_ms)
{
    usleep(u4_time_ms * 1000);
}

void ithread_usleep(UWORD32 u4_time_us)
{
    usleep(u4_time_us);
}

UWORD32 ithread_get_sem_struct_size(void)
{
    return(sizeof(sem_t));
}


WORD32 ithread_sem_init(void *sem,WORD32 pshared,UWORD32 value)
{
    return sem_init((sem_t *)sem,pshared,value);
}

WORD32 ithread_sem_post(void *sem)
{
    return sem_post((sem_t *)sem);
}


WORD32 ithread_sem_wait(void *sem)
{
    return sem_wait((sem_t *)sem);
}


WORD32 ithread_sem_destroy(void *sem)
{
    return sem_destroy((sem_t *)sem);
}

void ithread_set_name(CHAR *pc_thread_name)
{

#ifndef WIN32
#ifndef QNX
#ifndef IOS
    UNUSED(pc_thread_name);
//prctl(PR_SET_NAME, (unsigned long)pu1_thread_name, 0, 0, 0);
#endif
#endif
#endif

}
WORD32 ithread_set_affinity(WORD32 core_id)
{
#ifdef PTHREAD_AFFINITY
    cpu_set_t cpuset;
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    pthread_t cur_thread = pthread_self();

    if (core_id >= num_cores)
        return -1;

    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    return pthread_setaffinity_np(cur_thread, sizeof(cpu_set_t), &cpuset);

#elif SYSCALL_AFFINITY
    WORD32 i4_sys_res;
    UNUSED(core_id);

    pid_t pid = gettid();


    i4_sys_res = syscall(__NR_sched_setaffinity, pid, sizeof(i4_mask), &i4_mask);
    if (i4_sys_res)
    {
        //WORD32 err;
        //err = errno;
        //perror("Error in setaffinity syscall PERROR : ");
        //LOG_ERROR("Error in the syscall setaffinity: mask=0x%x err=0x%x", i4_mask, i4_sys_res);
        return -1;
    }
#else
    UNUSED(core_id);
#endif
    return 1;

}

WORD32 ithread_get_cond_struct_size(void)
{
    return (sizeof(pthread_cond_t));
}

WORD32 ithread_cond_init(void *cond)
{
    return pthread_cond_init((pthread_cond_t *)cond, NULL);
}

WORD32 ithread_cond_destroy(void *cond)
{
    return pthread_cond_destroy((pthread_cond_t *)cond);
}

WORD32 ithread_cond_wait(void *cond, void *mutex)
{
    return pthread_cond_wait((pthread_cond_t *)cond, (pthread_mutex_t *)mutex);
}

WORD32 ithread_cond_signal(void *cond)
{
    return pthread_cond_signal((pthread_cond_t *)cond);
}

UWORD32 ithread_get_cond_size(void)
{
    return sizeof(pthread_cond_t);
}

WORD32  ithread_cond_broadcast(void *cond)
{
    return pthread_cond_broadcast((pthread_cond_t *)cond);
}
