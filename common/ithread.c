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

/*
 * If the end target is bare metal, then there shall be no OS.
 * In this case, the functions ithread_* used inside the h264 encoder library to assist multicore
 * will not longer be functional. To resolve link issues, the functions are re-defined with no body.
 */
#ifndef BAREMETAL


#include "ithread.h"
#include <sys/types.h>


#define UNUSED(x) ((void)(x))

#ifndef X86_MSVC
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

#endif

#if defined(X86_MSVC) || defined (X86_MINGW)

#include <windows.h>
#define SEM_MAX_COUNT       100
#define SEM_INCREMENT_COUNT 1

UWORD32 ithread_get_handle_size(void)
{
    return (sizeof(HANDLE));
}

UWORD32 ithread_get_mutex_lock_size(void)
{
    return (sizeof(HANDLE));
}

WORD32 ithread_create(void *thread_handle, void *attribute, void *strt, void *argument)
{
    HANDLE *ppv_thread_handle;
    HANDLE thread_handle_value;

    UNUSED(attribute);

    if(0 == thread_handle)
        return -1;

    ppv_thread_handle = (HANDLE *)thread_handle;
    thread_handle_value = (void *)CreateThread
            (NULL,                             /* Attributes      */
            1024*128,                          /* Stack i4_size      */
            (LPTHREAD_START_ROUTINE)strt,      /* Thread function */
            argument,                          /* Parameters      */
            0,                                 /* Creation flags  */
            NULL);                             /* Thread ID       */
    *ppv_thread_handle = (HANDLE)thread_handle_value;

    return 0;
}

WORD32 ithread_join(void *thread_handle, void ** val_ptr)
{
    HANDLE *ppv_thread_handle;
    HANDLE thread_handle_value;

    UNUSED(val_ptr);

    if(0 == thread_handle)
        return -1;

    ppv_thread_handle = (HANDLE *)thread_handle;
    thread_handle_value = *ppv_thread_handle;

    if(WAIT_OBJECT_0 == WaitForSingleObject(thread_handle_value, INFINITE))
    {
        CloseHandle(thread_handle_value);
    }

    return 0;
}

void ithread_exit(void *thread_handle)
{
    HANDLE *ppv_thread_handle;
    HANDLE thread_handle_value;
    DWORD thread_exit_code;

    if(0 == thread_handle)
        return;

    ppv_thread_handle = (HANDLE *)thread_handle;
    thread_handle_value = *ppv_thread_handle;
    /* Get exit code for thread. If the return value is 0, means thread is busy */
    if( 0 != GetExitCodeThread(thread_handle_value, &thread_exit_code))
    {
        TerminateThread(thread_handle_value, thread_exit_code);
    }

    return;
}

WORD32 ithread_get_mutex_struct_size(void)
{
    return (sizeof(HANDLE));
}

WORD32 ithread_mutex_init(void *mutex)
{
    HANDLE *ppv_mutex_handle;
    HANDLE mutex_handle_value;

    if(0 == mutex)
        return -1;

    ppv_mutex_handle = (HANDLE *)mutex;
    mutex_handle_value = CreateSemaphore(NULL, 1, 1, NULL);
    *ppv_mutex_handle = mutex_handle_value;
    return 0;
}

WORD32 ithread_mutex_destroy(void *mutex)
{
    HANDLE *ppv_mutex_handle;
    HANDLE mutex_handle_value;

    if(0 == mutex)
        return -1;

    ppv_mutex_handle = (HANDLE *)mutex;
    mutex_handle_value = *ppv_mutex_handle;
    CloseHandle(mutex_handle_value);
    return 0;
}

WORD32 ithread_mutex_lock(void *mutex)
{
    HANDLE *ppv_mutex_handle;
    HANDLE mutex_handle_value;
    DWORD  result = 0;

    if(0 == mutex)
        return -1;

    ppv_mutex_handle = (HANDLE *)mutex;
    mutex_handle_value = *ppv_mutex_handle;
    result = WaitForSingleObject(mutex_handle_value, INFINITE);

    if(WAIT_OBJECT_0 == result)
        return 0;

    return 1;

}

WORD32 ithread_mutex_unlock(void *mutex)
{
    HANDLE *ppv_mutex_handle;
    HANDLE mutex_handle_value;
    DWORD  result = 0;

    if(0 == mutex)
        return -1;

    ppv_mutex_handle = (HANDLE *)mutex;
    mutex_handle_value = *ppv_mutex_handle;
    result = ReleaseSemaphore(mutex_handle_value, 1, NULL);

    if(0 == result)
        return -1;

    return 0;
}

void ithread_yield(void) { }

void ithread_usleep(UWORD32 u4_time_us)
{
    UWORD32 u4_time_ms = u4_time_us / 1000;
    Sleep(u4_time_ms);
}

void ithread_msleep(UWORD32 u4_time_ms)
{
    Sleep(u4_time_ms);
}

void ithread_sleep(UWORD32 u4_time)
{
    UWORD32 u4_time_ms = u4_time * 1000;
    Sleep(u4_time_ms);
}

UWORD32 ithread_get_sem_struct_size(void)
{
    return (sizeof(HANDLE));
}

WORD32 ithread_sem_init(void *sem,WORD32 pshared,UWORD32 value)
{
    HANDLE *sem_handle = (HANDLE *)sem;
    HANDLE sem_handle_value;

    if(0 == sem)
        return -1;

    sem_handle_value = CreateSemaphore(NULL,  /* Security Attribute*/
                         value,  /* Initial count     */
                        SEM_MAX_COUNT,/* Max value         */
                        NULL);        /* Name, not used    */
    *sem_handle = sem_handle_value;
    return 0;
}

WORD32 ithread_sem_post(void *sem)
{
    HANDLE *sem_handle = (HANDLE *)sem;
    HANDLE sem_handle_value;

    if(0 == sem)
        return -1;

    sem_handle_value = *sem_handle;

    /* Post on Semaphore by releasing the lock on mutex */
    if(ReleaseSemaphore(sem_handle_value, SEM_INCREMENT_COUNT, NULL))
        return 0;

    return -1;
}

WORD32 ithread_sem_wait(void *sem)
{
    DWORD          result = 0;
    HANDLE *sem_handle = (HANDLE *)sem;
    HANDLE sem_handle_value;

    if(0 == sem)
        return -1;

    sem_handle_value = *sem_handle;

    /* Wait on Semaphore object infinitly */
    result = WaitForSingleObject(sem_handle_value, INFINITE);

    /* If lock on semaphore is acquired, return SUCCESS */
    if(WAIT_OBJECT_0 == result)
        return 0;

    /* If call timeouts, return FAILURE */
    if(WAIT_TIMEOUT == result)
        return -1;

    return 0;
}

WORD32 ithread_sem_destroy(void *sem)
{
    HANDLE *sem_handle = (HANDLE *)sem;
    HANDLE sem_handle_value;

    if(0 == sem)
        return -1;

    sem_handle_value = *sem_handle;

    if(FALSE == CloseHandle(sem_handle_value) )
    {
        return -1;
    }
    return 0;
}

WORD32 ithread_set_affinity(WORD32 core_id)
{
        return 1;
}

void ithread_set_name(CHAR *pc_thread_name)
{
    return;
}

#else

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
    UNUSED(val_ptr);
    pthread_t *pthread_handle   = (pthread_t *)thread_handle;
    return pthread_join(*pthread_handle, NULL);
}

void ithread_exit(void *val_ptr)
{
    return pthread_exit(val_ptr);
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
#endif

#else

UWORD32 ithread_get_handle_size(void)
{
    return sizeof(int);
}

UWORD32 ithread_get_mutex_lock_size(void)
{
    return sizeof(int);
}

UWORD32 ithread_get_cond_size(void)
{
    return(sizeof(int));
}
WORD32 ithread_create(void *thread_handle, void *attribute, void *strt, void *argument)
{
    return 0;
}

WORD32 ithread_join(void *thread_handle, void ** val_ptr)
{
    return 0;
}

void ithread_exit(void *val_ptr)
{
    return;
}

WORD32 ithread_mutex_init(void *mutex)
{
    return 0;
}

WORD32 ithread_mutex_destroy(void *mutex)
{
    return 0;
}

WORD32 ithread_mutex_lock(void *mutex)
{
    return 0;
}

WORD32 ithread_mutex_unlock(void *mutex)
{
    return 0;
}

void ithread_yield(void)
{
    return;
}

void ithread_sleep(UWORD32 u4_time_in_us)
{
    return;
}

void ithread_usleep(UWORD32 u4_time_us)
{
    return;
}

UWORD32 ithread_get_sem_strcut_size(void)
{
    return(sizeof(int));
}


WORD32 ithread_sem_init(void *sem,WORD32 pshared,UWORD32 value)
{
    return 0;
}

WORD32 ithread_sem_post(void *sem)
{
    return 0;
}


WORD32 ithread_sem_wait(void *sem)
{
    return 0;
}

WORD32 ithread_sem_destroy(void *sem)
{
    return 0;
}

void ithread_set_name(UWORD8 *pu1_thread_name)
{
    return;
}

void ithread_condition_init(void *condition)
{
    return;
}

void ithread_condition_signal(void * condition)
{
    return;
}



void ithread_condition_wait(void *condition,void *mutex)
{
    return;
}

WORD32 ithread_set_affinity(WORD32 core_id)
{
    return 1;
}
#endif
