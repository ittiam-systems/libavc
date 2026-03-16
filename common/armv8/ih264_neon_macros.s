//******************************************************************************
//*
//* Copyright (C) 2015 The Android Open Source Project
//*
//* Licensed under the Apache License, Version 2.0 (the "License");
//* you may not use this file except in compliance with the License.
//* You may obtain a copy of the License at:
//*
//* http://www.apache.org/licenses/LICENSE-2.0
//*
//* Unless required by applicable law or agreed to in writing, software
//* distributed under the License is distributed on an "AS IS" BASIS,
//* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//* See the License for the specific language governing permissions and
//* limitations under the License.
//*
//*****************************************************************************
//* Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
//*/
//*******************************************************************************


.macro push_v_regs
    stp       d8, d9, [sp, #-16]!
    stp       d10, d11, [sp, #-16]!
    stp       d12, d13, [sp, #-16]!
    stp       d14, d15, [sp, #-16]!
.endm
.macro pop_v_regs
    ldp       d14, d15, [sp], #16
    ldp       d12, d13, [sp], #16
    ldp       d10, d11, [sp], #16
    ldp       d8, d9, [sp], #16
.endm

.macro swp reg1, reg2
    eor       \reg1, \reg1, \reg2
    eor       \reg2, \reg1, \reg2
    eor       \reg1, \reg1, \reg2
.endm

// --- Internal Security Dispatchers ---
// These expand to real instructions only if the compiler flags are present.

.macro BTI_ENABLE
#if defined(__ARM_FEATURE_BTI_DEFAULT)
    bti c
#endif
.endm

.macro PAC_ENTRY
#if defined(__ARM_FEATURE_PAC_DEFAULT)
    paciasp
#endif
.endm

.macro PAC_EXIT
#if defined(__ARM_FEATURE_PAC_DEFAULT)
    autiasp
#endif
.endm

// --- Main ENTRY and EXIT_FUNC Macros ---

.macro ENTRY name
    .p2align 2
\name:
    BTI_ENABLE
    PAC_ENTRY
.endm

.macro EXIT_FUNC
    PAC_EXIT
.endm

// --- GNU Property Note ---
// Signals BTI and PAC support to the Android linker.
#if defined(__linux__) && defined(__aarch64__)
    .pushsection .note.gnu.property, "a"  // Switch to Note section
    .p2align 3
    .word 4           // Name size
    .word 16          // Data size
    .word 5           // NT_GNU_PROPERTY_TYPE_0
    .asciz "GNU"      // Owner
    .word 0xc0000000  // GNU_PROPERTY_AARCH64_FEATURE_1_AND
    .word 4           // Data size
    .word 3           // Value: BTI (Bit 0) | PAC (Bit 1)
    .word 0           // Padding
    .popsection                           // Switch back to previous section
#endif
