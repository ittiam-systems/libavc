LOCAL_PATH := $(call my-dir)
include $(CLEAR_VARS)

LOCAL_CLANG := true
LOCAL_DETECT_INTEGER_OVERFLOWS := true

# encoder
include $(LOCAL_PATH)/encoder.mk

# decoder
include $(LOCAL_PATH)/decoder.mk
