list(APPEND SVC_DEC_APP_SRCS "${AVC_ROOT}/examples/svcdec/main.c")

libavc_add_executable(svcdec libsvcdec SOURCES ${SVC_DEC_APP_SRCS})
target_compile_definitions(svcdec PRIVATE PROFILE_ENABLE MD5_DISABLE)
