list(APPEND MVC_DEC_APP_SRCS "${AVC_ROOT}/examples/mvcdec/main.c")

libavc_add_executable(mvcdec libmvcdec SOURCES ${MVC_DEC_APP_SRCS})
target_compile_definitions(mvcdec PRIVATE PROFILE_ENABLE MD5_DISABLE)
