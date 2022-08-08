list(APPEND MVC_DEC_APP_SRCS "${AVC_ROOT}/test/mvcdec/main.c")

libavc_add_executable(mvcdec libmvcdec SOURCES ${MVC_DEC_APP_SRCS})
