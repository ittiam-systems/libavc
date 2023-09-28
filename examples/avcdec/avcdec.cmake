libavc_add_executable(avcdec libavcdec SOURCES ${AVC_ROOT}/examples/avcdec/main.c)
target_compile_definitions(avcdec PRIVATE PROFILE_ENABLE MD5_DISABLE)
