if(NOT "${SYSTEM_NAME}" STREQUAL "Darwin")
    libavc_add_fuzzer(avc_enc_fuzzer libavcenc SOURCES
                  ${AVC_ROOT}/fuzzer/avc_enc_fuzzer.cpp)
endif()