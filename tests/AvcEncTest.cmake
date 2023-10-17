include(ExternalProject)
ExternalProject_Add(googletest
    GIT_REPOSITORY https://android.googlesource.com/platform/external/googletest
    GIT_TAG main
    PREFIX ${AVC_ROOT}/third_party/build/googletest
    SOURCE_DIR ${AVC_ROOT}/third_party/googletest
    TMP_DIR ${AVC_ROOT}/third_party/build/googletest/tmp
    INSTALL_COMMAND ""
)

list(
  APPEND
  AVCENCTEST_SRCS
  "${AVC_ROOT}/tests/AvcEncTest.cpp")

libavc_add_executable(AvcEncTest libavcenc
    SOURCES ${AVCENCTEST_SRCS}
    INCLUDES "${AVC_ROOT}/third_party/googletest/googletest/include")

target_link_libraries(AvcEncTest
    ${AVC_ROOT}/third_party/build/googletest/src/googletest-build/lib/libgtest.a
    ${AVC_ROOT}/third_party/build/googletest/src/googletest-build/lib/libgtest_main.a)

add_dependencies(AvcEncTest googletest)
