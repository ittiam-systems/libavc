set(SYSTEM_NAME Linux)
set(SYSTEM_PROCESSOR aarch32)
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR arm)

# Modify these variables with paths to appropriate compilers that can produce
# armv7 targets
set(CMAKE_C_COMPILER arm-linux-gnueabihf-gcc)
set(CMAKE_CXX_COMPILER arm-linux-gnueabihf-g++)

# Build all binaries as static, so that they can be run using qemu
set(CMAKE_EXE_LINKER_FLAGS "-static")
