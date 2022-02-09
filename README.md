# LIBAVC
## Getting Started Document

# LibAVC build steps

Supports:
- aarch32/aarch64 on Linux.
- aarch32/aarch64 on Android.
- x86_32/x86_64 on Linux.

## Native Builds
Use the following commands for building on the target machine

```
$ cd external/libavc
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Cross-compiler based builds

### Building for x86_32 on a x86_64 Linux machine
```
$ cd external/libavc
$ mkdir build
$ cd build
$ CFLAGS="-m32" CXXFLAGS="-m32" LDFLAGS="-m32" cmake ..
$ make
```

### Building for aarch32/aarch64
Update 'CMAKE_C_COMPILER', 'CMAKE_CXX_COMPILER', 'CMAKE_C_COMPILER_AR', and
'CMAKE_CXX_COMPILER_AR' in CMAKE_TOOLCHAIN_FILE passed below

```
$ cd external/libavc
$ mkdir build
$ cd build
```

#### For aarch64
```
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/aarch64_toolchain.cmake
$ make
```

#### For aarch32
```
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/aarch32_toolchain.cmake
$ make
```
