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

### Building for android
NOTE: This assumes that you are building on a machine that has
 [Android NDK](https://developer.android.com/ndk/downloads).

```
$ cd external/libavc
$ mkdir build
$ cd build
```

#### Armv7 (32-bit)

    cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/android_toolchain.cmake\
        -DAVC_ANDROID_NDK_PATH=/opt/android-ndk-r26d/\
        -DANDROID_ABI=armeabi-v7a\
        -DANDROID_PLATFORM=android-23 ../
    make

#### Armv8 (64-bit)

    cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/android_toolchain.cmake\
        -DAVC_ANDROID_NDK_PATH=/opt/android-ndk-r26d/\
        -DANDROID_ABI=arm64-v8a\
        -DANDROID_PLATFORM=android-23 ../
    make