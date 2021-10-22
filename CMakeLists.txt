cmake_minimum_required(VERSION 3.19)
project(CURC LANGUAGES C CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GetVersion)
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CHECK_CUDA_EXE ${CMAKE_CURRENT_SOURCE_DIR}/env/check_cuda)
    set(CHECK_CUDA_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/env/check_cuda.cu)
    execute_process(COMMAND ${CMAKE_CUDA_COMPILER} -lcuda ${CHECK_CUDA_SOURCE} -o ${CHECK_CUDA_EXE})
    execute_process(COMMAND ${CHECK_CUDA_EXE}
            RESULT_VARIABLE CUDA_RETURN_CODE
            OUTPUT_VARIABLE ARCH)
    if(${CUDA_RETURN_CODE} EQUAL 0)
        message(STATUS "CUDA_ARCHITECTURE : ${ARCH}")
    else()
        message(FATAL_ERROR "CUDA_ARCHITECTURE require >= 60")
    endif()
else()
    message(FATAL_ERROR "No CUDA support")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.3.0)
        message(FATAL_ERROR "GCC version require >= 7.3.0")
    endif()
else()
    message(FATAL_ERROR "Only the GCC compiler is supported")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "-Ofast -m64 -march=native -fopenmp -DLIBBSC_OPENMP_SUPPORT")
set(CMAKE_CXX_FLAGS "-Ofast -m64 -march=native -fopenmp -DLIBBSC_OPENMP_SUPPORT")
set(CMAKE_CUDA_FLAGS "--expt-extended-lambda --compiler-bindir=${CMAKE_CXX_COMPILER} -Xcompiler -Ofast,-m64,-march=native,-fopenmp,-pthread -Xptxas -O3,--allow-expensive-optimizations=true")

find_package(CUDAToolkit 10.1 REQUIRED)
if(NOT DEFINED CUDAToolkit_VERSION)
    message(FATAL_ERROR "expected CUDAToolkit variable CUDAToolkit_VERSION not found")
endif()
message(STATUS "CUDAToolkit_VERSION: ${CUDAToolkit_VERSION}")
message(STATUS "CUDAToolkit_VERSION_MAJOR: ${CUDAToolkit_VERSION_MAJOR}")
message(STATUS "CUDAToolkit_VERSION_MINOR: ${CUDAToolkit_VERSION_MINOR}")
message(STATUS "CUDAToolkit_VERSION_PATCH: ${CUDAToolkit_VERSION_PATCH}")
message(STATUS "CUDAToolkit_BIN_DIR: ${CUDAToolkit_BIN_DIR}")
message(STATUS "CUDAToolkit_INCLUDE_DIRS: ${CUDAToolkit_INCLUDE_DIRS}")
message(STATUS "CUDAToolkit_LIBRARY_DIR: ${CUDAToolkit_LIBRARY_DIR}")
message(STATUS "CUDAToolkit_NVCC_EXECUTABLE ${CUDAToolkit_NVCC_EXECUTABLE}")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(SYSTEM third_party/)

set(LZMA_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/VarLenDNACoder/LzmaLib.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/VarLenDNACoder/VarLenDNACoder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/VarLenDNACoder/helper.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/LzmaEnc.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/LzmaDec.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/Alloc.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/LzFindMt.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/LzFind.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/Ppmd7.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/Ppmd7Enc.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/Ppmd7Dec.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lzma/Threads.c)

set(LIBBSC_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/bsc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/adler32/adler32.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/bwt/divsufsort/divsufsort.c
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/bwt/bwt.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/coder/qlfc/qlfc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/coder/qlfc/qlfc_model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/coder/coder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/filters/detectors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/filters/preprocessing.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/libbsc/libbsc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/lzp/lzp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libbsc/libbsc/platform/platform.cpp)

add_custom_target(
        FastLZMA2
        COMMAND make clean && make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third_party/fast-lzma2/
)

add_executable(curc src/main.cpp src/compress.cu src/decompress.cpp src/fast_lzma2_helper.cpp ${LZMA_FILES})
add_executable(bsc ${LIBBSC_SOURCE_FILES})
# add_executable(bench bench/benchmark.cpp)
# add_executable(sort_pe_read bench/sort_pe_read.cpp)

add_dependencies(curc FastLZMA2)
target_link_libraries(curc PRIVATE stdc++fs CUDA::toolkit Threads::Threads ${CMAKE_CURRENT_SOURCE_DIR}/third_party/fast-lzma2/libfast-lzma2.a)
target_link_libraries(bsc  PRIVATE Threads::Threads)
# target_link_libraries(bench PRIVATE stdc++fs Threads::Threads)
# target_link_libraries(sort_pe_read PRIVATE stdc++fs Threads::Threads)

set_target_properties(
        curc
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${ARCH}")
