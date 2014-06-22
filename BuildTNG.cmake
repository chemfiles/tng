set(TNG_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
file(RELATIVE_PATH TNG_ROOT_BINARY_DIR ${CMAKE_SOURCE_DIR} ${TNG_ROOT_SOURCE_DIR})
set(TNG_ROOT_BINARY_DIR ${CMAKE_BINARY_DIR}/${TNG_ROOT_BINARY_DIR})

function (TNG_GENERATE_VERSION_H)
    set(MAJOR_VERSION "1" CACHE STRING "" FORCE)
    set(MINOR_VERSION "6")
    set(VERSION_PATCH_LEVEL "0")
    set(TNG_IO_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${VERSION_PATCH_LEVEL}" CACHE STRING "" FORCE)
    set(API_VERSION "6")
    configure_file(${TNG_ROOT_SOURCE_DIR}/include/tng/version.h.in
                   ${TNG_ROOT_BINARY_DIR}/include/tng/version.h)
endfunction()

tng_generate_version_h()

include(TestBigEndian)
test_big_endian(TNG_INTEGER_BIG_ENDIAN)
include(CheckIncludeFile)
check_include_file(inttypes.h TNG_HAVE_INTTYPES_H)

macro(TNG_GET_SOURCE_LIST TNG_SOURCELIST TNG_COMPILEDEFS)
    include_directories(${TNG_ROOT_SOURCE_DIR}/include)
    include_directories(${TNG_ROOT_BINARY_DIR}/include)
    set(_tng_compression_sources bwlzh.c bwt.c coder.c dict.c fixpoint.c huffman.c huffmem.c lz77.c merge_sort.c mtf.c rle.c tng_compress.c vals16.c warnmalloc.c widemuldiv.c xtc2.c xtc3.c)
    set(_tng_io_sources tng_io.c md5.c)
    set(${TNG_SOURCELIST})
    set(${TNG_COMPILEDEFS})
    foreach(_file ${_tng_compression_sources})
        list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/compression/${_file})
    endforeach()
    foreach(_file ${_tng_io_sources})
        list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/lib/${_file})
    endforeach()
    if(TNG_BUILD_FORTRAN)
      list(APPEND ${TNG_SOURCELIST} ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io_fortran.c)
    endif()
    if (TNG_HAVE_INTTYPES_H)
        list(APPEND ${TNG_COMPILEDEFS} USE_STD_INTTYPES_H)
    endif()
endmacro()

macro(TNG_SET_SOURCE_PROPERTIES)
    set(_tng_with_zlib OFF)
    set(_curr_var)
    foreach (_arg ${ARGN})
        if (_arg STREQUAL "WITH_ZLIB")
            set(_curr_var with_zlib)
        elseif (_curr_var)
            set(_tng_${_curr_var} ${_arg})
            set(_curr_var "")
        else()
            message(FATAL_ERROR "Invalid argument ${_arg} to TNG_SET_SOURCE_PROPERTIES")
        endif()
    endforeach()
    if (_tng_with_zlib)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io.c
                     APPEND PROPERTY COMPILE_DEFINITIONS USE_ZLIB)
    endif()
    if (TNG_HAVE_INTTYPES_H)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/tng_io.c
                     APPEND PROPERTY COMPILE_DEFINITIONS USE_STD_INTTYPES_H)
    endif()
    if (TNG_INTEGER_BIG_ENDIAN)
        set_property(SOURCE ${TNG_ROOT_SOURCE_DIR}/src/lib/md5.c
                     APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
    endif()
endmacro()
