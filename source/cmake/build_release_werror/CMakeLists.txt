CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)
SET(LINK_RULES shared)

## Note: This build uses the CMAKE default compiler settings for your computer.
## Changing the "COMPILER" setting here will *not* change which compiler is used;
## it only changes the extension of the executable (e.g. .linuxgccdebug).
##
## To change the default compiler for CMAKE, set the $CC and $CXX shell environment variables.
## For example, in bash, "export CXX=clang++; export CC=clang;"
## Or "export CXX='distcc g++'; export CC='distcc gcc';"
##
## Alternatively, you can pass the -DCMAKE_CXX_COMPILER="<compiler name>" flag
## when you invoke cmake.
##
## IMPORTANT: These setting are only obeyed on the first invocation on a clean build.
## An incremental build will not change the compiler even if the $CC/$CXX settings change.
## Because of this - and unlike the scons system - the cmake system can't handle
## multiple compilers on the same target. If you need multiple simultaneous compiler
## builds, use additional build directories.


## Try to set the executable extension appropriately for people with a default compiler of clang.
## Assume that if you're not obviously running clang, you're running gcc.
## There's additional CMAKE flags that can tease other cases apart, but
## this /should/ work for most of Rosetta developers.

IF( ${CMAKE_CXX_COMPILER_ID} MATCHES ".*[Cc][Ll][Aa][Nn][Gg].*" )
    SET(COMPILER clang)
ELSE()
    SET(COMPILER gcc)
ENDIF()

MESSAGE( ">> CMAKE identifies C++ compiler as '${CMAKE_CXX_COMPILER_ID}', interpreting this as '${COMPILER}'" )
MESSAGE( ">> To change, set CXX and CC environment variables (or pass -DCMAKE_CXX_COMPILER) and do a clean rebuild.")
MESSAGE( ">>  current settings: CXX='$ENV{CXX}' CC='$ENV{CC}'" )

############ Build Rules

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)

#SET(COMPILE_FLAGS "${COMPILE_FLAGS} -std=c++11 -pipe -ffor-scope -O3 -ffast-math -fno-finite-math-only -funroll-loops -finline-functions")
ADD_DEFINITIONS(-DUNUSUAL_ALLOCATOR_DECLARATION)
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -stdlib=libstdc++ -pipe -ffast-math -fno-finite-math-only -O0 -Qunused-arguments -ftemplate-depth-256 -march=native -fPIC")
#ET(COMPILE_FLAGS "${COMPILE_FLAGS} -stdlib=libstdc++ -pipe -ffast-math -fno-finite-math-only -O3 -ffor-scope -funroll-loops -finline-functions -s")


## Inhibit non-error warning messages to keep ninja output clean.
## This is now turned off as it also inhibits the warnings as errors settings,
## which would mean a nasty surprise for people pushing a "clean" cmake build.
#SET(COMPILE_FLAGS "${COMPILE_FLAGS} -w")

## Other flags that sometimes occur with scons. Uncomment to add.
##
## More agressive inlining - this gets set for some platform, but causes compiler errors for others:
#
#SET(COMPILE_FLAGS "${COMPILE_FLAGS} -finline-limit=20000")
#
## More warnings - good for cleaning up code, but makes ninja output more noisy.
## (Remember to also comment out the -w above if uncommenting.)
## (no-unused-variable is here because asserts can leave variables unused.)

SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Wall -Wextra -pedantic -Wno-long-long -Wno-strict-aliasing")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Wno-unused-private-field -Wno-overloaded-virtual")

## Warnings as errors to match scons build
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Werror")

IF( ${COMPILER} STREQUAL "clang" )
    # Quiet warnings about non-clang compatable arguments (either gcc-specific ones, or newer-version ones.)
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Qunused-arguments")
    # Recent Macs can have issues with the default standard library - explicitly specify the GNU one.
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} -stdlib=libstdc++")
    ##Some clang builds may also need to set additional compile flags
    #SET(COMPILE_FLAGS "${COMPILE_FLAGS} -DUNUSUAL_ALLOCATOR_DECLARATION") # Works around issues with vector1 optimization
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
