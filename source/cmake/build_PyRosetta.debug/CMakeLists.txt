CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

##############################################################################
# This file only needs to specify which mode (i.e., release, debug),
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be
# implemented in ../build/build.settings.cmake.
##############################################################################

############ Build Type, Mode, Link Rules, Compiler, Extras
SET(MODE debug)
SET(EXTRAS zeromq)
SET(WITH_PYROSETTA ON CACHE INTERNAL "Add PyRosetta build options.")
SET(WITH_LIBSTDC++ OFF CACHE INTERNAL "Allow link against platform-native standard library.")
SET(LINK_RULES shared)

############ Build Rules
INCLUDE(../build/compiler.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)
INCLUDE(../build/build.settings.cmake)
INCLUDE(../build/build.post.cmake)
