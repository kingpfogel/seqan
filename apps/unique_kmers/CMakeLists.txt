# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/unqiue_kmers/CMakeLists.txt
#
# CMakeLists.txt file for unique_kmers.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_unique_kmers CXX)
message (STATUS "Configuring apps/unique_kmers")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib not found: building unique_kmers without gzipped input and BAM output.")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (STATUS "  BZip2 not found: building unique_kmers without bzipped input.")
endif (NOT BZIP2_FOUND)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "1.0.0")

# ----------------------------------------------------------------------------
# Disable Compiler Warnings
# ----------------------------------------------------------------------------

# Disables warning about conversion from 64 bit to 32 bit integrals.
# Note that there is a runtime exception thrown if the type is not sufficient.
# Unfortunately the VC does not detect it.
if (MSVC)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4244")  # When compiling on 32 bit windows.
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4267")  # When compiling on 64 bit windows.
endif (MSVC)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(unique_kmers unique_kmers.cpp
                            bits_matches.h
                            index_fm.h
                            misc_options_dis.h
                            misc_options.h
                            misc_types.h
                            store_seqs.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (unique_kmers ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (unique_kmers)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install unique_kmers in ${PREFIX}/bin directory
install (TARGETS unique_kmers
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/unique_kmers for SeqAn release builds.
# install (FILES LICENSE README.rst
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable unique_kmers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES}
                           unique_kmers CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:unique_kmers")
  set (CPACK_PACKAGE_NAME "unique_kmers")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "unique_kmers - Yet Another Read Aligner")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
  set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")

  seqan_configure_cpack_app (unique_kmers "unique_kmers")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:unique_kmers")
