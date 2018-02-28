# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
#
# CMakeLists.txt file for unique_kmers.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_yara CXX)
message (STATUS "Configuring unique_kmers")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  OpenMP 3.0 not found: building yara without multi-threading.")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib not found: building yara without gzipped input and BAM output.")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (STATUS "  BZip2 not found: building yara without bzipped input.")
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
include_directories (${SDSL_INCLUDE_DIRS})
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
                            store_seqs.h
                            misc_timer.h
                            misc_types.h
                            bits_matches.h
                            misc_options.h
                            misc_options_dis.h
                            kdx_filter.h
                            bloom_filter.h
                            index_fm.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (unique_kmers ${SEQAN_LIBRARIES})


# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------
