# ==============================================================================

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(oomph_convert VERSION 1.0.0)
enable_testing()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build in debugging mode by default
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# ---------------------------------[ PACKAGES ]---------------------------------

# include(FetchContent) FetchContent_Declare( google_benchmark GIT_REPOSITORY
# https://github.com/google/benchmark.git GIT_TAG v1.7.1 GIT_PROGRESS TRUE
# GIT_SHALLOW TRUE) FetchContent_MakeAvailable(google_benchmark)

# ------------------------------------------------------------------------------

# Define target and link to Boost libraries
add_executable(oomph-convert oomph-convert.cc)

# Symlink the data/ folder into the build/ directory when building
add_custom_command(
  TARGET oomph-convert COMMAND ln -sf "${CMAKE_CURRENT_LIST_DIR}/data"
                               "${CMAKE_CURRENT_BINARY_DIR}")

# ------------------------------------------------------------------------------

# # Define target and link to Boost libraries
# add_executable(oomph-convert-benchmark oomph-convert-benchmark.cc)
# target_link_libraries(oomph-convert-benchmark PRIVATE benchmark::benchmark)

# ------------------------------------------------------------------------------
