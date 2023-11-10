cmake_minimum_required(VERSION 3.15..3.26)

# Set a name and a version number for your project:
project(
  mcmt
  VERSION 0.0.7
  LANGUAGES CXX)

# Initialize some default paths
include(GNUInstallDirs)

# Define the minimum C++ standard that is required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Configuration options
option(BUILD_PYTHON "Enable building of Python bindings" ON)
option(BUILD_DOCS "Enable building of documentation" ON)


# add mcmt
add_subdirectory(differentiable_mcmt)
add_subdirectory(python)



# Add an alias target for use if this project is included as a subproject in
# another project
# add_library(mcmt::mcmt ALIAS mcmt_lib)

# install(
#   TARGETS mcmt_lib
#   EXPORT mcmt-config
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install(
#   EXPORT mcmt-config
#   NAMESPACE mcmt::
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mcmt)

# install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# # This prints a summary of found dependencies
# include(FeatureSummary)
# feature_summary(WHAT ALL)