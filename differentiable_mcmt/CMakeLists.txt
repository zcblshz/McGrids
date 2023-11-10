cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set(CMAKE_CUDA_ARCHITECTURES 60)

project(fast_MCMT LANGUAGES CXX CUDA)
 
SET(SRCS 
	fast_mcmt.cpp
	fast_mcmt.hpp
  kdtree.hpp
  nanoflann.hpp
  KDTreeVectorOfVectorsAdaptor.hpp
  )

add_library(fast_mcmt STATIC ${SRCS})
add_executable(fastCVT 	main.cpp sdfs.hpp)
 
if(MSVC)
  # #add ThirdPartyInclude
  # list( APPEND ThirdPartIncludePath "D:/Projects/geogram/geogram/src/lib")
  # list( APPEND ThirdPartIncludePath "D:/Projects/geogram/geogram/src/lib/geogram_gfx/third_party")
  # target_include_directories(fastCVT PUBLIC ${ThirdPartIncludePath})
  # #add ThirdPartylib
  # list( APPEND ThirdPartLibPath "D:/Projects/geogram/geogram/build/lib/Release/geogram.lib")
  # target_link_libraries(fastCVT PUBLIC ${ThirdPartLibPath})
  # target_link_libraries(fastCVT PUBLIC ${ThirdPartLibPath}

else()
  find_package(TBB REQUIRED)

  #add ThirdPartyInclude
  list( APPEND ThirdPartIncludePath "/home/daxuan/Documents/geogram/src/lib")
  list( APPEND ThirdPartIncludePath "/home/daxuan/Documents/geogram/src/lib/third_party")
  target_include_directories(fast_mcmt PUBLIC ${ThirdPartIncludePath})
  target_include_directories(fast_mcmt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  message(STATUS "ThirdPartIncludePath: ${ThirdPartIncludePath}")
  message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
  #add ThirdPartylib
  list( APPEND ThirdPartLibPath "/home/daxuan/Documents/geogram/build/Linux64-gcc-dynamic-Release/lib/libgeogram.so")
  target_link_libraries(fast_mcmt PUBLIC ${ThirdPartLibPath} TBB::tbb)
  target_link_libraries(fastCVT PUBLIC fast_mcmt)


endif()


set_property(TARGET fastCVT 
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

 