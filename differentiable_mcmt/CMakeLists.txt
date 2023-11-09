cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set(CMAKE_CUDA_ARCHITECTURES 60)

project(fastCVT LANGUAGES CXX CUDA)
 
SET(SRCS 
	main.cpp
	mcmt.cpp
	mcmt.hpp)
add_executable(fastCVT ${SRCS})
 
if(MSVC)
  #add ThirdPartyInclude
  list( APPEND ThirdPartIncludePath "D:/Projects/geogram/geogram/src/lib")
  list( APPEND ThirdPartIncludePath "D:/Projects/geogram/geogram/src/lib/geogram_gfx/third_party")
  target_include_directories(fastCVT PRIVATE ${ThirdPartIncludePath})
  #add ThirdPartylib
  list( APPEND ThirdPartLibPath "D:/Projects/geogram/geogram/build/lib/Release/geogram.lib")
  target_link_libraries(fastCVT PRIVATE ${ThirdPartLibPath})

else()
  find_package(TBB REQUIRED)

  #add ThirdPartyInclude
  list( APPEND ThirdPartIncludePath "/home/daxuan/Documents/geogram/src/lib")
  list( APPEND ThirdPartIncludePath "/home/daxuan/Documents/geogram/src/lib/third_party")
  target_include_directories(fastCVT PRIVATE ${ThirdPartIncludePath})
  #add ThirdPartylib
  list( APPEND ThirdPartLibPath "/home/daxuan/Documents/geogram/build/Linux64-gcc-dynamic-Release/lib/libgeogram.so")
  target_link_libraries(fastCVT PRIVATE ${ThirdPartLibPath} TBB::tbb)

endif()


set_property(TARGET fastCVT 
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

 