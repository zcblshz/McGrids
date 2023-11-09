cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set(CMAKE_CUDA_ARCHITECTURES 60)

project(fastCVT LANGUAGES CXX CUDA)
 
SET(SRCS 
	main.cpp
	fast_cvt3d.cpp
	fast_cvt3d.h)
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

endif()


set_property(TARGET fastCVT 
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

 