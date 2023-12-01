cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(fast_MCMT LANGUAGES CXX)

SET(SRCS 
	fast_mcmt.cpp
	fast_mcmt.hpp
  )

add_library(fast_mcmt STATIC ${SRCS})
add_executable(fastCVT 	main.cpp sdfs.hpp)

find_package(TBB REQUIRED)
find_package(CGAL REQUIRED)
# include fast_mcmt.hpp with fast_mcmt lib
target_include_directories(fast_mcmt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fast_mcmt PUBLIC TBB::tbb TBB::tbbmalloc CGAL::CGAL)
target_link_libraries(fastCVT PUBLIC fast_mcmt)

 