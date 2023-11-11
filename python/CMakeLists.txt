project(cmake_and_cuda LANGUAGES C CXX CUDA)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(mcmt_python mcmt.cpp)

if(MSVC)
  #add ThirdPartyInclude
  #find_package(TBB REQUIRED)

else()
  find_package(TBB REQUIRED)
  # target_include_directories(mcmt_python PUBLIC ${ThirdPartIncludePath} )
  # target_link_libraries(mcmt_python PUBLIC ${ThirdPartLibPath} )
  # target_link_libraries(mcmt_python PUBLIC ${ThirdPartLibPath} )
endif()


target_link_libraries(mcmt_python PUBLIC fast_mcmt)
set_target_properties(mcmt_python PROPERTIES OUTPUT_NAME mcmt)

install(TARGETS mcmt_python DESTINATION .)
