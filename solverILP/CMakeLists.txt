project(LifT_Solver)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# compiler options
add_definitions(-march=native)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(include)
#add_subdirectory(src)



find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
    add_definitions("-DWITH_GUROBI")

    add_library(MyCallback src/disjoint-paths/MyCallback.cxx)
    target_include_directories(MyCallback PUBLIC external/andres-graph/include)
    target_include_directories(MyCallback PUBLIC external/TCLAP/include)
    target_link_libraries(MyCallback ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})

    add_executable(run-disjoint-paths src/disjoint-paths/runDisjointPaths.cxx  ${headers} ${sources})
    target_include_directories(run-disjoint-paths PUBLIC external/andres-graph/include)
    target_include_directories(MyCallback PUBLIC external/TCLAP/include)
    target_link_libraries(run-disjoint-paths ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY} MyCallback)
endif()

if(COMPILER_SUPPORTS_CXX0X OR COMPILER_SUPPORTS_CXX11)
    include_directories(src)

endif()


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Some functionality will not be available.")
endif()




