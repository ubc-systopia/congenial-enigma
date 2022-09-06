cmake_minimum_required(VERSION 3.16)
project(graph_preprocess CXX)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_STANDARD 99)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -pedantic -O3")

FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG master
        )
FetchContent_MakeAvailable(fmt)

find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
find_package(igraph REQUIRED)
find_package(oneDPL REQUIRED)

SET(CMAKE_C_COMPILER icx)

add_executable(graph_preprocess main.cpp io.cpp io.h typedefs.h util.cpp util.h)
add_executable(slashburn slashburn.cpp typedefs.h order_slashburn.cpp order_slashburn.h util.cpp util.h io.cpp io.h)

if (Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(graph_preprocess LINK_PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(slashburn LINK_PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)
endif ()
