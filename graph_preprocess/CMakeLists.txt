cmake_minimum_required(VERSION 3.16)
project(graph_preprocess CXX)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_FLAGS "-fopenmp")

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -pedantic -O3")

FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG master
        )
FetchContent_MakeAvailable(fmt)

# boost <= 1.74 serialization breaks when reading binary archives of edgelists
# todo boost != 1.58 breaks error: static assertion failed: boost::atomic<T> requires T to be a trivially copyable type for union rabbit_order::aux::atomix<rabbit_order::aux::atom>
find_package(Boost 1.58.0 EXACT COMPONENTS system filesystem graph serialization REQUIRED)
find_package(igraph REQUIRED)
find_package(oneDPL REQUIRED)

SET(CMAKE_C_COMPILER icx)

add_executable(graph_preprocess main.cpp io.cpp io.h typedefs.h util.cpp util.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h Hilbert.cpp Hilbert.h)

add_executable(slashburn slashburn.cpp typedefs.h order_slashburn.cpp order_slashburn.h util.cpp util.h io.cpp io.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h sql.h Hilbert.cpp Hilbert.h)
set(ParSlashburnDir "../par_slashburn")

target_include_directories(slashburn PUBLIC ${ParSlashburnDir})

add_executable(cuthill_mckee cuthill_mckee.cpp typedefs.h util.cpp util.h io.cpp io.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h)

add_executable(pr_experiments pr_experiments.cpp pr_experiments.h typedefs.h util.cpp util.h io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h PageRank.cpp PageRank.h ${ParSlashburnDir}/pvector.h)

target_include_directories(pr_experiments PUBLIC ${ParSlashburnDir})

add_executable(sb_furhilbert sb_furhilbert.cpp sb_furhilbert.h typedefs.h util.cpp util.h io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h PageRank.cpp PageRank.h furhilbert.h QuadrantCalculator.cpp QuadrantCalculator.h furhilbert.cpp)

add_executable(convert_map_to_binary convert_map_to_binary.cpp io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h)

if (Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})

    target_link_libraries(graph_preprocess PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(slashburn PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(cuthill_mckee PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(pr_experiments PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(sb_furhilbert PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(convert_map_to_binary ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)
endif ()