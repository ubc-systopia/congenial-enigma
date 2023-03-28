cmake_minimum_required(VERSION 3.22.0)
project(graph_preprocess CXX)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_FLAGS "-fopenmp -DLIKWID_PERFMON -DTBB_ENABLE_IPO=OFF")

#set(CMAKE_CXX_COMPILER /opt/intel/oneapi/compiler/2023.0.0/linux/bin/icpx)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -pedantic -O3")

#-D CMAKE_CXX_COMPILER=/opt/intel/oneapi/compiler/latest/linux/bin-llvm/clang++

FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG master
        )
FetchContent_MakeAvailable(fmt)

##FetchContent_Declare(taskflow
##        GIT_REPOSITORY https://github.com/taskflow/taskflow.git
##        GIT_TAG master
##        )
##FetchContent_MakeAvailable(taskflow)
#
#FetchContent_Declare(work-stealing-queue
#        GIT_REPOSITORY https://github.com/taskflow/work-stealing-queue.git
#        GIT_TAG master
#        )
#FetchContent_MakeAvailable(work-stealing-queue)

FetchContent_Declare(Spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG master
        )
FetchContent_MakeAvailable(Spectra)


FetchContent_Declare(spray
        GIT_REPOSITORY https://github.com/atrostan/spray.git
        GIT_TAG master
        )

FetchContent_GetProperties(spray)
if (NOT spray_POPULATED)
    FetchContent_Populate(spray)
endif ()

message(spray_SOURCE_DIR="${spray_SOURCE_DIR}")

FetchContent_MakeAvailable(spray)

# boost <= 1.74 serialization breaks when reading binary archives of edgelists
# todo boost != 1.58 breaks error: static assertion failed: boost::atomic<T> requires T to be a trivially copyable type for union rabbit_order::aux::atomix<rabbit_order::aux::atom>
find_package(Boost 1.58.0 EXACT COMPONENTS system filesystem graph serialization REQUIRED)
find_package(igraph REQUIRED)
find_package(oneDPL  REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(IntelDPCPP REQUIRED)


#SET(CMAKE_C_COMPILER /opt/intel/oneapi/compiler/2023.0.0/linux/bin/icx) # todo use icx compiler when using intel dpl
set(ParSlashburnDir "../par_slashburn")
set(IPS4ODir "../par_slashburn/ips4o")
add_subdirectory(${IPS4ODir} ips4o)

add_executable(graph_preprocess main.cpp io.cpp io.h typedefs.h util.cpp util.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h Hilbert.cpp Hilbert.h ${ParSlashburnDir}/pvector.h)

target_include_directories(graph_preprocess PUBLIC ${ParSlashburnDir})


add_executable(slashburn slashburn.cpp typedefs.h order_slashburn.cpp order_slashburn.h util.cpp util.h io.cpp io.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h sql.h Hilbert.cpp Hilbert.h)
target_include_directories(slashburn PUBLIC ${ParSlashburnDir})

add_executable(cuthill_mckee cuthill_mckee.cpp typedefs.h util.cpp util.h io.cpp io.h rabbit_util.cpp rabbit_util.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h)
target_include_directories(cuthill_mckee PUBLIC ${ParSlashburnDir})

add_executable(pr_experiments pr_experiments.cpp pr_experiments.h typedefs.h util.cpp util.h io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h PageRank.cpp PageRank.h ${ParSlashburnDir}/pvector.h)
target_include_directories(pr_experiments PUBLIC ${ParSlashburnDir})

add_executable(sb_furhilbert sb_furhilbert.cpp sb_furhilbert.h typedefs.h util.cpp util.h io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h Hilbert.cpp Hilbert.h PageRank.cpp PageRank.h furhilbert.h QuadrantCalculator.cpp QuadrantCalculator.h furhilbert.cpp)

add_executable(convert_map_to_binary convert_map_to_binary.cpp io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h)

add_executable(compute_ccs compute_ccs.cpp io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h ${ParSlashburnDir}/cc.h ${ParSlashburnDir}/cc.cpp)
target_include_directories(compute_ccs PUBLIC ${ParSlashburnDir})

add_executable(stats stats.cpp io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h ${ParSlashburnDir}/cc.h ${ParSlashburnDir}/cc.cpp stats.h)
target_include_directories(stats PUBLIC ${ParSlashburnDir})

add_executable(par_hilbert_bench par_hilbert_bench.cpp io.cpp io.h sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h Hilbert.cpp Hilbert.h wsq.h)
target_include_directories(par_hilbert_bench PUBLIC ${ParSlashburnDir})

add_executable(hilburn_preproc hilburn_preproc.cpp hilburn_preproc.h io.cpp io.h util.h util.cpp sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h Hilbert.cpp Hilbert.h Quad.h CSR.cpp CSR.h)
target_include_directories(hilburn_preproc PUBLIC ${ParSlashburnDir} ${IPS4ODir})

add_executable(hilburn_pr hilburn_pr.cpp io.cpp io.h util.h util.cpp sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h Hilbert.cpp Hilbert.h Quad.h)
target_include_directories(hilburn_pr PUBLIC ${ParSlashburnDir} ${IPS4ODir})

add_executable(hilburn_stripe_preproc hilburn_stripe_preproc.cpp  io.cpp io.h util.h util.cpp sql.cpp sql.h rabbit_util.cpp rabbit_util.h ${ParSlashburnDir}/pvector.h Hilbert.cpp Hilbert.h Quad.h CSR.cpp CSR.h)
target_include_directories(hilburn_stripe_preproc PUBLIC ${ParSlashburnDir} ${IPS4ODir})

find_package(likwid REQUIRED)
#find_package(taskflow REQUIRED)
add_library(spray INTERFACE)

target_include_directories(spray INTERFACE ${spray_SOURCE_DIR}/include)

if (Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${likwid_INCLUDE_DIRS})

    target_link_libraries(graph_preprocess PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3 spray)

    target_link_libraries(slashburn PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(cuthill_mckee PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(pr_experiments PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(sb_furhilbert PUBLIC ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(convert_map_to_binary ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3)

    target_link_libraries(compute_ccs ${Boost_LIBRARIES} Eigen3::Eigen fmt::fmt igraph::igraph oneDPL sqlite3 Spectra)

    target_link_libraries(stats ${Boost_LIBRARIES} Eigen3::Eigen fmt::fmt igraph::igraph oneDPL sqlite3 Spectra)

    target_link_libraries(par_hilbert_bench ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3 likwid Taskflow pthread)

    target_link_libraries(hilburn_preproc ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3 likwid ips4o)
    target_link_libraries(hilburn_pr ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3 likwid ips4o)
    target_link_libraries(hilburn_stripe_preproc ${Boost_LIBRARIES} fmt::fmt igraph::igraph oneDPL sqlite3 likwid ips4o)

endif ()