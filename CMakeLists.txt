cmake_minimum_required (VERSION 3.16)
enable_testing()

project(daisy_tiramisu_backend VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

Include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.10.2
)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
    GIT_TAG        v3.9.1
)
set(JSON_BuildTests OFF CACHE INTERNAL "") 
FetchContent_MakeAvailable(nlohmann_json)

find_library(ISLLib NAMES libisl.so)
message(STATUS ${ISLLib})

add_library(isl STATIC IMPORTED)
set_target_properties(isl PROPERTIES IMPORTED_LOCATION ${ISLLib})

add_library(${PROJECT_NAME} STATIC)
target_sources(${PROJECT_NAME}
    PRIVATE
        src/core.cpp
        src/buffer.cpp
        src/function.cpp
        src/computation.cpp
        src/ast.cpp
        src/auto_scheduler.cpp
        src/dnn_accesses.cpp
        src/evaluator.cpp
        src/optimization_info.cpp
        src/sdfg_wrapper.cpp
        src/schedules_generator.cpp
        src/search_method.cpp
)
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> 
)
target_link_libraries(daisy_tiramisu_backend PRIVATE isl pybind11::embed nlohmann_json)

set(public_headers
    include/ast.h
    include/auto_scheduler.h
    include/buffer.h
    include/core.h
    include/computation.h
    include/function.h
    include/dnn_accesses.h
    include/evaluator.h
    include/function.h
    include/optimization_info.h
    include/sdfg_wrapper.h
    include/schedules_generator.h
    include/search_method.h
    include/utils.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

# install the target and create export-set
install(TARGETS ${PROJECT_NAME}
    EXPORT "${PROJECT_NAME}Targets"
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

add_subdirectory(test)
