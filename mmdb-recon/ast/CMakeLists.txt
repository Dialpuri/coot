cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(FindMMDBFunctions)

# Find LLVM and Clang
find_package(LLVM REQUIRED CONFIG)
find_package(Clang REQUIRED CONFIG)

# Add LLVM definitions
add_definitions(${LLVM_DEFINITIONS})

add_executable(find_mmdb_functions find_mmdb_functions.cpp)

target_include_directories(find_mmdb_functions PRIVATE
        ${CLANG_INCLUDE_DIRS}
        ${LLVM_INCLUDE_DIRS}
)

# Link against necessary clang and LLVM libraries
target_link_libraries(find_mmdb_functions PRIVATE
        clangTooling
        clangBasic
        clangAST
        clangFrontend
        clangParse
        clangLex
        clangDriver
        clangSerialization
        LLVMSupport
        LLVMCore
)

# Compiler flags
if(APPLE)
    target_compile_options(find_mmdb_functions PRIVATE -fno-rtti)
endif()