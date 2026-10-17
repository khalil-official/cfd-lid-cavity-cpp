cmake_minimum_required(VERSION 3.16)
project(LidCavityCPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LIDCAVITY_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

function(lidcavity_apply_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow
        )
    endif()
    if(LIDCAVITY_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

add_executable(lid_cavity src/lid_cavity.cpp)
lidcavity_apply_warnings(lid_cavity)

add_library(lidcavity_verification
    src/verification/operators.cpp
    src/verification/poisson.cpp
    src/verification/convergence.cpp
)
target_include_directories(lidcavity_verification PUBLIC include)
lidcavity_apply_warnings(lidcavity_verification)

enable_testing()

function(add_lidcavity_test target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE lidcavity_verification)
    lidcavity_apply_warnings(${target})
    add_test(NAME ${target} COMMAND ${target})
endfunction()

add_lidcavity_test(test_operators tests/test_operators.cpp)
add_lidcavity_test(test_poisson tests/test_poisson.cpp)
add_lidcavity_test(test_convergence tests/test_convergence.cpp)

add_test(
    NAME canonical_cavity_regression
    COMMAND lid_cavity
        --single
        --N 32
        --Re 100
        --scheme upwind
        --pressure RBSOR
        --strict
        --no-fields
)
set_tests_properties(canonical_cavity_regression PROPERTIES TIMEOUT 120)
