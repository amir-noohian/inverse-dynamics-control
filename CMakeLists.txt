cmake_minimum_required(VERSION 2.6)
project(inverse_dynamics_control)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} include)

## libbarrett
find_package(Barrett REQUIRED)
include_directories(${BARRETT_INCLUDE_DIRS})
link_directories(${BARRETT_LIBRARY_DIRS})
add_definitions(${BARRETT_DEFINITIONS})

## curses
find_package(Curses REQUIRED)
include_directories(${CURSES_INCLUDE_DIR})

## Boost
# add iostreams for gnuplot-iostream
find_package(Boost 1.46 COMPONENTS system iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# gnuplot-iostream header location (adjust if different)
include_directories(${CMAKE_SOURCE_DIR}/extern)

function(add_programs)
  foreach (prog ${ARGV})
    add_executable(${prog} ${prog}.cpp)
    target_link_libraries(${prog} ${BARRETT_LIBRARIES})
  endforeach()
endfunction()

add_programs(
  4dof_gravcomp_zeus
  4dof_dyncomp_zeus
  4dof_gravcomp_ares
  4dof_dyncomp_ares
  4dof_dyncompff_ares
  compare_torque_4dof
  compare_torque_7dof
  4dof_gravcomp_slax
  4dof_dyncomp_slax
)

# link Boost (iostreams) only to the plotting tool:
target_link_libraries(compare_torque_4dof ${Boost_LIBRARIES})
target_link_libraries(compare_torque_7dof ${Boost_LIBRARIES})
