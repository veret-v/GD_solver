cmake_minimum_required(VERSION 3.20)
project(cfd_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CFD_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
set(GMSH_DIR "" CACHE PATH "Path to Gmsh SDK root")

find_path(GMSH_INCLUDE_DIR
  NAMES gmsh.h
  HINTS ${GMSH_DIR} ENV GMSH_DIR
  PATH_SUFFIXES include Include
)

find_library(GMSH_LIBRARY
  NAMES gmsh libgmsh
  HINTS ${GMSH_DIR} ENV GMSH_DIR
  PATH_SUFFIXES lib Lib
)

if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
  message(FATAL_ERROR
    "Gmsh SDK not found.\n"
    "Set GMSH_DIR to the SDK root directory.\n"
    "Expected: <GMSH_DIR>/include/gmsh.h and <GMSH_DIR>/lib/gmsh(.lib/.a/.so).")
endif()

add_executable(cfd
  src/main.cpp
)

target_include_directories(cfd PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${GMSH_INCLUDE_DIR}
)

target_link_libraries(cfd PRIVATE ${GMSH_LIBRARY})

if(MSVC)
  target_compile_options(cfd PRIVATE /W4 /permissive-)
  if(CFD_WARNINGS_AS_ERRORS)
    target_compile_options(cfd PRIVATE /WX)
  endif()
else()
  target_compile_options(cfd PRIVATE -Wall -Wextra -Wpedantic)
  if(CFD_WARNINGS_AS_ERRORS)
    target_compile_options(cfd PRIVATE -Werror)
  endif()
endif()

if(WIN32)
  find_file(GMSH_DLL
    NAMES gmsh.dll libgmsh.dll
    HINTS ${GMSH_DIR} ENV GMSH_DIR
    PATH_SUFFIXES bin Bin lib Lib
  )
  if(GMSH_DLL)
    add_custom_command(TARGET cfd POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${GMSH_DLL}"
        "$<TARGET_FILE_DIR:cfd>"
    )
  endif()
endif()
