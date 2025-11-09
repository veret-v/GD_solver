cmake_minimum_required(VERSION 3.20)
project(GodunovSolver)  # ← Измените название проекта

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Добавьте все исходные файлы
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/parser.cpp
    src/point.cpp
    src/grid.cpp
    src/solver.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
)

# Уберите эту строку - она не нужна, если все файлы уже добавлены выше
# add_subdirectory("src")