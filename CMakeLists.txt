cmake_minimum_required(VERSION 3.14)
project(FactoryLayoutOptimizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =================================================================================
# 1. 查找 Qt5 Core 模块
# =================================================================================
find_package(Qt5 REQUIRED COMPONENTS Core)

# =================================================================================
# 2. 将 src 目录加入工程
# =================================================================================
#    这一步会去执行 src/CMakeLists.txt，并在子目录里创建名为 FactoryLayoutOptimizer 的可执行 target。
# =================================================================================
add_subdirectory(src)

# =================================================================================
# 3. 全局添加 src/ 目录到头文件搜索路径（可选，方便跨子目录引用）
# =================================================================================
include_directories(${PROJECT_SOURCE_DIR}/src)

# =================================================================================
# 4. 链接 Qt5::Core 到在 src 子目录里定义的 FactoryLayoutOptimizer target
#    （注意：因为 add_subdirectory(src) 已经在 src/CMakeLists.txt 里 add_executable 了，
#     所以此处可以直接引用 “FactoryLayoutOptimizer” 这个 target）
# =================================================================================
target_link_libraries(FactoryLayoutOptimizer PRIVATE Qt5::Core)
