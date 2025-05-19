cmake_minimum_required(VERSION 3.14)
# 这里可以不写 project(...)，CMake 会继承根目录的 project 名称和版本信息。
# 如果写了 project()，它会变成一个独立的子项目，但也不会影响 target 的创建。
project(FactoryLayoutOptimizerSrc LANGUAGES CXX)

# =================================================================================
# 1. 列出所有源文件和头文件（相对路径都从 src/ 开始）
# =================================================================================
set(SRC_FILES
    main.cpp                              # src/main.cpp
    ga/GeneticAlgorithmManager.cpp        # src/ga/GeneticAlgorithmManager.cpp
)

set(HEADER_FILES
    ga/GeneticAlgorithmManager.h          # src/ga/GeneticAlgorithmManager.h
    # 如果还有其他 .h 文件，也都列在这里

)

# =================================================================================
# 2. 定义可执行文件 FactoryLayoutOptimizer，包含所有 .cpp 和 .h
#    注意：头文件也要写在 add_executable 里，IDE 才能在项目树里把它们显示出来。
# =================================================================================
add_executable(FactoryLayoutOptimizer
    ${SRC_FILES}
    ${HEADER_FILES}
)

# =================================================================================
# 3. 将当前目录（src/）加入头文件搜索路径，方便 #include "ga/XXX.h"
#    这句话一定要在 add_executable() 之后，否则 CMake 会报“Cannot specify include directories for target …” 的错误。
# =================================================================================
target_include_directories(FactoryLayoutOptimizer
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# （注意：Qt5::Core 的链接在根目录的 CMakeLists.txt 里完成，这里不需要再写 target_link_libraries）
