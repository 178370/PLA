# src/CMakeLists.txt

# 查找 Qt5 Core 模块（用于 QCoreApplication、qInfo() 等控制台输出）
find_package(Qt5 REQUIRED COMPONENTS Core)

# 把当前目录（src/）加入头文件搜索路径，方便 #include "ga/GaOptimizer.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# 列出所有源文件（相对路径均从 src/ 开始）
set(SRC_FILES
    main.cpp              # src/main.cpp
    ga/GeneticAlgorithmManager.cpp    # src/ga/GaOptimizer.cpp
)

# 生成可执行文件
add_executable(FactoryLayoutOptimizer ${SRC_FILES})

# 链接 Qt5::Core 库
target_link_libraries(FactoryLayoutOptimizer PRIVATE Qt5::Core)
