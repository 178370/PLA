# src/CMakeLists.txt
cmake_minimum_required(VERSION 3.14)

# =================================================================================
# 1. 列出需要编译/链接的源文件和头文件
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
# =================================================================================
add_executable(FactoryLayoutOptimizer
    ${SRC_FILES}
    ${HEADER_FILES}
)

# =================================================================================
# 3. 将当前目录（src/）加入头文件搜索路径
# =================================================================================
target_include_directories(FactoryLayoutOptimizer
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# =================================================================================
# 4. 列出所有要复制的 CSV 文件名（它们和 main.cpp 同级）
# =================================================================================
set(CSV_FILES
    device_size.csv
    Pij.csv
    Qij.csv
    hij.csv
    hi0.csv
    de.csv
)

# =================================================================================
# 5. 确保 build/src/Debug 目录存在（configure_file 不会自动创建子目录）
#    ${CMAKE_CURRENT_BINARY_DIR} 通常是 build/src（因为当前 CMakeLists 在 src/ 下）
#    所以这里创建 build/src/Debug
# =================================================================================
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Debug")

# =================================================================================
# 6. 遍历 CSV_FILES，把它们从 src/ 复制到 build/src/Debug/
#    使用 configure_file(... COPYONLY) 保证每次 run cmake 时都刷新一次
# =================================================================================
foreach(_csv IN LISTS CSV_FILES)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${_csv}"         # 源：src/device_size.csv 等
        "${CMAKE_CURRENT_BINARY_DIR}/Debug/${_csv}"   # 目标：build/src/Debug/device_size.csv
        COPYONLY
    )
endforeach()
