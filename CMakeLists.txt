cmake_minimum_required(VERSION 3.10)
project(OptimizationProject)

find_library(M_LIB m)

set(CMAKE_CXX_STANDARD 17)

# Укажите исходные файлы основного проекта
set(SOURCE_FILES
    Global_method/Global_method.cpp
    Hill/HillProblem.cpp
    Shekel/ShekelProblem.cpp
    Hill/HillConstrainedProblem.cpp
    Shekel/ShekelConstrainedProblem.cpp
    Grishagin/GrishaginOption.cpp
    Grishagin/grishagin_function.cpp
    GKLS/gkls_random.cpp
    GKLS/GKLSOption.cpp
    GKLS/GKLSProblem.cpp
    GKLS/GKLSConstrainedProblem.cpp
    src/IOptProblem.cpp
    src/IOptProblemFamily.cpp
    src/IGeneralOptProblem.cpp
    src/IGeneralOptProblemFamily.cpp
    src/IConstrainedOptProblem.cpp
    src/IConstrainedOptProblemFamily.cpp
    src/OptSqConstrProblem.hpp  # Добавлено для включения файла с задачей
    direct_lib/DIRect.c
    direct_lib/DIRserial.c
    direct_lib/DIRsubrout.c
    direct_lib/direct_wrap.c
)

#add_subdirectory(${CMAKE_SOURCE_DIR}/Evolvents/output_system)

add_executable(OptimizationProject ${SOURCE_FILES})


target_include_directories(OptimizationProject PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Hill"
    "${CMAKE_CURRENT_SOURCE_DIR}/Shekel"
    "${CMAKE_CURRENT_SOURCE_DIR}/GKLS"
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/Grishagin"
    "${CMAKE_CURRENT_SOURCE_DIR}/direct_lib"  # Добавляем путь к direct.h
    #"${CMAKE_CURRENT_SOURCE_DIR}/Evolvents/output_system/include" # Добавлено для Evolvents
)

# Копирование DLL в папку с исполняемым файлом
#add_custom_command(TARGET OptimizationProject POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#    "${CMAKE_BINARY_DIR}/Evolvents/output_system/Release/output_system.dll"
#    $<TARGET_FILE_DIR:OptimizationProject>)
