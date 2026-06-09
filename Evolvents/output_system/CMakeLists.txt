cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME_STR output_system)
project(${PROJECT_NAME_STR})

# Укажите исходные файлы (только .cpp)
file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Укажите заголовочные файлы
file(GLOB HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

# Создайте динамическую библиотеку
add_library(${PROJECT_NAME_STR} SHARED ${SRC_FILES} ${HEADER_FILES})

# Укажите директорию с заголовочными файлами
target_include_directories(${PROJECT_NAME_STR} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Настройка экспорта символов для Windows
target_compile_definitions(${PROJECT_NAME_STR} PRIVATE OUTPUT_SYSTEM_EXPORTS)

# Укажите директорию для вывода библиотеки в зависимости от конфигурации
set_target_properties(${PROJECT_NAME_STR} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/Debug
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/Debug
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/Debug
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/Release
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/Release
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/Release
)

# Определения для разных уровней вывода
if(${EXAMIN_OUTPUT_LEVEL} MATCHES "L1")
    target_compile_definitions(${PROJECT_NAME_STR} PUBLIC -DEXAMIN_OUTPUT_L1)
elseif(${EXAMIN_OUTPUT_LEVEL} MATCHES "L2")
    target_compile_definitions(${PROJECT_NAME_STR} PUBLIC -DEXAMIN_OUTPUT_L2)
endif()

# Определения для отладочной сборки
if(CMAKE_BUILD_TYPE MATCHES "DEBUG|Debug")
    target_compile_definitions(${PROJECT_NAME_STR} PUBLIC -D_DEBUG)
endif()

# Настройки для Microsoft компилятора
if (IS_MICROSOFT_COMPILER)
    set_property(TARGET ${PROJECT_NAME_STR} PROPERTY FOLDER "3rd-party")
endif()
