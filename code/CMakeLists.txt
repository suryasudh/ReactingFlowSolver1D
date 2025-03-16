cmake_minimum_required(VERSION 3.10) # Specify the minimum CMake version required

project(JsonParserExample) # Set the name of the project

# Specify the C++ standard to use (optional, but recommended)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a library target for the json_parser code
add_library(json_parser json_parser.cpp json_parser.h)
# 'json_parser' is the name of the library target
# 'json_parser.cpp' is the source file for the library
# 'json_parser.h' is the header file for the library

# Create an executable target for the main program
add_executable(main001 main01.cpp)
# 'my_program' is the name of the executable
# 'main.cpp' is the source file for the executable

# Link the executable with the json_parser library
target_link_libraries(main001 json_parser)
# This tells CMake to link the 'my_program' executable with the 'json_parser' library

# Optional: Install targets (if you want to install the executable and library)
# install(TARGETS my_program DESTINATION bin)
# install(TARGETS json_parser DESTINATION lib)
# install(FILES json_parser.h DESTINATION include)