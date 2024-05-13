
file(GLOB SOURCE_V0 *.cpp)
add_executable(n_body ${SOURCE_V0})
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(n_body PUBLIC "C:\\Program Files (x86)\\Intel\\oneAPI\\advisor\\latest\\include")
