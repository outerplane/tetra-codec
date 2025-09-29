cmake_minimum_required(VERSION 2.8)

project(tetra-codec)

set(CMAKE_BUILD_TYPE Release)

add_library(tetra-codec SHARED tetra-codec.c tetra-codec-impl.c)

IF(WIN32)
set_target_properties(tetra-codec PROPERTIES PREFIX "")
ENDIF()

install(TARGETS tetra-codec DESTINATION lib64)
install(FILES tetra-codec.h DESTINATION include)
install(FILES tetra-codec.pc DESTINATION lib64/pkgconfig)
