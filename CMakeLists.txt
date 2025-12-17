cmake_minimum_required(VERSION 3.20)

#------------------------------------------------------------------------------

project(tetra-codec)

include(GNUInstallDirs)

set(CMAKE_BUILD_TYPE Release)

#------------------------------------------------------------------------------

if (NOT DEFINED TETRA_CODEC_LIBRARY_TYPE)
set(TETRA_CODEC_LIBRARY_TYPE SHARED)
endif ()

add_library(
	tetra-codec ${TETRA_CODEC_LIBRARY_TYPE}
	source/tetra-codec.c
	source/tetra-codec-impl.c
)

target_include_directories(
	tetra-codec PUBLIC
	include
)

if (TETRA_CODEC_ENABLE_LTO)
set_property(
	TARGET tetra-codec
	PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
)
endif ()

set_target_properties(
	tetra-codec PROPERTIES
	PUBLIC_HEADER "include/tetra-codec.h"
)

install(
	TARGETS tetra-codec
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

#------------------------------------------------------------------------------

set(TETRA_CODEC_PC "${CMAKE_CURRENT_BINARY_DIR}/tetra-codec.pc")

configure_file("pkgconfig/tetra-codec.pc.in" "${TETRA_CODEC_PC}" @ONLY)

install(
	FILES "${TETRA_CODEC_PC}"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

#------------------------------------------------------------------------------
