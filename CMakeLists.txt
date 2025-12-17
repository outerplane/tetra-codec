cmake_minimum_required(VERSION 3.20)

#------------------------------------------------------------------------------

project(tetra-codec)

include(GNUInstallDirs)

set(CMAKE_BUILD_TYPE Release)

#------------------------------------------------------------------------------

add_library(
	tetra-codec SHARED
	tetra-codec.c
	tetra-codec-impl.c
)

set_target_properties(
	tetra-codec PROPERTIES
	PUBLIC_HEADER "tetra-codec.h"
)

install(
	TARGETS tetra-codec
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

#------------------------------------------------------------------------------

set(TETRA_CODEC_PC "${CMAKE_CURRENT_BINARY_DIR}/tetra-codec.pc")

configure_file("tetra-codec.pc.in" "${TETRA_CODEC_PC}" @ONLY)

install(
	FILES "${TETRA_CODEC_PC}"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

#------------------------------------------------------------------------------
