vcpkg_from_github(
	OUT_SOURCE_PATH SOURCE_PATH
	REPO outerplane/tetra-codec
	REF d933efed57e5588ba1b606b1e38bf9b93af24666
	SHA512 f2703a498224b953dfa0a1a55cff0089c1fce534af620c5fbd92d2142c7449177f5eef3ba4fadea9bd86c0f1651884ad8dc4e7a47e16e8daf705454c68ef8cec
	HEAD_REF master
)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})
vcpkg_cmake_install()
vcpkg_fixup_pkgconfig()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/README.md")
