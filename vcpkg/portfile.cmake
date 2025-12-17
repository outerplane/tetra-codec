vcpkg_from_github(
	OUT_SOURCE_PATH SOURCE_PATH
	REPO outerplane/tetra-codec
	REF 2604d923fccefe109e3138b8f3d873876a870105
	SHA512 a8812bd0a470be61d6d1c18aef38be89291181d9de3e4b1457c6026bc0652bb5c70c88211305eb89d9ba9a6ada4fd223bdf632b91fe5fed001bfbb3af5643439
	HEAD_REF master
)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})
vcpkg_cmake_install()
vcpkg_fixup_pkgconfig()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/README.md")
