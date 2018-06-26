find_package(PkgConfig)

pkg_check_modules (LIMESUITE LimeSuite)

find_path(LIMESUITE_INCLUDE_DIR
	NAMES LimeSuite.h
	PATHS ${LIMESUITE_INCLUDE_DIRS}
  	/usr/include/lime
  	/usr/local/include/lime
)

find_library(LIMESUITE_LIBRARIES 
	NAMES LimeSuite
	PATHS ${LIMESUITE_LIBRARY_DIRS}
  	/usr/lib
  	/usr/local/lib
)

find_package_handle_standard_args(LimeSuite DEFAULT_MSG LIMESUITE_LIBRARIES LIMESUITE_INCLUDE_DIR)

mark_as_advanced(LIMESUITE_LIBRARIES LIMESUITE_INCLUDE_DIR)

