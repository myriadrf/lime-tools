find_package(PkgConfig)

pkg_check_modules (FFTW3 "fftw3 >= 3.0")

find_path(FFTW3_INCLUDE_DIR
            NAMES fftw3.h
            PATHS ${FFTW3_INCLUDE_DIRS} 
		  /usr/local/include 
                  /usr/include )

find_library(FFTW3_LIBRARIES
            NAMES fftw3
            PATHS ${FFTW3_LIBRARY_DIRS}
	          /usr/local/lib
                  /usr/lib)

find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARIES)
