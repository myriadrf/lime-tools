find_package(PkgConfig REQUIRED)
pkg_search_module(TURBOFEC turbofec)
if (TURBOFEC_FOUND)
message(STATUS "TURBOFEC LIBRARY: " ${TURBOFEC_LIBRARIES})
message(STATUS "TURBOFEC INCLUDE DIR: " ${TURBOFEC_INCLUDE_DIR})
endif()

find_package(Sensors)
if (SENSORS_FOUND)
message(STATUS "SENSORS LIBRARY: " ${SENSORS_LIBRARY})
message(STATUS "SENSORS INCLUDE DIR: " ${SENSORS_INCLUDE_DIR})
endif()

include(CMakeDependentOption)
cmake_dependent_option(ENABLE_LTESCAN "Enable LTEscan" ON "TURBOFEC_FOUND;SENSORS_FOUND" OFF)

if (ENABLE_LTESCAN)
    set(LTESCAN_SOURCES
	src/CellItem.cpp
	src/earfcn.c
	src/LimeScanLTE.cpp
	src/pbch.cpp
	src/SDRlib.cpp
	src/CellList.cpp
	src/evm.cpp
	src/LTEscan.cpp
	src/PschSsch.cpp
	src/turbofecLib.c
    )
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
			${TURBOFEC_INCLUDE_DIR}
			)

    add_executable(LTEscan ${LTESCAN_SOURCES})

    target_link_libraries(LTEscan  
			${LIMESUITE_LIBRARIES}
			${FFTW3_LIBRARIES}
			${TURBOFEC_LIBRARIES}
			${SENSORS_LIBRARY}
			)

    install(TARGETS LTEscan DESTINATION bin)
endif()









