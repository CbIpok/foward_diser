cmake_minimum_required(VERSION 3.8) 
project(TsunamiRunner)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

add_subdirectory(NetCDFViewer)
add_subdirectory(export_cdf)

set(SOURCES
    TsunamiRunner.cpp
    DataArray.cpp
    TsunamiProcessor.cpp
    MacCormackProcessor.cpp
    FPGAProcessor.cpp
    DataProcessors.cpp
)

set(HEADERS
    TsunamiRunner.hpp
    DataArray.hpp
    TsunamiProcessor.hpp
    MacCormackProcessor.hpp
    FPGAProcessor.hpp
    DataProcessors.hpp
    TsunamiReg.hpp
    json.hpp
)

add_executable(TsunamiRunner ${SOURCES} ${HEADERS})
target_link_libraries(TsunamiRunner PUBLIC OpenMP::OpenMP_CXX netcdf netcdf_c++4)
target_compile_options(TsunamiRunner PUBLIC "-mavx2")
install(TARGETS TsunamiRunner RUNTIME DESTINATION bin)
