cmake_minimum_required(VERSION 3.8)

project(NetCDFExporter)

set(SOURCES
    NetCDFExporter.cpp
)

set(HEADERS
)

add_executable(netcdf-exporter ${SOURCES} ${HEADERS})
target_link_libraries(netcdf-exporter PUBLIC netcdf netcdf_c++4)
