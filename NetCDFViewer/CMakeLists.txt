set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 20)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(netcdf-viewer
    main.cpp
    MainWindow.h
    MainWindow.cpp
    MainWindow.ui
    NetCDFWidget.h
    NetCDFWidget.cpp
    LayerModel.h
    LayerModel.cpp
    Shaders.qrc
)
find_package(Threads REQUIRED)
target_link_libraries(netcdf-viewer Qt5::Widgets netcdf netcdf_c++4 ${CMAKE_THREAD_LIBS_INIT})
