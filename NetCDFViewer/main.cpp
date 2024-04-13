#include <QApplication>
#include <QSurfaceFormat>
#include <netcdf>

#include "MainWindow.h"

int main(int argc, char** argv)
{
    QSurfaceFormat format;
    format.setMajorVersion(3);
    format.setMinorVersion(2);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setOption(QSurfaceFormat::DebugContext);
    QSurfaceFormat::setDefaultFormat(format);

    qRegisterMetaType<netCDF::NcVar>("netCDF::NcVar");

    QApplication app(argc, argv);

    MainWindow w;
    w.show();

    if (app.arguments().length() == 2)
        w.openNetCDF(app.arguments()[1]);

    return app.exec();
}
