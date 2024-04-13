#pragma once

#include <memory>
#include <QMainWindow>
#include <netcdf>

#include "LayerModel.h"
#include "ui_MainWindow.h"
#include <algorithm>

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);

public Q_SLOTS:
    void openNetCDF(const QString& fn);

private Q_SLOTS:
    void addVar();
    void removeVar();
    void moveVar(bool up);
    void updateValues(float x, float y);
    void exportPointData(float x, float y);
    void takeScreenShot();

private:
    std::unique_ptr<Ui::MainWindow> m_ui;
    std::shared_ptr<netCDF::NcFile> m_data;
    std::multimap<std::string, netCDF::NcDim> m_dimensions;
    std::multimap<std::string, netCDF::NcVar> m_vars;
    std::shared_ptr<LayerModel> m_model;
    float m_x = 0.0, m_y = 0.0;
    std::map<std::string, QLabel*> m_varViews;
};
