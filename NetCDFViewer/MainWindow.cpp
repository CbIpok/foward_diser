#include "MainWindow.h"

#include <QFileDialog>
#include <QInputDialog>
#include <QProgressDialog>
#include <cmath>
#include <future>

namespace {
const int POINT_EXPORT_CHUNK = 100;
}

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , m_ui(new Ui::MainWindow)
    , m_model(std::make_shared<LayerModel>())
{
    m_ui->setupUi(this);

    m_ui->varList->setModel(m_model.get());
    m_ui->varList->setItemDelegateForColumn(
        (int)LayerModel::Column::Shader, m_model->createDelegate(LayerModel::Column::Shader));

    m_ui->viewer->setModel(m_model);

    connect(m_ui->actionOpen, &QAction::triggered, [this] { openNetCDF(QString()); });
    connect(m_ui->varAdd, &QToolButton::clicked, this, &MainWindow::addVar);
    connect(m_ui->varRemove, &QToolButton::clicked, this, &MainWindow::removeVar);
    connect(m_ui->varUp, &QToolButton::clicked, [this] { moveVar(true); });
    connect(m_ui->varDown, &QToolButton::clicked, [this] { moveVar(false); });
    connect(m_ui->timeSlider, &QSlider::valueChanged, m_ui->viewer, &NetCDFWidget::updateTime);
    connect(m_ui->viewer, &NetCDFWidget::coordUpdated, this, &MainWindow::updateValues);
    connect(m_ui->viewer, &NetCDFWidget::pointExportRequested, this, &MainWindow::exportPointData);
    connect(
        m_ui->timeSlider, &QSlider::valueChanged, [this](int value) { updateValues(m_x, m_y); });
    connect(m_ui->displayLinearFilter, &QCheckBox::toggled, [this](bool value) {
        m_ui->viewer->setFiltering(value ? QOpenGLTexture::Linear : QOpenGLTexture::Nearest);
    });
    connect(m_ui->screenShot, &QPushButton::clicked, this, &MainWindow::takeScreenShot);

    connect(m_ui->dimTime, (void(QComboBox::*)(int)) & QComboBox::currentIndexChanged, [this] {
        auto tdim
            = m_dimensions.find(m_ui->dimTime->currentData().toString().toStdString())->second;
        m_ui->timeSlider->setMinimum(0);
        m_ui->timeSlider->setMaximum(tdim.getSize() - 1);
    });

    resize(1280, 720);
}

void MainWindow::openNetCDF(const QString& fn)
{
    QString file = fn;
    if (file.isEmpty()) {
        file = QFileDialog::getOpenFileName(
            this, tr("Open NetCDF file"), QString(), tr("NetCDF file (*.nc)"));
        if (file.isEmpty())
            return;
    }

    m_data = std::make_shared<netCDF::NcFile>(file.toStdString(), netCDF::NcFile::read);
    m_dimensions = m_data->getDims();
    m_vars = m_data->getVars();

    m_ui->dimHorizontal->clear();
    m_ui->dimVertical->clear();
    m_ui->dimTime->clear();

    for (auto dim : m_dimensions) {
        m_ui->dimHorizontal->addItem(
            QString::fromStdString(dim.first), QString::fromStdString(dim.first));
        m_ui->dimVertical->addItem(
            QString::fromStdString(dim.first), QString::fromStdString(dim.first));
        m_ui->dimTime->addItem(
            QString::fromStdString(dim.first), QString::fromStdString(dim.first));
    }
    m_ui->dimVertical->setCurrentIndex(0);
    m_ui->dimHorizontal->setCurrentIndex(1);
    m_ui->dimTime->setCurrentIndex(2);

    auto tdim = m_dimensions.find(m_ui->dimTime->currentData().toString().toStdString())->second;
    m_ui->timeSlider->setMinimum(0);
    m_ui->timeSlider->setMaximum(tdim.getSize() - 1);
    m_ui->timeSlider->setSingleStep(1);
    m_ui->timeSlider->setValue(0);

    m_ui->varSelect->clear();
    QFormLayout* layout = qobject_cast<QFormLayout*>(m_ui->valuesGroup->layout());
    for (auto var : m_vars) {
        m_ui->varSelect->addItem(
            QString::fromStdString(var.first), QString::fromStdString(var.first));

        QLabel* label = new QLabel();
        m_varViews[var.second.getName()] = label;
        layout->addRow(QString("%1:").arg(QString::fromStdString(var.second.getName())), label);
    }
}

void MainWindow::addVar()
{
    std::string name = m_ui->varSelect->currentData().toString().toStdString();
    netCDF::NcVar var = m_data->getVar(name);
    int xdim = -1, ydim = -1, tdim = -1;

    for (int dim = 0; dim < var.getDimCount(); dim++) {
        std::string name = var.getDim(dim).getName();
        if (name == m_ui->dimHorizontal->currentData().toString().toStdString())
            xdim = dim;
        else if (name == m_ui->dimVertical->currentData().toString().toStdString())
            ydim = dim;
        else if (name == m_ui->dimTime->currentData().toString().toStdString())
            tdim = dim;
    }

    m_model->addLayer(var, xdim, ydim, tdim);
}

void MainWindow::removeVar()
{
    QModelIndexList list = m_ui->varList->selectionModel()->selectedRows();
    for (QModelIndex idx : list)
        m_model->removeRow(idx.row());
}

void MainWindow::moveVar(bool up)
{
    QModelIndexList list = m_ui->varList->selectionModel()->selectedRows();
    for (QModelIndex idx : list)
        m_model->moveLayer(idx, up);
}

void MainWindow::updateValues(float x, float y)
{
    if (m_dimensions.empty() || m_vars.empty())
        return;

    m_x = x;
    m_y = y;

    int ix = (int)round(x);
    int iy = (int)round(y);
    int it = m_ui->timeSlider->value();

    netCDF::NcDim xdim
        = m_dimensions.find(m_ui->dimHorizontal->currentData().toString().toStdString())->second;
    netCDF::NcDim ydim
        = m_dimensions.find(m_ui->dimVertical->currentData().toString().toStdString())->second;
    netCDF::NcDim tdim
        = m_dimensions.find(m_ui->dimTime->currentData().toString().toStdString())->second;

    ix = std::clamp(ix, 0, (int)(xdim.getSize() - 1));
    iy = std::clamp(iy, 0, (int)(ydim.getSize() - 1));

    m_ui->valueIndex->setText(tr("(%1;%2;%3)").arg(ix).arg(iy).arg(it));

    for (auto var : m_vars) {
        if (m_varViews.find(var.first) == m_varViews.end())
            return;

        QLabel* view = m_varViews.at(var.first);
        netCDF::NcVar variable = var.second;

        std::vector<size_t> idx(variable.getDimCount(), 0);

        for (int i = 0; i < variable.getDimCount(); i++) {
            if (variable.getDim(i) == xdim)
                idx[i] = (size_t)ix;
            else if (variable.getDim(i) == ydim)
                idx[i] = (size_t)iy;
            else if (variable.getDim(i) == tdim)
                idx[i] = (size_t)it;
        }

        float value;
        variable.getVar(idx, &value);
        view->setText(QString("%1").arg(value));
    }
}

void MainWindow::exportPointData(float x, float y)
{
    if (m_dimensions.empty() || m_vars.empty())
        return;

    int ix = (int)round(x);
    int iy = (int)round(y);

    netCDF::NcDim xdim
        = m_dimensions.find(m_ui->dimHorizontal->currentData().toString().toStdString())->second;
    netCDF::NcDim ydim
        = m_dimensions.find(m_ui->dimVertical->currentData().toString().toStdString())->second;
    netCDF::NcDim tdim
        = m_dimensions.find(m_ui->dimTime->currentData().toString().toStdString())->second;

    ix = std::clamp(ix, 0, (int)(xdim.getSize() - 1));
    iy = std::clamp(iy, 0, (int)(ydim.getSize() - 1));

    QStringList varTexts;
    for (auto kv : m_vars)
        varTexts.append(QString::fromStdString(kv.first));

    std::string vname
        = QInputDialog::getItem(this, tr("Select variable"), tr("Variable"), varTexts, 0, false)
              .toStdString();

    if (vname.empty())
        return;

    netCDF::NcVar variable = m_vars.find(vname)->second;
    netCDF::NcVar tVariable = m_vars.find(tdim.getName())->second;

    std::vector<size_t> idx(variable.getDimCount(), 0);
    std::vector<size_t> size(variable.getDimCount(), 0);
    int ti = 0;

    for (int i = 0; i < variable.getDimCount(); i++) {
        if (variable.getDim(i) == xdim) {
            idx[i] = (size_t)ix;
            size[i] = 1;
        } else if (variable.getDim(i) == ydim) {
            idx[i] = (size_t)iy;
            size[i] = 1;
        } else if (variable.getDim(i) == tdim) {
            idx[i] = 0;
            size[i] = (size_t)tdim.getSize();
            ti = i;
        }
    }

    QString file = QFileDialog::getSaveFileName(
        this, tr("Save point data"), QString(), tr("Text file (*.txt)"));

    if (file.isEmpty())
        return;

    QFile dataFile(file);
    dataFile.open(QIODevice::WriteOnly);
    if (!dataFile.isOpen())
        return;

    QTextStream ts(&dataFile);

    double* data = new double[POINT_EXPORT_CHUNK];

    QProgressDialog progress(tr("Exporting data"), tr("Cancel"), 0, tdim.getSize() - 1, this);

    for (int i = 0; i < tdim.getSize(); i += POINT_EXPORT_CHUNK) {
        if (progress.wasCanceled())
            break;
        progress.setValue(i);
        qApp->processEvents();

        size[ti] = std::min((int)tdim.getSize() - i, POINT_EXPORT_CHUNK);
        idx[ti] = i;
        variable.getVar(idx, size, data);

        for (int i = 0; i < size[ti]; i++) {
            double t;
            tVariable.getVar({ idx[ti] + i }, &t);
            ts << QString("%1 %2\n").arg(t, 0, 'f').arg(data[i], 0, 'f');
        }
    }
    dataFile.close();
    progress.setValue(tdim.getSize() - 1);

    delete[] data;
}

void MainWindow::takeScreenShot()
{
    QString file = QFileDialog::getSaveFileName(
        this, tr("Save screenshot"), QString(), tr("PNG file (*.png);;JPEG file (*.jpg)"));

    if (file.isEmpty())
        return;

    if (QFileInfo(file).suffix().isEmpty())
        file += ".png";

    QImage img = m_ui->viewer->screenShot();

    auto result = std::async([img, file] { img.save(file); });
}
