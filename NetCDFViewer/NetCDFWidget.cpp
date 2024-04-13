#include "NetCDFWidget.h"
#include "LayerModel.h"

#include <QMenu>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLDebugLogger>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QVector3D>
#include <QWheelEvent>
#include <cmath>
#include <vector>

static const std::map<std::string, std::tuple<std::string, std::string>> SHADERS = {
    { "null", { ":/null.vert", ":/null.frag" } },
    { "jet", { ":/null.vert", ":/jet.frag" } },
    { "geo", { ":/null.vert", ":/geo.frag" } },
    { "water", { ":/null.vert", ":/water.frag" } },
    { "depth", { ":/null.vert", ":/depth.frag" } },
};

NetCDFWidget::NetCDFWidget(QWidget* parent)
    : QOpenGLWidget(parent)
{
    setMouseTracking(true);
}

void NetCDFWidget::setModel(std::shared_ptr<LayerModel> model)
{
    m_model = model;

    connect(m_model.get(), &QAbstractItemModel::rowsInserted, this, &NetCDFWidget::addLayer);
    connect(
        m_model.get(), &QAbstractItemModel::rowsAboutToBeRemoved, this, &NetCDFWidget::removeLayer);
    connect(m_model.get(), &QAbstractItemModel::dataChanged, [this]() { update(); });
    connect(m_model.get(), &QAbstractItemModel::rowsMoved, [this]() { update(); });
}

void NetCDFWidget::addLayer(const QModelIndex& parent, int first, int last)
{
    makeCurrent();
    auto gl = context()->functions();

    for (int i = first; i <= last; i++) {
        auto& layer = m_model->layer(m_model->index(i, 0, parent));
        LayerData d;

        netCDF::NcDim xd = layer.var.getDim(layer.xdim);
        netCDF::NcDim yd = layer.var.getDim(layer.ydim);

        d.textureBuffer = QOpenGLBuffer(QOpenGLBuffer::PixelUnpackBuffer);
        d.textureBuffer.create();
        d.textureBuffer.bind();
        d.textureBuffer.allocate(xd.getSize() * yd.getSize() * sizeof(float));
        d.textureBuffer.setUsagePattern(QOpenGLBuffer::StreamDraw);
        d.textureBuffer.release();
        d.rect = QRectF(0.0, 0.0, xd.getSize(), yd.getSize());

        d.texture = std::make_unique<QOpenGLTexture>(QOpenGLTexture::Target2D);
        d.texture->setSize(xd.getSize(), yd.getSize());
        d.texture->setFormat(QOpenGLTexture::R32F);
        d.texture->allocateStorage();
        d.texture->setMinMagFilters(m_filter, m_filter);

        std::vector<std::array<float, 5>> vertices {
            { -0.5f, -0.5f, 0.0, 0.0, 0.0 },
            { (float)xd.getSize() - 0.5f, 0.0, 0.0, 1.0, 0.0 },
            { 0.0, (float)yd.getSize() - 0.5f, 0.0, 0.0, 1.0 },
            { (float)xd.getSize() - 0.5f, (float)yd.getSize() - 0.5f, 0.0, 1.0, 1.0 },
        };

        d.vao = std::make_unique<QOpenGLVertexArrayObject>();
        d.vao->create();
        d.vao->bind();

        d.vertexBuffer.create();
        d.vertexBuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
        d.vertexBuffer.bind();
        d.vertexBuffer.allocate(vertices.data(), vertices.size() * sizeof(float) * 5);

        gl->glEnableVertexAttribArray(0);
        gl->glEnableVertexAttribArray(1);
        gl->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 5, (void*)0);
        gl->glVertexAttribPointer(
            1, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 5, (void*)(sizeof(float) * 3));

        d.vertexBuffer.release();
        d.vao->release();

        m_layers[layer.id] = std::move(d);

        updateLayerT(layer.id, 0);
    }
}

void NetCDFWidget::removeLayer(const QModelIndex& parent, int first, int last)
{
    auto gl = context()->functions();

    for (int i = first; i <= last; i++) {
        auto& layer = m_model->layer(m_model->index(i, 0, parent));
        m_layers.erase(layer.id);
    }
    update();
}

void NetCDFWidget::setFiltering(QOpenGLTexture::Filter filter)
{
    m_filter = filter;
    for (auto& d : m_layers)
        d.second.texture->setMinMagFilters(filter, filter);

    update();
}

QImage NetCDFWidget::screenShot()
{
    makeCurrent();
    auto gl = context()->functions();

    QRectF rect;
    for (auto& l : m_layers)
        rect = rect.united(l.second.rect);

    QSize size = rect.size().toSize();

    QOpenGLFramebufferObject fb(
        rect.size().toSize(), QOpenGLFramebufferObject::NoAttachment, GL_TEXTURE_2D, GL_RGBA);
    fb.bind();
    gl->glViewport(0, 0, size.width(), size.height());

    QMatrix4x4 transform = QMatrix4x4();
    transform.ortho(rect.left(), rect.right(), rect.bottom(), rect.top(), -1000.0f, 1000.0f);

    for (auto& shader : m_shaders) {
        shader.second->bind();
        shader.second->setUniformValue("transform", transform);
        shader.second->release();
    }

    paintGL();

    QImage img = fb.toImage();
    fb.release();

    recalcMatrices();
    update();

    return img;
}

void NetCDFWidget::updateTime(int t)
{
    for (auto& l : m_layers)
        updateLayerT(l.first, t);
}

void NetCDFWidget::initializeGL()
{
    QOpenGLWidget::initializeGL();
    auto gl = context()->functions();

    QOpenGLDebugLogger* logger = new QOpenGLDebugLogger(this);
    logger->initialize();
    logger->startLogging();

    connect(logger, &QOpenGLDebugLogger::messageLogged,
        [](const QOpenGLDebugMessage& debugMessage) { qDebug() << debugMessage; });

    for (auto sources : SHADERS) {
        std::unique_ptr<QOpenGLShaderProgram> shader = std::make_unique<QOpenGLShaderProgram>();
        shader->create();
        shader->addShaderFromSourceFile(
            QOpenGLShader::Vertex, QString::fromStdString(std::get<0>(sources.second)));
        shader->addShaderFromSourceFile(
            QOpenGLShader::Fragment, QString::fromStdString(std::get<1>(sources.second)));
        shader->bindAttributeLocation("inPosition", 0);
        shader->bindAttributeLocation("inTexCoord", 1);
        shader->link();

        shader->bind();
        shader->setUniformValue("tex", (GLuint)0);
        shader->release();

        if (shader->isLinked())
            m_shaders[sources.first] = std::move(shader);
    }

    gl->glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    gl->glEnable(GL_BLEND);

    recalcMatrices();
}

void NetCDFWidget::paintGL()
{
    QOpenGLWidget::paintGL();
    auto gl = context()->functions();

    gl->glClear(GL_COLOR_BUFFER_BIT);

    if (!m_model)
        return;

    auto order = m_model->order();
    for (auto id : order) {
        auto& layer = m_model->layer(id);
        auto& data = m_layers[layer.id];
        auto& shader = m_shaders[layer.shader];

        data.texture->bind();
        shader->bind();
        shader->setUniformValue("min", layer.min);
        shader->setUniformValue("max", layer.max);
        data.vao->bind();
        gl->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        data.vao->release();
        shader->bind();
        data.texture->release();
    }
}

void NetCDFWidget::resizeGL(int w, int h)
{
    QOpenGLWidget::resizeGL(w, h);

    m_projectionTransform = QMatrix4x4();
    m_projectionTransform.ortho(0.0f, w, h, 0.0f, -1000.0f, 1000.0f);

    recalcMatrices();
}

void NetCDFWidget::mouseMoveEvent(QMouseEvent* ev)
{
    QVector3D viewPortCoord = { (float)ev->pos().x() * 2.0f / width() - 1.0f,
        1.0f - (float)ev->pos().y() * 2.0f / height(), 0.0f };

    if (ev->buttons() & Qt::LeftButton) {
        QVector3D delta = (m_mouseTransform * viewPortCoord) - (m_mouseTransform * m_lastMousePos);
        m_translate += delta;
        recalcMatrices();
        update();
    }

    QVector3D pos = m_mouseTransform * viewPortCoord;
    emit(coordUpdated(pos.x(), pos.y()));

    m_lastMousePos = viewPortCoord;
}

void NetCDFWidget::wheelEvent(QWheelEvent* e)
{
    float amount = e->angleDelta().y() / 120.0;
    QVector3D viewPortCoord;//{ (float)e->position().x() * 2.0f / width() - 1.0f,
//        1.0f - (float)e->position().y() * 2.0f / height(), 0.0f };

    QVector3D center = m_mouseTransform * viewPortCoord;
    m_scale *= pow(1.25, amount);
    recalcMatrices();

    QVector3D newCenter = m_mouseTransform * viewPortCoord;
    m_translate += newCenter - center;
    recalcMatrices();

    update();
}

void NetCDFWidget::contextMenuEvent(QContextMenuEvent* e)
{
    QVector3D viewPortCoord = { (float)e->pos().x() * 2.0f / width() - 1.0f,
        1.0f - (float)e->pos().y() * 2.0f / height(), 0.0f };
    QVector3D pos = m_mouseTransform * viewPortCoord;

    QMenu menu(this);
    menu.addAction(
        "Export data at point", [this, pos] { emit(pointExportRequested(pos.x(), pos.y())); });
    menu.exec(e->globalPos());
}

void NetCDFWidget::recalcMatrices()
{
    makeCurrent();

    m_worldTransform = QMatrix4x4();
    m_worldTransform.scale(m_scale);
    m_worldTransform.translate(m_translate);

    QMatrix4x4 mvp = m_projectionTransform * m_worldTransform;
    m_mouseTransform = mvp.inverted();
    for (auto& shader : m_shaders) {
        shader.second->bind();
        shader.second->setUniformValue("transform", mvp);
        shader.second->release();
    }
}

void NetCDFWidget::updateLayerT(int id, int t)
{
    LayerModel::Layer& layer = m_model->layer(id);
    auto& d = m_layers.at(id);

    if ((layer.tdim < 0) && (t != 0))
        return;

    d.textureBuffer.bind();

    float* data = (float*)d.textureBuffer.map(QOpenGLBuffer::WriteOnly);

    std::vector<size_t> start(layer.var.getDimCount(), 0);
    std::vector<size_t> count(layer.var.getDimCount(), 1);

    if (layer.tdim >= 0)
        start[layer.tdim] = t;

    count[layer.xdim] = layer.var.getDim(layer.xdim).getSize();
    count[layer.ydim] = layer.var.getDim(layer.ydim).getSize();

    layer.var.getVar(start, count, data);
    d.textureBuffer.unmap();

    d.texture->setData(QOpenGLTexture::Red, QOpenGLTexture::Float32, (const void*)nullptr,
        (QOpenGLPixelTransferOptions*)nullptr);

    d.textureBuffer.release();

    update();
}
