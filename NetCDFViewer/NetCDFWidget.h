#pragma once

#include <memory>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <memory>
#include <netcdf>

class LayerModel;

class NetCDFWidget : public QOpenGLWidget {
    Q_OBJECT
public:
    NetCDFWidget(QWidget* parent = nullptr);

    void setModel(std::shared_ptr<LayerModel> model);
    QImage screenShot();

Q_SIGNALS:
    void coordUpdated(float x, float y);
    void pointExportRequested(float x, float y);

public Q_SLOTS:
    void updateTime(int t);
    void setFiltering(QOpenGLTexture::Filter filter);

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;
    void mouseMoveEvent(QMouseEvent* ev) override;
    void wheelEvent(QWheelEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* event) override;

private Q_SLOTS:
    void addLayer(const QModelIndex& parent, int first, int last);
    void removeLayer(const QModelIndex& parent, int first, int last);

private:
    void recalcMatrices();

    struct LayerData {
        std::unique_ptr<QOpenGLVertexArrayObject> vao;
        QOpenGLBuffer vertexBuffer;
        QOpenGLBuffer textureBuffer;
        std::unique_ptr<QOpenGLTexture> texture;
        QRectF rect;
    };
    void updateLayerT(int id, int t);

    std::map<std::string, std::unique_ptr<QOpenGLShaderProgram>> m_shaders;
    QMatrix4x4 m_projectionTransform, m_worldTransform, m_mouseTransform;
    std::map<int, LayerData> m_layers;
    std::shared_ptr<LayerModel> m_model;
    float m_scale = 1.0f;
    QVector3D m_translate;
    QVector3D m_lastMousePos;
    QOpenGLTexture::Filter m_filter = QOpenGLTexture::Linear;
};
