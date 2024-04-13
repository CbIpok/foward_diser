#pragma once

#include <memory>
#include <QAbstractTableModel>
#include <QStyledItemDelegate>
#include <netcdf>

class LayerModel : public QAbstractTableModel {
    Q_OBJECT

public:
    struct Layer {
        int id;
        netCDF::NcVar var;
        int xdim = -1, ydim = -1, tdim = -1;
        std::string shader = "null";
        float min = 0.0f, max = 1.0f;
    };
    enum class Column { Name = 0, Shader, Min, Max, Count };

    LayerModel(QObject* parent = nullptr);

    Qt::ItemFlags flags(const QModelIndex& index) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    int rowCount(const QModelIndex& parent) const override;
    int columnCount(const QModelIndex& parent) const override;
    bool setData(const QModelIndex& index, const QVariant& value, int role) override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
    bool removeRows(int row, int count, const QModelIndex& parent) override;

    void addLayer(const netCDF::NcVar& variable, int xdim, int ydim, int tdim = -1);

    QAbstractItemDelegate* createDelegate(Column col);

    std::map<int, Layer>& layers();
    const std::map<int, Layer>& layers() const;
    const Layer& layer(int id) const;
    Layer& layer(int id);
    const Layer& layer(const QModelIndex& index) const;
    Layer& layer(const QModelIndex& index);

    std::vector<int> order() const;

    void moveLayer(const QModelIndex& index, bool up);

private:
    class Delegate;

    std::map<int, Layer> m_layers;
    std::vector<int> m_order;
};

class LayerModel::Delegate : public QStyledItemDelegate {
    Q_OBJECT

public:
    Delegate(LayerModel::Column column);

    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option,
        const QModelIndex& index) const override;
    void setEditorData(QWidget* editor, const QModelIndex& index) const override;
    void setModelData(
        QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const override;
    // void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const
    // QModelIndex &index) const override;

private:
    LayerModel::Column m_column;
};
