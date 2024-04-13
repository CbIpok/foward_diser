#include "LayerModel.h"

#include <QComboBox>

static const std::vector<std::string> SHADERS = { { "null", "jet", "geo", "water", "depth" } };

LayerModel::LayerModel(QObject* parent)
    : QAbstractTableModel(parent)
{
}

Qt::ItemFlags LayerModel::flags(const QModelIndex& index) const
{
    Qt::ItemFlags f = Qt::NoItemFlags;

    switch ((Column)index.column()) {
    case Column::Shader:
    case Column::Min:
    case Column::Max:
        f |= Qt::ItemIsEditable;
        break;
    }

    return f /*| Qt::ItemIsDropEnabled | Qt::ItemIsDragEnabled*/ | Qt::ItemIsEnabled
        | Qt::ItemIsSelectable;
}

QVariant LayerModel::data(const QModelIndex& index, int role) const
{
    if ((index.row() >= m_order.size()) || (index.row() < 0))
        return QVariant();

    auto& var = m_layers.at(m_order[index.row()]);

    if (role == Qt::DisplayRole) {
        switch ((Column)index.column()) {
        case Column::Name:
            return QString::fromStdString(var.var.getName());
        case Column::Shader:
            return QString::fromStdString(var.shader);
        case Column::Min:
            return var.min;
        case Column::Max:
            return var.max;
        default:
            return QVariant();
        }
    }

    return QVariant();
}

bool LayerModel::setData(const QModelIndex& index, const QVariant& value, int role)
{
    if ((index.row() >= m_order.size()) || (index.row() < 0))
        return false;

    auto& var = m_layers.at(m_order[index.row()]);
    Column col = (Column)index.column();

    if (col == Column::Shader)
        var.shader = value.toString().toStdString();
    else if (col == Column::Min)
        var.min = value.toFloat();
    else if (col == Column::Max)
        var.max = value.toFloat();
    else
        return false;

    emit(dataChanged(index, index));

    return true;
}

QVariant LayerModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role != Qt::DisplayRole)
        return QVariant();

    if (orientation == Qt::Vertical) {
        if ((section >= m_order.size()) || (section < 0))
            return QVariant();

        return m_layers.at(m_order.at(section)).id;
    } else {
        switch ((Column)section) {
        case Column::Name:
            return tr("Variable");
        case Column::Shader:
            return tr("Shader");
        case Column::Min:
            return tr("Min");
        case Column::Max:
            return tr("Max");
        default:
            return QVariant();
        }
    }
    return QVariant();
}

int LayerModel::rowCount(const QModelIndex& parent) const { return m_layers.size(); }

int LayerModel::columnCount(const QModelIndex& parent) const { return (int)Column::Count; }

bool LayerModel::removeRows(int row, int count, const QModelIndex& parent)
{
    beginRemoveRows(QModelIndex(), row, row + count - 1);

    auto it = m_order.begin() + row;
    for (int i = 0; i < count; i++) {
        auto& layer = m_layers.at(*it);
        it = m_order.erase(it);
        m_layers.erase(layer.id);
    }

    endRemoveRows();

    return true;
}

void LayerModel::addLayer(const netCDF::NcVar& variable, int xdim, int ydim, int tdim)
{
    beginInsertRows(QModelIndex(), m_layers.size(), m_layers.size());
    Layer v;

    static std::atomic_int id = 0;

    v.id = ++id;
    v.var = variable;
    v.xdim = xdim;
    v.ydim = ydim;
    v.tdim = tdim;
    v.min = 0.0;
    v.max = 1.0;
    v.shader = "null";

    m_layers[v.id] = v;
    m_order.push_back(v.id);

    endInsertRows();
}

std::map<int, LayerModel::Layer>& LayerModel::layers() { return m_layers; }

const std::map<int, LayerModel::Layer>& LayerModel::layers() const { return m_layers; }

const LayerModel::Layer& LayerModel::layer(int id) const { return m_layers.at(id); }

LayerModel::Layer& LayerModel::layer(int id) { return m_layers[id]; }

const LayerModel::Layer& LayerModel::layer(const QModelIndex& index) const
{
    int id = m_order.at(index.row());
    return m_layers.at(id);
}

LayerModel::Layer& LayerModel::layer(const QModelIndex& index)
{
    int id = m_order.at(index.row());
    return m_layers[id];
}

std::vector<int> LayerModel::order() const { return m_order; }

QAbstractItemDelegate* LayerModel::createDelegate(Column col)
{
    if (col == Column::Shader) {
        return new Delegate(col);
    }
    return nullptr;
}

void LayerModel::moveLayer(const QModelIndex& index, bool up)
{
    if ((up && (index.row() <= 0)) || (!up && (index.row() >= m_order.size() - 1)))
        return;

    beginMoveRows(
        QModelIndex(), index.row(), index.row(), QModelIndex(), index.row() + (up ? -1 : 2));
    std::iter_swap(m_order.begin() + index.row(), m_order.begin() + index.row() + (up ? -1 : 1));
    endMoveRows();
}

LayerModel::Delegate::Delegate(LayerModel::Column column)
    : m_column(column)
{
}

QWidget* LayerModel::Delegate::createEditor(
    QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
    QComboBox* box = new QComboBox(parent);
    for (auto shader : SHADERS) {
        auto name = QString::fromStdString(shader);
        box->addItem(name, name);
        if (name == index.model()->data(index).toString())
            box->setCurrentIndex(box->count() - 1);
    }
    return box;
}

void LayerModel::Delegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{
    QComboBox* box = qobject_cast<QComboBox*>(editor);
    QString value = index.model()->data(index, Qt::DisplayRole).toString();
    box->setCurrentText(value);
}

void LayerModel::Delegate::setModelData(
    QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{
    QComboBox* box = qobject_cast<QComboBox*>(editor);
    QString value = box->currentData().toString();
    model->setData(index, value);
}
