#include "DataArray.hpp"
#include <cmath>

namespace tsunami {

template <class F>
DataArray<F>::DataArray()
    : m_width(0)
    , m_height(0)
{
}

template <class F>
DataArrayPtr<F> DataArray<F>::create(int width, int height, double xStep, double yStep, Type type)
{
    DataArrayPtr<F> b = std::shared_ptr<DataArray<F>>(new DataArray<F>());

    b->m_type = type;

    b->m_width = width;
    b->m_height = height;

    b->m_data = std::vector<F>(b->m_width * b->m_height);

    b->m_xStep = xStep;
    b->m_yStep = yStep;

    b->m_xCoord = std::vector<F>(b->m_width);
    for (int i = 0; i < width; i++)
        b->m_xCoord[i] = i * xStep;

    b->m_yCoord = std::vector<F>(b->m_height);
    for (int i = 0; i < height; i++)
        b->m_yCoord[i] = i * yStep;

    return b;
}

template <class F> DataArrayPtr<F> DataArray<F>::fromMOST(std::istream& stream)
{
    int width, height;
//    char s[512];
//    stream.getline(s,512);
    stream >> width >> height;

    DataArrayPtr<F> b = create(width, height, 0.0, 0.0, Spheric);

    for (int i = 0; i < width; i++)
        stream >> b->xCoord(i);

    for (int i = 0; i < height; i++)
        stream >> b->yCoord(i);

    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            stream >> b->data(j, i);

    return b;
}

template <class F> void DataArray<F>::toASCII(std::ostream& stream) const
{
    stream << m_width << " " << m_height << std::endl;

    for (int i = 0; i < m_width; i++)
        stream << xCoord(i) << std::endl;

    for (int i = 0; i < m_height; i++)
        stream << yCoord(i) << std::endl;

    for (int i = 0; i < m_height; i++) {
        for (int j = 0; j < m_width; j++) {
            stream << data(j, i) << " ";
        }
        stream << std::endl;
    }
}

template <class F> DataArrayPtr<F> DataArray<F>::clone(int decimate) const
{
    if (decimate < 1)
        decimate = 1;

    int nwidth = m_width / decimate;
    int nheight = m_height / decimate;

    DataArrayPtr<F> result
        = create(nwidth, nheight, m_xStep * decimate, m_yStep * decimate, m_type);

    for (size_t i = 0; i < nwidth; i++)
        result->m_xCoord[i] = m_xCoord[i * decimate];

    for (size_t i = 0; i < nheight; i++)
        result->m_yCoord[i] = m_yCoord[i * decimate];

    for (size_t y = 0; y < nheight; y++) {
        for (size_t x = 0; x < nwidth; x++)
            result->m_data[y * nwidth + x] = m_data[y * decimate * m_width + x * decimate];
    }

    return result;
}

template <class F> DataArrayPtr<F> DataArray<F>::region(int x, int y, int width, int height) const
{
    DataArrayPtr<F> result = create(width, height, m_xStep, m_yStep, m_type);

    for (size_t i = 0; i < width; i++)
        result->m_xCoord[i] = m_xCoord[x + i];

    for (size_t i = 0; i < height; i++)
        result->m_yCoord[i] = m_yCoord[y + i];

    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++)
            result->m_data[i * width + j] = m_data[(y + i) * m_width + x + j];
    }

    return result;
}

template <class F> DataArrayPtr<F> DataArray<F>::scale(double ratio) const
{
    int nwidth = m_width * ratio;
    int nheight = m_height * ratio;

    DataArrayPtr<F> result = create(nwidth, nheight, m_xStep / ratio, m_yStep / ratio, m_type);

    for (size_t i = 0; i < nwidth; i++) {
        double x = i / ratio;
        int x1 = (int)floor(x), x2 = x1 + 1;
        F Q1, Q2;

        if (x1 >= m_width)
            Q1 = m_xCoord[m_width - 1] + (m_xCoord[m_width - 1] - m_xCoord[m_width - 2]);
        else
            Q1 = m_xCoord[x1];

        if (x2 >= m_width)
            Q2 = m_xCoord[m_width - 1] + (m_xCoord[m_width - 1] - m_xCoord[m_width - 2]);
        else
            Q2 = m_xCoord[x2];

        result->m_xCoord[i] = (x2 - x) * Q1 + (x - x1) * Q2;
    };

    for (size_t i = 0; i < nheight; i++) {
        double y = i / ratio;
        int y1 = (int)floor(y), y2 = y1 + 1;
        F Q1, Q2;

        if (y1 >= m_height)
            Q1 = m_yCoord[m_height - 1] + (m_yCoord[m_height - 1] - m_yCoord[m_height - 2]);
        else
            Q1 = m_yCoord[y1];

        if (y2 >= m_height)
            Q2 = m_yCoord[m_height - 1] + (m_yCoord[m_height - 1] - m_yCoord[m_height - 2]);
        else
            Q2 = m_yCoord[y2];

        result->m_yCoord[i] = (y2 - y) * Q1 + (y - y1) * Q2;
    };

    for (size_t i = 0; i < nheight; i++) {
        double y = i / ratio;
        int y1 = (int)floor(y), y2 = y1 + 1;

        for (size_t j = 0; j < nwidth; j++) {
            double x = j / ratio;
            int x1 = (int)floor(x), x2 = x1 + 1;
            F Q11, Q12, Q21, Q22;

            if (x1 >= m_width) {
                if (y1 >= m_height)
                    Q11 = m_data[(m_height - 1) * m_width + m_width - 1];
                else
                    Q11 = m_data[y1 * m_width + m_width - 1];
            } else {
                if (y1 >= m_height)
                    Q11 = m_data[(m_height - 1) * m_width + x1];
                else
                    Q11 = m_data[y1 * m_width + x1];
            }

            if (x2 >= m_width) {
                if (y1 >= m_height)
                    Q21 = m_data[(m_height - 1) * m_width + m_width - 1];
                else
                    Q21 = m_data[y1 * m_width + m_width - 1];
            } else {
                if (y1 >= m_height)
                    Q21 = m_data[(m_height - 1) * m_width + x2];
                else
                    Q21 = m_data[y1 * m_width + x2];
            }

            if (x1 >= m_width) {
                if (y2 >= m_height)
                    Q12 = m_data[(m_height - 1) * m_width + m_width - 1];
                else
                    Q12 = m_data[y2 * m_width + m_width - 1];
            } else {
                if (y2 >= m_height)
                    Q12 = m_data[(m_height - 1) * m_width + x1];
                else
                    Q12 = m_data[y2 * m_width + x1];
            }

            if (x2 >= m_width) {
                if (y2 >= m_height)
                    Q22 = m_data[(m_height - 1) * m_width + m_width - 1];
                else
                    Q22 = m_data[y2 * m_width + m_width - 1];
            } else {
                if (y2 >= m_height)
                    Q22 = m_data[(m_height - 1) * m_width + x2];
                else
                    Q22 = m_data[y2 * m_width + x2];
            }

            result->m_data[i * nwidth + j] = (x2 - x) * (y2 - y) * Q11 + (x - x1) * (y2 - y) * Q21
                + (x2 - x) * (y - y1) * Q12 + (x - x1) * (y - y1) * Q22;
        }
    }
    return result;
}

template <class F> DataArrayPtr<F> DataArray<F>::fromBinary(std::istream& stream)
{
    return std::shared_ptr<DataArray>(new DataArray());
}

template <class F>
DataArrayPtr<F> DataArray<F>::fromRawASCII(
    std::istream& stream, int width, int height, double xStep, double yStep)
{
    DataArrayPtr<F> b = create(width, height, xStep, yStep, Orthogonal);

    for (int i = 0; i < width; i++)
        b->m_xCoord[i] = i * xStep;

    for (int i = 0; i < height; i++)
        b->m_yCoord[i] = i * yStep;

    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            stream >> b->data(j, i);

    return b;
}

template <class F> void DataArray<F>::toRawASCII(std::ostream& stream) const
{
    for (int i = 0; i < m_height; i++) {
        for (int j = 0; j < m_width; j++)
            stream << m_data[i * m_width + j] << " ";
        stream << std::endl;
    }
}

template <class F> typename DataArray<F>::Type DataArray<F>::type() const { return m_type; }

template <class F> int DataArray<F>::width() const { return m_width; }

template <class F> int DataArray<F>::height() const { return m_height; }

template <class F> const F& DataArray<F>::xCoord(int x) const { return m_xCoord[x]; }

template <class F> F& DataArray<F>::xCoord(int x) { return m_xCoord[x]; }

template <class F> const F* DataArray<F>::xCoord() const { return m_xCoord.data(); }

template <class F> F* DataArray<F>::xCoord() { return m_xCoord.data(); }

template <class F> const F& DataArray<F>::yCoord(int y) const { return m_yCoord[y]; }

template <class F> F& DataArray<F>::yCoord(int y) { return m_yCoord[y]; }

template <class F> const F* DataArray<F>::yCoord() const { return m_yCoord.data(); }

template <class F> F* DataArray<F>::yCoord() { return m_yCoord.data(); }

template <class F> const F* DataArray<F>::data() const { return m_data.data(); }

template <class F> F* DataArray<F>::data() { return m_data.data(); }

template <class F> F& DataArray<F>::xStep() { return m_xStep; }

template <class F> const F& DataArray<F>::xStep() const { return m_xStep; }

template <class F> F& DataArray<F>::yStep() { return m_yStep; }

template <class F> const F& DataArray<F>::yStep() const { return m_yStep; }

template class DataArray<float>;
template class DataArray<double>;
}
