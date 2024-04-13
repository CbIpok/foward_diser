#pragma once

#include <iostream>
#include <memory>
#include <vector>

namespace tsunami {
template <class F> class DataArray;

template <class F> using DataArrayPtr = std::shared_ptr<DataArray<F>>;

template <class F> class DataArray {
public:
    enum Type { Spheric, Orthogonal };

    Type type() const;

    int width() const;
    int height() const;

    const F& xCoord(int x) const;
    F& xCoord(int x);
    const F* xCoord() const;
    F* xCoord();

    const F& yCoord(int y) const;
    F& yCoord(int y);
    const F* yCoord() const;
    F* yCoord();

    inline const F& data(int x, int y) const { return m_data[y * m_width + x]; }
    inline F& data(int x, int y) { return m_data[y * m_width + x]; }

    const F* data() const;
    F* data();

    const F& xStep() const;
    F& xStep();
    const F& yStep() const;
    F& yStep();

    void toRawASCII(std::ostream& stream) const;
    void toASCII(std::ostream& stream) const;

    DataArrayPtr<F> clone(int decimate = 1) const;
    DataArrayPtr<F> region(int x, int y, int width, int height) const;
    DataArrayPtr<F> scale(double ratio) const;

    static DataArrayPtr<F> create(int width, int height, double xStep, double yStep, Type type);
    static DataArrayPtr<F> fromMOST(std::istream& stream);
    static DataArrayPtr<F> fromBinary(std::istream& stream);
    static DataArrayPtr<F> fromRawASCII(
        std::istream& stream, int width, int height, double xStep, double yStep);

private:
    DataArray();
    DataArray(const DataArray<F>&) = delete;
    DataArray<F>& operator=(const DataArray<F>&) = delete;

    Type m_type;
    int m_width, m_height;
    F m_xStep, m_yStep;
    std::vector<F> m_xCoord, m_yCoord;
    std::vector<F> m_data;
};

}
