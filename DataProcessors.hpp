#pragma once

#include "DataArray.hpp"
#include "json.hpp"

using json = nlohmann::json;

namespace tsunami {
template <class F> class DataProcessor {
public:
    DataProcessor(bool inplace)
        : m_inPlace(inplace)
    {
    }

    bool inPlace() const { return m_inPlace; }

    virtual void operator()(
        int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy) {};

    static std::shared_ptr<DataProcessor<F>> fromJSON(const json& node);

protected:
    bool m_inPlace;
};

template <class F> class BorderSinWaveProcessor : public DataProcessor<F> {
public:
    BorderSinWaveProcessor(F amplitude, F period, bool top, bool bottom, bool left, bool right);
    void operator()(
        int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy) override;

private:
    F m_amplitude, m_period;
    bool m_top, m_bottom, m_left, m_right;
};

template <class F> class OffshoreSmoothWave : public DataProcessor<F> {
public:
    OffshoreSmoothWave(int period, F grnd, F bound, F coeff);
    void operator()(
        int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy) override;

private:
    int m_period;
    F m_grnd, m_bound, m_coeff;
};

template <class F> class LimitWaveOnDeep : public DataProcessor<F> {
public:
    LimitWaveOnDeep(int period, F deep, F grnd);
    void operator()(
        int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy) override;

private:
    int m_period;
    F m_deep, m_grnd;
};

template <class F> class RemoveSingleIslands : public DataProcessor<F> {
public:
    RemoveSingleIslands(F grnd);
    void operator()(
        int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy) override;

private:
    F m_grnd;
};

}
