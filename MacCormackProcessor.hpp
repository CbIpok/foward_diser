#pragma once

#include "TsunamiProcessor.hpp"

namespace tsunami {
template <class F> class MacCormackProcessor : public TsunamiProcessor<F> {
public:
    void prepare(DataArrayPtr<F> bathymetry) override;
    void addDeformation(DataArrayPtr<F> deformation) override;
    void init(DataArrayPtr<F> e, DataArrayPtr<F> vx, DataArrayPtr<F> vy) override;
    void parameters(const json& params) override;
    void step(F dT, int& number) override;
    void finalize() override;

    DataArrayPtr<F> height() const override;

private:
    F *u, *v, *e, *D, *dx, *dy;
    F *uh, *vh, *eh;
    F *uf, *vf, *ef;
    F* cosv;
    F grnd;
    int nx, ny;
    DataArrayPtr<F> m_outHeight;
    mutable bool m_outDataValid;
};
}
