#pragma once

#include "DataArray.hpp"
#include "DataProcessors.hpp"
#include "json.hpp"

using json = nlohmann::json;

namespace tsunami {
template <class F> class TsunamiProcessor {
public:
    virtual void prepare(DataArrayPtr<F> bathymetry);
    virtual void addDeformation(DataArrayPtr<F> deformation);
    virtual void init(DataArrayPtr<F> e, DataArrayPtr<F> vx, DataArrayPtr<F> vy);
    virtual void parameters(const json& params);
    virtual void step(F dT, int& number) = 0;
    virtual void finalize();

    virtual DataArrayPtr<F> bathymetry() const;
    virtual DataArrayPtr<F> height() const = 0;

    void addPreprocessor(std::shared_ptr<DataProcessor<F>> proc);
    void addPostprocessor(std::shared_ptr<DataProcessor<F>> proc);

protected:
    std::vector<std::shared_ptr<DataProcessor<F>>> m_preProcessors;
    std::vector<std::shared_ptr<DataProcessor<F>>> m_postProcessors;

private:
    DataArrayPtr<F> m_bathymetry;
};
}
