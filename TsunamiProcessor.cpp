#include "TsunamiProcessor.hpp"

namespace tsunami {
template <class F> void TsunamiProcessor<F>::prepare(DataArrayPtr<F> bathymetry)
{
    m_bathymetry = bathymetry;
}

template <class F> void TsunamiProcessor<F>::addDeformation(DataArrayPtr<F> deformation) { }

template <class F>
void TsunamiProcessor<F>::init(DataArrayPtr<F> e, DataArrayPtr<F> vx, DataArrayPtr<F> vy)
{
}

template <class F> void TsunamiProcessor<F>::finalize() { }

template <class F> void TsunamiProcessor<F>::parameters(const json& params) { }

template <class F> DataArrayPtr<F> TsunamiProcessor<F>::bathymetry() const { return m_bathymetry; }

template <class F> void TsunamiProcessor<F>::addPreprocessor(std::shared_ptr<DataProcessor<F>> proc)
{
    m_preProcessors.push_back(proc);
}

template <class F>
void TsunamiProcessor<F>::addPostprocessor(std::shared_ptr<DataProcessor<F>> proc)
{
    m_postProcessors.push_back(proc);
}

template class TsunamiProcessor<float>;
template class TsunamiProcessor<double>;
}
