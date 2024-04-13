#pragma once

#include "TsunamiProcessor.hpp"

namespace tsunami {
class FPGAProcessor : public TsunamiProcessor<float> {
public:
    struct Point {
        float e;
        float u;
        float v;
        float D;
    };
    enum Type { ZCU106, VC709 };

    void prepare(DataArrayPtr<float> bathymetry) override;
    void addDeformation(DataArrayPtr<float> deformation) override;
    void init(DataArrayPtr<float> e, DataArrayPtr<float> vx, DataArrayPtr<float> vy) override;
    void parameters(const json& params) override;
    void step(float dT, int& number) override;
    void finalize() override;

    DataArrayPtr<float> height() const override;

private:
    template <class T> inline void writeReg(const T& reg) { m_csr[T::Address()] = reg.pack(); }

    template <class T> inline T readReg()
    {
        T reg;
        reg.unpack(m_csr[T::Address()]);
        return reg;
    }

    void uploadData(uint64_t dst, void* src, size_t size) const;
    void downloadData(void* dst, uint64_t src, size_t size) const;
    void waitProcessor();

    Point* p;
    float *dx, dy;
    float grnd;
    int nx, ny;
    std::string m_devicePath;
    DataArrayPtr<float> m_outHeight;
    mutable bool m_outDataValid;

    Type m_type;
    int m_fd;
    uint32_t* m_csr;
    float* m_memory;
};
}
