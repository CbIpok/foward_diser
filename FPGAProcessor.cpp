#include "FPGAProcessor.hpp"
#include "TsunamiReg.hpp"

#include <chrono>
#include <cmath>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <thread>
#include <unistd.h>

#include "tsunami-ioctl.h"

#define GA 9.81

using namespace forward::Tsunami;

namespace tsunami {
const size_t MEMORY_SIZE = 2ULL * 1024 * 1024 * 1024;
const size_t DX_OFFSET = 2ULL * 1024 * 1024 * 1024 - 8192 * 4;

void FPGAProcessor::prepare(DataArrayPtr<float> bathymetry)
{
    TsunamiProcessor::prepare(bathymetry);

    nx = bathymetry->width();
    ny = bathymetry->height();

    if (m_type == ZCU106) {
        m_fd = open(m_devicePath.c_str(), O_RDWR | O_SYNC);
        m_csr = (uint32_t*)mmap(0, 4096, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_SYNC, m_fd, 4096);
        m_memory = (float*)mmap(0, MEMORY_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, m_fd, 0);
    } else if (m_type == VC709) {
        size_t size = nx * ny * sizeof(Point) + ny * sizeof(float);
        size = ((size - 1) / 4096 + 1) * 4096;

        m_fd = open(m_devicePath.c_str(), O_RDWR | O_SYNC);
        m_csr = (uint32_t*)mmap(0, 4096, PROT_READ | PROT_WRITE, MAP_SHARED, m_fd, 0);
        void* memptr;
        posix_memalign(&memptr, 4096, size);
        m_memory = (float*)memptr;

        m_csr[0x04] = 0x00000001;
    }

    TsunamiControl ctl = readReg<TsunamiControl>();
    std::cout << ctl.to_string() << std::endl;

    p = (Point*)m_memory;
    dx = (float*)m_memory + nx * ny * sizeof(Point) / 4;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            p[ij].D = bathymetry->data(i, j);
            p[ij].e = 0.0;
            p[ij].u = 0.0;
            p[ij].v = 0.0;
        }
    }

    if (bathymetry->type() == DataArray<float>::Spheric) {
        for (int j = 0; j < ny; j++)
            dx[j] = 1.0
                / (fabs(bathymetry->xCoord(1) - bathymetry->xCoord(0))
                    * cos(M_PI / 180.0 * bathymetry->yCoord(j)) * 111320.0);
        dy = 1.0 / (fabs(bathymetry->xCoord(1) - bathymetry->xCoord(0)) * 111320.0);
    } else {
        for (int j = 0; j < ny; j++)
            dx[j] = 1.0 / bathymetry->xStep();

        dy = 1.0 / bathymetry->yStep();
    }
    m_outHeight = DataArray<float>::create(
        nx, ny, bathymetry->xStep(), bathymetry->yStep(), bathymetry->type());
    m_outDataValid = false;
}

void FPGAProcessor::addDeformation(DataArrayPtr<float> deformation)
{
    TsunamiProcessor::addDeformation(deformation);
}

void FPGAProcessor::parameters(const json& params)
{
    grnd = params["ground"].get<float>();
    m_devicePath = params["device"];
    std::string type = params["fpga-type"].get<std::string>();
    if (type == "ZCU106")
        m_type = ZCU106;
    else if (type == "VC709")
        m_type = VC709;
    else
        m_type = ZCU106;
}

void FPGAProcessor::init(DataArrayPtr<float> h, DataArrayPtr<float> vx, DataArrayPtr<float> vy)
{
    TsunamiProcessor::init(h, vx, vy);

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            if (vx)
                p[ij].u = vx->data(i, j);
            else
                p[ij].u = 0.0f;
            if (vy)
                p[ij].v = vy->data(i, j);
            else
                p[ij].v = 0.0f;
            if (h)
                p[ij].e = h->data(i, j);
            else
                p[ij].e = 0.0f;
        }
    }
    TsunamiControl ctl;
    ctl.nx = nx;
    ctl.ny = ny;
    ctl.loadDx = 1;
    ctl.start = 0;
    writeReg(ctl);

    TsunamiDY dyr;
    dyr.dy = *((uint32_t*)&dy);
    writeReg(dyr);

    uploadData(0x00000000, p, nx * ny * sizeof(Point));
    uploadData(DX_OFFSET, dx, ny * sizeof(float));
}

void FPGAProcessor::finalize()
{
    if (m_type == ZCU106) {
        munmap(m_memory, MEMORY_SIZE);
        munmap(m_csr, 4096);
        close(m_fd);
    } else if (m_type == VC709) {
        munmap(m_csr, 4096);
        free(m_memory);
        close(m_fd);
    }

    TsunamiProcessor::finalize();
}

void FPGAProcessor::step(float dT, int& number)
{
    m_outDataValid = false;

    TsunamiDT dt;
    dt.dt = *((uint32_t*)&dT);
    writeReg(dt);

    TsunamiControl ctl = readReg<TsunamiControl>();
    ctl.start = 1;
    ctl.loadDx = (number == 0) ? 1 : 0;
    writeReg(ctl);

    waitProcessor();

    if (m_type == ZCU106)
        number += 3;
    else if (m_type == VC709)
        number += 7;
}

static void outputData(FPGAProcessor::Point* data, DataArrayPtr<float> out)
{
    for (int i = 0; i < out->height(); i++)
        for (int j = 0; j < out->width(); j++)
            out->data(j, i) = data[i * out->width() + j].e;
}

DataArrayPtr<float> FPGAProcessor::height() const
{
    if (!m_outDataValid) {
        downloadData(p, 0x00000000, nx * ny * sizeof(Point));
        outputData(p, m_outHeight);
        m_outDataValid = true;
    }
    return m_outHeight;
}

void FPGAProcessor::uploadData(uint64_t dst, void* src, size_t size) const
{
    if (m_type == ZCU106) {
        struct tsunamic_dma_xfer xfer;
        xfer.dst_address = dst;
        xfer.src_address = (uint8_t*)src - (uint8_t*)m_memory;
        xfer.size = size;

        auto start = std::chrono::high_resolution_clock::now();

        ioctl(m_fd, tsunamic_start_dma_tx, &xfer);
        ioctl(m_fd, tsunamic_wait_dma, &xfer.dma_handle);

        auto us = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start)
                      .count();

        std::cout << "Data uploaded in " << us << "us, " << size * 1000000ULL / us / 1024 / 1024
                  << "MiB/s" << std::endl;
    }
}

void FPGAProcessor::downloadData(void* dst, uint64_t src, size_t size) const
{
    if (m_type == ZCU106) {
        struct tsunamic_dma_xfer xfer;
        xfer.dst_address = (uint8_t*)dst - (uint8_t*)m_memory;
        xfer.src_address = src;
        xfer.size = size;

        auto start = std::chrono::high_resolution_clock::now();

        ioctl(m_fd, tsunamic_start_dma_rx, &xfer);
        ioctl(m_fd, tsunamic_wait_dma, &xfer.dma_handle);

        auto us = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start)
                      .count();

        std::cout << "Data downloaded in " << us << "us, " << size * 1000000ULL / us / 1024 / 1024
                  << "MiB/s" << std::endl;
    }
}

void FPGAProcessor::waitProcessor()
{
    uint32_t irq = 0;

    auto start = std::chrono::high_resolution_clock::now();
    if (m_type == ZCU106) {
        ioctl(m_fd, tsunamic_wait_interrupt, &irq);
    } else if (m_type == VC709) {
        TsunamiControl ctl;
        while (!ctl.start) {
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            ctl = readReg<TsunamiControl>();
        }
    }
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start)
                  .count();
    std::cout << "Proc done in " << us << "us, " << nx * ny * 1000000ULL / us << "Hz" << std::endl;
}
}
