#include "DataArray.hpp"
#include "DataProcessors.hpp"
#include "FPGAProcessor.hpp"
#include "MacCormackProcessor.hpp"
#include "TsunamiProcessor.hpp"
#include "json.hpp"

#include <chrono>
#include <csignal>
#include <fstream>
#include <future>

#include <netcdf>
#include <omp.h>

using json = nlohmann::json;
using namespace tsunami;

bool running = false;
void shutdown(int signal) { running = false; }

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Usage: TsunamiRunner <input>" << std::endl;
        return 0;
    }
    std::ifstream inf(argv[1]);

    if (!inf.is_open()) {
        std::cout << "Cannot read " << argv[1] << std::endl;
        return -1;
    }

    json input;
    inf >> input;

    if (!input.contains("bathymetry")) {
        std::cout << "No bathymetry" << std::endl;
        return -1;
    }

    DataArrayPtr<float> bath;
    int nx = 0;
    int ny = 0;

    json jsonBath = input["bathymetry"];
    if (jsonBath["type"] == "raw") {
        std::ifstream bf(jsonBath["file"]);
        nx = jsonBath["width"].get<int>();
        ny = jsonBath["height"].get<int>();
        double xStep = jsonBath["x-step"].get<double>();
        double yStep = jsonBath["y-step"].get<double>();
        bath = DataArray<float>::fromRawASCII(bf, nx, ny, xStep, yStep);
    } else if (jsonBath["type"] == "most") {
        std::ifstream bf(jsonBath["file"]);
        bath = DataArray<float>::fromMOST(bf);
        nx = bath->width();
        ny = bath->height();
    } else if (jsonBath["type"] == "flat") {
        nx = jsonBath["width"].get<int>();
        ny = jsonBath["height"].get<int>();
        double xStep = jsonBath["x-step"].get<double>();
        double yStep = jsonBath["y-step"].get<double>();
        double deep = jsonBath["deep"].get<double>();

        bath = DataArray<float>::create(nx, ny, xStep, yStep, DataArray<float>::Orthogonal);

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                bath->data(i, j) = deep;
    } else if (jsonBath["type"] == "slope") {
        nx = jsonBath["width"].get<int>();
        ny = jsonBath["height"].get<int>();
        double xStep = jsonBath["x-step"].get<double>();
        double yStep = jsonBath["y-step"].get<double>();
        double min_deep = jsonBath["min-deep"].get<double>();
        double max_deep = jsonBath["max-deep"].get<double>();

        bath = DataArray<float>::create(nx, ny, xStep, yStep, DataArray<float>::Orthogonal);

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                bath->data(i, j) = min_deep + (double)j / (ny - 1) * (max_deep - min_deep);
    }

    json run = input["run"];
    DataArrayPtr<float> e, vx, vy;
    if (run.contains("init-wave")) {
        json wave = run["init-wave"];
        if (wave["type"] == "most") {
            std::ifstream wf(wave["file"]);
            e = DataArray<float>::fromMOST(wf);
            if (e) {
                for (int i = 0; i < e->height(); i++) {
                    for (int j = 0; j < e->width(); j++) {
                        e->data(j, i) = -e->data(j, i);
                    }
                }
            }
        } else if (wave["type"] == "raw") {
            std::ifstream wf(wave["file"]);
            e = DataArray<float>::fromRawASCII(wf, nx, ny, 0.0, 0.0);
        } else if (wave["type"] == "gauss") {
            e = DataArray<float>::create(nx, ny, bath->xStep(), bath->yStep(), bath->type());
            int centerX = wave["center"][0].get<int>();
            int centerY = wave["center"][1].get<int>();
            int width = wave["size"][0].get<int>();
            int height = wave["size"][1].get<int>();
            double amplitude = wave["amplitude"].get<double>();
            double angle = 0.0;
            if (wave.contains("rotation"))
                angle = wave["rotation"].get<double>() * M_PI / 180.0;

            for (int i = 0; i < bath->height(); i++) {
                for (int j = 0; j < bath->width(); j++) {
                    double rx = j - centerX, ry = i - centerY;
                    double x = (rx * cos(angle) - ry * sin(angle)) / width;
                    double y = (rx * sin(angle) + ry * cos(angle)) / height;
                    double v = x * x + y * y;
                    if (v > 10)
                        continue;

                    e->data(j, i) = amplitude * exp(-v);
                }
            }
        }
    }

    if (run.contains("threads"))
        omp_set_num_threads(run["threads"].get<int>());

    if (run.contains("region")) {
        auto region = run["region"];
        int x = region[0].get<int>();
        int y = region[1].get<int>();
        int w = region[2].get<int>();
        int h = region[3].get<int>();
        bath = bath->region(x, y, w, h);
        e = e->region(x, y, w, h);
        nx = bath->width();
        ny = bath->height();
    }

    if (run.contains("scale")) {
        double scale = run["scale"].get<double>();
        bath = bath->scale(scale);
        e = e->scale(scale);
        nx = bath->width();
        ny = bath->height();
    }

    std::string name = run["name"].get<std::string>();
    double dT = run["timestep"].get<double>();
    int steps = run["steps"].get<int>();
    int saveInterval = 1;
    if (run.contains("save-interval"))
        saveInterval = run["save-interval"].get<int>();
    int saveStride = 1;
    if (run.contains("save-stride"))
        saveStride = run["save-stride"].get<int>();
    bool compression = true;
    if (run.contains("save-compression"))
        compression = run["save-compression"].get<bool>();

    std::unique_ptr<TsunamiProcessor<float>> proc = nullptr;

    if (run.contains("processor")) {
        if (run["processor"] == "maccormack")
            proc = std::unique_ptr<MacCormackProcessor<float>>(new MacCormackProcessor<float>());
        else if (run["processor"] == "fpga")
            proc = std::unique_ptr<FPGAProcessor>(new FPGAProcessor());
    } else {
        proc = std::unique_ptr<MacCormackProcessor<float>>(new MacCormackProcessor<float>());
    }

    for (auto pre : run["preprocess"]) {
        proc->addPreprocessor(DataProcessor<float>::fromJSON(pre));
    }
    for (auto post : run["postprocess"]) {
        proc->addPostprocessor(DataProcessor<float>::fromJSON(post));
    }

    proc->parameters(run);
    proc->prepare(bath);
    proc->init(e, vx, vy);

    std::string outDir = run.value("output-dir", "");

    if (!outDir.empty())
        outDir += '/';

    netCDF::NcFile nc(outDir + run["name"].get<std::string>() + "_h.nc", netCDF::NcFile::replace);

    std::string xname = "lon";
    std::string yname = "lat";

    if (bath->type() == DataArray<float>::Orthogonal) {
        xname = "x";
        yname = "y";
    }

    int ncnx = nx / saveStride, ncny = ny / saveStride;
    netCDF::NcDim ncX = nc.addDim(xname, ncnx);
    netCDF::NcDim ncY = nc.addDim(yname, ncny);
    netCDF::NcDim ncTime = nc.addDim("time");

    netCDF::NcVar ncXVar = nc.addVar(xname, netCDF::ncFloat, ncX);
    ncXVar.putVar(bath->xCoord());
    if (bath->type() == DataArray<float>::Orthogonal)
        ncXVar.putAtt("units", "meter");
    else
        ncXVar.putAtt("units", "degrees_east");

    netCDF::NcVar ncYVar = nc.addVar(yname, netCDF::ncFloat, ncY);
    ncYVar.putVar(bath->yCoord());
    if (bath->type() == DataArray<float>::Orthogonal)
        ncYVar.putAtt("units", "meter");
    else
        ncYVar.putAtt("units", "degrees_north");

    netCDF::NcVar ncTVar = nc.addVar("time", netCDF::ncFloat, ncTime);
    ncTVar.putAtt("units", "second");

    netCDF::NcVar ncHeight = nc.addVar("height", netCDF::ncFloat, { ncTime, ncY, ncX });
    std::vector<size_t> chunking = { 1, (size_t)ncny, (size_t)ncnx };
    ncHeight.setChunking(netCDF::NcVar::nc_CHUNKED, chunking);
    if (compression)
        ncHeight.setCompression(true, true, 5);

    netCDF::NcVar ncMaxHeight = nc.addVar("max_height", netCDF::ncFloat, { ncY, ncX });
    if (compression)
        ncMaxHeight.setCompression(true, true, 5);

    netCDF::NcVar bathymetry = nc.addVar("bathymetry", netCDF::ncFloat, { ncY, ncX });
    if (compression)
        bathymetry.setCompression(true, true, 5);
    bathymetry.putVar({ 0, 0 }, { (size_t)ncny, (size_t)ncnx }, bath->clone(saveStride)->data());

    DataArrayPtr<float> maxE
        = DataArray<float>::create(nx, ny, bath->xStep(), bath->yStep(), bath->type());

    DataArrayPtr<float> height = proc->height();
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            maxE->data(i, j) = height->data(i, j);

    signal(SIGINT, &shutdown);
    signal(SIGTERM, &shutdown);

    running = true;
    std::future<void> saveFuture;
    std::future<void> maxFuture;
    int step;
    for (step = 0; (step < steps) && running; step++) {
        auto start = std::chrono::high_resolution_clock::now();
        if ((step % saveInterval) == 0) {
            if (saveFuture.valid())
                saveFuture.wait();

            DataArrayPtr<float> height = proc->height();
            DataArrayPtr<float> data = height->clone(saveStride);
            saveFuture = std::async(std::launch::async | std::launch::deferred,
                [&nc, &ncHeight, &ncTVar, data, step, dT, ncnx, ncny, saveInterval] {
                    ncTVar.putVar({ (size_t)step / saveInterval }, step * dT);
                    ncHeight.putVar({ (size_t)step / saveInterval, 0, 0 },
                        { 1, (size_t)ncny, (size_t)ncnx }, data->data());
                    nc.sync();
                });
        }

        if (maxFuture.valid())
            maxFuture.wait();

        height = proc->height();
        maxFuture = std::async(std::launch::async | std::launch::deferred, [=] {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    if (maxE->data(i, j) < height->data(i, j))
                        maxE->data(i, j) = height->data(i, j);
                }
            }
        });

        proc->step(dT, step);

        auto time = std::chrono::high_resolution_clock::now() - start;
        std::cout << "Step " << step << "/" << steps << ": "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << "ms"
                  << std::endl;
    }

    if (saveFuture.valid())
        saveFuture.wait();

    if (maxFuture.valid())
        maxFuture.wait();

    ncTVar.putVar({ (size_t)step / saveInterval }, step * dT);
    ncHeight.putVar({ (size_t)step / saveInterval, 0, 0 }, { 1, (size_t)ncny, (size_t)ncnx },
        proc->height()->clone(saveStride)->data());
    ncMaxHeight.putVar({ 0, 0 }, { (size_t)ncny, (size_t)ncnx }, maxE->clone(saveStride)->data());

    std::ofstream maxf(outDir + run["name"].get<std::string>() + "_h_max.txt");
    maxE->toRawASCII(maxf);

    proc->finalize();

    nc.sync();

    return 0;
}
