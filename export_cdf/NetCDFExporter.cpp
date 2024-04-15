#include <netcdf>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>

using namespace netCDF;

int main(int argc, char **argv)
{

    std::cout<<argc<<std::endl;

    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <netcdf_file> <variable> <output_file> [timestep]" << std::endl;
        return 0;
    }

    NcFile nc(argv[1], netCDF::NcFile::read);
    NcDim x = nc.getDim("lon");
    NcDim y = nc.getDim("lat");
    NcDim time = nc.getDim("time");


    if (x.isNull() || y.isNull()) {
        x = nc.getDim("x");
        y = nc.getDim("y");
    }
    if (x.isNull() || y.isNull() || time.isNull()) {
        std::cerr << "Wrong NetCDF file - invalid dimensions" << std::endl;
        return -1;
    }

    NcVar var = nc.getVar(argv[2]);
    if (var.isNull()) {
        std::cerr << "No such var" << std::endl;
        return -1;
    }


    double *data = new double[x.getSize() * y.getSize()];
    if (var.getDimCount() == 3) {
        if (argc < 5) {
            std::cerr << "Timestep required" << std::endl;
            return -1;
        }
        int time = std::atoi(argv[4]);
        var.getVar({(size_t)time, (size_t)0, (size_t)0}, {1, (size_t)y.getSize(), (size_t)x.getSize()}, data);
    } else if (var.getDimCount() == 2) {
        var.getVar({(size_t)0, (size_t)0}, {(size_t)y.getSize(), (size_t)x.getSize()}, data);
    } else {
        std::cerr << "Invalid variable" << std::endl;
        return -1;
    }

    unsigned long numTimesteps = 4500;
    int px_start = std::stoi(argv[5]);
    int px_stop = std::stoi(argv[6]);
    int py_start = std::stoi(argv[7]);
    int py_stop = std::stoi(argv[8]);
    std::string type = argv[9];
    std::vector<std::vector<double>> mariograms(35);
    for (int time_c = 0;time_c < time.getSize();time_c++)
    {
        std::cout << time_c << std::endl;
        var.getVar({(size_t)time_c, (size_t)0, (size_t)0}, {1, (size_t)y.getSize(), (size_t)x.getSize()}, data);
        int index = 0;
        for (unsigned long py = py_start;py < py_stop;py += 150)
        {
            for (unsigned long px = px_start;px < px_stop;px += 100)
            {
            {
                mariograms[index++].push_back(data[px + py*x.getSize()]);
            }
        }
    }
    }


    int index = 0;
    for (unsigned long py = py_start;py < py_stop;py += 150)
    {

        for (unsigned long px = px_start;px < px_stop;px += 100)
        {
            auto mar =  mariograms[index];
            std::ofstream outmf(type + "_" + std::to_string(px) + "_" + std::to_string(py));
            for (int time_c = 0;time_c < time.getSize();time_c++)
            {
                outmf << std::setprecision(std::numeric_limits<double>::digits10 + 1) << mar[time_c] << std::endl;
                if (mar[time_c] != 0)
                {
                    std::cout << index << " " << mar[time_c] << std::endl;
                }

            }
            index++;
        }


    }




    //    for (unsigned long px = px_start;px < px_stop;px += 100)
    //    {
    //        for (unsigned long py = py_start;py < px_stop;py += 150)
    //        {
    //            std::cout << px << " " << py << std::endl;
    //            var.getVar({0, py, px}, {time.getSize(), 1, 1}, data);
    //            std::ofstream outmf(type +"_" + std::to_string((px)) + "_" + std::to_string(py));
    //            for (int i = 0;i < time.getSize();i++)
    //            {
    //                outmf << std::setprecision(std::numeric_limits<double>::digits10 + 1) << data[i] << " ";
    //                outmf << std::endl;
    //            }
    //        }
    //    }



    // std::ofstream outf(argv[3]);
    // for (int i = 0; i < y.getSize(); i++)
    //  {
    //for (int j = 0; j < x.getSize(); j++)
    //    outf << std::setprecision(std::numeric_limits<double>::digits10 + 1) << data[i * x.getSize() + j] << " ";
    //  outf << std::endl;
    //}



    delete [] data;

    return 0;
}
