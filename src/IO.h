#ifndef IO_H
#define IO_H

#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Config.h"
#include "Grid.h"

/**
    Read initial conditions from file

    TODO: some exception checking for proper file format?
 */
template <typename T>
std::pair<VelocityField<T>, PressureField<T>>
read_from_file(std::string fname) {
    std::ifstream ifile;
    ifile.open(fname);
    std::string line;
    // read first line
    std::getline(ifile, line, ',');
    size_t width = std::stoi(line);
    std::getline(ifile, line, ',');
    size_t height = std::stoi(line);
    std::getline(ifile, line, ',');
    size_t depth = std::stoi(line);
    std::getline(ifile, line, '\n');
    T resolution = std::stod(line);

    auto grid = std::make_shared<Grid<T>>(width, height, depth, resolution);
    PressureField<T> pressures(grid);
    VelocityField<T> velocities(grid);

    for (size_t z = 0; z < grid->getDepth(); ++z) {
        std::getline(ifile, line, '\n');
        assert(line == "");
        for (size_t y = 0; y < grid->getHeight(); ++y) {
            std::getline(ifile, line, '\n');
            std::istringstream line_stream(line);
            std::string val;
            for (size_t x = 0; x < grid->getWidth(); ++x) {
                std::getline(line_stream, val, ',');
                T p = std::stod(val);
                pressures.setPressure(x, y, z, p);
                std::getline(line_stream, val, ',');
                T u = std::stod(val);
                velocities.setLeftU(x, y, z, u);
                std::getline(line_stream, val, ',');
                T v = std::stod(val);
                velocities.setTopV(x, y, z, v);
                std::getline(line_stream, val, ',');
                T w = std::stod(val);
                velocities.setFrontW(x, y, z, w);
            }
        }
    }
    ifile.close();

    return std::make_pair(velocities, pressures);
}

/* Documentation of raw output file:
 * First line: grid dimensions followed by resolution
 * Followed By:
 * Comma seperated list of values
 * For each field:pressure, x,y,z velocities
 * each line lists all field in x direction,
 * the lines are written in y direction
 * each z entry is seperated by blank lines between "y blocks"
 */

template <typename T>
void write_to_raw_file(const VelocityField<T>& velocities,
                       const PressureField<T>& pressures, std::string fname) {
    assert(velocities.getDepth() == pressures.getDepth());
    assert(velocities.getHeight() == pressures.getHeight());
    assert(velocities.getWidth() == pressures.getWidth());

    std::ofstream ofile;
    ofile.open(fname);

    ofile << velocities.getWidth() << "," << velocities.getHeight() << ","
          << velocities.getDepth() << "," << velocities.getCellSize() << "\n\n";
    for (size_t z = 0; z < velocities.getDepth(); ++z) {
        for (size_t y = 0; y < velocities.getHeight(); ++y) {
            for (size_t x = 0; x < velocities.getWidth(); ++x) {
                ofile << pressures.getPressure(x, y, z) << ","
                      << velocities.getLeftU(x, y, z) << ","
                      << velocities.getTopV(x, y, z) << ","
                      << velocities.getFrontW(x, y, z) << ",";
            }
            ofile << "\n";
        }
        ofile << "\n";
    }
    ofile.close();
}

template <typename T>
void write_to_csv_file(const VelocityField<T>& velocities,
                       const PressureField<T>& pressures, std::string fname) {

    assert(velocities.getDepth() == pressures.getDepth());
    assert(velocities.getHeight() == pressures.getHeight());
    assert(velocities.getWidth() == pressures.getWidth());

    auto ds = velocities.getCellSize();

    std::ofstream ofile;
    ofile.open(fname);

    const char delimiter = ',';

    ofile << "x" << delimiter << "y" << delimiter << "z" << delimiter << "u"
          << delimiter << "v" << delimiter << "w" << delimiter << "p\n";

    for (size_t z = 0; z < velocities.getDepth(); ++z) {
        for (size_t y = 0; y < velocities.getHeight(); ++y) {
            for (size_t x = 0; x < velocities.getWidth(); ++x) {
                // Interpolate velocities at cell center
                auto vCenter = velocities.trilerp(
                    {(x + 0.5) * ds, (y + 0.5) * ds, (z + 0.5) * ds});
                ofile << x << delimiter << y << delimiter << z << delimiter
                      << vCenter.x << delimiter << vCenter.y << delimiter
                      << vCenter.z << delimiter
                      << pressures.getPressure(x, y, z) << "\n";
            }
        }
    }
    ofile.close();
}

/**
    Write state (velocity field, pressure field) to file
 */
template <typename T>
void write_to_file(const VelocityField<T>& velocities,
                   const PressureField<T>& pressures, std::string prefix,
                   unsigned step, FileFormat format) {
    switch (format) {
    case FileFormat::CSV: {
        std::string filename = prefix + ".csv." + std::to_string(step);
        write_to_csv_file(velocities, pressures, filename);
        break;
    }
    case FileFormat::RAW: {
        std::string filename = prefix + ".txt." + std::to_string(step);
        write_to_raw_file(velocities, pressures, filename);
        break;
    }
    default:
        break;
    }
}

#endif // IO_H
