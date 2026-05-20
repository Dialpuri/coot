#pragma once
#include <vector>
#include <gemmi/unitcell.hpp>

namespace coot {
namespace minimol {

struct molecule {
    gemmi::UnitCell unit_cell;
    bool have_cell = false;
    std::vector<float> mmdb_cell;

    std::vector<float> get_cell() const {
        std::vector<float> c;
        if (have_cell) {
            c = std::vector<float>(6);
            c[0] = unit_cell.a;
            c[1] = unit_cell.b;
            c[2] = unit_cell.c;
            c[3] = unit_cell.alpha;
            c[4] = unit_cell.beta;
            c[5] = unit_cell.gamma;
        }
        return c;
    }

    void set_cell(float a[6]) {
        mmdb_cell = std::vector<float>(6);
        for (int i = 0; i < 6; i++)
            mmdb_cell[i] = a[i];
        unit_cell.a = a[0];
        unit_cell.b = a[1];
        unit_cell.c = a[2];
        unit_cell.alpha = a[3];
        unit_cell.beta = a[4];
        unit_cell.gamma = a[5];
        have_cell = true;
    }
};

} // namespace minimol
} // namespace coot