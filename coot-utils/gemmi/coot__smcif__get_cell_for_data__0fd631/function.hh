#pragma once

#include <clipper/core/cell.h>
#include <clipper/core/clipper_util.h>
#include <gemmi/cifdoc.hpp>
#include <string>

namespace coot {
namespace smcif {

clipper::Cell get_cell_for_data_gemmi(const gemmi::cif::Block& data) {
    clipper::Cell cell;

    double a = 0, b = 0, c = 0;
    double alpha = 0, beta = 0, gamma = 0;
    bool ok = true;

    const std::string* val_a = data.find_value("_cell_length_a");
    if (!val_a || val_a->empty()) {
        ok = false;
    } else {
        try { a = std::stod(*val_a); } catch (...) { ok = false; }
    }

    if (ok) {
        const std::string* val_b = data.find_value("_cell_length_b");
        if (!val_b || val_b->empty()) {
            ok = false;
        } else {
            try { b = std::stod(*val_b); } catch (...) { ok = false; }
        }
    }

    if (ok) {
        const std::string* val_c = data.find_value("_cell_length_c");
        if (!val_c || val_c->empty()) {
            ok = false;
        } else {
            try { c = std::stod(*val_c); } catch (...) { ok = false; }
        }
    }

    if (ok) {
        const std::string* val_alpha = data.find_value("_cell_angle_alpha");
        if (!val_alpha || val_alpha->empty()) {
            ok = false;
        } else {
            try { alpha = std::stod(*val_alpha); } catch (...) { ok = false; }
        }
    }

    if (ok) {
        const std::string* val_beta = data.find_value("_cell_angle_beta");
        if (!val_beta || val_beta->empty()) {
            ok = false;
        } else {
            try { beta = std::stod(*val_beta); } catch (...) { ok = false; }
        }
    }

    if (ok) {
        const std::string* val_gamma = data.find_value("_cell_angle_gamma");
        if (!val_gamma || val_gamma->empty()) {
            ok = false;
        } else {
            try { gamma = std::stod(*val_gamma); } catch (...) { ok = false; }
        }
    }

    if (ok) {
        clipper::Cell_descr cell_descr(a, b, c,
                                       clipper::Util::d2rad(alpha),
                                       clipper::Util::d2rad(beta),
                                       clipper::Util::d2rad(gamma));
        cell = clipper::Cell(cell_descr);
    }
    return cell;
}

} // namespace smcif
} // namespace coot