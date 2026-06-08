#pragma once

#include <clipper/core/xmap.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/cell.h>
#include <clipper/core/coords.h>
#include <clipper/contrib/edcalc.h>
#include <gemmi/model.hpp>
#include <cmath>
#include <iostream>
#include <vector>

namespace coot::util {

clipper::Xmap<float> calc_atom_map_gemmi(
    std::vector<gemmi::Atom const*> atoms,
    const clipper::Cell& cell,
    const clipper::Spacegroup& space_group,
    const clipper::Grid_sampling& sampling) {

    clipper::Xmap<float> xmap;
    xmap.init(space_group, cell, sampling);

    std::vector<clipper::Atom> l;
    float rescale_b_u = 1.0f / (8.0f * static_cast<float>(M_PI) * static_cast<float>(M_PI));

    for (gemmi::Atom const* at : atoms) {
        // gemmi has no terminator atoms (isTer equivalent), so process all
        clipper::Coord_orth pt(at->pos.x, at->pos.y, at->pos.z);
        std::string ele(at->element.name());
        clipper::Atom cat;
        cat.set_element(ele);
        cat.set_coord_orth(pt);
        float u_iso = at->b_iso * rescale_b_u;
        if (u_iso < 0.1f) u_iso = 0.1f;  // B < ~8.0: EDcalc_iso produces NaN from 0/0
        cat.set_u_iso(u_iso);
        cat.set_occupancy(at->occ);
        l.push_back(cat);
    }

    try {
        clipper::Atom_list al(l);
        clipper::EDcalc_iso<float> e;
        e(xmap, al);

        // EDcalc_iso can produce NaN for certain atom configurations.
        // Clean them to avoid corrupting downstream FFTs.
        unsigned int nan_count = 0;
        clipper::Xmap<float>::Map_reference_index inx;
        for (inx = xmap.first(); !inx.last(); inx.next()) {
            if (std::isnan(xmap[inx])) {
                clipper::Coord_grid grid_coords = inx.coord();
                clipper::Coord_orth orth_coords = inx.coord_orth();
                std::cout << "WARNING:: NaN value at " << grid_coords.format()
                          << " coord-orth " << orth_coords.x() << " " << orth_coords.y() << " " << orth_coords.z()
                          << std::endl;
            }
            if (std::isinf(xmap[inx])) {
                clipper::Coord_grid grid_coords = inx.coord();
                std::cout << "WARNING:: Inf value at " << grid_coords.format() << std::endl;
            }
            if (std::isnan(xmap[inx]) || std::isinf(xmap[inx])) {
                xmap[inx] = 0.0f;
                nan_count++;
            }
        }
        if (nan_count > 0)
            std::cout << "WARNING:: calc_atom_map() cleaned " << nan_count
                      << " NaN/inf values from map" << std::endl;
    }
    catch (const clipper::Message_generic &e) {
        std::cout << "ERROR:: some sort of clipper map problem" << std::endl;
        std::cout << e.text() << std::endl;
    }
    return xmap;
}

} // namespace coot::util
