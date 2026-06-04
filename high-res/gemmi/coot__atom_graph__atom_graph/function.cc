#include "function.hh"
#include <iostream>

void coot::atom_graph::atom_graph_gemmi(gemmi::Structure* st,
                                        const std::vector<std::vector<coot::node_info>> &connection_indices,
                                        const std::vector<gemmi::Vec3> &coords,
                                        coot::atom_graph &ag) {
    if (connection_indices.size() != coords.size()) {
        std::cout << "ERROR:: atom_graph initialization error\n";
    } else {
        int ic = connection_indices.size();
        ag.atom_info.resize(ic);
        ag.nodes = connection_indices;

        // Create a dummy atom
        coot::minimol::atom dummy_atom;
        dummy_atom.name = " DUM";
        dummy_atom.element = " C";
        dummy_atom.altLoc = "";
        dummy_atom.pos = gemmi::Position(0.0, 0.0, 0.0);
        dummy_atom.occupancy = -1;
        dummy_atom.temperature_factor = 30.0;
        dummy_atom.int_user_data = -1;

        ag.atoms.resize(ic, dummy_atom);

        for (unsigned int iat = 0; iat < coords.size(); iat++) {
            coot::minimol::atom new_atom;
            new_atom.name = " DUM";
            new_atom.element = " C";
            new_atom.altLoc = "";
            new_atom.pos = gemmi::Position(coords[iat].x, coords[iat].y, coords[iat].z);
            new_atom.occupancy = -1;
            new_atom.temperature_factor = 30.0;
            new_atom.int_user_data = -1;
            ag.atoms[iat] = new_atom;
        }

        // assign cell and space group spgr
        ag.cell.resize(6);
        gemmi::UnitCell uc = st->cell;
        ag.cell[0] = uc.a;
        ag.cell[1] = uc.b;
        ag.cell[2] = uc.c;
        ag.cell[3] = uc.alpha;
        ag.cell[4] = uc.beta;
        ag.cell[5] = uc.gamma;

        ag.spgr = st->spacegroup_hm;
    }
}