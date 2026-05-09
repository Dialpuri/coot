#pragma once

#include <string>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>

namespace coot {

class molecule_t {
public:
    // Ported function that uses gemmi for residue lookup
    // Takes a gemmi::Structure as parameter since molecule_t doesn't inherently store one
    std::string get_svg_for_2d_ligand_environment_view_gemmi(
        const std::string &residue_cid,
        void *geom,  // coot::protein_geometry* - passed through to SVG generator
        bool add_key,
        const gemmi::Structure &st) const;
};

} // namespace coot