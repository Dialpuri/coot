#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

class Bond_lines_container {
public:
    float b_factor_scale = 1.0f;

    // Stores computed b-factor fractions for non-HETATM atoms
    // Key = "chain_name/res_seqnum/atom_name"
    std::map<std::string, float> b_factor_fractions_;

    int set_b_factor_colours_gemmi(gemmi::Structure& st);

    // Get the b-factor fraction for a given atom
    bool get_b_factor_fraction(const gemmi::Chain& chain,
                               const gemmi::Residue& res,
                               const gemmi::Atom& atom,
                               float& out_val) const;
};