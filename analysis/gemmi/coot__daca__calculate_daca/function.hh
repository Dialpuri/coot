#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <clipper/core/coords.h>


// Ported helper forward declarations (defined in included headers, included in function.cc)
namespace coot {

clipper::Coord_orth co_gemmi(const gemmi::Atom* at);

bool is_main_chain_p_gemmi(const gemmi::Atom& atom, const gemmi::Residue& residue);

std::vector<std::vector<const gemmi::Atom*>> get_daca_fragments_gemmi(const gemmi::Residue* reference_residue_p);

std::pair<bool, clipper::RTop_orth> get_frag_to_reference_rtop_gemmi(
    const std::string& res_name,
    const unsigned int& frag_idx,
    const std::vector<const gemmi::Atom*>& fragment_atoms,
    const std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>>& reference_fragments);

bool atom_is_close_to_a_residue_atom_gemmi(const gemmi::Atom& at, const gemmi::Residue& residue);

bool atom_is_neighbour_mainchain_gemmi(const gemmi::CRA& atom_cra, const gemmi::CRA& ref_cra);

} // namespace coot

namespace coot {

namespace daca {

// DACA calculation mode
enum mode_t { REFERENCE, ANALYSIS };

// Quantizes a Coord_orth position into integer grid indices
struct box_index_t {
    int idx_x, idx_y, idx_z;

    box_index_t(const clipper::Coord_orth& pos) {
        idx_x = static_cast<int>(std::floor(pos.x()));
        idx_y = static_cast<int>(std::floor(pos.y()));
        idx_z = static_cast<int>(std::floor(pos.z()));
    }
};

// Key for the daca_data reference map
struct daca_key_t {
    std::string res_name;
    bool helical_flag;
    unsigned int frag_idx;
    box_index_t box_idx;
    std::string atom_type;

    bool operator<(const daca_key_t& o) const {
        if (res_name != o.res_name) return res_name < o.res_name;
        if (helical_flag != o.helical_flag) return helical_flag < o.helical_flag;
        if (frag_idx != o.frag_idx) return frag_idx < o.frag_idx;
        if (box_idx.idx_x != o.box_idx.idx_x) return box_idx.idx_x < o.box_idx.idx_x;
        if (box_idx.idx_y != o.box_idx.idx_y) return box_idx.idx_y < o.box_idx.idx_y;
        if (box_idx.idx_z != o.box_idx.idx_z) return box_idx.idx_z < o.box_idx.idx_z;
        return atom_type < o.atom_type;
    }
};

// Helpers (no port yet — defined in function.cc)
void add_to_box(
    coot::daca::mode_t mode,
    const std::string& res_name,
    bool helical_flag,
    unsigned int frag_idx,
    const box_index_t& box_index,
    const std::string& atom_type);

int get_reference_counts(
    const std::string& res_name,
    bool helical_flag,
    unsigned int frag_idx,
    const box_index_t& box_index,
    const std::string& atom_type);

// Print helper
inline std::string residue_spec_t_gemmi(const gemmi::CRA& cra) {
    std::ostringstream oss;
    if (!cra.residue) return "???";
    oss << "chain=" << cra.chain->name
        << " seq=" << cra.residue->seqid.num.value
        << " res=" << cra.residue->name;
    return oss.str();
}

// Print helper
inline std::string atom_spec_t_gemmi(const gemmi::Atom& atom) {
    return atom.name;
}

// Main calculation function
int calculate_daca_gemmi(
    const gemmi::CRA& reference_cra,
    const std::vector<std::pair<gemmi::CRA, std::string>>& typed_atoms,
    coot::daca::mode_t mode);

} // namespace daca

} // namespace coot