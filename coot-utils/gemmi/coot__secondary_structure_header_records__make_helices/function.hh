#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <string>
#include <vector>

namespace coot {
namespace secondary_structure_header_records {

class helix_info_t {
public:
    gemmi::CRA start_res;
    gemmi::CRA end_res;
    unsigned int length;
    helix_info_t(gemmi::CRA r1, gemmi::CRA r2, unsigned int l)
        : start_res(r1), end_res(r2), length(l) {}
};

inline std::string helix_index_to_helix_id(int idx) {
    return std::to_string(idx);
}

inline void make_helices_gemmi(gemmi::Structure& structure,
                               const std::vector<helix_info_t>& helices_in) {
    for (std::size_t i = 0; i < helices_in.size(); i++) {
        gemmi::Helix h;

        const gemmi::CRA& s_cra = helices_in[i].start_res;
        const gemmi::CRA& e_cra = helices_in[i].end_res;

        // Start residue — fill AtomAddress
        h.start.chain_name = s_cra.chain->name;
        h.start.res_id     = *s_cra.residue;

        // End residue
        h.end.chain_name   = e_cra.chain->name;
        h.end.res_id       = *e_cra.residue;

        // Helix class (MMDB defaults to 1, which maps to Halpha)
        h.set_helix_class_as_int(1);

        h.length = static_cast<int>(helices_in[i].length);

        structure.helices.push_back(h);
    }
}

} // namespace secondary_structure_header_records
} // namespace coot