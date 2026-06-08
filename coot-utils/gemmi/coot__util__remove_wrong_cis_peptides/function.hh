#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <vector>
#include <string>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__cis_peptides_info_from_coords/gemmi/function.hh"

namespace coot {
namespace util {

inline void remove_wrong_cis_peptides_gemmi(gemmi::Structure &st) {
    std::vector<coot::util::cis_peptide_info_t> v_coords =
        coot::util::cis_peptides_info_from_coords_gemmi(st);

    std::vector<gemmi::CisPep> good_cis_peptides;

    for (auto &cispep : st.cispeps) {
        std::string c1 = cispep.partner_c.chain_name;
        int r1 = cispep.partner_c.res_id.seqid.num.value;
        std::string c2 = cispep.partner_n.chain_name;
        int r2 = cispep.partner_n.res_id.seqid.num.value;

        bool ifound = false;
        for (auto &vc : v_coords) {
            if (vc.chain_id_1 == c1 &&
                vc.resno_1 == r1 &&
                vc.chain_id_2 == c2 &&
                vc.resno_2 == r2) {
                ifound = true;
                break;
            }
        }
        if (ifound) {
            good_cis_peptides.push_back(cispep);
        }
    }

    st.cispeps = std::move(good_cis_peptides);
}

} // namespace util
} // namespace coot
