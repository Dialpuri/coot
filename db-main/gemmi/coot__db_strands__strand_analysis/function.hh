#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__db_strands__orient_strand_on_z/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__db_strands__apply_rtop_to_strand/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__create_mmdbmanager_from_res_selection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__db_strands__trim_to_mainchain/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__molecule__molecule__a16eeb/gemmi/function.hh"

namespace coot {
namespace db_strands {

inline std::vector<coot::minimol::molecule>
strand_analysis_gemmi(gemmi::Structure& st, int strand_length) {

    std::vector<coot::minimol::molecule> rv;
    auto& model = st.models[0];

    std::cout << "has " << st.sheets.size() << " sheets" << std::endl;

    for (const auto& sheet : st.sheets) {
        int nstrands = static_cast<int>(sheet.strands.size());
        std::cout << "   Sheet has " << nstrands << " strands " << std::endl;

        for (const auto& strand : sheet.strands) {

            std::cout << "      strand " << strand.name
                      << " from " << strand.start.chain_name
                      << " " << strand.start.res_id.seqid.num.value
                      << " " << strand.start.res_id.seqid.icode
                      << " to "
                      << strand.end.chain_name
                      << " " << strand.end.res_id.seqid.num.value
                      << " " << strand.end.res_id.seqid.icode
                      << std::endl;

            int start_seq = strand.start.res_id.seqid.num.value;
            int end_seq   = strand.end.res_id.seqid.num.value;
            int lo = std::min(start_seq, end_seq);
            int hi = std::max(start_seq, end_seq);

            std::vector<gemmi::Residue*> residues;
            for (auto& chain : model.chains) {
                if (chain.name != strand.start.chain_name) continue;
                for (auto& res : chain.residues) {
                    int seq = res.seqid.num.value;
                    if (seq >= lo && seq <= hi) {
                        residues.push_back(&res);
                    }
                }
            }

            if (int(residues.size()) == strand_length) {

                // Orient strand on Z
                std::vector<const gemmi::Residue*> const_residues;
                for (auto* r : residues) const_residues.push_back(r);
                std::pair<bool, clipper::RTop_orth> ori =
                    coot::db_strands::orient_strand_on_z_gemmi(const_residues);

                if (ori.first) {

                    // Create a new structure from the residue selection
                    std::pair<gemmi::Structure, int> p =
                        coot::util::create_mmdbmanager_from_res_selection_gemmi(
                            st, residues,
                            0, 0, "", strand.start.chain_name, 0);

                    if (p.second) {

                        // Collect ALL atom indices in the new structure
                        std::vector<size_t> indices;
                        size_t idx = 0;
                        for (auto& ch : p.first.models[0].chains) {
                            for (auto& res : ch.residues) {
                                for (auto& atom : res.atoms) {
                                    indices.push_back(idx++);
                                }
                            }
                        }

                        // Apply the rotation/translation
                        coot::db_strands::apply_rtop_to_strand_gemmi(
                            p.first, indices, ori.second);

                        // Trim to mainchain
                        coot::db_strands::trim_to_mainchain_gemmi(p.first);

                        // Create molecule
                        coot::minimol::molecule m =
                            coot::minimol::molecule_gemmi(&p.first);
                        rv.push_back(m);
                    }
                }
            }
        }
    }

    return rv;
}

} // namespace db_strands
} // namespace coot