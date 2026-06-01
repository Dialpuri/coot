#include "function.hh"
#include <gemmi/model.hpp>
#include <algorithm>
#include <sstream>

namespace coot {

// Helper: trim whitespace from a string
static std::string trim_ws(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

std::vector<start_pos_match_result>
add_extra_start_pos_restraints_gemmi(
    gemmi::Model& model,
    const extra_restraints_t& extra_restraints)
{
    std::vector<start_pos_match_result> results;

    for (unsigned int i = 0; i < extra_restraints.start_pos_restraints.size(); i++) {
        const auto& spec = extra_restraints.start_pos_restraints[i].atom_1;

        // Find matching chain
        gemmi::Chain* matched_chain = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == spec.chain_id) {
                matched_chain = &chain;
                break;
            }
        }
        if (!matched_chain) continue;

        // Normalize insertion code: MMDB uses "" for "no icode", gemmi uses ' '
        std::string norm_ins = spec.ins_code.empty() ? std::string(" ") : spec.ins_code;

        // Find matching residue
        gemmi::Residue* matched_residue = nullptr;
        for (auto& res : matched_chain->residues) {
            if (res.seqid.num.value == spec.res_no) {
                // Compare insertion codes with normalization
                std::string res_icode(1, res.seqid.icode);
                if (norm_ins == res_icode) {
                    matched_residue = &res;
                    break;
                }
            }
        }
        if (!matched_residue) continue;

        // Find matching atom by name and alt_conf
        gemmi::Atom* matched_atom = nullptr;
        for (auto& atom : matched_residue->atoms) {
            std::string atom_name_trimmed = trim_ws(atom.name);
            if (atom_name_trimmed == spec.atom_name) {
                // Check alt_conf (altloc in gemmi)
                // MMDB uses "" for no altloc; gemmi uses '\0' for no altloc
                char alt = atom.altloc;
                std::string alt_str;
                if (alt != '\0') alt_str = std::string(1, alt);

                std::string spec_alt = spec.alt_conf.empty() ? std::string() : spec.alt_conf;

                if (alt_str == spec_alt) {
                    matched_atom = &atom;
                    break;
                }
            }
        }

        if (matched_atom) {
            start_pos_match_result result;
            result.cra = {matched_chain, matched_residue, matched_atom};
            result.restraint_index = static_cast<int>(i);
            result.esd = extra_restraints.start_pos_restraints[i].esd;
            results.push_back(result);
        }
    }

    return results;
}

} // namespace coot