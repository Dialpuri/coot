#pragma once

#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace coot {
namespace util {

std::vector<std::pair<coot::atom_spec_t, std::string>>
gln_asn_b_factor_outliers_gemmi(const gemmi::Model& model) {
    std::vector<std::pair<coot::atom_spec_t, std::string>> v;

    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            std::string res_name = residue.name;
            if (res_name != "ASN" && res_name != "GLN") continue;

            float b_sum = 0.0f;
            float b_sum_sq = 0.0f;
            const gemmi::Atom* oatom = nullptr;
            const gemmi::Atom* natom = nullptr;
            const gemmi::Atom* go_to_atom = nullptr;
            int n_residue_atoms = 0;

            for (const auto& atom : residue.atoms) {
                // Original: if (altloc == "") skip
                // gemmi: altloc is '\0' when no alternate location, also accept ' '
                if (atom.altloc != '\0' && atom.altloc != ' ') continue;

                std::string aname = atom.name;

                bool is_special = false;
                if (res_name == "GLN") {
                    if (aname == " OE1") { oatom = &atom; is_special = true; }
                    else if (aname == " NE2") { natom = &atom; is_special = true; }
                } else if (res_name == "ASN") {
                    if (aname == "OD1") { oatom = &atom; is_special = true; }
                    else if (aname == "ND2") { natom = &atom; is_special = true; }
                }

                if (!is_special) {
                    b_sum += atom.b_iso;
                    b_sum_sq += atom.b_iso * atom.b_iso;
                    n_residue_atoms++;
                }

                // gemmi trims trailing spaces: " CD " → " CD", " CG " → " CG"
                if (res_name == "GLN" && aname == " CD") go_to_atom = &atom;
                if (res_name == "ASN" && aname == " CG") go_to_atom = &atom;
            }

            if (oatom && natom && n_residue_atoms > 2) {
                float mean = b_sum / static_cast<float>(n_residue_atoms);
                float var = b_sum_sq / static_cast<float>(n_residue_atoms) - mean * mean;
                if (var < 0) var = 0;
                float std_dev = static_cast<float>(std::sqrt(static_cast<double>(var)));
                float diff = (oatom->b_iso - natom->b_iso) / 2.0f;

                if (diff > 0.0f && std_dev > 0.0f) {
                    float z = static_cast<float>(std::fabs(static_cast<double>(diff / std_dev)));

                    if (go_to_atom && z > 2.25f) {
                        // Normalize insertion code: gemmi ' ' → MMDB-style ""
                        std::string ins_code;
                        if (residue.seqid.icode != ' ')
                            ins_code = std::string(1, residue.seqid.icode);

                        coot::atom_spec_t as(
                            chain.name,
                            residue.seqid.num.value,
                            ins_code,
                            go_to_atom->name,
                            ""
                        );
                        as.float_user_data = z;

                        char buf[64];
                        std::string button_label = "Z score = ";
                        std::snprintf(buf, sizeof(buf), "%.2f", z);
                        button_label += buf;
                        button_label += "   ";
                        button_label += chain.name;
                        button_label += " ";
                        std::snprintf(buf, sizeof(buf), "%d", residue.seqid.num.value);
                        button_label += buf;
                        button_label += " ";
                        button_label += residue.name;

                        v.push_back({as, button_label});
                    }
                }
            }
        }
    }

    if (v.size() > 1) {
        std::sort(v.begin(), v.end(), coot::compare_atom_specs_user_float_in_pair);
    }

    return v;
}

} // namespace util
} // namespace coot