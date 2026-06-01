#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "utils/coot-utils.hh"

// gemmi-compatible peptide flip atom triplet
struct flip_atom_triplet_t_gemmi {
    gemmi::Atom* CA_this;
    gemmi::Atom* O_this;
    gemmi::Atom* CA_next;
    gemmi::Residue* res_this;
    gemmi::Residue* res_next;
    std::string bond_type;

    flip_atom_triplet_t_gemmi(gemmi::Atom* ca, gemmi::Atom* o,
                              gemmi::Atom* ca_next,
                              gemmi::Residue* r1, gemmi::Residue* r2)
        : CA_this(ca), O_this(o), CA_next(ca_next),
          res_this(r1), res_next(r2), bond_type("PEPTIDE") {}
};

namespace coot {
namespace pepflip_using_difference_map {

std::vector<flip_atom_triplet_t_gemmi>
get_peptide_atom_triplets_gemmi(const gemmi::Structure& st) {
    std::vector<flip_atom_triplet_t_gemmi> rps;

    if (st.models.empty()) return rps;

    const gemmi::Model& model = st.models[0];

    for (const gemmi::Chain& chain : model.chains) {
        if (chain.residues.size() < 2) continue;
        for (size_t ires = 0; ires + 1 < chain.residues.size(); ++ires) {
            const gemmi::Residue& res_this = chain.residues[ires];
            const gemmi::Residue& res_next = chain.residues[ires + 1];

            std::string res_name_this = res_this.name;
            std::string res_name_next = res_next.name;
            int res_no_this = res_this.seqid.num.value;
            int res_no_next = res_next.seqid.num.value;

            if (res_no_next != (res_no_this + 1)) continue;
            if (!coot::util::is_standard_amino_acid_name(res_name_this)) continue;
            if (!coot::util::is_standard_amino_acid_name(res_name_next)) continue;

            const gemmi::Atom* O_this = nullptr;
            const gemmi::Atom* CA_this = nullptr;
            const gemmi::Atom* CA_next = nullptr;

            // gemmi atom names are unpadded (e.g. "O", "CA") not padded (" O  ", " CA ")
            for (const gemmi::Atom& atom : res_this.atoms) {
                if (atom.name == "O" && atom.altloc == '\0') {
                    O_this = &atom;
                }
                if (atom.name == "CA" && atom.altloc == '\0') {
                    CA_this = &atom;
                }
            }

            if (CA_this && O_this) {
                for (const gemmi::Atom& atom : res_next.atoms) {
                    if (atom.name == "CA" && atom.altloc == '\0') {
                        CA_next = &atom;
                        break;
                    }
                }

                if (CA_next) {
                    rps.push_back(flip_atom_triplet_t_gemmi(
                        const_cast<gemmi::Atom*>(CA_this),
                        const_cast<gemmi::Atom*>(O_this),
                        const_cast<gemmi::Atom*>(CA_next),
                        const_cast<gemmi::Residue*>(&res_this),
                        const_cast<gemmi::Residue*>(&res_next)));
                }
            }
        }
    }

    return rps;
}

}  // namespace pepflip_using_difference_map
}  // namespace coot