#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__cis_trans_conversion__ddbeb8/gemmi/function.hh"

int coot::molecule_t::cis_trans_conversion_gemmi(const std::string &atom_cid,
                                                 gemmi::Structure &st,
                                                 gemmi::Model *standard_residues_model) {

    // is_valid_model_molecule: check if the structure has models
    if (st.models.empty()) return 0;

    gemmi::Model *mol = &st.models[0];
    int status = 0;
    bool is_N_flag = false;

    gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi(atom_cid, st);
    std::string atom_name = at->name;

    if (atom_name == "N") is_N_flag = true;

    if (at) {
        // Build CRA from the atom — we need to find chain and residue
        gemmi::Residue *res = nullptr;
        gemmi::Chain *chain = nullptr;
        for (gemmi::Model &m : st.models) {
            for (gemmi::Chain &c : m.chains) {
                for (gemmi::Residue &r : c.residues) {
                    for (gemmi::Atom &a : r.atoms) {
                        if (&a == at) {
                            res = &r;
                            chain = &c;
                            break;
                        }
                    }
                    if (res) break;
                }
                if (res) break;
            }
            if (res) break;
        }

        gemmi::CRA cra{chain, res, at};
        status = coot::util::cis_trans_conversion_gemmi(cra, is_N_flag, mol, standard_residues_model);
    }

    return status;
}