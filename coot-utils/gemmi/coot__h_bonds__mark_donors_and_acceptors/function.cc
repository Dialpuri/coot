#include "function.hh"

#include <sstream>

namespace coot {
namespace h_bonds {

/* ------------------------------------------------------------------ */
/*  Helper: trim whitespace from gemmi atom names                      */
/* ------------------------------------------------------------------ */
static std::string trim_atom_name(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    size_t end   = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

/* ------------------------------------------------------------------ */
/*  select_atoms_by_path                                              */
/* ------------------------------------------------------------------ */

std::vector<gemmi::CRA>
select_atoms_by_path(gemmi::Model& model, const std::string& path) {
    std::vector<gemmi::CRA> result;

    // Expected path forms (from test): "//A/10" or "//A"
    // Strip leading '/' to get a clean token stream.
    std::string trimmed = path;
    while (!trimmed.empty() && trimmed.front() == '/') trimmed.erase(0, 1);

    std::stringstream ss(trimmed);
    std::vector<std::string> tokens;
    std::string tok;
    while (std::getline(ss, tok, '/')) {
        if (!tok.empty()) tokens.push_back(tok);
    }

    // tokens[0] = chain_id, optionally tokens[1] = residue_seqnum
    if (tokens.empty()) return result;

    std::string chain_id  = tokens[0];
    std::string seq_token = (tokens.size() > 1) ? tokens[1] : "";

    for (gemmi::Chain& chain : model.chains) {
        if (!chain_id.empty() && chain.name != chain_id) continue;

        for (gemmi::Residue& res : chain.residues) {
            if (!seq_token.empty()) {
                int target_seq = 0;
                try { target_seq = std::stoi(seq_token); } catch (...) { continue; }
                if (res.seqid.num.value != target_seq) continue;
            }

            for (gemmi::Atom& atom : res.atoms) {
                result.push_back(gemmi::CRA{&chain, &res, &atom});
            }
        }
    }

    return result;
}

/* ------------------------------------------------------------------ */
/*  mark_donors_and_acceptors_gemmi                                   */
/* ------------------------------------------------------------------ */

std::map<gemmi::Atom*, int>
mark_donors_and_acceptors_gemmi(
    const std::vector<gemmi::CRA>& sel_1_atoms,
    const std::vector<gemmi::CRA>& sel_2_atoms,
    const coot::protein_geometry& geom,
    int imol) {

    std::map<gemmi::Atom*, int> results;

    auto process = [&](const std::vector<gemmi::CRA>& cbras) {
        for (const gemmi::CRA& cra : cbras) {
            gemmi::Atom* atom  = cra.atom;
            gemmi::Residue* res = cra.residue;
            if (!atom || !res) continue;

            // gemmi atom names come padded (e.g. " N  "). Trim, then expand
            // to the 4-char dictionary format via coot::atom_id_mmdb_expand.
            std::string raw_name = trim_atom_name(atom->name);
            std::string name     = coot::atom_id_mmdb_expand(raw_name);
            std::string res_name = res->name;

            int h_bond_type = geom.get_h_bond_type(name, res_name, imol);
            results[atom] = h_bond_type;
        }
    };

    process(sel_1_atoms);

    // Process sel_2 only if it's a different vector (not the same pointer)
    if (&sel_1_atoms != &sel_2_atoms) {
        process(sel_2_atoms);
    }

    return results;
}

} // namespace h_bonds
} // namespace coot