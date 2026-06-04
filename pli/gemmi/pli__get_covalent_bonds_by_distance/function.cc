#include "function.hh"

#include <gemmi/model.hpp>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Pad an atom name to 4 chars the PDB way: leading space, then name, then trailing spaces
// e.g. "N"  -> " N  "
//      "CA" -> " CA "
static inline std::string pad_atom_name(const std::string& name) {
    std::string padded = " " + name;
    while (padded.size() < 4)
        padded += ' ';
    return padded.substr(0, 4);
}

namespace pli {

std::vector<pli::fle_ligand_bond_t>
get_covalent_bonds_by_distance_gemmi(const gemmi::Model& model,
                                     const coot::residue_spec_t& ligand_spec) {
    std::vector<fle_ligand_bond_t> v;

    // Normalize insertion code: MMDB uses "" for no icode, gemmi uses ' '
    auto norm_icode = [](const std::string& ic) {
        return ic.empty() ? std::string(" ") : ic;
    };
    std::string norm_lig_icode = norm_icode(ligand_spec.ins_code);

    // Find ligand residue and collect its atoms
    const gemmi::Chain* lig_chain = nullptr;
    const gemmi::Residue* lig_residue = nullptr;
    std::vector<const gemmi::Atom*> lig_atoms;

    for (const auto& chain : model.chains) {
        if (chain.name != ligand_spec.chain_id)
            continue;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value != ligand_spec.res_no)
                continue;
            if (std::string(1, res.seqid.icode) != norm_lig_icode)
                continue;
            lig_chain = &chain;
            lig_residue = &res;
            for (const auto& atom : res.atoms) {
                lig_atoms.push_back(&atom);
            }
            break;
        }
        if (lig_residue) break;
    }

    if (!lig_residue || lig_atoms.empty())
        return v;

    // Build a set of ligand atom pointers for fast lookup
    std::set<const gemmi::Atom*> lig_atom_set(lig_atoms.begin(), lig_atoms.end());

    // Residue identity for dedup
    struct ResKey {
        std::string chain_id;
        int res_no;
        char icode;
        bool operator<(const ResKey& o) const {
            if (chain_id != o.chain_id) return chain_id < o.chain_id;
            if (res_no != o.res_no) return res_no < o.res_no;
            return icode < o.icode;
        }
    };
    ResKey lig_key{ligand_spec.chain_id, ligand_spec.res_no,
                   norm_lig_icode.empty() ? ' ' : norm_lig_icode[0]};

    std::set<std::pair<ResKey, ResKey>> contacting_pairs_set;

    double max_dist_sq = 2.3 * 2.3;
    double cno_max_dist_sq = 1.8 * 1.8;

    for (const gemmi::Atom* lig_atom : lig_atoms) {
        // Skip hydrogen atoms
        if (lig_atom->element.is_hydrogen())
            continue;

        std::string ele_1 = lig_atom->element.name();

        // Search all other atoms in the model
        for (const auto& chain : model.chains) {
            for (const auto& res : chain.residues) {
                // Skip same residue
                if (&res == lig_residue)
                    continue;

                for (const auto& other_atom : res.atoms) {
                    // Skip if other atom is in ligand
                    if (lig_atom_set.count(&other_atom))
                        continue;

                    // Alternative location check
                    if (lig_atom->altloc != '\0' && other_atom.altloc != '\0') {
                        if (lig_atom->altloc != other_atom.altloc)
                            continue;
                    }

                    // Skip hydrogen in other
                    if (other_atom.element.is_hydrogen())
                        continue;

                    // Compute distance squared
                    double d = (lig_atom->pos - other_atom.pos).length_sq();

                    // Element-based distance constraint
                    std::string ele_2 = other_atom.element.name();

                    double dist_for_bond_sq = max_dist_sq;
                    if (((ele_1 == "C") || (ele_1 == "N") || (ele_1 == "O")) &&
                        ((ele_2 == "C") || (ele_2 == "N") || (ele_2 == "O"))) {
                        dist_for_bond_sq = cno_max_dist_sq;
                    }

                    if (d < dist_for_bond_sq) {
                        ResKey other_key{chain.name, res.seqid.num.value, res.seqid.icode};

                        // Dedup by residue pair (ordered)
                        auto pair1 = std::make_pair(lig_key, other_key);
                        auto pair2 = std::make_pair(other_key, lig_key);
                        if (contacting_pairs_set.count(pair1) || contacting_pairs_set.count(pair2))
                            continue;
                        contacting_pairs_set.insert(pair1);

                        // Build atom specs
                        std::string lig_icode_str = (lig_key.icode == ' ' ? "" : std::string(1, lig_key.icode));
                        std::string other_icode_str = (other_key.icode == ' ' ? "" : std::string(1, other_key.icode));
                        std::string lig_alt = (lig_atom->altloc != '\0' ? std::string(1, lig_atom->altloc) : "");
                        std::string other_alt = (other_atom.altloc != '\0' ? std::string(1, other_atom.altloc) : "");

                        coot::atom_spec_t lig_atom_spec(
                            lig_key.chain_id, lig_key.res_no, lig_icode_str,
                            pad_atom_name(lig_atom->name), lig_alt);
                        coot::atom_spec_t other_atom_spec(
                            other_key.chain_id, other_key.res_no, other_icode_str,
                            pad_atom_name(other_atom.name), other_alt);

                        fle_ligand_bond_t bond(lig_atom_spec, other_atom_spec,
                                               fle_ligand_bond_t::BOND_COVALENT, d, false);
                        v.push_back(bond);
                    }
                }
            }
        }
    }

    return v;
}

} // namespace pli