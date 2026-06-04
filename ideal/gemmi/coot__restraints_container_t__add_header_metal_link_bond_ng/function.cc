#include "function.hh"

namespace coot {

coot::atom_spec_t atom_spec_from_gemmi(const gemmi::Atom& atom,
                                       const gemmi::Residue& res,
                                       const gemmi::Chain& chain,
                                       int model_number)
{
    coot::atom_spec_t spec;
    spec.model_number = model_number;
    spec.chain_id = chain.name;
    spec.res_no = res.seqid.num.value;
    spec.ins_code = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
    spec.atom_name = atom.name;
    spec.alt_conf = (atom.altloc == '\0') ? "" : std::string(1, atom.altloc);
    spec.int_user_data = -1;
    spec.float_user_data = -1;
    return spec;
}

bool add_header_metal_link_bond_ng_gemmi(
    const std::vector<atom_entry>& atoms,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::vector<simple_restraint>& restraints_vec,
    const coot::atom_spec_t& atom_spec_1,
    const coot::atom_spec_t& atom_spec_2,
    double dist)
{
    int index_1 = -1;
    int index_2 = -1;

    for (int i = 0; i < static_cast<int>(atoms.size()); i++) {
        gemmi::Atom* at = atoms[i].atom;
        int rn_at = atoms[i].residue->seqid.num.value;
        if (rn_at == atom_spec_1.res_no) {
            coot::atom_spec_t s = atom_spec_from_gemmi(*at, *atoms[i].residue, *atoms[i].chain);
            if (s == atom_spec_1) {
                index_1 = i;
                continue;
            }
        }
        if (rn_at == atom_spec_2.res_no) {
            coot::atom_spec_t s = atom_spec_from_gemmi(*at, *atoms[i].residue, *atoms[i].chain);
            if (s == atom_spec_2) {
                index_2 = i;
            }
        }
        if ((index_1 != -1) && (index_2 != -1))
            break;
    }

    if (index_1 != -1) {
        if (index_2 != -1) {
            bonded_atom_indices[index_1].insert(index_2);
            bonded_atom_indices[index_2].insert(index_1);
            std::vector<bool> fixed_flags(4, false);
            fixed_flags[index_1] = true;
            fixed_flags[index_2] = true;
            restraints_vec.push_back(simple_restraint(
                BOND_RESTRAINT, index_1, index_2, fixed_flags, dist, 0.1, 1.2));
            return true;
        }
    }
    return false;
}

} // namespace coot