#include "function.hh"

namespace lig_build {

// Inline the get_bond_type logic (private member of molfile_molecule_t)
static bond_t::bond_type_t get_bond_type_static(const std::string &bond_restraint_type) {
    bond_t::bond_type_t bond_type = bond_t::BOND_UNDEFINED;
    if (bond_restraint_type == "single")
        bond_type = bond_t::SINGLE_BOND;
    if (bond_restraint_type == "double")
        bond_type = bond_t::DOUBLE_BOND;
    if (bond_restraint_type == "triple")
        bond_type = bond_t::TRIPLE_BOND;
    if (bond_restraint_type == "aromatic")
        bond_type = bond_t::AROMATIC_BOND;
    return bond_type;
}

// Convert a gemmi atom name to the 4-char PDB format
// PDB format: prepend 1 space, then left-justify in 4 chars
// "N" -> " N  ", "CA" -> " CA ", "CB" -> " CB "
static std::string gemmi_name_to_pdb_4c(const std::string& name) {
    if (name.empty()) return "    ";
    std::string p = ' ' + name;
    while (p.size() < 4) p += ' ';
    return p;
}

molfile_molecule_t molfile_molecule_t_gemmi(
    gemmi::Residue* residue_p,
    const coot::dictionary_residue_restraints_t &restraints)
{
    molfile_molecule_t mol;

    molfile_atom_t blank_atom(clipper::Coord_orth(0, 0, 0), "", "");
    mol.atoms.push_back(blank_atom);

    // Build atoms from residue
    for (const gemmi::Atom& atom : residue_p->atoms) {
        clipper::Coord_orth pos(atom.pos.x, atom.pos.y, atom.pos.z);
        molfile_atom_t a(pos, atom.element.name(), atom.name);
        mol.atoms.push_back(a);
    }

    // Build atom map: padded PDB name -> list of (index, altloc)
    std::map<std::string, std::vector<std::pair<int, char>>> atom_map;
    int n_residue_atoms = static_cast<int>(residue_p->atoms.size());
    for (int iat = 0; iat < n_residue_atoms; iat++) {
        const gemmi::Atom& atom = residue_p->atoms[iat];
        std::string key = gemmi_name_to_pdb_4c(atom.name);
        atom_map[key].push_back({iat, atom.altloc});
    }

    // Build bonds from restraints
    for (unsigned int ibond = 0; ibond < restraints.bond_restraint.size(); ibond++) {
        const coot::dict_bond_restraint_t &bond_restraint = restraints.bond_restraint[ibond];
        std::string atom_name_1 = bond_restraint.atom_id_1_4c();
        std::string atom_name_2 = bond_restraint.atom_id_2_4c();

        auto it_1_atom_map = atom_map.find(atom_name_1);
        auto it_2_atom_map = atom_map.find(atom_name_2);

        if (it_1_atom_map != atom_map.end()) {
            if (it_2_atom_map != atom_map.end()) {
                const std::vector<std::pair<int, char>> &v_1 = it_1_atom_map->second;
                const std::vector<std::pair<int, char>> &v_2 = it_2_atom_map->second;

                for (unsigned int iat_1 = 0; iat_1 < v_1.size(); iat_1++) {
                    char alt_conf_1 = v_1[iat_1].second;
                    if (alt_conf_1 == '\0' || alt_conf_1 == ' ') {
                        for (unsigned int iat_2 = 0; iat_2 < v_2.size(); iat_2++) {
                            char alt_conf_2 = v_2[iat_2].second;
                            if (alt_conf_2 == '\0' || alt_conf_2 == ' ') {
                                int idx_1 = v_1[iat_1].first;
                                int idx_2 = v_2[iat_2].first;

                                bond_t::bond_type_t bond_type = get_bond_type_static(bond_restraint.type());

                                molfile_bond_t bond(idx_1 + 1, idx_2 + 1, bond_type);
                                mol.bonds.push_back(bond);
                            }
                        }
                    }
                }
            }
        }
    }

    return mol;
}

} // namespace lig_build