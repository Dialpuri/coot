#pragma once

#include <string>
#include <vector>
#include <utility>

namespace coot { namespace protein_geometry { constexpr int IMOL_ENC_AUTO = 0; } }

namespace coot {

// --- Minimal definitions mirroring coot types needed by mon_lib_add_bond ---

struct dict_bond_restraint_t {
    enum aromaticity_t { NON_AROMATIC, AROMATIC };
    enum bond_length_type_t { UNKNOWN, NUCLEAR_POSITION, NUCLEAR_POSITION_PLUS };

    std::string atom_id_1;
    std::string atom_id_2;
    std::string type;
    double value_dist;
    double value_dist_esd;
    double value_dist_nuclear;
    double value_dist_nuclear_esd;
    bool value_dist_nuclear_was_set;
    aromaticity_t arom_in;
    bond_length_type_t type_in;

    dict_bond_restraint_t(
        const std::string& a1,
        const std::string& a2,
        const std::string& t,
        double vd,
        double vde,
        double vdn,
        double vdne,
        bool vdnws,
        aromaticity_t ar,
        bond_length_type_t blt)
        : atom_id_1(a1), atom_id_2(a2), type(t),
          value_dist(vd), value_dist_esd(vde),
          value_dist_nuclear(vdn), value_dist_nuclear_esd(vdne),
          value_dist_nuclear_was_set(vdnws), arom_in(ar), type_in(blt) {}
};

struct dict_chem_comp_t {
    std::string comp_id;
};

struct dictionary_residue_restraints_t {
    dict_chem_comp_t residue_info;
    std::vector<dict_bond_restraint_t> bond_restraint;

    dictionary_residue_restraints_t(const std::string& comp_id, int /*read_number*/)
        : residue_info(), bond_restraint()
    {
        residue_info.comp_id = comp_id;
    }
};

// --- gemmi port of mon_lib_add_bond ---

inline void mon_lib_add_bond_gemmi(
    const std::string& comp_id,
    int imol_enc,
    const std::string& atom_id_1,
    const std::string& atom_id_2,
    const std::string& type,
    double value_dist,
    double value_dist_esd,
    double value_dist_nuclear,
    double value_dist_nuclear_esd,
    dict_bond_restraint_t::aromaticity_t arom_in,
    dict_bond_restraint_t::bond_length_type_t type_in,
    std::vector<std::pair<int, dictionary_residue_restraints_t>>& dict_res_restraints,
    int& read_number)
{
    bool value_dist_nuclear_was_set = false;
    if (value_dist_nuclear_esd > 0.0)
        value_dist_nuclear_was_set = true;

    bool ifound = false;

    for (unsigned int i = 0; i < dict_res_restraints.size(); i++) {
        if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
            if (dict_res_restraints[i].first == imol_enc) {
                ifound = true;
                dict_bond_restraint_t bond_restr(
                    atom_id_1, atom_id_2, type,
                    value_dist, value_dist_esd,
                    value_dist_nuclear, value_dist_nuclear_esd,
                    value_dist_nuclear_was_set, arom_in, type_in);
                dict_res_restraints[i].second.bond_restraint.push_back(bond_restr);
            }
        }
    }

    if (!ifound) {
        dictionary_residue_restraints_t rest(comp_id, read_number);
        std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
        dict_res_restraints.push_back(p);
        dict_bond_restraint_t bond_restr(
            atom_id_1, atom_id_2, type,
            value_dist, value_dist_esd,
            value_dist_nuclear, value_dist_nuclear_esd,
            value_dist_nuclear_was_set, arom_in, type_in);
        dict_res_restraints[dict_res_restraints.size() - 1].second.bond_restraint.push_back(bond_restr);
    }
}

} // namespace coot
