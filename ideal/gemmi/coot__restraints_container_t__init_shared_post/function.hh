#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <vector>
#include <set>
#include <string>
#include <utility>
#include <iostream>

namespace coot {

struct atom_spec_t {
    std::string chain_id;
    std::string atom_name;
    char alt_conf;
    int res_no;
    std::string ins_code;

    atom_spec_t()
        : chain_id(), atom_name(), alt_conf('\0'), res_no(0), ins_code() {}

    atom_spec_t(const gemmi::Atom& atom, const gemmi::Residue& res, const std::string& chain_id_str) {
        chain_id = chain_id_str;
        atom_name = atom.name;
        alt_conf = atom.altloc;
        res_no = res.seqid.num.value;
        ins_code = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
    }

    bool matches_spec(const gemmi::Atom* atom_p,
                      const gemmi::Residue* res_p,
                      const std::string& chain_id_str) const {
        if (!atom_p || !res_p) return false;
        std::string atom_name_str(atom_p->name);
        if (atom_name != atom_name_str) return false;
        if (alt_conf != '\0' && alt_conf != atom_p->altloc) return false;
        if (res_no != res_p->seqid.num.value) return false;
        std::string query_ic = ins_code;
        if (query_ic.empty()) query_ic = " ";
        std::string res_ic(1, res_p->seqid.icode);
        if (query_ic != res_ic) return false;
        if (!chain_id.empty() && chain_id != chain_id_str) return false;
        return true;
    }
};

class restraints_container_t {
public:
    struct residue_identifier {
        int seqnum;
        char icode;
        bool operator<(const residue_identifier& other) const {
            if (seqnum != other.seqnum) return seqnum < other.seqnum;
            return icode < other.icode;
        }
    };

    std::vector<gemmi::Atom*> atom;
    std::vector<gemmi::Residue*> residue;
    std::vector<std::string> chain_ids;
    size_t n_atoms;

    std::vector<double> initial_position_params_vec;
    std::vector<std::vector<int>> bonded_atom_indices;
    std::vector<bool> use_map_gradient_for_atom;
    std::set<int> fixed_atom_indices;

    bool from_residue_vector;
    bool do_hydrogen_atom_refinement;
    bool model_has_hydrogen_atoms;
    bool cryo_em_mode;
    bool do_neutron_refinement;
    int istart_res;
    int iend_res;
    int udd_bond_angle;
    int udd_atom_index_handle;

    std::set<restraints_container_t::residue_identifier> residues_vec_moving_set;

    std::vector<float> atom_z_occ_weight;

    restraints_container_t()
        : n_atoms(0), from_residue_vector(false),
          do_hydrogen_atom_refinement(false),
          model_has_hydrogen_atoms(false), cryo_em_mode(false),
          do_neutron_refinement(false), istart_res(0), iend_res(0),
          udd_bond_angle(-1), udd_atom_index_handle(-1) {}

    bool is_hydrogen(const gemmi::Atom* a) const {
        if (!a) return false;
        return a->element.is_hydrogen();
    }

    bool is_a_moving_residue_p_gemmi(const gemmi::Residue* r) const {
        if (!r) return false;
        residue_identifier id;
        id.seqnum = r->seqid.num.value;
        id.icode = r->seqid.icode;
        return residues_vec_moving_set.find(id) != residues_vec_moving_set.end();
    }

    void set_has_hydrogen_atoms_state_gemmi() {
        model_has_hydrogen_atoms = false;
        for (size_t i = 0; i < n_atoms; i++) {
            if (atom[i] && atom[i]->element.is_hydrogen()) {
                model_has_hydrogen_atoms = true;
                break;
            }
        }
    }

    void set_z_occ_weights_gemmi();

    void init_shared_post_gemmi(const std::vector<coot::atom_spec_t>& fixed_atom_specs);
};

} // namespace coot
