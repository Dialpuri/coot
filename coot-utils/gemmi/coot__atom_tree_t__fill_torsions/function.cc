#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_name_map/gemmi/function.hh"

#include <algorithm>

// ── Inline copy of fill_atom_vertex_vec_gemmi (from existing port) ──
// Not auto-linked by build system, so we must include the body here.

namespace coot {

static bool fill_atom_vertex_vec_gemmi(
    const dictionary_residue_restraints_t& rest,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, map_index_t>& name_to_index,
    std::vector<atom_vertex>& atom_vertex_vec)
{
    bool retval = false;
    bool found_start = false;
    int rest_tree_start_index = -1;

    for (unsigned int i = 0; i < rest.tree.size(); i++) {
        if (rest.tree[i].connect_type == "START") {
            found_start = true;
            rest_tree_start_index = i;
            break;
        }
    }

    if (found_start) {
        int n_residue_atoms = static_cast<int>(res->atoms.size());
        atom_vertex_vec.resize(n_residue_atoms);

        for (unsigned int itree = 0; itree < rest.tree.size(); itree++) {
            map_index_t atom_id_index = name_to_index[rest.tree[itree].atom_id];
            if (atom_id_index.is_assigned()) {
                retval = true;
                int idx = atom_id_index.index();

                map_index_t atom_back_index = name_to_index[rest.tree[itree].atom_back];
                if (rest.tree[itree].atom_back != "") {
                    if (atom_back_index.is_assigned()) {
                        if (atom_back_index.index() < n_residue_atoms) {
                            atom_vertex_vec[idx].backward.push_back(atom_back_index.index());
                            int back_idx = atom_back_index.index();
                            const std::vector<int>& fwd = atom_vertex_vec[back_idx].forward;
                            bool already = (std::find(fwd.begin(), fwd.end(), idx) != fwd.end());
                            if (!already) {
                                atom_vertex_vec[back_idx].forward.push_back(idx);
                            }
                        }
                    }
                }

                map_index_t atom_forward_index = name_to_index[rest.tree[itree].atom_forward];
                if (atom_forward_index.is_assigned()) {
                    if (atom_forward_index.index() < n_residue_atoms) {
                        int fwd_idx = atom_forward_index.index();
                        const std::vector<int>& fwd = atom_vertex_vec[idx].forward;
                        bool already = (std::find(fwd.begin(), fwd.end(), fwd_idx) != fwd.end());
                        if (!already) {
                            atom_vertex_vec[idx].forward.push_back(fwd_idx);
                        }
                    }
                }

                atom_vertex_vec[idx].connection_type = atom_vertex::STANDARD;
                if (rest.tree[itree].connect_type == "START")
                    atom_vertex_vec[idx].connection_type = atom_vertex::START;
                if (rest.tree[itree].connect_type == "END")
                    atom_vertex_vec[idx].connection_type = atom_vertex::END;
            }
        }
    }

    return retval;
}

void construct_internal_gemmi(
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf,
    std::vector<std::pair<int, int>> &bonds_out,
    std::map<std::string, map_index_t> &name_to_index_out,
    std::vector<atom_vertex> &atom_vertex_vec_out)
{
    name_to_index_out = coot::fill_name_map_gemmi(*res, altconf);

    fill_atom_vertex_vec_gemmi(
        rest, res, altconf,
        name_to_index_out,
        atom_vertex_vec_out);

    // Fill bonds_out from rest.bond_restraint using name_to_index
    for (size_t ib = 0; ib < rest.bond_restraint.size(); ib++) {
        const dict_bond_restraint_t &br = rest.bond_restraint[ib];
        auto it1 = name_to_index_out.find(br.atom_id_1());
        auto it2 = name_to_index_out.find(br.atom_id_2());
        if (it1 != name_to_index_out.end() && it2 != name_to_index_out.end()) {
            bonds_out.push_back(
                std::make_pair(it1->second.index(), it2->second.index()));
        }
    }
}

} // namespace coot

// ── fill_torsions_gemmi ──

namespace coot {

bool fill_torsions_gemmi(
    std::vector<coot::atom_vertex> &atom_vertex_vec,
    const dictionary_residue_restraints_t &rest,
    gemmi::Residue *res,
    const std::string &altconf)
{
    bool r = false;
    int n_torsions_inserted = 0;

    if (rest.torsion_restraint.size() > 0) {
        std::vector<coot::atom_index_quad> quads;
        for (unsigned int itr = 0; itr < rest.torsion_restraint.size(); itr++) {
            coot::dict_torsion_restraint_t tr = rest.torsion_restraint[itr];
            std::pair<bool, coot::atom_index_quad> quad_pair =
                coot::get_atom_index_quad_gemmi(tr, *res, altconf);
            if (quad_pair.first) {
                quads.push_back(quad_pair.second);
            }
        }

        // Now we have a set of atom index quads, put them in the
        // atom_vertex_vec, at the position of the second atom in the
        // torsion (the position of the first atom in the rotation vector).

        for (unsigned int iquad = 0; iquad < quads.size(); iquad++) {
            bool inserted = false;
            int n_atom_vertex = static_cast<int>(atom_vertex_vec.size());
            for (int iv = 0; iv < n_atom_vertex; iv++) {
                if (iv == quads[iquad].index2) {
                    for (unsigned int ifo = 0; ifo < atom_vertex_vec[iv].forward.size(); ifo++) {
                        if (atom_vertex_vec[iv].forward[ifo] == quads[iquad].index3) {
                            // now check that the forward atom of this
                            // forward atom is index4
                            int this_forward = atom_vertex_vec[iv].forward[ifo];
                            for (unsigned int ifo2 = 0; ifo2 < atom_vertex_vec[this_forward].forward.size(); ifo2++) {
                                if (atom_vertex_vec[this_forward].forward[ifo2] == quads[iquad].index4) {
                                    atom_vertex_vec[iv].torsion_quad.first = true;
                                    atom_vertex_vec[iv].torsion_quad.second = quads[iquad];
                                    r = true;
                                    inserted = true;
                                    n_torsions_inserted++;
                                }
                            }
                        }
                    }
                }
                if (inserted)
                    break;
            }
        }
    }

    return r;
}

} // namespace coot
