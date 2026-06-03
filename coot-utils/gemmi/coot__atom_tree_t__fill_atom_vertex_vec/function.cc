#include "function.hh"
#include "coot-utils/atom-tree.hh"
#include <algorithm>

// --- atom_tree_t_gemmi: constructor port (from existing verified port) ---

namespace coot {

void atom_tree_t_gemmi(
    const std::vector<std::vector<int> >& contact_indices,
    int base_atom_index,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, int>& name_to_index_out,
    std::vector<atom_vertex>& atom_vertex_vec_out)
{
    name_to_index_out.clear();
    atom_vertex_vec_out.clear();

    std::string altloc;
    if (!altconf.empty()) {
        altloc = altconf;
    }

    for (int i = 0; i < static_cast<int>(res->atoms.size()); i++) {
        gemmi::Atom& atom = res->atoms[i];
        if (!altconf.empty() && atom.altloc != 0 && atom.altloc != altconf[0]) {
            continue;
        }

        std::string atom_name = atom.name;
        // Trim spaces for comparison
        std::string trimmed = atom_name;
        size_t start = trimmed.find_first_not_of(" ");
        if (start != std::string::npos) {
            size_t end = trimmed.find_last_not_of(" ");
            trimmed = trimmed.substr(start, end - start + 1);
        } else {
            trimmed = "";
        }

        atom_vertex av;
        atom_vertex_vec_out.push_back(av);
        name_to_index_out[atom_name] = static_cast<int>(atom_vertex_vec_out.size()) - 1;
    }
}

// --- fill_atom_vertex_vec_gemmi ---

bool fill_atom_vertex_vec_gemmi(
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
                            // add_unique_forward_atom logic inline:
                            // Add idx to atom_vertex_vec[atom_back_index.index()].forward if not already there
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
                        // add_unique_forward_atom(idx, atom_forward_index.index())
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

} // namespace coot