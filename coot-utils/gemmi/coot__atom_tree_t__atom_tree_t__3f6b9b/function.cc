#include "function.hh"
#include <stdexcept>

namespace coot {

void atom_tree_t_gemmi(
    const std::vector<std::vector<int>>& contact_indices,
    int base_atom_index,
    gemmi::Residue* res,
    const std::string& altconf,
    std::map<std::string, int>& name_to_index,
    std::vector<atom_vertex>& atom_vertex_vec) {

    if (!res) {
        throw std::runtime_error("null residue in alternate atom_tree_t constructor");
    }

    // Clear outputs
    name_to_index.clear();
    atom_vertex_vec.clear();

    // fill_name_map — inline the logic:
    // iterate over residue atoms, map (atom_name + altconf) -> index
    for (size_t i = 0; i < res->atoms.size(); i++) {
        std::string atom_name = res->atoms[i].name;
        std::string key = atom_name + altconf;
        name_to_index[key] = static_cast<int>(i);
    }

    // fill_atom_vertex_vec_using_contacts — use the verified gemmi port
    coot::fill_atom_vertex_vec_using_contacts_gemmi(
        contact_indices, base_atom_index, res->atoms, atom_vertex_vec
    );
}

} // namespace coot