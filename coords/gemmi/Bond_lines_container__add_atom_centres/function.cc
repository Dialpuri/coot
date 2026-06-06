#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {

void add_atom_centres_gemmi(
    std::vector<gemmi::Vec3>& atom_centres,
    std::vector<int>& atom_centres_colour,
    const gemmi::Structure& st,
    int atom_colour_type,
    int model_number,
    bool do_bonds_to_hydrogens,
    int n_selected_atoms_hint,
    std::map<std::string, int>* atom_colour_map_p) {

    atom_centres.clear();
    atom_centres_colour.clear();

    int idx = 0;
    for (size_t mi = 0; mi < st.models.size(); mi++) {
        const auto& model = st.models[mi];

        // Filter by model number if specified (model_number != 0 means specific model)
        // model_number is 1-based in coot (MMDB), models in gemmi are 0-based
        if (model_number != 0) {
            if (static_cast<int>(mi) != model_number - 1) {
                continue;
            }
        }

        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                for (const auto& atom : residue.atoms) {
                    bool is_H_flag = atom.element.is_hydrogen();

                    // Filter by hydrogen bonds setting
                    if (do_bonds_to_hydrogens || !is_H_flag) {
                        gemmi::Vec3 pos(atom.pos.x, atom.pos.y, atom.pos.z);
                        atom_centres.push_back(pos);

                        // Simplified colour: just use atom_colour_type
                        int colour = atom_colour_type;
                        atom_centres_colour.push_back(colour);
                    }

                    idx++;
                }
            }
        }
    }
}

}  // namespace coot