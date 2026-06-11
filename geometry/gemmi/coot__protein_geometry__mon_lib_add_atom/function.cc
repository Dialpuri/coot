#include "function.hh"

namespace coot {

void mon_lib_add_atom_gemmi(
    std::vector<std::pair<int, dictionary_residue_restraints_t>>& dict_res_restraints,
    int read_number,
    const std::string &comp_id,
    int imol_enc,
    const std::string &atom_id,
    const std::string &atom_id_4c,
    const std::string &type_symbol,
    const std::string &type_energy,
    const std::pair<bool, double> &partial_charge,
    const std::pair<bool, int> &formal_charge,
    dict_atom::aromaticity_t arom_in,
    const std::pair<bool, clipper::Coord_orth> &model_pos,
    const std::pair<bool, clipper::Coord_orth> &model_pos_ideal)
{
    bool debug = false;

    if (debug) {
        std::cout << "   mon_lib_add_atom  " << comp_id << " atom-id:" << atom_id << ": ::"
                  << atom_id_4c << ": " << type_symbol << " " << type_energy << " ( "
                  << partial_charge.first << "," << partial_charge.second << ")";
        std::cout << " model-pos: " << model_pos.first << " ";
        if (model_pos.first)
            std::cout << "( "
                      << model_pos.second.x() << " "
                      << model_pos.second.y() << " "
                      << model_pos.second.z() << " ) ";
        std::cout << "model-pos-ideal: " << model_pos_ideal.first << " ";
        if (model_pos_ideal.first)
            std::cout << "( "
                      << model_pos_ideal.second.x() << " "
                      << model_pos_ideal.second.y() << " "
                      << model_pos_ideal.second.z() << " ) ";
        std::cout << std::endl;
    }

    dict_atom at_info(atom_id, atom_id_4c, type_symbol, type_energy, partial_charge);
    at_info.aromaticity = arom_in;
    at_info.formal_charge = formal_charge;

    if (debug) {
        std::cout << "   mon_lib_add_atom model_pos       " << model_pos.first << " "
                  << model_pos.second.format() << std::endl;
        std::cout << "   mon_lib_add_atom model_pos_ideal " << model_pos_ideal.first << " "
                  << model_pos_ideal.second.format() << std::endl;
    }

    if (model_pos.first)
        at_info.add_pos(dict_atom::REAL_MODEL_POS, model_pos);
    if (model_pos_ideal.first)
        at_info.add_pos(dict_atom::IDEAL_MODEL_POS, model_pos_ideal);

    bool ifound = false;
    int this_index = -1; // unset

    for (unsigned int i = 0; i < dict_res_restraints.size(); i++) {
        if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
            if (dict_res_restraints[i].first == imol_enc) {
                if (dict_res_restraints[i].second.read_number == read_number) {
                    ifound = true;
                    this_index = static_cast<int>(i);
                    dict_res_restraints[i].second.atom_info.push_back(at_info);
                    break;
                } else {
                    std::cout << "######## trash the old one " << comp_id << std::endl;
                    dict_res_restraints[i].second.clear_dictionary_residue();
                }
            }
        }
    }

    if (!ifound) {
        dictionary_residue_restraints_t rest(comp_id, read_number);
        std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
        dict_res_restraints.push_back(p);
        this_index = static_cast<int>(dict_res_restraints.size()) - 1;
        dict_res_restraints[this_index].second.atom_info.push_back(at_info);
    }

    if (debug) {
        std::cout << "   dictionary for " << dict_res_restraints[this_index].second.residue_info.comp_id
                  << " now contains " << dict_res_restraints[this_index].second.atom_info.size()
                  << " atoms" << std::endl;
        for (unsigned int i = 0; i < dict_res_restraints[this_index].second.atom_info.size(); i++) {
        }
    }
}

} // namespace coot
