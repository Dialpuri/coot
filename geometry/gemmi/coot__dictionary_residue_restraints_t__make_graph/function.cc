#include "function.hh"

namespace coot {

static bool is_hydrogen_elem(const std::string &ele) {
    return (ele == " H" || ele == "H");
}

static bool is_hydrogen_by_name(const std::string &atom_name_4c, const std::vector<dict_atom> &atom_info) {
    for (const auto &ai : atom_info) {
        if (ai.atom_id_4c == atom_name_4c) {
            return is_hydrogen_elem(ai.type_symbol);
        }
    }
    return false;
}

std::unique_ptr<dict_graph_t> make_graph_gemmi(const dictionary_residue_restraints_t &drr, bool use_hydrogens) {
    std::map<std::string, unsigned int> name_map;

    dict_graph_t *graph = new dict_graph_t;
    int i_atom = 0;

    for (unsigned int iat = 0; iat < drr.atom_info.size(); iat++) {
        std::string ele = drr.atom_info[iat].type_symbol;
        if (use_hydrogens || !is_hydrogen_elem(ele)) {
            std::string name = drr.atom_info[iat].atom_id_4c;
            dict_graph_vertex_t v;
            v.element = ele;
            v.name    = name;
            graph->vertices.push_back(v);
            name_map[name] = i_atom;
            i_atom++;
        }
    }

    for (unsigned int ib = 0; ib < drr.bond_restraint.size(); ib++) {
        const dict_bond_restraint_t &br = drr.bond_restraint[ib];
        int mmdb_bond_type = br.mmdb_bond_type();
        if (mmdb_bond_type != -1) {
            std::map<std::string, unsigned int>::const_iterator it_1;
            std::map<std::string, unsigned int>::const_iterator it_2;
            it_1 = name_map.find(br.atom_id_1_4c());
            it_2 = name_map.find(br.atom_id_2_4c());

            if (it_1 == name_map.end()) {
                // would warn in original
            } else {
                if (it_2 == name_map.end()) {
                    // would warn in original
                } else {
                    if (use_hydrogens || (!is_hydrogen_by_name(br.atom_id_1_4c(), drr.atom_info) &&
                                          !is_hydrogen_by_name(br.atom_id_2_4c(), drr.atom_info))) {
                        dict_graph_edge_t e;
                        e.v1        = static_cast<unsigned int>(it_1->second + 1);  // 1-indexed
                        e.v2        = static_cast<unsigned int>(it_2->second + 1);
                        e.bond_type = mmdb_bond_type;
                        graph->edges.push_back(e);
                    }
                }
            }
        }
    }

    return std::unique_ptr<dict_graph_t>(graph);
}

} // namespace coot