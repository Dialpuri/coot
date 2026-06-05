#pragma once

#include <gemmi/model.hpp>
#include <string>
#include "coords/Bond_lines.hh"
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__atom_colour/gemmi/function.hh"

// Helper class that inherits from Bond_lines_container to expose protected/private members
class Bond_lines_container_gemmi : public Bond_lines_container {
public:
    using Bond_lines_container::bonds;

    // Expose protected addBond
    void addBond_exposed(int colour, const coot::Cartesian &first, const coot::Cartesian &second,
                         graphics_line_t::cylinder_class_t cc, int model_number,
                         int udd_atom_index_handle, int atom_index_2,
                         bool is_dashed, bool is_highlight)
    {
        Bond_lines_container::addBond(colour, first, second, cc, model_number,
                                      udd_atom_index_handle, atom_index_2, is_dashed, is_highlight);
    }

    // Count total bonds across all bond groups
    int total_bond_count() const {
        int total = 0;
        for (const auto& bl : bonds) total += bl.size();
        return total;
    }
};

namespace coot {

// Struct holding gemmi-compatible link bond information (replaces mmdb::Link)
struct gemmi_link_info_t {
    std::string chain1;
    int seq_num1;
    std::string res_name1;
    std::string atom_name1;
    char icode1 = ' ';

    std::string chain2;
    int seq_num2;
    std::string res_name2;
    std::string atom_name2;
    char icode2 = ' ';
};

// Find atom in gemmi::Model by identifiers, return gemmi::CRA
inline gemmi::CRA find_atom_cra_in_model(
    gemmi::Model& model,
    const std::string& chain_name,
    int seq_num,
    const std::string& res_name,
    const std::string& atom_name,
    char icode = ' ')
{
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_name) continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value != seq_num) continue;
            if (res.name != res_name) continue;
            if (res.seqid.icode != icode) continue;
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.altloc != ' ') continue;
                std::string trimmed = atom.name;
                trimmed.erase(trimmed.find_last_not_of(' ') + 1);
                if (trimmed == atom_name) {
                    return gemmi::CRA{&chain, &res, &atom};
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

// gemmi port of Bond_lines_container::add_link_bond
// Finds two atoms by identifiers in the gemmi model and adds a dashed bond
// between them (mimicking add_link_bond_templ logic from the original).
inline void add_link_bond_gemmi(
    Bond_lines_container_gemmi& bonds,
    gemmi::Model& model,
    int udd_atom_index_handle,
    int udd_user_defined_atom_colour_index_handle,
    int atom_colour_type,
    const gemmi_link_info_t& link)
{
    gemmi::CRA cra_1 = find_atom_cra_in_model(
        model, link.chain1, link.seq_num1, link.res_name1, link.atom_name1, link.icode1);
    gemmi::CRA cra_2 = find_atom_cra_in_model(
        model, link.chain2, link.seq_num2, link.res_name2, link.atom_name2, link.icode2);

    if (!cra_1.atom || !cra_2.atom) {
        return;
    }

    // Atom indices (use -1 when UDD not available, matching original behavior)
    int atom_index_1 = -1;
    int atom_index_2 = -1;

    coot::Cartesian pos_1(cra_1.atom->pos.x, cra_1.atom->pos.y, cra_1.atom->pos.z);
    coot::Cartesian pos_2(cra_2.atom->pos.x, cra_2.atom->pos.y, cra_2.atom->pos.z);

    std::string ele_1 = cra_1.atom->element.name();
    std::string ele_2 = cra_2.atom->element.name();

    int model_number = 0;  // default model number

    if (ele_1 == ele_2) {
        // Same element: single dashed bond with one colour
        int col = atom_colour_gemmi(cra_1, atom_colour_type,
                                     udd_user_defined_atom_colour_index_handle, nullptr);
        bonds.addBond_exposed(col, pos_1, pos_2, graphics_line_t::SINGLE,
                               model_number, atom_index_1, atom_index_2, false, true);
    } else {
        // Different elements: two half bonds with different colours
        coot::Cartesian bond_mid_point = pos_1.mid_point(pos_2);

        int col1 = atom_colour_gemmi(cra_1, atom_colour_type,
                                      udd_user_defined_atom_colour_index_handle, nullptr);
        bonds.addBond_exposed(col1, pos_1, bond_mid_point, graphics_line_t::SINGLE,
                               model_number, atom_index_1, atom_index_2, false, true);

        int col2 = atom_colour_gemmi(cra_2, atom_colour_type,
                                      udd_user_defined_atom_colour_index_handle, nullptr);
        bonds.addBond_exposed(col2, bond_mid_point, pos_2, graphics_line_t::SINGLE,
                               model_number, atom_index_1, atom_index_2, false, true);
    }
}

} // namespace coot