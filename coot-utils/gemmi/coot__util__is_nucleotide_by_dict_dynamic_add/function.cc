#include "function.hh"

namespace coot {
namespace util {

bool is_nucleotide_by_dict_dynamic_add_gemmi(
    gemmi::Residue* residue_p,
    coot::protein_geometry* geom_p
) {
    bool is_nuc = false;
    std::string residue_name = residue_p->name;

    std::pair<bool, dictionary_residue_restraints_t> p =
        geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);

    if (p.first) {
        if (p.second.residue_info.group == "RNA" ||
            p.second.residue_info.group == "DNA") {
            is_nuc = true;
        }
    } else {
        int read_number = 40;
        int status = geom_p->try_dynamic_add(residue_name, read_number);
        if (status != 0) {
            p = geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
            if (p.first) {
                if (p.second.residue_info.group == "RNA" ||
                    p.second.residue_info.group == "DNA") {
                    is_nuc = true;
                }
            }
        }
    }
    return is_nuc;
}

} // namespace util
} // namespace coot