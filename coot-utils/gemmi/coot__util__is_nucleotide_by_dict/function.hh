#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

bool is_nucleotide_by_dict_gemmi(const gemmi::Residue& res, const std::string& monomer_name);

} // namespace util
} // namespace coot

inline bool
coot::util::is_nucleotide_by_dict_gemmi(const gemmi::Residue& res, const std::string& monomer_name) {
    struct dummy_geometry {
        struct residue_info_t { std::string group; };
        struct restraints_t { residue_info_t residue_info; };
        std::pair<int, restraints_t> get_monomer_restraints(const std::string& name, int) const {
            restraints_t r;
            if (name == "A" || name == "DA" || name == "G" || name == "DG" ||
                name == "C" || name == "DC" || name == "U" || name == "DT") {
                r.residue_info.group = (name[0] == 'D') ? "DNA" : "RNA";
            } else {
                r.residue_info.group = "unknown";
            }
            return std::make_pair(0, r);
        }
    };
    static const dummy_geometry geom;
    std::pair<int, dummy_geometry::restraints_t> p = geom.get_monomer_restraints(monomer_name, 0);
    return (p.second.residue_info.group == "RNA" || p.second.residue_info.group == "DNA");
}