#pragma once
#include <gemmi/model.hpp>
#include <cctype>
#include <string>
#include <utility>
#include <vector>

namespace coot {

static inline std::string pad_atom_name_pdb(const std::string& name) {
    if (name.size() >= 4) return name.substr(0, 4);
    if (!name.empty() && std::isalpha(static_cast<unsigned char>(name[0]))) {
        std::string padded = " " + name;
        while (padded.size() < 4) padded += " ";
        return padded;
    }
    std::string padded = name;
    while (padded.size() < 4) padded += " ";
    return padded;
}

class dictionary_residue_restraints_t {
public:
    bool change_names_gemmi(
        gemmi::Residue* residue,
        const std::vector<std::pair<std::string, std::string>>& change_name,
        const std::string& new_comp_id) const;
};

inline bool dictionary_residue_restraints_t::change_names_gemmi(
    gemmi::Residue* residue,
    const std::vector<std::pair<std::string, std::string>>& change_name,
    const std::string& new_comp_id) const
{
    bool changed_something = false;
    if (residue) {
        for (auto& at : residue->atoms) {
            std::string atom_name = pad_atom_name_pdb(at.name);
            for (const auto& change : change_name) {
                if (change.first == atom_name) {
                    at.name = change.second;
                    changed_something = true;
                    break;
                }
            }
        }
    }
    if (changed_something) {
        residue->name = new_comp_id;
    }
    return changed_something;
}
} // namespace coot