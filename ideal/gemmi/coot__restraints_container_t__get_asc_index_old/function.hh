#pragma once
#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>

namespace coot {

class restraints_container_t_gemmi {
public:
    restraints_container_t_gemmi(
        const std::vector<std::pair<bool, gemmi::Residue*>>& residues,
        gemmi::Structure* st);

    int get_asc_index_old_gemmi(
        const std::string& at_name, int resno, const std::string& chain_id) const;

    size_t get_residues_vec_size() const { return residues_.size(); }
    const std::vector<gemmi::Vec3>& get_atom_vec() const { return atom_; }

private:
    std::vector<std::pair<bool, gemmi::Residue*>> residues_;
    std::vector<gemmi::Vec3> atom_;
    gemmi::Structure* st_;
};

} // namespace coot
