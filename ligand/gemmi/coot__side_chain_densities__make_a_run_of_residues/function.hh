#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class side_chain_densities {
public:
    std::vector<gemmi::CRA> make_a_run_of_residues_gemmi(const gemmi::Structure& st,
                                                         const std::string& chain_id,
                                                         int resno_start, int resno_end) const;
};

}