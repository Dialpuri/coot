#pragma once

#include <gemmi/model.hpp>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__glyco_tree_t__get_level/gemmi/function.hh"

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

struct residue_id_t_gemmi {
    int level = 0;
    int prime_flag = 0;
    std::string res_type;
    std::string link_type;
    std::string parent_res_type;
    residue_spec_t parent_spec;

    residue_id_t_gemmi() = default;

    residue_id_t_gemmi(int lvl, int pf, const std::string& rt,
                       const std::string& lt, const std::string& prt,
                       residue_spec_t ps)
        : level(lvl), prime_flag(pf), res_type(rt), link_type(lt),
          parent_res_type(prt), parent_spec(std::move(ps)) {}
};

residue_id_t_gemmi get_id_gemmi(const tree<linked_residue_t_gemmi>& glyco_tree,
                                gemmi::Residue* residue_p);

}  // namespace coot