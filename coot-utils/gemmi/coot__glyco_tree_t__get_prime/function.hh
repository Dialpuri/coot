#pragma once

#include "coot-utils/tree.hh"
#include <gemmi/model.hpp>
#include <string>

namespace coot {

struct linked_residue_t_gemmi {
  gemmi::Residue* residue = nullptr;
  std::string residue_name;
  std::string link_type;

  linked_residue_t_gemmi() = default;
  linked_residue_t_gemmi(gemmi::Residue* r, const std::string& rn, const std::string& lt)
    : residue(r), residue_name(rn), link_type(lt) {}
};

// gemmi port of glyco_tree_t::get_prime — free function because C++ does not
// allow re-opening coot::glyco_tree_t (which is mmdb-based).
//
// Returns:
//   0 = UNSET     (residue not in tree or not connected to BMA)
//   1 = PRIME     (connected via ALPHA1-6 to BMA)
//   2 = NON_PRIME (connected via ALPHA1-3 to BMA)
int get_prime_gemmi(const tree<linked_residue_t_gemmi>& glyco_tree,
                    gemmi::Residue* residue_p) {
    int result = 0; // UNSET

    for (auto it = glyco_tree.begin(); it != glyco_tree.end(); ++it) {
        if (it->residue == residue_p) {
            auto this_one = it;
            bool has_parent = true;

            while (has_parent) {
                if (!this_one.node->parent) {
                    has_parent = false;
                } else {
                    if (this_one.node->parent->data.residue_name == "BMA") {
                        if (this_one->link_type == "ALPHA1-3")
                            result = 2; // NON_PRIME
                        if (this_one->link_type == "ALPHA1-6")
                            result = 1; // PRIME
                    }
                    this_one = this_one.node->parent;
                }
            }
        }
    }
    return result;
}

} // namespace coot