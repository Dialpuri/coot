#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_pair_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot { namespace reduce {

void add_2_sp3_hydrogens_gemmi(const std::string &H_at_name_1,
                               const std::string &H_at_name_2,
                               const std::string &at_name_1,
                               const std::string &at_name_2,
                               const std::string &at_name_3,
                               double bond_length,
                               double angle_between_Hs,
                               gemmi::Residue &residue,
                               bool choose_only_farthest_position);

}} // namespace coot::reduce