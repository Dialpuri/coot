#pragma once
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_his_ring_H__9cb67e/gemmi/function.hh"

namespace coot {
namespace reduce {

inline void
add_his_ring_C_Hs_gemmi(gemmi::Residue& residue, double bl_arom) {
   add_his_ring_H_gemmi(" HD2", " CG ", " CD2", "NE2", bl_arom, residue);
   add_his_ring_H_gemmi(" HE1", " ND1", " CE1", "NE2", bl_arom, residue);
}

} // namespace reduce
} // namespace coot
