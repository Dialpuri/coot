#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__crankshaft_set__crankshaft_set/gemmi/function.hh"

namespace coot {

class triple_crankshaft_set {
public:
    crankshaft_set cs[3];
    std::vector<std::string> residue_types;

    triple_crankshaft_set(gemmi::Residue *res_0,
                          gemmi::Residue *res_1,
                          gemmi::Residue *res_2,
                          gemmi::Residue *res_3,
                          gemmi::Residue *res_4,
                          gemmi::Residue *res_5,
                          const std::vector<std::string> &residue_types_in)
        : cs{ crankshaft_set(*res_0, *res_1, *res_2, *res_3),
              crankshaft_set(*res_1, *res_2, *res_3, *res_4),
              crankshaft_set(*res_2, *res_3, *res_4, *res_5) } {
        residue_types = residue_types_in;
    }

    const crankshaft_set &operator[](unsigned int i) const { return cs[i]; }
    crankshaft_set &operator[](unsigned int i) { return cs[i]; }
    const std::string &residue_type(unsigned int idx) const { return residue_types[idx]; }
};

} // namespace coot