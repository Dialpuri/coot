#pragma once

#include "mini-mol/mini-mol.hh"

namespace coot {
namespace minimol {

// Static factory method replacing the MMDB constructor
// Original: residue(int i, const std::string &resname) : ins_code(""), name(resname) { seqnum = i; }
inline residue residue_gemmi(int i, const std::string& resname) {
    residue r;
    r.seqnum = i;
    r.name = resname;
    r.ins_code = "";
    return r;
}

} // namespace minimol
} // namespace coot
