#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <ostream>

namespace coot {
namespace molecule_t {
    inline unsigned int get_number_of_atoms_gemmi(const gemmi::Structure& st) {
        if (st.models.empty()) return 0;
        return static_cast<unsigned int>(gemmi::count_atom_sites(st) + 1);
    }
}
}

/// Gemmi port of operator<< for MyMolecule.
/// Outputs the structure name and total atom count to the stream.
std::ostream& operator_ltlt_gemmi(std::ostream& o, const gemmi::Structure& st);

inline std::ostream& operator_ltlt_gemmi(std::ostream& o, const gemmi::Structure& st)
{
    o << "Original name:" << st.name << "\n"
      << "nAtoms:" << coot::molecule_t::get_number_of_atoms_gemmi(st);
    return o;
}