#include "function.hh"

#include <string>
#include <iostream>

namespace coot {
namespace protein_geometry {

void add_chem_mod_atom_gemmi(const ChemModAtomLoop& loop,
                             std::map<std::string, chem_mod>& mods) {
    for (int j = 0; j < loop.GetLoopLength(); ++j) {
        int ierr_tot = 0;
        int ierr;

        std::string mod_id;
        std::string function;
        std::string atom_id;
        std::string new_atom_id;
        std::string new_type_symbol;
        std::string new_type_energy;
        double new_partial_charge = 0.0;

        const char* s = nullptr;

        s = loop.GetString("mod_id", j, ierr);
        ierr_tot += ierr;
        if (ierr) std::cout << "   oops getting mod_id " << std::endl;
        if (s) mod_id = s;

        s = loop.GetString("function", j, ierr);
        if (ierr) std::cout << "   oops getting function " << std::endl;
        ierr_tot += ierr;
        if (s) function = s;

        s = loop.GetString("atom_id", j, ierr);
        if (ierr) std::cout << "   oops getting atom_id " << std::endl;
        ierr_tot += ierr;
        if (s) atom_id = s;

        s = loop.GetString("new_atom_id", j, ierr);
        if (ierr) std::cout << "   oops getting new_atom_id " << std::endl;
        ierr_tot += ierr;
        if (s) new_atom_id = s;

        s = loop.GetString("new_type_symbol", j, ierr);
        if (ierr) std::cout << "   oops getting new_type_symbol " << std::endl;
        ierr_tot += ierr;
        if (s) new_type_symbol = s;

        s = loop.GetString("new_type_energy", j, ierr);
        if (ierr) std::cout << "   oops getting new_type_energy " << std::endl;
        ierr_tot += ierr;
        if (s) new_type_energy = s;

        ierr = loop.GetReal(new_partial_charge, "new_partial_charge", j);
        if (ierr) {
            new_partial_charge = 0.0; // dummy value
        }

        if (ierr_tot == 0) {
            chem_mod_atom cma(function, atom_id, new_atom_id,
                              new_type_symbol, new_type_energy,
                              new_partial_charge);
            mods[mod_id].add_mod_atom(cma);
        } else {
            std::cout << "oops in add_chem_mod_atom ierr_tot = "
                      << ierr_tot << std::endl;
            std::cout << "   mod_id: \"" << mod_id
                      << "\"    function: \"" << function << "\" atom_id: \""
                      << atom_id_mmdb_expand(atom_id) << "\" new_atom_id: \""
                      << atom_id_mmdb_expand(new_atom_id) << "\" new_type_symbol: \""
                      << new_type_symbol << "\" new_type_energy: \""
                      << new_type_energy << "\" new_partial_charge: \""
                      << new_partial_charge << "\"" << std::endl;
        }
    }
}

} // namespace protein_geometry
} // namespace coot