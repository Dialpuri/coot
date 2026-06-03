#include "function.hh"

namespace coot {

void add_chem_mod_tree_gemmi(const ChemModTreeLoop& loop,
                             protein_geometry& geom)
{
    for (int j = 0; j < loop.GetLoopLength(); j++) {
        int ierr_tot = 0;

        std::string mod_id;
        std::string function;
        std::string atom_id;
        std::string atom_back;
        std::string back_type;
        std::string atom_forward;
        std::string connect_type;

        int ierr = 0;
        const char* s = nullptr;

        s = loop.GetString("mod_id", j, ierr);
        ierr_tot += ierr;
        if (s) mod_id = s;

        s = loop.GetString("function", j, ierr);
        ierr_tot += ierr;
        if (s) function = s;

        s = loop.GetString("atom_id", j, ierr);
        ierr_tot += ierr;
        if (s) atom_id = s;

        s = loop.GetString("atom_back", j, ierr);
        ierr_tot += ierr;
        if (s) atom_back = s;

        s = loop.GetString("back_type", j, ierr);
        ierr_tot += ierr;
        if (s) back_type = s;

        s = loop.GetString("atom_forward", j, ierr);
        ierr_tot += ierr;
        if (s) atom_forward = s;

        s = loop.GetString("connect_type", j, ierr);
        ierr_tot += ierr;
        if (s) connect_type = s;

        if (ierr_tot == 0) {
            chem_mod_tree cmt(
                function,
                atom_id_mmdb_expand(atom_id),
                atom_id_mmdb_expand(atom_back),
                back_type,
                atom_id_mmdb_expand(atom_forward),
                connect_type);
            geom.mods[mod_id].add_mod_tree(cmt);
        }
    }
}

} // namespace coot