#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

class Bond_lines_container {
public:
    bool draw_these_atom_contacts_gemmi(gemmi::CRA cra_this,
                                        gemmi::CRA cra_env,
                                        coot::protein_geometry* protein_geom) {
        bool draw_flag = true;

        // Access parent info via CRA pointers
        gemmi::Chain* ch_this = cra_this.chain;
        gemmi::Chain* ch_env  = cra_env.chain;

        gemmi::Residue* this_residue = cra_this.residue;
        gemmi::Residue* env_residue  = cra_env.residue;

        if (ch_this != ch_env) {
            return true;
        } else {
            if (this_residue == env_residue) {
                return false;
            } else {
                if (std::abs(this_residue->seqid.num.value - env_residue->seqid.num.value) > 1) {
                    return true;
                } else {
                    // OK, we have neighbouring residues in the same chain
                    std::string this_res_type = this_residue->name;
                    std::string env_residue_res_type = env_residue->name;
                    if (!protein_geom->linkable_residue_types_p(this_res_type, env_residue_res_type)) {
                        return true;
                    } else {
                        // In gemmi, atom.name is unpadded (e.g. "N" not " N  ")
                        std::string this_atom_name = cra_this.atom->name;
                        std::string env_atom_name  = cra_env.atom->name;

                        // PDBv3 FIXME
                        if (this_atom_name == "N")
                            if (env_atom_name == "CA")
                                draw_flag = false;

                        if ((this_atom_name == "N") || (this_atom_name == "CA") ||
                            (this_atom_name == "C")  || (this_atom_name == "O")  ||
                            (this_atom_name == "H"))
                            if ((env_atom_name == "N") || (env_atom_name == "CA") ||
                                (env_atom_name == "C")  || (env_atom_name == "O")  ||
                                (env_atom_name == "H"))
                                draw_flag = false;

                        if ((this_atom_name == "O3'") || (this_atom_name == "C3'") ||
                            (this_atom_name == "P")   || (this_atom_name == "OP1") ||
                            (this_atom_name == "OP2") || (this_atom_name == "O5'") ||
                            (this_atom_name == "C5'"))
                            if ((env_atom_name == "O3'") || (env_atom_name == "C3'") ||
                                (env_atom_name == "P")   || (env_atom_name == "OP1") ||
                                (env_atom_name == "OP2") || (env_atom_name == "O5'") ||
                                (env_atom_name == "C5'"))
                                draw_flag = false;
                    }
                }
            }
        }
        return draw_flag;
    }
};