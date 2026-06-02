#include "function.hh"
#include <gemmi/elem.hpp>
#include "utils/coot-utils.hh"

namespace coot {

gemmi::Residue* get_residue_gemmi(
    protein_geometry& geom,
    const std::string& comp_id,
    int imol_enc,
    bool idealised_flag,
    bool try_autoload_if_needed,
    float b_factor)
{
    gemmi::Residue* residue_p = nullptr;

    bool r = geom.have_dictionary_for_residue_type(comp_id, imol_enc, try_autoload_if_needed);

    if (r) {
        for (unsigned int i = 0; i < geom.size(); i++) {
            const std::pair<int, dictionary_residue_restraints_t>& entry = geom[i];
            const dictionary_residue_restraints_t& rest = entry.second;
            if (rest.residue_info.comp_id == comp_id) {
                if (entry.first == imol_enc) {
                    // Build the residue from dictionary restraints — same logic
                    // as dictionary_residue_restraints_t::GetResidue but using
                    // gemmi types instead of mmdb.

                    std::vector<gemmi::Atom> atoms;
                    bool make_hetatoms = !coot::util::is_standard_residue_name(rest.residue_info.comp_id);

                    for (unsigned int iat = 0; iat < rest.atom_info.size(); iat++) {
                        gemmi::Position p(0, 0, 0);
                        bool flag_and_have_coords = false;

                        if (idealised_flag && rest.atom_info[iat].pdbx_model_Cartn_ideal.first) {
                            const clipper::Coord_orth& c = rest.atom_info[iat].pdbx_model_Cartn_ideal.second;
                            p = gemmi::Position(c.x(), c.y(), c.z());
                            flag_and_have_coords = true;
                        }

                        if (!flag_and_have_coords) {
                            if (rest.atom_info[iat].model_Cartn.first) {
                                const clipper::Coord_orth& c = rest.atom_info[iat].model_Cartn.second;
                                p = gemmi::Position(c.x(), c.y(), c.z());
                                flag_and_have_coords = true;
                            }
                        }

                        if (flag_and_have_coords) {
                            gemmi::Atom atom;
                            atom.name = rest.atom_info[iat].atom_id_4c;
                            atom.pos = p;
                            atom.occ = 1.0;
                            atom.b_iso = b_factor;
                            atom.element = gemmi::Element(rest.atom_info[iat].type_symbol);
                            // (hetero flag handled by residue name; standard residues
                            // like ALA/GLY will have make_hetatoms==false)
                            atoms.push_back(atom);
                        }
                    }

                    if (!atoms.empty()) {
                        residue_p = new gemmi::Residue();
                        residue_p->name = rest.residue_info.comp_id;
                        for (gemmi::Atom& atom : atoms) {
                            residue_p->atoms.push_back(std::move(atom));
                        }
                        break;
                    }
                }
            }
        }
    }

    return residue_p;
}

} // namespace coot