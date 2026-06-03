#include "function.hh"

namespace coot {

std::pair<bool, bool>
atom_overlaps_container_t::is_h_bond_H_and_acceptor_gemmi(
    const gemmi::Atom& ligand_atom,
    const gemmi::Residue& ligand_res,
    const gemmi::Atom& env_atom,
    const gemmi::Residue& env_res
) {
    bool status = false;
    bool H_on_ligand = false;

    coot::hb_t hb_1 = get_h_bond_type_gemmi(ligand_atom, ligand_res);
    coot::hb_t hb_2 = get_h_bond_type_gemmi(env_atom, env_res);

    if (hb_1 != coot::hb_t::HB_NEITHER && hb_2 != coot::hb_t::HB_NEITHER) {

        if (hb_1 == coot::hb_t::HB_HYDROGEN) {
            if (hb_2 == coot::hb_t::HB_ACCEPTOR || hb_2 == coot::hb_t::HB_BOTH) {
                status = true;
                H_on_ligand = true;
            }
        }

        if (hb_1 == coot::hb_t::HB_ACCEPTOR || hb_1 == coot::hb_t::HB_BOTH) {
            if (hb_2 == coot::hb_t::HB_HYDROGEN) {
                status = true;
                H_on_ligand = false;
            }
        }
    }

    if (status == false) {
        // allow HOH to H-bond
        std::string resname_1 = ligand_res.name;
        std::string resname_2 = env_res.name;
        if (resname_1 == "HOH")
            if (hb_2 == coot::hb_t::HB_ACCEPTOR || hb_2 == coot::hb_t::HB_DONOR || hb_2 == coot::hb_t::HB_BOTH || hb_2 == coot::hb_t::HB_HYDROGEN)
                status = true;
        if (resname_2 == "HOH")
            if (hb_1 == coot::hb_t::HB_ACCEPTOR || hb_1 == coot::hb_t::HB_DONOR || hb_1 == coot::hb_t::HB_BOTH || hb_1 == coot::hb_t::HB_HYDROGEN)
                status = true;
    }

    return std::pair<bool, bool>(status, H_on_ligand);
}

} // namespace coot