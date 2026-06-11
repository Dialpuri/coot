#include "function.hh"

namespace coot {

std::vector<coot::torsion_atom_quad_gemmi>
torsionable_bonds_monomer_internal_quads_gemmi(
    const gemmi::Residue &residue,
    std::vector<gemmi::Atom*> atom_selection,
    bool include_pyranose_ring_torsions_flag,
    coot::protein_geometry *geom_p) {

    std::vector<coot::torsion_atom_quad_gemmi> quads;
    bool hydrogen_torsions = false;
    std::string rn = residue.name;
    std::vector<coot::dict_torsion_restraint_t> tors_restraints =
        geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);

    bool is_pyranose = false;
    std::string group = coot::get_group_gemmi(*geom_p, residue);
    if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
        group == "D-SACCHARIDE" || group == "SACCHARIDE")
        is_pyranose = true;

    std::vector<std::string> residue_alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);

    int n_selected_atoms = static_cast<int>(atom_selection.size());

    for (unsigned int itor = 0; itor < tors_restraints.size(); itor++) {
        if (!tors_restraints[itor].is_const()) {
            std::string tor_atom_name[5];
            std::vector<gemmi::Atom*> ats(5, static_cast<gemmi::Atom*>(nullptr));

            tor_atom_name[1] = tors_restraints[itor].atom_id_1_4c();
            tor_atom_name[2] = tors_restraints[itor].atom_id_2_4c();
            tor_atom_name[3] = tors_restraints[itor].atom_id_3_4c();
            tor_atom_name[4] = tors_restraints[itor].atom_id_4_4c();

            if ((include_pyranose_ring_torsions_flag == 1) ||
                (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
                (!is_pyranose)) {
                for (unsigned int ialt = 0; ialt < residue_alt_confs.size(); ialt++) {
                    for (int iat = 0; iat < n_selected_atoms; iat++) {
                        // Use pad_atom_name to get PDB-format 4-char name for matching
                        std::string atom_name = coot::pad_atom_name(atom_selection[iat]->name,
                                                                    atom_selection[iat]->element.name());
                        char alt_char = atom_selection[iat]->altloc;
                        std::string alt_conf = (alt_char == ' ' || alt_char == '\0') ? "" : std::string(1, alt_char);
                        if (alt_conf == residue_alt_confs[ialt]) {
                            for (unsigned int jtor = 1; jtor < 5; jtor++) {
                                if (atom_name == tor_atom_name[jtor])
                                    ats[jtor] = atom_selection[iat];
                            }
                        }
                    }
                    if (ats[1] && ats[2] && ats[3] && ats[4]) {
                        coot::torsion_atom_quad_gemmi q(ats[1], ats[2], ats[3], ats[4],
                                                       tors_restraints[itor].angle(),
                                                       tors_restraints[itor].esd(),
                                                       tors_restraints[itor].periodicity());
                        q.name = tors_restraints[itor].id();
                        q.residue_name = rn;
                        quads.push_back(q);
                    }
                }
            }
        }
    }

    return quads;
}

} // namespace coot
