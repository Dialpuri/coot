#include "function.hh"

namespace {
    std::string trim(const std::string &s) {
        size_t start = s.find_first_not_of(" ");
        size_t end = s.find_last_not_of(" ");
        if (start != std::string::npos)
            return s.substr(start, end - start + 1);
        return s;
    }
} // namespace

bool
Bond_lines_container::invert_deloc_bond_displacement_vector_gemmi(
    const clipper::Coord_orth &vect,
    int iat_1,
    int iat_2,
    const std::vector<const gemmi::Atom*> &residue_atoms,
    int n_atoms,
    const std::vector<coot::dict_bond_restraint_t> &bond_restraints) const
{
    bool r = false;

    // Use atom_id_mmdb_expand to match the _4c() comparison format
    std::string atom_name_iat = coot::atom_id_mmdb_expand(trim(residue_atoms[iat_1]->name));
    std::string atom_name_jat = coot::atom_id_mmdb_expand(trim(residue_atoms[iat_2]->name));

    std::map<std::string, int> atom_name_map;
    for (int iat = 0; iat < n_atoms; iat++) {
        atom_name_map[coot::atom_id_mmdb_expand(trim(residue_atoms[iat]->name))] = iat;
    }

    for (unsigned int ib = 0; ib < bond_restraints.size(); ib++) {
        if (bond_restraints[ib].atom_id_1_4c() == atom_name_iat) {
            if (bond_restraints[ib].atom_id_2_4c() != atom_name_jat) {
                if (bond_restraints[ib].type() == "deloc") {
                    clipper::Coord_orth pt_1(residue_atoms[iat_1]->pos.x,
                                             residue_atoms[iat_1]->pos.y,
                                             residue_atoms[iat_1]->pos.z);
                    std::map<std::string, int>::const_iterator it;
                    it = atom_name_map.find(bond_restraints[ib].atom_id_2_4c());
                    if (it != atom_name_map.end()) {
                        clipper::Coord_orth pt_2(residue_atoms[it->second]->pos.x,
                                                 residue_atoms[it->second]->pos.y,
                                                 residue_atoms[it->second]->pos.z);
                        clipper::Coord_orth diff = pt_2 - pt_1;
                        double d = clipper::Coord_orth::dot(vect, diff);
                        if (d < 0) {
                            r = true;
                            break;
                        }
                    }
                }
            }
        }

        // same again, restraints ordered differently
        if (bond_restraints[ib].atom_id_2_4c() == atom_name_iat) {
            if (bond_restraints[ib].atom_id_1_4c() != atom_name_jat) {
                if (bond_restraints[ib].type() == "deloc") {
                    clipper::Coord_orth pt_1(residue_atoms[iat_1]->pos.x,
                                             residue_atoms[iat_1]->pos.y,
                                             residue_atoms[iat_1]->pos.z);
                    std::map<std::string, int>::const_iterator it;
                    it = atom_name_map.find(bond_restraints[ib].atom_id_1_4c());
                    if (it != atom_name_map.end()) {
                        clipper::Coord_orth pt_2(residue_atoms[it->second]->pos.x,
                                                 residue_atoms[it->second]->pos.y,
                                                 residue_atoms[it->second]->pos.z);
                        clipper::Coord_orth diff = pt_2 - pt_1;
                        double d = clipper::Coord_orth::dot(vect, diff);
                        if (d < 0) {
                            r = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    return r;
}