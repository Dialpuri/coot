#include "function.hh"

#include <string>
#include <iostream>

namespace coot {
namespace rotamer {

// Pad atom name to MMDB GetAtomName style: 1 space + name, right-padded to 4
static std::string pad4(const std::string& s) {
    std::string r = " " + s;
    while (r.size() < 4) r += ' ';
    return r;
}

std::vector<gemmi::Atom*> ordered_residue_atoms_gemmi(const gemmi::Residue& residue) {
    std::vector<gemmi::Atom*> atom_vec;
    std::vector<gemmi::Atom*> store_vec;

    std::string residue_name(residue.name);
    std::vector<std::string> atom_order;

    if (residue_name == "VAL") {
        atom_order.push_back(" CG1");
        atom_order.push_back(" CG2");
    }
    if (residue_name == "ASP") {
        atom_order.push_back(" OD1");
        atom_order.push_back(" OD2");
    }
    if (residue_name == "ASN") {
        atom_order.push_back(" OD1");
        atom_order.push_back(" ND2");
    }
    if (residue_name == "GLN") {
        atom_order.push_back(" OE1");
        atom_order.push_back(" NE2");
    }
    if (residue_name == "GLU") {
        atom_order.push_back(" OE1");
        atom_order.push_back(" OE2");
    }
    if (residue_name == "PHE") {
        atom_order.push_back(" CD1");
        atom_order.push_back(" CD2");
    }
    if (residue_name == "TYR") {
        atom_order.push_back(" CD1");
        atom_order.push_back(" CD2");
    }
    if (residue_name == "TRP") {
        atom_order.push_back(" CD1");
        atom_order.push_back(" CD2");
    }
    if (residue_name == "THR") {
        atom_order.push_back(" OG1");
        atom_order.push_back(" CG2");
    }
    if (residue_name == "LEU") {
        atom_order.push_back(" CD1");
        atom_order.push_back(" CD2");
    }
    if (residue_name == "ILE") {
        atom_order.push_back(" CG1");
        atom_order.push_back(" CG2");
    }
    if (residue_name == "PRO") {
        atom_order.push_back(" N  ");
        atom_order.push_back(" CA ");
        atom_order.push_back(" CB ");
        atom_order.push_back(" CG ");
        atom_order.push_back(" CD ");
        atom_order.push_back(" C  ");
        atom_order.push_back(" O  ");
    }

    const std::vector<gemmi::Atom>& residue_atoms = residue.atoms;
    int n_residue_atoms = static_cast<int>(residue_atoms.size());

    if (atom_order.size() == 0) {
        for (int i = 0; i < n_residue_atoms; i++) {
            atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
        }
    } else {
        if (n_residue_atoms == 0) {
            std::cout << " something broken in atom residue selection in "
                      << "ordered_residue_atoms, got 0 atoms" << std::endl;
        } else {
            int sensitive_count = 0;
            short int needs_reorder = 0;

            if (atom_order.size() == 2) {
                for (int i = 0; i < n_residue_atoms; i++) {
                    std::string atom_name = pad4(residue_atoms[i].name);

                    if (atom_name == atom_order[1]) {
                        sensitive_count++;
                    } else {
                        if (sensitive_count && atom_name == atom_order[0]) {
                            needs_reorder = 1;
                        }
                    }
                }

                if (needs_reorder) {
                    int passed_sensitive_atoms = 0;
                    for (int i = 0; i < n_residue_atoms; i++) {
                        std::string atom_name = pad4(residue_atoms[i].name);
                        if (passed_sensitive_atoms < sensitive_count) {
                            if (atom_name == atom_order[0]) {
                                passed_sensitive_atoms++;
                                atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                            } else {
                                if (atom_name == atom_order[1]) {
                                    store_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                                } else {
                                    atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                                }
                            }
                        } else {
                            atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                        }
                        if (passed_sensitive_atoms == sensitive_count) {
                            for (unsigned int j = 0; j < store_vec.size(); j++) {
                                atom_vec.push_back(store_vec[j]);
                            }
                            passed_sensitive_atoms++;
                        }
                    }
                } else {
                    for (int i = 0; i < n_residue_atoms; i++)
                        atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                }
            } else {
                std::vector<gemmi::Atom*> other_atoms;
                for (unsigned int ian = 0; ian < atom_order.size(); ian++) {
                    for (int i = 0; i < n_residue_atoms; i++) {
                        std::string atom_name = pad4(residue_atoms[i].name);
                        if (atom_name == atom_order[ian]) {
                            atom_vec.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                        }
                    }
                }

                for (int i = 0; i < n_residue_atoms; i++) {
                    bool found = false;
                    for (unsigned int in_atom = 0; in_atom < atom_vec.size(); in_atom++) {
                        if (atom_vec[in_atom] == &residue_atoms[i]) {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        other_atoms.push_back(const_cast<gemmi::Atom*>(&residue_atoms[i]));
                }
                for (unsigned int i = 0; i < other_atoms.size(); i++)
                    atom_vec.push_back(other_atoms[i]);
            }
        }
    }

    return atom_vec;
}

} // namespace rotamer
} // namespace coot