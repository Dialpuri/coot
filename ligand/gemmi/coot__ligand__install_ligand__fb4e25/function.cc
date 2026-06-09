#include "function.hh"
#include <iostream>

namespace coot {

void install_ligand_gemmi(ligand& lig, std::string ligand_pdb_filename) {
    std::cout << "Reading ligand pdb file: " << ligand_pdb_filename << std::endl;
    coot::minimol::molecule mol;
    coot::minimol::read_file_gemmi(mol, ligand_pdb_filename);
    lig.install_ligand(mol);
}

} // namespace coot
