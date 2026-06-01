#include "function.hh"
#include <gemmi/elem.hpp>
#include <iostream>
#include <map>

namespace coot {

std::unique_ptr<gemmi::Residue> residue_from_rdkit_mol_gemmi(
    const RDKit::ROMol &mol, int conf_id, const std::string &new_comp_id) {

  std::cout << "DEBUG:: residue_from_rdkit_mol_gemmi()::::::::::::::::::::::::::::::: --- start ---" << std::endl;
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  std::unique_ptr<gemmi::Residue> r = nullptr;
  std::vector<gemmi::Atom> atoms;
  unsigned int n_atoms = mol.getNumAtoms();
  std::map<std::string, unsigned int> ele_count;

  if (n_atoms > 0) {
    // Use the specified conf_id, or fall back to 0 if invalid
    const RDKit::Conformer *conf = nullptr;
    if (mol.getNumConformers() > 0) {
      conf = &mol.getConformer(conf_id);
    }

    for (unsigned int i = 0; i < n_atoms; i++) {
      const RDKit::Atom *rat = mol.getAtomWithIdx(i);
      int atomic_number = rat->getAtomicNum();
      std::string ele = tbl->getElementSymbol(atomic_number);
      ele_count[ele]++;

      try {
        std::string name;
        rat->getProp("name", name); // They must have been made by caller

        // Pad name to 4 characters (PDB style)
        if (name.size() == 1) name = " " + name + "  ";
        if (name.size() == 2) name = " " + name + " ";
        if (name.size() == 3) name = " " + name;

        gemmi::Atom at;
        at.name = name;
        at.element = gemmi::Element(atomic_number);

        RDGeom::Point3D p = conf->getAtomPos(i);
        at.pos = gemmi::Position(p.x, p.y, p.z);
        at.occ = 1.0;
        at.b_iso = 30.0;

        std::cout << "DEBUG:: residue atom " << i << " \"" << name << "\" at " << p.x << " " << p.y << " " << p.z << std::endl;
        atoms.push_back(at);
      } catch (const std::runtime_error &rte) {
        std::cout << "DEBUG:: rte " << rte.what() << std::endl;
      }
    }

    if (!atoms.empty()) {
      r = std::make_unique<gemmi::Residue>();
      r->name = new_comp_id;
      // Copy atoms into the residue
      r->atoms.insert(r->atoms.end(), atoms.begin(), atoms.end());
    }
  }

  std::cout << "DEBUG:: residue_from_rdkit_mol_gemmi()::::::::::::::::::::::::::::::: returning " << r.get() << std::endl;
  return r;
}

} // namespace coot