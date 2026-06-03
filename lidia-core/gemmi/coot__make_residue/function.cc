#include "function.hh"
#include <gemmi/elem.hpp>
#include "lidia-core/rdkit-interface.hh"

namespace coot {

gemmi::Residue* make_residue_gemmi(const RDKit::ROMol& rdkm, int iconf, const std::string& res_name) {
    // replace this function by making a residue directly instead of via a molfile.
    // If there are no atom names, make them from the element and atom number

    gemmi::Residue* residue_p = nullptr;
    lig_build::molfile_molecule_t mol = coot::make_molfile_molecule(rdkm, iconf);

    // now convert mol to a gemmi::Residue *
    if (mol.atoms.size()) {
        residue_p = new gemmi::Residue;
        residue_p->name = res_name;
        residue_p->seqid.num.value = 1;
        residue_p->seqid.icode = ' ';

        for (unsigned int iat = 0; iat < mol.atoms.size(); iat++) {
            gemmi::Atom at;
            at.name = mol.atoms[iat].name; // overridden hopefully
            at.element = gemmi::Element(mol.atoms[iat].element);
            at.pos = gemmi::Position(
                mol.atoms[iat].atom_position.x(),
                mol.atoms[iat].atom_position.y(),
                mol.atoms[iat].atom_position.z());
            at.occ = 1.0;
            at.b_iso = 30.0;
            // Set hetatm flag (is_hetatm)
            at.flag = 'H'; // HETATM
            residue_p->atoms.push_back(at);
        }
    }

    return residue_p;
}

} // namespace coot