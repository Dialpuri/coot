#pragma once

#include <gemmi/model.hpp>

// gemmi port of fix_element_name_lengths
// In MMDB, element names were stored as mutable char arrays that needed padding.
// In gemmi, element is stored as gemmi::El (enum) — formatting for PDB output
// is handled automatically by gemmi's PDB writer. This function performs the
// same traversal but the "fix" is a no-op since gemmi doesn't store element
// names as mutable strings.
void fix_element_name_lengths_gemmi(gemmi::Structure& st) {
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom& atom : residue.atoms) {
                    // In gemmi, element is an enum; element.name() returns the
                    // canonical name (e.g. "N", "Mg"). No padding needed —
                    // PDB writer handles formatting automatically.
                    (void)atom.element.name();
                }
            }
        }
    }
}