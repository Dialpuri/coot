#include "function.hh"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

namespace coot {
namespace util {

namespace {

// Helper to get atom name string from gemmi atom
std::string atom_name(const gemmi::Atom& atom) {
    return atom.name;
}

// Helper to check if residue is purine or pyrimidine
std::pair<bool, bool> classify_residue(const gemmi::Residue& res) {
    const std::string& name = res.name;
    bool is_purine = false;
    bool is_pyrimidine = false;
    
    if (name == "A" || name == "G" || name == "DA" || name == "DG" ||
        name == "Ar" || name == "Ad" || name == "Gr" || name == "Gd") {
        is_purine = true;
        is_pyrimidine = false;
    } else if (name == "C" || name == "T" || name == "U" || 
               name == "DC" || name == "DT" || name == "DU" ||
               name == "Cr" || name == "Td" || name == "Ur" || name == "Ud") {
        is_purine = false;
        is_pyrimidine = true;
    }
    
    return {is_purine, is_pyrimidine};
}

// Build atom name vectors for standard nucleotides
std::vector<std::string> make_atom_vector(const std::vector<const char*>& names) {
    std::vector<std::string> result;
    for (const char* name : names) {
        result.push_back(std::string(name));
    }
    return result;
}

} // anonymous namespace

std::pair<bool, clipper::RTop_orth>
nucleotide_to_nucleotide_gemmi(const gemmi::Residue& residue,
                               const gemmi::Residue& std_base,
                               bool use_old_style_naming) {
    bool good_rtop_flag = false;
    clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
    clipper::Coord_orth pt_dum(0,0,0);
    clipper::RTop_orth rtop(m_dum, pt_dum);
    
    // Atom name lists for standard nucleotides
    std::vector<const char*> adenine_names = {
        " N9 ", " C8 ", " N7 ", " C5 ", " C4 ", 
        " N1 ", " C2 ", " N3 ", " C6 ", " N6 "
    };
    std::vector<const char*> guanine_names = {
        " N9 ", " C8 ", " N7 ", " C5 ", " C4 ", 
        " N1 ", " C2 ", " N3 ", " C6 ", " O6 ", " N2 "
    };
    std::vector<const char*> thymine_names = {
        " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ",
        " O2 ", " O4 ", " C5M"
    };
    std::vector<const char*> cytosine_names = {
        " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ",
        " O2 ", " N4 "
    };
    std::vector<const char*> uracil_names = {
        " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ",
        " O2 ", " O4 "
    };
    std::vector<const char*> purine_names = {
        " N9 ", " C4 ", " C5 ", " N7 ", " C8 "
    };
    std::vector<const char*> pyrimidine_names = {
        " N1 ", " C2 ", " N3 ", " C5 ", " C6 ", " C4 "
    };
    
    std::vector<std::string> adenine = make_atom_vector(adenine_names);
    std::vector<std::string> guanine = make_atom_vector(guanine_names);
    std::vector<std::string> thymine = make_atom_vector(thymine_names);
    std::vector<std::string> cytosine = make_atom_vector(cytosine_names);
    std::vector<std::string> uracil = make_atom_vector(uracil_names);
    std::vector<std::string> purine = make_atom_vector(purine_names);
    std::vector<std::string> pyrimidine = make_atom_vector(pyrimidine_names);
    
    // Classify residues
    auto mol_class = classify_residue(residue);
    auto std_class = classify_residue(std_base);
    
    bool mol_is_purine = mol_class.first;
    bool mol_is_pyrimidine = mol_class.second;
    bool std_is_purine = std_class.first;
    bool std_is_pyrimidine = std_class.second;
    
    if (!mol_is_purine && !mol_is_pyrimidine) {
        return {false, rtop};
    }
    if (!std_is_purine && !std_is_pyrimidine) {
        return {false, rtop};
    }
    
    int n_match_atoms = 5;
    if (mol_is_pyrimidine && std_is_pyrimidine)
        n_match_atoms = 6;
    
    std::vector<std::string> moving_name_vector;
    std::vector<std::string> refrce_name_vector;
    
    if (std_is_purine)
        moving_name_vector = purine;
    else
        moving_name_vector = pyrimidine;
    
    if (mol_is_purine)
        refrce_name_vector = purine;
    else
        refrce_name_vector = pyrimidine;
    
    // Build coordinate vectors
    std::vector<clipper::Coord_orth> refrce_atom_positions;
    std::vector<clipper::Coord_orth> moving_atom_positions;
    
    // Match reference atoms
    for (int j = 0; j < n_match_atoms; j++) {
        for (const auto& atom : residue.atoms) {
            std::string atom_name = atom.name;
            if (refrce_name_vector[j] == atom_name) {
                refrce_atom_positions.push_back(
                    clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
            }
        }
    }
    
    // Match moving atoms
    for (int j = 0; j < n_match_atoms; j++) {
        for (const auto& atom : std_base.atoms) {
            std::string atom_name = atom.name;
            if (moving_name_vector[j] == atom_name) {
                moving_atom_positions.push_back(
                    clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
            }
        }
    }
    
    if (static_cast<int>(refrce_atom_positions.size()) != n_match_atoms) {
        return {false, rtop};
    }
    
    if (static_cast<int>(moving_atom_positions.size()) != n_match_atoms) {
        return {false, rtop};
    }
    
    // Add sugar/phosphate atoms for better alignment
    std::vector<std::string> const_nuc_atoms;
    std::string prime_char = "'";
    if (use_old_style_naming)
        prime_char = "*";
    const_nuc_atoms.push_back(" C1" + prime_char);
    const_nuc_atoms.push_back(" C2" + prime_char);
    const_nuc_atoms.push_back(" C3" + prime_char);
    const_nuc_atoms.push_back(" C4" + prime_char);
    const_nuc_atoms.push_back(" C5" + prime_char);
    const_nuc_atoms.push_back(" O3" + prime_char);
    const_nuc_atoms.push_back(" O4" + prime_char);
    const_nuc_atoms.push_back(" O5" + prime_char);
    const_nuc_atoms.push_back(" P " + prime_char);
    
    for (const auto& sugar_atom : const_nuc_atoms) {
        for (const auto& atom1 : std_base.atoms) {
            if (atom1.name == sugar_atom) {
                for (const auto& atom2 : residue.atoms) {
                    if (atom2.name == sugar_atom) {
                        // Check altloc match
                        if (atom1.altloc == atom2.altloc) {
                            moving_atom_positions.push_back(
                                clipper::Coord_orth(atom1.pos.x, atom1.pos.y, atom1.pos.z));
                            refrce_atom_positions.push_back(
                                clipper::Coord_orth(atom2.pos.x, atom2.pos.y, atom2.pos.z));
                        }
                    }
                }
            }
        }
    }
    
    // Compute the rotation/translation
    if (moving_atom_positions.size() >= 3 && refrce_atom_positions.size() >= 3) {
        rtop = clipper::RTop_orth(moving_atom_positions, refrce_atom_positions);
        good_rtop_flag = true;
    }
    
    return {good_rtop_flag, rtop};
}

} // namespace util
} // namespace coot