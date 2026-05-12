#include "function.hh"
#include <string>
#include <vector>

namespace coot {

void pucker_analysis_info_t::assign_base_atom_coords_gemmi(gemmi::Residue* residue) {
    std::vector<std::string> cytidine_base_names;
    std::vector<std::string> uracil_base_names;
    std::vector<std::string> adenine_base_names;
    std::vector<std::string> guanine_base_names;
    std::vector<std::string> thymine_base_names;

    cytidine_base_names.push_back(" N1 ");
    cytidine_base_names.push_back(" C2 ");
    cytidine_base_names.push_back(" N3 ");
    cytidine_base_names.push_back(" C4 ");
    cytidine_base_names.push_back(" C5 ");
    cytidine_base_names.push_back(" C6 ");
    cytidine_base_names.push_back(" O2 ");
    cytidine_base_names.push_back(" N4 ");

    uracil_base_names.push_back(" N1 ");
    uracil_base_names.push_back(" C2 ");
    uracil_base_names.push_back(" N3 ");
    uracil_base_names.push_back(" C4 ");
    uracil_base_names.push_back(" C5 ");
    uracil_base_names.push_back(" C6 ");
    uracil_base_names.push_back(" O2 ");
    uracil_base_names.push_back(" O4 ");

    adenine_base_names.push_back(" N9 ");
    adenine_base_names.push_back(" C8 ");
    adenine_base_names.push_back(" N7 ");
    adenine_base_names.push_back(" C5 ");
    adenine_base_names.push_back(" C4 ");
    adenine_base_names.push_back(" N1 ");
    adenine_base_names.push_back(" C2 ");
    adenine_base_names.push_back(" N3 ");
    adenine_base_names.push_back(" C6 ");
    adenine_base_names.push_back(" N6 ");

    guanine_base_names.push_back(" N9 ");
    guanine_base_names.push_back(" C8 ");
    guanine_base_names.push_back(" N7 ");
    guanine_base_names.push_back(" C5 ");
    guanine_base_names.push_back(" C4 ");
    guanine_base_names.push_back(" N1 ");
    guanine_base_names.push_back(" C2 ");
    guanine_base_names.push_back(" N3 ");
    guanine_base_names.push_back(" C6 ");
    guanine_base_names.push_back(" O6 ");
    guanine_base_names.push_back(" N2 ");

    thymine_base_names.push_back(" N1 ");
    thymine_base_names.push_back(" C2 ");
    thymine_base_names.push_back(" N3 ");
    thymine_base_names.push_back(" C4 ");
    thymine_base_names.push_back(" C5 ");
    thymine_base_names.push_back(" C6 ");
    thymine_base_names.push_back(" O2 ");
    thymine_base_names.push_back(" O4 ");
    thymine_base_names.push_back(" C5M");

    base_atoms_coords.clear();
    N1_or_9 = nullptr;
    C1_prime = nullptr;

    // Assign N1_or_9 and C1_prime
    for (gemmi::Atom& atom : residue->atoms) {
        std::string alt_name(1, atom.altloc);
        if (alt_name == altconf) {
            if (atom.name == " N1 ")
                N1_or_9 = &atom;
            if (atom.name == " N9 ")
                N1_or_9 = &atom;
            if (atom.name == " C1*")
                C1_prime = &atom;
            if (atom.name == " C1'")
                C1_prime = &atom;
        }
    }

    // Fill base_names according to residue type/name
    std::vector<std::string> base_names;
    std::string residue_name = residue->name;

    // current names
    if (residue_name == "C") base_names = cytidine_base_names;
    if (residue_name == "U") base_names = uracil_base_names;
    if (residue_name == "A") base_names = adenine_base_names;
    if (residue_name == "G") base_names = guanine_base_names;
    // old names
    if (residue_name == "Cr") base_names = cytidine_base_names;
    if (residue_name == "Ur") base_names = uracil_base_names;
    if (residue_name == "Ar") base_names = adenine_base_names;
    if (residue_name == "Gr") base_names = guanine_base_names;
    // modern (3.x) RNA base names
    if (residue_name == "CYT") base_names = cytidine_base_names;
    if (residue_name == "URA") base_names = uracil_base_names;
    if (residue_name == "ADE") base_names = adenine_base_names;
    if (residue_name == "GUA") base_names = guanine_base_names;

    if (base_names.size() > 0) {
        for (gemmi::Atom& atom : residue->atoms) {
            for (const std::string& base_name : base_names) {
                if (base_name == atom.name) {
                    base_atoms_coords.push_back(atom.pos);
                }
            }
        }
    }
}

} // namespace coot