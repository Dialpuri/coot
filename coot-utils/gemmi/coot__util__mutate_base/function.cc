#include "function.hh"
#include <gemmi/elem.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__residue_atoms_segid/gemmi/function.hh"

std::string coot::convert_base_name(const std::string &std_base_name, bool use_old_style_naming) {
    if (use_old_style_naming) {
        return std_base_name;
    } else {
        if (std_base_name == "Cd")  return "DC";
        if (std_base_name == "Ad")  return "DA";
        if (std_base_name == "Gd")  return "DG";
        if (std_base_name == "Td")  return "DT";
        if (std_base_name == "Cr")  return "C";
        if (std_base_name == "Ar")  return "A";
        if (std_base_name == "Gr")  return "G";
        if (std_base_name == "Ur")  return "U";
        if (std_base_name == "Tr")  return "T";
    }
    return std_base_name;
}

namespace coot {
namespace util {

void mutate_base_gemmi(
  gemmi::Residue &residue,
  const gemmi::Residue &std_base,
  bool use_old_style_naming,
  bool print_match_stats_flag,
  float b_factor)
{
    bool debug = false;

    // Unpadded atom name lists for each base type (gemmi stores names without padding)
    auto adenine = std::vector<std::string>{
        "N9","C8","N7","C5","C4",
        "N1","C2","N3","C6","N6"};

    auto guanine = std::vector<std::string>{
        "N9","C8","N7","C5","C4",
        "N1","C2","N3","C6","O6","N2"};

    auto thymine = std::vector<std::string>{
        "N1","C2","N3","C4","C5","C6",
        "O2","O4","C7","C5M"};

    auto cytosine = std::vector<std::string>{
        "N1","C2","N3","C4","C5","C6",
        "O2","N4"};

    auto uracil = std::vector<std::string>{
        "N1","C2","N3","C4","C5","C6",
        "O2","O4"};

    // Match-order atom lists (unpadded)
    auto purine = std::vector<std::string>{
        "N9","C4","C5","N7","C8"};

    auto pyrimidine = std::vector<std::string>{
        "N1","C2","N3","C5","C6","C4"};

    std::string old_seg_id_for_residue_atoms;
    bool use_old_seg_id = false;
    try {
        old_seg_id_for_residue_atoms = coot::residue_atoms_segid_gemmi(residue);
        use_old_seg_id = true;
    } catch (const std::runtime_error &) {
    }

    short int mol_base_is_pyrimidine = -1;
    short int mol_base_is_purine     = -1;
    short int std_base_is_pyrimidine = -1;
    short int std_base_is_purine     = -1;

    std::string mol_base_name = residue.name;
    std::string std_base_name = std_base.name;

    if (mol_base_name == "A" || mol_base_name == "DA" ||
        mol_base_name == "G" || mol_base_name == "DG") {
        mol_base_is_purine = 1;
        mol_base_is_pyrimidine = 0;
    }
    if (mol_base_name == "C" || mol_base_name == "DC" ||
        mol_base_name == "U" || mol_base_name == "DU" ||
        mol_base_name == "T" || mol_base_name == "DT") {
        mol_base_is_pyrimidine = 1;
        mol_base_is_purine = 0;
    }

    if (std_base_name == "A" || std_base_name == "DA" ||
        std_base_name == "G" || std_base_name == "DG") {
        std_base_is_purine = 1;
        std_base_is_pyrimidine = 0;
    }
    if (std_base_name == "C" || std_base_name == "DC" ||
        std_base_name == "T" || std_base_name == "DT" ||
        std_base_name == "U" || std_base_name == "DU") {
        std_base_is_pyrimidine = 1;
        std_base_is_purine = 0;
    }

    if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) ||
        (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1)) {
        std::cout << "ERROR:: mutate_base() unassigned type "
                  << "mol_base_is_pyrimidine:" << " " << mol_base_is_pyrimidine
                  << " mol_base_is_purine: " << " " << mol_base_is_purine
                  << " std_base_is_pyrimidine: " << " " << std_base_is_pyrimidine
                  << " std_base_is_purine: " << " " << std_base_is_purine
                  << " " << mol_base_name << " " << std_base_name << std::endl;
        return;
    }

    int n_match_atoms = 5;
    if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
        n_match_atoms = 6;

    const std::vector<std::string> &moving_name_vector =
        (std_base_is_purine) ? purine : pyrimidine;
    const std::vector<std::string> &refrce_name_vector =
        (mol_base_is_purine) ? purine : pyrimidine;

    const auto &mol_base_atoms = residue.atoms;
    const auto &std_base_atoms = std_base.atoms;

    int n_mol_base_atoms = static_cast<int>(mol_base_atoms.size());
    int n_std_base_atoms = static_cast<int>(std_base_atoms.size());

    std::vector<clipper::Coord_orth> refrce_atom_positions;
    std::vector<clipper::Coord_orth> moving_atom_positions;

    if (debug) {
        for (const auto &n : refrce_name_vector)
            std::cout << "ref base search atom :" << n << ":" << std::endl;
        for (const auto &n : moving_name_vector)
            std::cout << "mov base search atom :" << n << ":" << std::endl;
    }

    for (int j = 0; j < n_match_atoms; j++) {
        for (int i = 0; i < n_mol_base_atoms; i++) {
            std::string atom_name = mol_base_atoms[i].name;
            if (refrce_name_vector[j] == atom_name) {
                refrce_atom_positions.push_back(
                    clipper::Coord_orth(mol_base_atoms[i].pos.x,
                                        mol_base_atoms[i].pos.y,
                                        mol_base_atoms[i].pos.z));
                if (debug)
                    std::cout << "Found " << atom_name << " in reference " << std::endl;
            }
        }
    }

    for (int j = 0; j < n_match_atoms; j++) {
        for (int i = 0; i < n_std_base_atoms; i++) {
            std::string atom_name = std_base_atoms[i].name;
            if (moving_name_vector[j] == atom_name) {
                moving_atom_positions.push_back(
                    clipper::Coord_orth(std_base_atoms[i].pos.x,
                                        std_base_atoms[i].pos.y,
                                        std_base_atoms[i].pos.z));
                if (debug)
                    std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
            }
        }
    }

    if (static_cast<int>(refrce_atom_positions.size()) != n_match_atoms) {
        std::cout << "ERROR:: wrong number of reference atoms found! "
                  << refrce_atom_positions.size() << std::endl;
        return;
    }

    if (static_cast<int>(moving_atom_positions.size()) != n_match_atoms) {
        std::cout << "ERROR:: wrong number of moving atoms found! "
                  << moving_atom_positions.size() << std::endl;
        return;
    }

    clipper::RTop_orth rtop(moving_atom_positions, refrce_atom_positions);

    double sum_dist = 0.0;
    double sum_dist2 = 0.0;
    double mind  =  999999999.9;
    double maxd = -999999999.9;
    for (unsigned int i = 0; i < refrce_atom_positions.size(); i++) {
        double d = clipper::Coord_orth::length(
            refrce_atom_positions[i],
            clipper::Coord_orth(moving_atom_positions[i].transform(rtop)));
        sum_dist  += d;
        sum_dist2 += d * d;
        if (d > maxd) maxd = d;
        if (d < mind) mind = d;
    }
    double mean = sum_dist / static_cast<double>(moving_atom_positions.size());
    double var  = sum_dist2 / static_cast<double>(moving_atom_positions.size());

    std::cout << "INFO:: " << moving_atom_positions.size() << " matched atoms had: \n"
              << "   mean devi: " << mean << "\n"
              << "    rms devi: " << std::sqrt(var) << "\n"
              << "    max devi: " << maxd << "\n"
              << "    min devi: " << mind << std::endl;

    std::vector<std::string> mol_base_atom_names;
    if (mol_base_name == "A" || mol_base_name == "DA")      mol_base_atom_names = adenine;
    else if (mol_base_name == "G" || mol_base_name == "DG") mol_base_atom_names = guanine;
    else if (mol_base_name == "C" || mol_base_name == "DC") mol_base_atom_names = cytosine;
    else if (mol_base_name == "T" || mol_base_name == "DT") mol_base_atom_names = thymine;
    else if (mol_base_name == "U" || mol_base_name == "DU") mol_base_atom_names = uracil;

    if (mol_base_atom_names.empty()) {
        std::cout << "ERROR:: muate_base(): ";
        std::cout << "failed to find mol_base_name for mol_base_atom_names\n";
        return;
    }

    std::vector<std::string> std_base_atom_names;
    if (std_base_name == "A" || std_base_name == "DA")      std_base_atom_names = adenine;
    else if (std_base_name == "G" || std_base_name == "DG") std_base_atom_names = guanine;
    else if (std_base_name == "C" || std_base_name == "DC") std_base_atom_names = cytosine;
    else if (std_base_name == "T" || std_base_name == "DT") std_base_atom_names = thymine;
    else if (std_base_name == "U" || std_base_name == "DU") std_base_atom_names = uracil;

    if (std_base_atom_names.empty()) {
        std::cout << "ERROR:: muate_base(): ";
        std::cout << "failed to find std_base_name for std_base_atom_names\n";
        return;
    }

    // Delete atoms from the residue that match mol_base_atom_names
    std::vector<bool> to_delete(residue.atoms.size(), false);
    for (const auto &aname : mol_base_atom_names) {
        for (size_t i = 0; i < residue.atoms.size(); i++) {
            if (aname == residue.atoms[i].name) {
                if (debug)
                    std::cout << ".... Deleting Atom " << residue.atoms[i].name
                              << " i = " << i << std::endl;
                to_delete[i] = true;
                break;
            }
        }
    }

    // gemmi does not store TER atoms; no TER handling needed

    // Remove marked atoms (iterate backwards to keep indices valid)
    for (int i = static_cast<int>(residue.atoms.size()) - 1; i >= 0; i--) {
        if (to_delete[i]) {
            residue.atoms.erase(residue.atoms.begin() + i);
        }
    }

    // Add atoms from std_base, transformed
    for (const auto &std_aname : std_base_atom_names) {
        bool found = false;
        for (int i = 0; i < n_std_base_atoms; i++) {
            if (std_aname == std_base_atoms[i].name) {
                clipper::Coord_orth p(std_base_atoms[i].pos.x,
                                      std_base_atoms[i].pos.y,
                                      std_base_atoms[i].pos.z);
                clipper::Coord_orth pt = p.transform(rtop);
                std::string ele = std_aname.substr(0, 2);

                if (debug)
                    std::cout << ".... Adding Atom " << std_base_atoms[i].name << std::endl;

                gemmi::Atom new_atom;
                new_atom.pos = gemmi::Position(pt.x(), pt.y(), pt.z());
                new_atom.occ = 1.0;
                new_atom.b_iso = b_factor;

                std::string new_atom_name = std_base_atoms[i].name;
                if (std_base_name == "DT") {
                    if (new_atom_name == "C5M") {
                        if (!use_old_style_naming) {
                            new_atom_name = "C7";
                        }
                    }
                }
                new_atom.name = new_atom_name;
                new_atom.element = gemmi::find_element(ele.c_str());
                new_atom.altloc = '\0';

                residue.atoms.push_back(new_atom);
                found = true;
                break;
            }
        }
    }

    std::string new_base_name = coot::convert_base_name(std_base_name, use_old_style_naming);
    residue.name = new_base_name;
}

} // util
} // coot