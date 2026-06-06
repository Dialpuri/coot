#include "function.hh"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>

// Include coot headers FIRST (defines coot::atom_index_quad without torsion_gemmi)
#include "geometry/protein-geometry.hh"

// Then include gemmi ports that don't conflict on atom_index_quad
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__nomenclature__test_and_fix_PHE_TYR_nomenclature_errors/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_inverted_chiral_atom_p/gemmi/function.hh"

// --- Inline the ASP_GLU port to avoid atom_index_quad header conflict ---
// We need atom_index_quad with torsion_gemmi, so we define a local version
struct local_atom_index_quad {
    int index1;
    int index2;
    int index3;
    int index4;
    local_atom_index_quad(int i1, int i2, int i3, int i4)
        : index1(i1), index2(i2), index3(i3), index4(i4) {}

    double torsion_gemmi(const gemmi::Residue& residue) const {
        if (index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0)
            throw std::runtime_error("missing atoms");
        if (index1 >= int(residue.atoms.size()) || index2 >= int(residue.atoms.size()) ||
            index3 >= int(residue.atoms.size()) || index4 >= int(residue.atoms.size()))
            throw std::runtime_error("index out of range");
        const gemmi::Vec3& a1 = residue.atoms[index1].pos;
        const gemmi::Vec3& a2 = residue.atoms[index2].pos;
        const gemmi::Vec3& a3 = residue.atoms[index3].pos;
        const gemmi::Vec3& a4 = residue.atoms[index4].pos;
        gemmi::Vec3 b1 = a2 - a1;
        gemmi::Vec3 b2 = a3 - a2;
        gemmi::Vec3 b3 = a4 - a3;
        gemmi::Vec3 n1 = b1.cross(b2);
        gemmi::Vec3 n2 = b2.cross(b3);
        gemmi::Vec3 p  = n1.cross(b2);
        double x = n1.dot(n2);
        double y = p.dot(n2);
        return std::atan2(y, x) * 180.0 / 3.14159265358979323846;
    }
};

static int test_and_fix_ASP_GLU_inline(gemmi::Residue& residue, bool apply_swap_if_found) {
    int iswapped = 0;
    int n_atoms = int(residue.atoms.size());
    std::string residue_name = residue.name;

    if (residue_name == "ASP" || residue_name == "GLU") {
        // gemmi atom names are unpadded (e.g. "OD1" not " OD1")
        std::string test_atom_name = "OD1";
        if (residue_name == "GLU")
            test_atom_name = "OE1";

        std::vector<char> alt_conf_list;
        for (int i = 0; i < n_atoms; i++) {
            if (residue.atoms[i].name == test_atom_name) {
                alt_conf_list.push_back(residue.atoms[i].altloc);
            }
        }

        for (unsigned int ialtconf = 0; ialtconf < alt_conf_list.size(); ialtconf++) {
            local_atom_index_quad quad(-1, -1, -1, -1);

            for (int i = 0; i < n_atoms; i++) {
                if (residue.atoms[i].altloc == alt_conf_list[ialtconf]) {
                    if (residue_name == "ASP") {
                        if (residue.atoms[i].name == "CA") quad.index1 = i;
                        if (residue.atoms[i].name == "CB") quad.index2 = i;
                        if (residue.atoms[i].name == "CG") quad.index3 = i;
                        if (residue.atoms[i].name == "OD1") quad.index4 = i;
                    }
                    if (residue_name == "GLU") {
                        if (residue.atoms[i].name == "CB") quad.index1 = i;
                        if (residue.atoms[i].name == "CG") quad.index2 = i;
                        if (residue.atoms[i].name == "CD") quad.index3 = i;
                        if (residue.atoms[i].name == "OE1") quad.index4 = i;
                    }
                }
            }

            try {
                double torsion = quad.torsion_gemmi(residue);
                bool torsion_is_good = false;
                if (torsion >= -90.0) {
                    if (torsion <= 90.0) {
                        torsion_is_good = true;
                    }
                }
                if (!torsion_is_good) {
                    std::string swap_name_1 = "OD1";
                    std::string swap_name_2 = "OD2";
                    if (residue_name == "GLU") {
                        swap_name_1 = "OE1";
                        swap_name_2 = "OE2";
                    }

                    int at_1_idx = -1;
                    int at_2_idx = -1;
                    for (int i = 0; i < n_atoms; i++) {
                        if (residue.atoms[i].altloc == alt_conf_list[ialtconf]) {
                            if (residue.atoms[i].name == swap_name_1) at_1_idx = i;
                            if (residue.atoms[i].name == swap_name_2) at_2_idx = i;
                        }
                    }
                    if (at_1_idx >= 0 && at_2_idx >= 0) {
                        if (apply_swap_if_found) {
                            std::swap(residue.atoms[at_1_idx].name,
                                      residue.atoms[at_2_idx].name);
                        }
                        iswapped = 1;
                    }
                }
            } catch (const std::runtime_error &rte) {
                std::cout << "WARNING:: missing atoms " << rte.what() << std::endl;
            }
        }
    }

    return iswapped;
}
// --- End inline ASP_GLU ---

// Helper: convert alt_conf string to char for gemmi comparison
static char alt_conf_to_char(const std::string& s) {
    if (s.empty()) return '\0';
    return s[0];
}

std::vector<gemmi::CRA> coot::nomenclature::fix_and_swap_maybe_gemmi(
    gemmi::Structure& st,
    const coot::protein_geometry* geom_p,
    bool apply_swaps)
{
    std::vector<gemmi::CRA> vr;

    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::string residue_name = residue.name;

                gemmi::CRA cra{&chain, &residue, nullptr};

                if ((residue_name == "PHE") || (residue_name == "TYR")) {
                    int isw = coot::nomenclature::test_and_fix_PHE_TYR_nomenclature_errors_gemmi(cra, apply_swaps);
                    if (isw) {
                        if (apply_swaps)
                            std::cout << "INFO:: (result) " << residue_name << " swapped atoms in "
                                      << residue.name << std::endl;
                        vr.push_back(cra);
                    }
                }

                if ((residue_name == "ASP") || (residue_name == "GLU")) {
                    int isw = test_and_fix_ASP_GLU_inline(residue, apply_swaps);
                    if (isw) {
                        if (apply_swaps)
                            std::cout << "INFO:: (result) " << residue_name << " swapped atoms in "
                                      << residue.name << std::endl;
                        vr.push_back(cra);
                    }
                }

                if (residue_name == "THR") {
                    int imol = 0;
                    std::vector<coot::dict_chiral_restraint_t> chiral_restraints =
                        geom_p->get_monomer_chiral_volumes(residue_name, imol);
                    coot::dict_chiral_restraint_t chiral_restraint;
                    for (unsigned int irestr = 0; irestr < chiral_restraints.size(); irestr++) {
                        chiral_restraint = chiral_restraints[irestr];

                        if (chiral_restraint.atom_id_c_4c() == " CB ") {
                            std::vector<std::pair<short int, coot::atom_spec_t>> c =
                                coot::is_inverted_chiral_atom_p_gemmi(
                                    residue,
                                    chiral_restraint.atom_id_c_4c(),
                                    chiral_restraint.atom_id_1_4c(),
                                    chiral_restraint.atom_id_2_4c(),
                                    chiral_restraint.atom_id_3_4c(),
                                    chiral_restraint.volume_sign);

                            for (unsigned int ibad = 0; ibad < c.size(); ibad++) {
                                if (c[ibad].first) {
                                    std::cout << "INFO:: found bad THR chiral atom: "
                                              << chain.name << " "
                                              << residue.seqid.num.value << " "
                                              << residue.seqid.icode << " "
                                              << chiral_restraint.atom_id_c_4c() << " "
                                              << c[ibad].second.alt_conf << std::endl;

                                    std::string alt_conf_bad = c[ibad].second.alt_conf;
                                    char alt_char = alt_conf_to_char(alt_conf_bad);
                                    gemmi::Atom* og1 = nullptr;
                                    gemmi::Atom* cg2 = nullptr;
                                    for (gemmi::Atom& atom : residue.atoms) {
                                        if (atom.name == " OG1" && atom.altloc == alt_char)
                                            og1 = &atom;
                                        if (atom.name == " CG2" && atom.altloc == alt_char)
                                            cg2 = &atom;
                                    }
                                    if (og1 && cg2) {
                                        if (apply_swaps) {
                                            og1->name = " CG2";
                                            cg2->name = " OG1";
                                            std::cout << "        CG2 and OG1 atoms swapped\n";
                                            std::cout << "INFO:: swapped atoms in "
                                                      << residue.name << " "
                                                      << residue.name << std::endl;
                                        }
                                        vr.push_back(cra);
                                    } else {
                                        std::cout << "ERROR:: Bizarre missing atom scenario "
                                                  << "in fix_nomenclature_errors\n";
                                    }
                                }
                            }
                        }
                    }
                }

                if ((residue_name == "LEU") || (residue_name == "VAL")) {
                    int volume_sign = -1;
                    coot::dict_chiral_restraint_t synthetic_restraint;

                    if (residue_name == "VAL")
                        synthetic_restraint =
                            coot::dict_chiral_restraint_t(residue_name,
                                                          " CB ", " CA ", " CG1", " CG2",
                                                          volume_sign);
                    if (residue_name == "LEU")
                        synthetic_restraint =
                            coot::dict_chiral_restraint_t(residue_name,
                                                          " CG ", " CB ", " CD1", " CD2",
                                                          volume_sign);

                    std::vector<std::pair<short int, coot::atom_spec_t>> c =
                        coot::is_inverted_chiral_atom_p_gemmi(
                            residue,
                            synthetic_restraint.atom_id_c_4c(),
                            synthetic_restraint.atom_id_1_4c(),
                            synthetic_restraint.atom_id_2_4c(),
                            synthetic_restraint.atom_id_3_4c(),
                            synthetic_restraint.volume_sign);

                    for (unsigned int ibad = 0; ibad < c.size(); ibad++) {
                        if (c[ibad].first) {
                            std::cout << "INFO:: found bad " << residue_name
                                      << " chiral atom: "
                                      << chain.name << " "
                                      << residue.seqid.num.value << " "
                                      << residue.seqid.icode << " "
                                      << synthetic_restraint.atom_id_c_4c() << " "
                                      << c[ibad].second.alt_conf << std::endl;

                            std::string alt_conf_bad = c[ibad].second.alt_conf;
                            char alt_char = alt_conf_to_char(alt_conf_bad);
                            std::string target_atom_1 = " CG1";
                            std::string target_atom_2 = " CG2";
                            if (residue_name == "LEU") {
                                target_atom_1 = " CD1";
                                target_atom_2 = " CD2";
                            }

                            gemmi::Atom* cg1 = nullptr;
                            gemmi::Atom* cg2 = nullptr;
                            for (gemmi::Atom& atom : residue.atoms) {
                                if (atom.name == target_atom_1 && atom.altloc == alt_char)
                                    cg1 = &atom;
                                if (atom.name == target_atom_2 && atom.altloc == alt_char)
                                    cg2 = &atom;
                            }
                            if (cg1 && cg2) {
                                if (apply_swaps) {
                                    cg1->name = target_atom_2;
                                    cg2->name = target_atom_1;
                                    std::cout << "        " << target_atom_1 << " and "
                                              << target_atom_2 << " atoms swapped\n";
                                    std::cout << "INFO:: swapped atoms in "
                                              << residue.name << " "
                                              << residue.name << std::endl;
                                }
                                vr.push_back(cra);
                            } else {
                                std::cout << "ERROR:: Bizarre missing atom scenario "
                                          << "in fix_nomenclature_errors, residue type:"
                                          << residue_name << "\n";
                            }
                        }
                    }
                }
            }
        }
    }

    if (apply_swaps)
        std::cout << "INFO:: " << vr.size() << " residues had their atoms swapped\n";

    return vr;
}