#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <vector>
#include <string>

namespace coot {
class protein_geometry;
}

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

namespace coot {
namespace nomenclature {

namespace detail {

int test_and_fix_PHE_TYR_nomenclature_errors_gemmi(gemmi::CRA cra, bool apply_swap) {
    if (!cra.residue) return 0;
    std::string rname = cra.residue->name;
    if (rname != "PHE" && rname != "TYR") return 0;

    int iswapped = 0;

    // Collect unique alt confs from CD1 atoms
    std::vector<char> alt_confs;
    for (auto &atom : cra.residue->atoms) {
        if (atom.name == "CD1") {
            bool found = false;
            for (char a : alt_confs) if (a == atom.altloc) { found = true; break; }
            if (!found) alt_confs.push_back(atom.altloc);
        }
    }
    if (alt_confs.empty()) return 0;

    for (char alt : alt_confs) {
        gemmi::Atom *CA=nullptr,*CB=nullptr,*CG=nullptr,*CD1=nullptr,*CD2=nullptr;
        for (auto &atom : cra.residue->atoms) {
            if (atom.altloc == alt) {
                if (atom.name=="CA")  CA  = &atom;
                if (atom.name=="CB")  CB  = &atom;
                if (atom.name=="CG")  CG  = &atom;
                if (atom.name=="CD1") CD1 = &atom;
                if (atom.name=="CD2") CD2 = &atom;
            }
        }
        if (!CA||!CB||!CG) {
            for (auto &atom : cra.residue->atoms) {
                if (atom.altloc=='\0') {
                    if (!CA  && atom.name=="CA")  CA  = &atom;
                    if (!CB  && atom.name=="CB")  CB  = &atom;
                    if (!CG  && atom.name=="CG")  CG  = &atom;
                    if (!CD2 && atom.name=="CD2") CD2 = &atom;
                }
            }
        }
        if (!CA||!CB||!CG||!CD1) continue;
        if (!CD2) continue;
        double tors = gemmi::calculate_dihedral(CA->pos, CB->pos, CG->pos, CD1->pos)*180.0/M_PI;
        if (tors < -90.0 || tors > 90.0) {
            double tors2 = gemmi::calculate_dihedral(CA->pos, CB->pos, CG->pos, CD2->pos)*180.0/M_PI;
            if (tors2 >= -90.0 && tors2 <= 90.0) {
                gemmi::Atom *CE1=nullptr,*CE2=nullptr,*HD1=nullptr,*HD2=nullptr,*HE1=nullptr,*HE2=nullptr;
                for (auto &atom : cra.residue->atoms) {
                    if (atom.altloc == alt) {
                        if (atom.name=="CE1") CE1 = &atom;
                        if (atom.name=="CE2") CE2 = &atom;
                        if (atom.name=="HD1") HD1 = &atom;
                        if (atom.name=="HD2") HD2 = &atom;
                        if (atom.name=="HE1") HE1 = &atom;
                        if (atom.name=="HE2") HE2 = &atom;
                    }
                }
                if (apply_swap) {
                    CD1->pos = CD2->pos;
                    CD1->name = CD2->name;
                }
                if (CE1 && CE2) {
                    if (apply_swap) {
                        CE1->pos = CE2->pos;
                        CE1->name = CE2->name;
                    }
                    iswapped++;
                } else {
                    iswapped++;
                }
                if (HD1 && HD2 && apply_swap) HD1->name = HD2->name;
                if (HE1 && HE2 && apply_swap) HE1->name = HE2->name;
            }
        }
    }
    return iswapped;
}

int test_and_fix_ASP_GLU_nomenclature_errors_gemmi(gemmi::Residue &residue, bool apply_swap) {
    int iswapped = 0;
    std::string rname = residue.name;
    if (rname != "ASP" && rname != "GLU") return 0;
    std::string test_atom = (rname=="GLU") ? "OE1" : "OD1";
    int n = int(residue.atoms.size());

    std::vector<char> alt_confs;
    for (int i=0;i<n;i++) {
        if (residue.atoms[i].name == test_atom) {
            bool found=false;
            for (char a:alt_confs) if (a==residue.atoms[i].altloc) {found=true;break;}
            if (!found) alt_confs.push_back(residue.atoms[i].altloc);
        }
    }

    for (char alt : alt_confs) {
        int idx_a=-1,idx_b=-1,idx_c=-1,idx_d=-1;
        if (rname=="ASP") {
            for (int i=0;i<n;i++) {
                if (residue.atoms[i].altloc==alt) {
                    if (residue.atoms[i].name=="CA")  idx_a=i;
                    if (residue.atoms[i].name=="CB")  idx_b=i;
                    if (residue.atoms[i].name=="CG")  idx_c=i;
                    if (residue.atoms[i].name=="OD1") idx_d=i;
                }
            }
        } else {
            for (int i=0;i<n;i++) {
                if (residue.atoms[i].altloc==alt) {
                    if (residue.atoms[i].name=="CB")  idx_a=i;
                    if (residue.atoms[i].name=="CG")  idx_b=i;
                    if (residue.atoms[i].name=="CD")  idx_c=i;
                    if (residue.atoms[i].name=="OE1") idx_d=i;
                }
            }
        }
        if (idx_a<0||idx_b<0||idx_c<0||idx_d<0) continue;
        double tors = gemmi::calculate_dihedral(
            residue.atoms[idx_a].pos, residue.atoms[idx_b].pos,
            residue.atoms[idx_c].pos, residue.atoms[idx_d].pos
        )*180.0/M_PI;

        std::string swap1 = (rname=="GLU")?"OE1":"OD1";
        std::string swap2 = (rname=="GLU")?"OE2":"OD2";
        if (tors < -90.0 || tors > 90.0) {
            int i1=-1,i2=-1;
            for (int i=0;i<n;i++) {
                if (residue.atoms[i].altloc==alt) {
                    if (residue.atoms[i].name==swap1) i1=i;
                    if (residue.atoms[i].name==swap2) i2=i;
                }
            }
            if (i1>=0 && i2>=0) {
                if (apply_swap) {
                    residue.atoms[i1].name = residue.atoms[i2].name;
                    residue.atoms[i1].pos = residue.atoms[i2].pos;
                }
                iswapped = 1;
            }
        }
    }
    return iswapped;
}

} // namespace detail

std::vector<gemmi::CRA> list_gemmi(gemmi::Model &m, const coot::protein_geometry *geom_p) {
    (void)geom_p; // unused in list (fix_and_swap_maybe with apply_swap=false doesn't use geom)
    std::vector<gemmi::CRA> swapped_residues;
    for (auto &chain : m.chains) {
        for (auto &residue : chain.residues) {
            gemmi::CRA cra{&chain, &residue, nullptr};
            if (detail::test_and_fix_PHE_TYR_nomenclature_errors_gemmi(cra, false))
                swapped_residues.push_back(cra);
            if (detail::test_and_fix_ASP_GLU_nomenclature_errors_gemmi(residue, false))
                swapped_residues.push_back(cra);
        }
    }
    return swapped_residues;
}

} // namespace nomenclature
} // namespace coot