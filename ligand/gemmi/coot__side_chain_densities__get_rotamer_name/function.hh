#pragma once
#include <string>
#include <cmath>
#include <cctype>
#include <gemmi/model.hpp>

namespace coot {
namespace side_chain_densities {

// Helper: find atom by name in residue (gemmi stores names like "N","CA","CB","CG1")
static const gemmi::Atom* find_atom(const gemmi::Residue& res, const std::string& name) {
    for (const gemmi::Atom& a : res.atoms) {
        if (a.name == name) return &a;
    }
    return nullptr;
}

// Helper: compute torsion angle from 4 atom positions (degrees)
static double compute_torsion(const gemmi::Vec3& p1, const gemmi::Vec3& p2,
                               const gemmi::Vec3& p3, const gemmi::Vec3& p4) {
    gemmi::Vec3 b1 = p2 - p1;
    gemmi::Vec3 b2 = p3 - p2;
    gemmi::Vec3 b3 = p4 - p3;

    gemmi::Vec3 n1 = b1.cross(b2);
    gemmi::Vec3 n2 = b2.cross(b3);

    double len1 = std::sqrt(n1.length_sq());
    double len2 = std::sqrt(n2.length_sq());

    if (len1 < 1e-6 || len2 < 1e-6) return 0.0;

    gemmi::Vec3 u1 = n1 * (1.0 / len1);
    gemmi::Vec3 u2 = n2 * (1.0 / len2);

    double x = -(u1.dot(u2));
    gemmi::Vec3 b2n = b2 * (1.0 / std::sqrt(b2.length_sq()));
    double y = u1.cross(u2).dot(b2n);

    return std::atan2(y, x) * 180.0 / M_PI;
}

// Helper: strip whitespace
static std::string remove_whitespace(const std::string& s) {
    std::string result;
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
            result += c;
        }
    }
    return result;
}

// Determine rotamer name for residues with chi angles
static std::string determine_rotamer(const gemmi::Residue& res) {
    std::string rname = res.name;

    if (rname == "GLY" || rname == "ALA") return "none";

    auto classify_chi = [](double chi) -> char {
        while (chi > 180.0) chi -= 360.0;
        while (chi < -180.0) chi += 360.0;
        // t: roughly [-120, 120], p: roughly (120, 180], m: roughly [-180, -120)
        if (chi > -120.0 && chi <= 120.0) return 't';
        if (chi > 120.0) return 'p';
        return 'm';
    };

    // VAL: chi1=N-CA-CB-CG1, chi2=CA-CB-CG1-CG2
    if (rname == "VAL") {
        const gemmi::Atom* n   = find_atom(res, "N");
        const gemmi::Atom* ca  = find_atom(res, "CA");
        const gemmi::Atom* cb  = find_atom(res, "CB");
        const gemmi::Atom* cg1 = find_atom(res, "CG1");
        const gemmi::Atom* cg2 = find_atom(res, "CG2");

        if (!n || !ca || !cb || !cg1 || !cg2) return "none";

        double chi1 = compute_torsion(n->pos, ca->pos, cb->pos, cg1->pos);
        double chi2 = compute_torsion(ca->pos, cb->pos, cg1->pos, cg2->pos);

        // VAL reports chi1 only (single letter)
        return std::string(1, classify_chi(chi1));
    }

    // ILE: chi1=N-CA-CB-CG1, chi2=CA-CB-CG1-CG2, chi3=CB-CG1-CG2-CD1
    if (rname == "ILE") {
        const gemmi::Atom* n   = find_atom(res, "N");
        const gemmi::Atom* ca  = find_atom(res, "CA");
        const gemmi::Atom* cb  = find_atom(res, "CB");
        const gemmi::Atom* cg1 = find_atom(res, "CG1");
        const gemmi::Atom* cg2 = find_atom(res, "CG2");
        const gemmi::Atom* cd1 = find_atom(res, "CD1");

        if (!n || !ca || !cb || !cg1 || !cg2 || !cd1) return "none";

        double chi1 = compute_torsion(n->pos, ca->pos, cb->pos, cg1->pos);
        double chi2 = compute_torsion(ca->pos, cb->pos, cg1->pos, cg2->pos);
        double chi3 = compute_torsion(cb->pos, cg1->pos, cg2->pos, cd1->pos);

        // ILE reports chi1+chi2 (two letters)
        return std::string(1, classify_chi(chi1)) + classify_chi(chi2);
    }

    return "none";
}

// Ported from MMDB version that took mmdb::Residue*
inline std::string get_rotamer_name_gemmi(const gemmi::CRA& cra) {
    if (!cra.residue) return "none";
    std::string rname = determine_rotamer(*cra.residue);
    return remove_whitespace(rname);
}

} // namespace side_chain_densities
} // namespace coot