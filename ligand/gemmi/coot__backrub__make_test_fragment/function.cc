#include "function.hh"

namespace coot {
namespace backrub {

// Trim atom name to match MMDB-style space-padded names
static std::string trim_atom(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

coot::minimol::fragment
make_test_fragment_gemmi(
    const gemmi::Residue* this_residue,
    const gemmi::Residue* prev_residue,
    const gemmi::Residue* next_residue,
    const std::string& chain_id,
    clipper::Coord_orth ca_prev,
    clipper::Coord_orth ca_next,
    double rotation_angle)
{
    coot::minimol::fragment f(chain_id);

    // Atoms to include in prev/next residues (trimmed)
    std::vector<std::string> prev_res_atoms;
    prev_res_atoms.push_back("C");
    prev_res_atoms.push_back("O");
    std::vector<std::string> next_res_atoms;
    next_res_atoms.push_back("N");
    next_res_atoms.push_back("H");

    // Build residues using the static from_gemmi method
    coot::minimol::residue this_res =
        coot::minimol::residue::from_gemmi(this_residue, std::vector<std::string>());
    coot::minimol::residue prev_res =
        coot::minimol::residue::from_gemmi(prev_residue, prev_res_atoms);
    coot::minimol::residue next_res =
        coot::minimol::residue::from_gemmi(next_residue, next_res_atoms);

    // addresidue() fails when adding residues with the same residue number
    try {
        f.addresidue(prev_res, 0);
        f.addresidue(this_res, 0);
        f.addresidue(next_res, 0);
    }
    catch (const std::runtime_error &rte) {
        std::cout << "ERROR:: make_test_fragment() " << rte.what() << std::endl;
    }

    // Now rotate fragment around the ca_prev -> ca_next vector
    clipper::Coord_orth dir = ca_next - ca_prev;

    for (int ires = f.min_res_no(); ires <= f.max_residue_number(); ires++) {
        for (unsigned int iat = 0; iat < f[ires].n_atoms(); iat++) {
            clipper::Coord_orth pt(f[ires][iat].pos);
            double ra = M_PI * rotation_angle / 180.0;
            clipper::Coord_orth pt_new =
                coot::util::rotate_around_vector(dir, pt, ca_prev, ra);
            f[ires][iat].pos = pt_new;
        }
    }

    // rotate_individual_peptides_back_best — no gemmi port, stub: do nothing

    if (f.n_filled_residues() != 3) {
        std::string mess = "  Failed to get 3 residues with atoms in test fragment. Got ";
        mess += coot::util::int_to_string(f.n_filled_residues());
        throw std::runtime_error(mess);
    }

    return f;
}

} // namespace backrub

// Implement from_gemmi
namespace minimol {

residue residue::from_gemmi(const gemmi::Residue* res, const std::vector<std::string>& include_atoms) {
    if (!res) {
        return residue(0);
    }
    residue r(res->seqid.num.value, res->name);
    r.ins_code = res->seqid.icode;

    for (const auto& atom : res->atoms) {
        std::string aname = coot::backrub::trim_atom(atom.name);
        // If include_atoms is non-empty, only include matching atoms
        if (!include_atoms.empty()) {
            bool found = false;
            for (const auto& inc : include_atoms) {
                if (aname == inc) {
                    found = true;
                    break;
                }
            }
            if (!found) continue;
        }
        r.atoms.push_back(
            coot::minimol::atom(
                aname,
                atom.pos.x, atom.pos.y, atom.pos.z
            )
        );
    }
    return r;
}

} // namespace minimol
} // namespace coot