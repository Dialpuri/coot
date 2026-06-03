#include "function.hh"

namespace coot {
namespace backrub {

void setup_this_and_prev_next_ca_positions_gemmi(
    const gemmi::Residue* this_res,
    const gemmi::Residue* prev_res,
    const gemmi::Residue* next_res,
    const std::string& alt_conf,
    clipper::Coord_orth& ca_this,
    clipper::Coord_orth& ca_prev,
    clipper::Coord_orth& ca_next
) {
    if (!this_res) {
        throw std::runtime_error(" Null this residue ");
    }
    if (!prev_res) {
        throw std::runtime_error(" Null previous residue ");
    }
    if (!next_res) {
        throw std::runtime_error(" Null next residue ");
    }

    // Convert alt_conf string to char for gemmi comparison
    // In gemmi, altloc is a char; empty string means '\0' (no alt loc)
    char alt_char = alt_conf.empty() ? '\0' : alt_conf[0];

    // Find CA in this_res
    bool found = false;
    for (const gemmi::Atom& atom : this_res->atoms) {
        if (atom.name == "CA" && atom.altloc == alt_char) {
            found = true;
            ca_this = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            break;
        }
    }
    if (!found) {
        std::string mess(" CA atom of this residue in alt conf \"");
        mess += alt_conf;
        mess += "\" not found";
        throw std::runtime_error(mess);
    }

    // Find CA in prev_res
    found = false;
    for (const gemmi::Atom& atom : prev_res->atoms) {
        if (atom.name == "CA" && atom.altloc == alt_char) {
            found = true;
            ca_prev = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            break;
        }
    }
    if (!found) {
        std::string mess(" CA atom of previous residue in alt conf \"");
        mess += alt_conf;
        mess += "\" not found";
        throw std::runtime_error(mess);
    }

    // Find CA in next_res
    found = false;
    for (const gemmi::Atom& atom : next_res->atoms) {
        if (atom.name == "CA" && atom.altloc == alt_char) {
            found = true;
            ca_next = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
            break;
        }
    }
    if (!found) {
        std::string mess(" CA atom of next residue in alt conf \"");
        mess += alt_conf;
        mess += "\" not found";
        throw std::runtime_error(mess);
    }
}

} // namespace backrub
} // namespace coot