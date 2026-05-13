#include "function.hh"

namespace coot { namespace util {

std::map<std::string, clipper::RTop_orth>
get_ori_to_this_res_gemmi(const gemmi::Residue& residue) {
    std::map<std::string, clipper::RTop_orth> orientations;

    std::map<char, std::vector<const gemmi::Atom*> > atoms;
    
    // Collect atoms grouped by altloc
    for (const gemmi::Atom& atom : residue.atoms) {
        char alt_conf = atom.altloc;
        if (alt_conf == '\0') alt_conf = ' ';  // gemmi uses ' ' for no altloc
        if (atoms.find(alt_conf) == atoms.end()) {
            atoms[alt_conf].resize(3, nullptr);
        }
    }
    
    // Map mainchain atoms to indices
    for (const gemmi::Atom& atom : residue.atoms) {
        char alt_conf = atom.altloc;
        if (alt_conf == '\0') alt_conf = ' ';  // gemmi uses ' ' for no altloc
        
        std::string atom_name = atom.name;
        int name_index = -1;
        
        // Check common formats for atom names
        if (atom_name == " N  " || atom_name == "N" || atom_name == " N ") name_index = 0;
        if (atom_name == " CA " || atom_name == "CA" || atom_name == " CA") name_index = 1;
        if (atom_name == " C  " || atom_name == "C" || atom_name == " C ") name_index = 2;
        
        if (name_index != -1) {
            atoms[alt_conf][name_index] = &atom;
        }
    }
    
    // Build orientation matrices for altlocs with all three atoms
    for (auto& pair : atoms) {
        char alt_conf = pair.first;
        std::vector<const gemmi::Atom*>& atom_vec = pair.second;
        
        if (atom_vec[0] && atom_vec[1] && atom_vec[2]) {
            const gemmi::Atom* n_atom = atom_vec[0];
            const gemmi::Atom* ca_atom = atom_vec[1];
            const gemmi::Atom* c_atom = atom_vec[2];
            
            clipper::Coord_orth n(n_atom->pos.x, n_atom->pos.y, n_atom->pos.z);
            clipper::Coord_orth ca(ca_atom->pos.x, ca_atom->pos.y, ca_atom->pos.z);
            clipper::Coord_orth c(c_atom->pos.x, c_atom->pos.y, c_atom->pos.z);
            
            clipper::Coord_orth can_unit = clipper::Coord_orth((n - ca).unit());
            clipper::Coord_orth cac_unit = clipper::Coord_orth((c - ca).unit());
            
            clipper::Coord_orth bisector((can_unit + cac_unit).unit());
            clipper::Coord_orth diff_unit((can_unit - cac_unit).unit());
            
            clipper::Coord_orth cross_prod(clipper::Coord_orth::cross(diff_unit, bisector));
            clipper::Coord_orth cpu(cross_prod.unit());
            
            clipper::Mat33<double> m(bisector.x(),  bisector.y(),  bisector.z(),
                                     cpu.x(),       cpu.y(),       cpu.z(),
                                     diff_unit.x(), diff_unit.y(), diff_unit.z());
            
            clipper::RTop_orth rtop(m.transpose(), ca);
            
            // Convert gemmi's ' ' (space) to "" to match MMDB's empty string
            std::string alt_str;
            if (alt_conf == ' ') {
                alt_str = "";
            } else {
                alt_str = std::string(1, alt_conf);
            }
            orientations[alt_str] = rtop;
        }
    }
    
    return orientations;
}

}}