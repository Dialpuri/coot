#include "function.hh"
#include <cmath>

namespace {

std::string trim_atom_name(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

} // anonymous namespace

coot::restraints_container_t_gemmi::restraints_container_t_gemmi(
    const std::vector<std::pair<bool, gemmi::Residue*>>& residues,
    gemmi::Structure* st)
  : residues_(residues), st_(st) {
    for (const auto& [unused, res] : residues_) {
        if (res) {
            for (const gemmi::Atom& atom : res->atoms) {
                atom_.push_back(atom.pos);
            }
        }
    }
}

int coot::restraints_container_t_gemmi::get_asc_index_old_gemmi(
    const std::string& at_name, int resno, const std::string& chain_id) const {

    // Find target atom in gemmi structure
    gemmi::Vec3 target_pos;
    bool found = false;

    for (const gemmi::Model& model : st_->models) {
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (const gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != resno) continue;
                for (const gemmi::Atom& atom : res.atoms) {
                    if (trim_atom_name(atom.name) == at_name) {
                        target_pos = atom.pos;
                        found = true;
                        break;
                    }
                }
            }
            if (found) break;
        }
        if (found) break;
    }

    if (!found) return -1;

    // Linear search through atom_ for matching coordinates (with tolerance)
    double tol_sq = 0.001;
    for (size_t i = 0; i < atom_.size(); i++) {
        double dx = static_cast<double>(atom_[i].x) - target_pos.x;
        double dy = static_cast<double>(atom_[i].y) - target_pos.y;
        double dz = static_cast<double>(atom_[i].z) - target_pos.z;
        if (dx * dx + dy * dy + dz * dz < tol_sq) {
            return static_cast<int>(i);
        }
    }

    return -1;
}
