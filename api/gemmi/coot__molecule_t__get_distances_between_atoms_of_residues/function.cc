#include "function.hh"
#include <cmath>
#include <string>

namespace coot {

namespace {

static std::string pad_atom_name(const std::string& name) {
    if (name.size() >= 4) return name.substr(0, 4);
    if (name.size() == 1) return " " + name + "  ";
    // 2 or 3 chars: left-justified, padded with trailing spaces
    std::string result = name;
    result.resize(4, ' ');
    return result;
}

static std::string extract_chain(const std::string& cid) {
    // Format: "//CHAIN/RESNO"
    size_t start = 2; // skip leading "//"
    size_t end = cid.find('/', start);
    return cid.substr(start, end - start);
}

} // anonymous namespace

std::vector<atom_distance_t> get_distances_between_atoms_of_residues_gemmi(
    const std::string &cid_res_1,
    const std::string &cid_res_2,
    float dist_max,
    gemmi::Model &model) 
{
    std::vector<atom_distance_t> v;
    
    gemmi::Residue *residue_1 = cid_to_residue_gemmi(cid_res_1, model);
    gemmi::Residue *residue_2 = cid_to_residue_gemmi(cid_res_2, model);
    
    if (residue_1 && residue_2) {
        std::string chain_1 = extract_chain(cid_res_1);
        std::string chain_2 = extract_chain(cid_res_2);
        int res_no_1 = residue_1->seqid.num.value;
        int res_no_2 = residue_2->seqid.num.value;
        
        for (const gemmi::Atom &at_1 : residue_1->atoms) {
            for (const gemmi::Atom &at_2 : residue_2->atoms) {
                double dd = (at_2.pos.x - at_1.pos.x) * (at_2.pos.x - at_1.pos.x) +
                            (at_2.pos.y - at_1.pos.y) * (at_2.pos.y - at_1.pos.y) +
                            (at_2.pos.z - at_1.pos.z) * (at_2.pos.z - at_1.pos.z);
                double d = std::sqrt(dd);
                if (d < dist_max) {
                    atom_spec_t spec_1;
                    spec_1.chain_id  = chain_1;
                    spec_1.res_no    = res_no_1;
                    spec_1.atom_name = pad_atom_name(at_1.name);
                    
                    atom_spec_t spec_2;
                    spec_2.chain_id  = chain_2;
                    spec_2.res_no    = res_no_2;
                    spec_2.atom_name = pad_atom_name(at_2.name);
                    
                    atom_distance_t ad(spec_1, spec_2, static_cast<float>(d));
                    v.push_back(ad);
                }
            }
        }
    }
    
    return v;
}

} // namespace coot