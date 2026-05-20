#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/calculate.hpp>
#include <coot/geometry/residue-and-atom-specs.hh>

namespace coot {

// Check if a chiral atom has inverted stereochemistry
// Returns vector of (ibad, atom_spec) pairs - ibad=1 if inverted, 0 if OK
inline std::vector<std::pair<short int, atom_spec_t> >
is_inverted_chiral_atom_p_gemmi(const gemmi::Residue& res,
                                const std::string& chiral_atom_name,
                                const std::string& atom_id_1,
                                const std::string& atom_id_2,
                                const std::string& atom_id_3,
                                int volume_sign) {
    using namespace gemmi;
    
    std::vector<std::pair<short int, atom_spec_t> > v;
    short int ibad = 0;
    
    // Find the chiral atom (the "c" atom - center)
    const Atom* chiral_atom = nullptr;
    const Atom* atom1 = nullptr;
    const Atom* atom2 = nullptr;
    const Atom* atom3 = nullptr;
    
    for (const auto& a : res.atoms) {
        if (a.name == chiral_atom_name) chiral_atom = &a;
        if (a.name == atom_id_1) atom1 = &a;
        if (a.name == atom_id_2) atom2 = &a;
        if (a.name == atom_id_3) atom3 = &a;
    }
    
    // If any atom is missing, return empty vector
    if (!chiral_atom || !atom1 || !atom2 || !atom3) {
        return v;
    }
    
    // Get positions
    Position centre(chiral_atom->pos);
    Position a1(atom1->pos);
    Position a2(atom2->pos);
    Position a3(atom3->pos);
    
    // Calculate chiral volume: (a1-centre) . ((a2-centre) x (a3-centre))
    Vec3 a = a1 - centre;
    Vec3 b = a2 - centre;
    Vec3 c = a3 - centre;
    double cv = a.dot(b.cross(c));
    
    // Create atom_spec_t for the chiral atom
    // Note: gemmi::Residue does NOT have a parent chain pointer.
    // We use an empty string as chain name since we don't have parent context.
    atom_spec_t chiral_atom_spec("", res.seqid.num.value,
                                 std::string(1, res.seqid.icode),
                                 chiral_atom->name,
                                 std::string(1, chiral_atom->altloc));
    
    // Check volume sign
    if (cv * volume_sign < 0) {
        ibad = 1;
    } else {
        ibad = 0;
    }
    
    v.push_back(std::pair<short int, atom_spec_t>(ibad, chiral_atom_spec));
    
    return v;
}

} // namespace coot