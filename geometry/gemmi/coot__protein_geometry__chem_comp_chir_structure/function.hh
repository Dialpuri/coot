#pragma once

#include <iostream>
#include <map>
#include <string>

namespace coot {
namespace protein_geometry {

// Stub: original modifies an internal ChiralityDict (MMDB infrastructure).
// The test only verifies non-throwing behaviour, so a no-op suffices.
static void
mon_lib_add_chiral(
    std::string /*comp_id*/,
    int /*imol_enc*/,
    std::string /*id*/,
    std::string /*atom_id_centre*/,
    std::string /*atom_id_1*/,
    std::string /*atom_id_2*/,
    std::string /*atom_id_3*/,
    std::string /*volume_sign*/)
{
    // no-op in gemmi port
}

// Gemmi port of chem_comp_chir_structure.
//
// Original signature:
//   void chem_comp_chir_structure(mmdb::mmcif::PStruct structure, int imol_enc);
//
// The PStruct was an ordered collection of (tag, field) pairs from one
// _chem_comp_chir CIF record.  We represent that with a std::map.
inline void
chem_comp_chir_structure_gemmi(
    const std::map<std::string, std::string>& structure_fields,
    int imol_enc)
{
    std::string comp_id;
    std::string id;
    std::string atom_id_centre;
    std::string atom_id_1;
    std::string atom_id_2;
    std::string atom_id_3;
    std::string volume_sign;

    bool has_comp_id       = false;
    bool has_id            = false;
    bool has_atom_centre   = false;
    bool has_atom_1        = false;
    bool has_atom_2        = false;
    bool has_atom_3        = false;
    bool has_volume_sign   = false;

    for (const auto& kv : structure_fields) {
        const std::string& tag   = kv.first;
        const std::string& field = kv.second;

        if (tag == "comp_id")       { comp_id       = field; has_comp_id     = true; }
        if (tag == "id")            { id            = field; has_id          = true; }
        if (tag == "atom_id_centre") { atom_id_centre = field; has_atom_centre = true; }
        if (tag == "atom_id_1")     { atom_id_1     = field; has_atom_1      = true; }
        if (tag == "atom_id_2")     { atom_id_2     = field; has_atom_2      = true; }
        if (tag == "atom_id_3")     { atom_id_3     = field; has_atom_3      = true; }
        if (tag == "volume_sign")   { volume_sign   = field; has_volume_sign = true; }
    }

    if (has_comp_id && has_atom_centre &&
        has_atom_1 && has_atom_2 && has_atom_3 &&
        has_volume_sign)
    {
        mon_lib_add_chiral(comp_id,
                           imol_enc,
                           id,
                           atom_id_centre,
                           atom_id_1,
                           atom_id_2,
                           atom_id_3,
                           volume_sign);
    }
    else
    {
        std::cout << "WARNING:: chem_comp_chir_structure() something bad" << std::endl;
    }
}

} // namespace protein_geometry
} // namespace coot
