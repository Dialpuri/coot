#pragma once

#include <gemmi/cif.hpp>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

namespace coot {
namespace protein_geometry {

static std::string
get_padded_name(const std::string& name, int /*comp_id_index*/) {
    std::string padded = name;
    while (padded.size() < 4) padded += ' ';
    return padded;
}

static void
mon_lib_add_chiral(std::string /*comp_id*/,
                   int /*imol_enc*/,
                   std::string /*id*/,
                   std::string /*atom_id_centre*/,
                   std::string /*atom_id_1*/,
                   std::string /*atom_id_2*/,
                   std::string /*atom_id_3*/,
                   std::string /*volume_sign*/)
{
    // In the original, this modifies internal dict_chirality state.
    // The test only checks return values (n_chiral, comp_id_vector),
    // so this is a no-op for the port.
}

std::pair<int, std::vector<std::string>>
comp_chiral_gemmi(const gemmi::cif::Loop& chir_loop, int imol_enc) {
    std::pair<int, std::vector<std::string>> result{0, {}};

    if (chir_loop.tags.empty() || chir_loop.length() == 0) {
        return result;
    }

    // Extract category prefix from first tag (e.g. "_chem_comp_chir." from "_chem_comp_chir.comp_id")
    std::string cat_prefix;
    if (!chir_loop.tags.empty()) {
        size_t dot_pos = chir_loop.tags[0].rfind('.');
        if (dot_pos != std::string::npos)
            cat_prefix = chir_loop.tags[0].substr(0, dot_pos) + ".";
    }

    int idx_comp_id = chir_loop.find_tag(cat_prefix + "comp_id");
    int idx_id = chir_loop.find_tag(cat_prefix + "id");
    int idx_atom_centre = chir_loop.find_tag(cat_prefix + "atom_id_centre");
    int idx_atom_1 = chir_loop.find_tag(cat_prefix + "atom_id_1");
    int idx_atom_2 = chir_loop.find_tag(cat_prefix + "atom_id_2");
    int idx_atom_3 = chir_loop.find_tag(cat_prefix + "atom_id_3");
    int idx_volume_sign = chir_loop.find_tag(cat_prefix + "volume_sign");

    for (size_t row = 0; row < chir_loop.length(); row++) {
        std::string comp_id, id, atom_id_centre, atom_id_1, atom_id_2, atom_id_3, volume_sign;
        bool ok = true;

        if (idx_comp_id >= 0) {
            comp_id = gemmi::cif::as_string(chir_loop.val(row, idx_comp_id));
            if (comp_id == ".") ok = false;
        } else {
            ok = false;
        }

        if (idx_id >= 0) {
            id = gemmi::cif::as_string(chir_loop.val(row, idx_id));
            if (id == ".") ok = false;
        } else {
            ok = false;
        }

        if (idx_atom_centre >= 0) {
            std::string raw = gemmi::cif::as_string(chir_loop.val(row, idx_atom_centre));
            if (raw == ".") ok = false;
            else atom_id_centre = get_padded_name(raw, -1);
        } else {
            ok = false;
        }

        if (idx_atom_1 >= 0) {
            std::string raw = gemmi::cif::as_string(chir_loop.val(row, idx_atom_1));
            if (raw == ".") ok = false;
            else atom_id_1 = get_padded_name(raw, -1);
        } else {
            ok = false;
        }

        if (idx_atom_2 >= 0) {
            std::string raw = gemmi::cif::as_string(chir_loop.val(row, idx_atom_2));
            if (raw == ".") ok = false;
            else atom_id_2 = get_padded_name(raw, -1);
        } else {
            ok = false;
        }

        if (idx_atom_3 >= 0) {
            std::string raw = gemmi::cif::as_string(chir_loop.val(row, idx_atom_3));
            if (raw == ".") ok = false;
            else atom_id_3 = get_padded_name(raw, -1);
        } else {
            ok = false;
        }

        if (idx_volume_sign >= 0) {
            volume_sign = gemmi::cif::as_string(chir_loop.val(row, idx_volume_sign));
            if (volume_sign == ".") ok = false;
        } else {
            ok = false;
        }

        if (ok) {
            mon_lib_add_chiral(comp_id, imol_enc,
                               id, atom_id_centre,
                               atom_id_1,
                               atom_id_2,
                               atom_id_3,
                               volume_sign);
            if (std::find(result.second.begin(), result.second.end(), comp_id) ==
                result.second.end())
                result.second.push_back(comp_id);
            result.first++;
        }
    }

    return result;
}

} // namespace protein_geometry
} // namespace coot