#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>

// Gemmi replacement for mmdb::LinkR
struct gemmi_linkr_t {
    std::string chainID1;
    int seqNum1 = 0;
    std::string atName1;
    std::string insCode1;
    std::string chainID2;
    int seqNum2 = 0;
    std::string atName2;
    std::string insCode2;
};

// Find atom in gemmi::Model matching the given LinkR-side info.
// Returns pointer to atom, or nullptr on failure.
// Sets *out_cra on success.
inline const gemmi::Atom* find_atom_from_linkr_templ_gemmi(
    gemmi::Model *model_p,
    int udd_atom_index_handle,
    int udd_user_defined_atom_colour_index_handle,
    const std::string &chain_id,
    int seq_num,
    const std::string &atom_name,
    const std::string &ins_code,
    gemmi::CRA *out_cra)
{
    (void)udd_atom_index_handle;
    (void)udd_user_defined_atom_colour_index_handle;
    if (!model_p) return nullptr;

    // Normalize insertion codes: empty string (MMDB) and ' ' (gemmi) both mean "none"
    auto norm = [](const std::string& s) { return s.empty() ? std::string(" ") : s; };
    std::string norm_ic = norm(ins_code);

    for (gemmi::Chain &chain : model_p->chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue &res : chain.residues) {
            if (res.seqid.num.value != seq_num) continue;
            std::string res_ic(1, res.seqid.icode);
            if (norm_ic != norm(res_ic)) continue;
            for (const gemmi::Atom &at : res.atoms) {
                // Compare trimmed atom names (gemmi stores padded names like " CA ")
                std::string at_name_trimmed = at.name;
                at_name_trimmed.erase(at_name_trimmed.find_last_not_of(" \t") + 1);
                at_name_trimmed.erase(0, at_name_trimmed.find_first_not_of(" \t"));
                std::string query_trimmed = atom_name;
                query_trimmed.erase(query_trimmed.find_last_not_of(" \t") + 1);
                query_trimmed.erase(0, query_trimmed.find_first_not_of(" \t"));
                if (at_name_trimmed == query_trimmed) {
                    if (out_cra) {
                        *out_cra = gemmi::CRA{&chain, &res, const_cast<gemmi::Atom*>(&at)};
                    }
                    return &at;
                }
            }
        }
    }
    return nullptr;
}

// Bond_line struct (simplified from coot)
struct Bond_line {
    float x1, y1, z1;
    float x2, y2, z2;
    int colour;
    int thickness;
};

// Forward declare Bond_lines_container
class Bond_lines_container {
public:
    std::vector<std::vector<Bond_line> > bonds;

    Bond_lines_container() {}
    explicit Bond_lines_container(int n_colour_types) {
        bonds.resize(n_colour_types);
    }

    // Gemmi port of add_link_bond
    void add_link_bond_gemmi(gemmi::Model *model_p,
                             int udd_atom_index_handle,
                             int udd_user_defined_atom_colour_index_handle,
                             int atom_colour_type,
                             const gemmi_linkr_t *linkr);
};

// Inline definition of add_link_bond_gemmi
// NOTE: This preserves the original copy-vs-reference bug from add_link_bond_templ
// where bonds_bonds_coloured is a local copy that gets modified but never written back.
inline void Bond_lines_container::add_link_bond_gemmi(
    gemmi::Model *model_p,
    int udd_atom_index_handle,
    int udd_user_defined_atom_colour_index_handle,
    int atom_colour_type,
    const gemmi_linkr_t *linkr)
{
    if (!model_p || !linkr) return;

    gemmi::CRA cra1{};
    const gemmi::Atom *at1 = find_atom_from_linkr_templ_gemmi(
        model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle,
        linkr->chainID1, linkr->seqNum1, linkr->atName1, linkr->insCode1,
        &cra1);
    if (!at1) return;

    gemmi::CRA cra2{};
    const gemmi::Atom *at2 = find_atom_from_linkr_templ_gemmi(
        model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle,
        linkr->chainID2, linkr->seqNum2, linkr->atName2, linkr->insCode2,
        &cra2);
    if (!at2) return;

    // Ensure bonds vector has enough space
    if ((int)bonds.size() <= atom_colour_type) {
        bonds.resize(atom_colour_type + 1);
    }

    // Compute colour using atom colour type
    int col = 0;
    if (atom_colour_type == 5) {
        if (at1->element == gemmi::Element("C")) col = 0;
        else if (at1->element == gemmi::Element("N")) col = 3;
        else if (at1->element == gemmi::Element("O")) col = 2;
        else if (at1->element == gemmi::Element("S")) col = 1;
        else col = 6;
    } else {
        col = atom_colour_type;
    }

    // BUG PRESERVATION: copy the vector, modify the copy, discard it.
    // This mirrors the original add_link_bond_templ which does:
    //   std::vector<Bond_line> bonds_bonds_coloured = this->bonds_bonds[colour_type];
    // and then push_back onto the local copy.
    std::vector<Bond_line> bonds_copy = bonds[atom_colour_type];
    Bond_line bl;
    bl.x1 = (float)at1->pos.x;
    bl.y1 = (float)at1->pos.y;
    bl.z1 = (float)at1->pos.z;
    bl.x2 = (float)at2->pos.x;
    bl.y2 = (float)at2->pos.y;
    bl.z2 = (float)at2->pos.z;
    bl.colour = col;
    bl.thickness = 2;
    bonds_copy.push_back(bl);
    // bonds_copy is discarded here — the actual bonds[atom_colour_type] is untouched
}

// Helper: count total bonds
inline size_t count_bonds_gemmi(const Bond_lines_container &bonds) {
    size_t total = 0;
    for (const auto &bl : bonds.bonds) {
        total += bl.size();
    }
    return total;
}