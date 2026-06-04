#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_main_chain_p/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_hydrophobic_atom/gemmi/function.hh"

// Bond colour type constants (mirror coot::bond_colour_type enum values)
static const int BC_MOLECULE = 1;
static const int BC_CHAIN = 2;
static const int BC_HYDROPHOBIC = 3;
static const int BC_SEC_STRUCT = 4;
static const int BC_ATOM_TYPE = 5;
static const int BC_CHAIN_C_ONLY = 6;
static const int BC_DISULFIDE = 7;
static const int BC_OCCUPANCY = 8;
static const int BC_B_FACTOR = 9;
static const int BC_RAINBOW = 10;
static const int BC_USER_DEFINED = 11;
static const int BC_CHAIN_GOODSELL = 12;

// Bond colour constants (mirror coot::bond_colours enum values)
static const int BC_COLOR_CARBON = 0;
static const int BC_COLOR_YELLOW = 1;
static const int BC_COLOR_RED = 2;
static const int BC_COLOR_BLUE = 3;
static const int BC_COLOR_GREEN = 4;
static const int BC_COLOR_MAGENTA = 5;
static const int BC_COLOR_GREY = 6;
static const int BC_COLOR_ORANGE = 7;
static const int BC_COLOR_CYAN = 8;
static const int BC_COLOR_HYDROGEN_GREY = 9;
static const int BC_COLOR_DARK_BROWN = 10;
static const int BC_COLOR_DARK_GREEN = 11;
static const int BC_COLOR_DARK_ORANGE = 12;
static const int BC_COLOR_DEUTERIUM_PINK = 13;
static const int BC_COLOR_DARK_VIOLET = 14;
static const int BC_COLOR_VIOLET = 15;
static const int BC_COLOR_BORON_PINK = 1002;

// Inline translation of coot::util::is_standard_amino_acid_name
static inline bool is_standard_amino_acid_name(const std::string &residue_name) {
    static const char* names[] = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
        "HIS","ILE","LEU","LYS","MET","MSE","PHE","PRO",
        "SER","THR","TRP","TYR","VAL"
    };
    for (auto* n : names) {
        if (residue_name == n) return true;
    }
    return false;
}

// Forward declaration for atom colour map (opaque pointer, no coot include needed)
struct my_atom_colour_map_t {
    virtual int index_for_chain(const std::string&) = 0;
    virtual ~my_atom_colour_map_t() = default;
};

// gemmi port of Bond_lines_container::atom_colour
// bond_colour_type: use BC_* constants defined above
// atom_colour_map_p: nullable pointer for chain-based colouring; nullptr in most cases
int atom_colour_gemmi(gemmi::CRA cra,
                      int bond_colour_type,
                      int udd_user_defined_atom_colour_index_handle,
                      my_atom_colour_map_t *atom_colour_map_p,
                      float b_factor_scale = 1.0f);

inline int atom_colour_gemmi(gemmi::CRA cra,
                             int bond_colour_type,
                             int udd_user_defined_atom_colour_index_handle,
                             my_atom_colour_map_t *atom_colour_map_p,
                             float b_factor_scale) {

    if (!cra.atom || !cra.residue) return 20;

    const gemmi::Atom& at = *cra.atom;
    const gemmi::Residue& res = *cra.residue;
    const gemmi::Element& elem = at.element;
    std::string element = elem.name();

    int col = 0;

    // UDD data — gemmi has no user-defined data, skip entirely.

    if (bond_colour_type == BC_MOLECULE) return col;

    if (bond_colour_type == BC_CHAIN) {
        if (atom_colour_map_p) {
            std::string chain_id;
            if (cra.chain) chain_id = cra.chain->name;
            col = atom_colour_map_p->index_for_chain(chain_id);
        }
    } else {

        if (bond_colour_type == BC_CHAIN_GOODSELL) {
            if (atom_colour_map_p) {
                std::string chain_id;
                if (cra.chain) chain_id = cra.chain->name;
                int col_idx = atom_colour_map_p->index_for_chain(chain_id);
                col = 2 * col_idx;
                if (element != "C")
                    col += 1;
            }
        } else {

            if (bond_colour_type == BC_HYDROPHOBIC) {
                std::string res_name = res.name;
                if (is_standard_amino_acid_name(res_name)) {
                    std::string atom_name = at.name;
                    if (coot::is_main_chain_p_gemmi(at, res)) {
                        col = 50;
                    } else {
                        if (coot::is_hydrophobic_atom_gemmi(res_name, atom_name))
                            col = 1;
                        else
                            col = 2;
                    }
                }
            }

            if (bond_colour_type == BC_SEC_STRUCT) {
                // gemmi::Residue has no SSE field; default to col=0 (None)
                col = 0;
            } else {

                if (bond_colour_type == BC_ATOM_TYPE) {
                    if (element == "C") return BC_COLOR_CARBON;
                    if (element == "N") return BC_COLOR_BLUE;
                    if (element == "O") return BC_COLOR_RED;
                    if (element == "S") return BC_COLOR_YELLOW;
                    if (elem.is_hydrogen()) {
                        if (element == "D") return BC_COLOR_DEUTERIUM_PINK;
                        return BC_COLOR_HYDROGEN_GREY;
                    }
                    if (element == "P") return BC_COLOR_ORANGE;
                    if (element == "F") return BC_COLOR_GREEN;
                    if (element == "Cl") return BC_COLOR_GREEN;
                    if (element == "Br") return BC_COLOR_DARK_BROWN;
                    if (element == "I") return BC_COLOR_DARK_VIOLET;
                    if (element == "B") return BC_COLOR_BORON_PINK;
                    if (element == "Mg" || element == "Be" || element == "Ca" || element == "Sr" || element == "Ba") return BC_COLOR_DARK_GREEN;
                    if (element == "Fe") return BC_COLOR_DARK_ORANGE;
                    if (element == "Li" || element == "Na" || element == "K" || element == "Rb" || element == "Cs" || element == "Fr") return BC_COLOR_VIOLET;
                    return BC_COLOR_GREY;
                } else {

                    if (bond_colour_type == BC_CHAIN_C_ONLY) {
                        if (element == "C") {
                            if (atom_colour_map_p) {
                                std::string chain_id;
                                if (cra.chain) chain_id = cra.chain->name;
                                return atom_colour_map_p->index_for_chain(chain_id);
                            }
                            return col;
                        }
                        if (element == "N") return BC_COLOR_BLUE;
                        if (element == "O") return BC_COLOR_RED;
                        if (element == "S") return BC_COLOR_YELLOW;
                        if (element == "P") return BC_COLOR_ORANGE;
                        if (elem.is_hydrogen()) {
                            if (element == "D") return BC_COLOR_DEUTERIUM_PINK;
                            return BC_COLOR_HYDROGEN_GREY;
                        }
                        if (element == "F") return BC_COLOR_GREEN;
                        if (element == "Cl") return BC_COLOR_GREEN;
                        if (element == "Br") return BC_COLOR_DARK_BROWN;
                        if (element == "I") return BC_COLOR_DARK_VIOLET;
                        if (element == "B") return BC_COLOR_BORON_PINK;
                        if (element == "Mg" || element == "Be" || element == "Ca" || element == "Sr" || element == "Ba") return BC_COLOR_DARK_GREEN;
                        if (element == "Fe") return BC_COLOR_DARK_ORANGE;
                        if (element == "Li" || element == "Na" || element == "K" || element == "Rb" || element == "Cs" || element == "Fr") return BC_COLOR_VIOLET;
                        return BC_COLOR_GREY;
                    } else {

                        if (bond_colour_type == BC_DISULFIDE) return BC_COLOR_YELLOW;

                        if (bond_colour_type == BC_OCCUPANCY) {
                            if (at.occ > 0.95) return BC_COLOR_BLUE;
                            if (at.occ < 0.05) return BC_COLOR_RED;
                            if (at.occ > 0.7) return BC_COLOR_CYAN;
                            if (at.occ > 0.45) return BC_COLOR_GREEN;
                            if (at.occ > 0.25) return BC_COLOR_YELLOW;
                            return BC_COLOR_ORANGE;
                        }

                        if (bond_colour_type == BC_B_FACTOR) {
                            float scaled_b = at.b_iso * b_factor_scale;
                            float f = scaled_b / 100.0;
                            if (f > 0.999f) f = 0.999f;
                            if (f < 0.0f)   f = 0.0f;
                            col = static_cast<int>(f * 45);
                        } else {
                            col = 20;
                        }
                    }
                }
            }
        }
    }

    return col;
}