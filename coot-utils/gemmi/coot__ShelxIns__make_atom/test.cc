#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to create a clipper::Cell from parameters
static clipper::Cell make_cell(const std::vector<double> &params) {
    clipper::Cell_descr cell_d(params[0], params[1], params[2],
                               clipper::Util::d2rad(params[3]),
                               clipper::Util::d2rad(params[4]),
                               clipper::Util::d2rad(params[5]));
    clipper::Cell cell;
    cell.init(cell_d);
    return cell;
}

// Helper function to create a shelx card
static coot::shelx_card_info_t make_card(const std::string &card_line,
                                         const std::vector<std::string> &words) {
    coot::shelx_card_info_t card;
    card.card = card_line;
    card.words = words;
    return card;
}

// Helper function to get padded element name like MMDB's GetElementName()
static std::string get_padded_element_name(const gemmi::Atom& atom) {
    // MMDB's GetElementName returns 2-char padded string
    // gemmi's element.name() returns unpadded like "C", "N", "O"
    // Pad it to match MMDB format
    std::string elem = atom.element.name();
    if (elem.length() == 1) {
        return " " + elem;
    } else {
        return elem;
    }
}

TEST(OracleTest, ShelxInsMakeAtom) {
    // Setup: Create a ShelxIns object
    coot::ShelxIns shelx_ins;
    
    // Populate sfac vector - this is required for make_atom_element to work
    // SFAC typically contains element symbols or types (forrtl: 1-based indexing)
    shelx_ins.sfac.push_back("C");  // index 0 -> sfac_index=1
    shelx_ins.sfac.push_back("N");  // index 1 -> sfac_index=2
    shelx_ins.sfac.push_back("O");  // index 2 -> sfac_index=3
    
    // Create a cell
    std::vector<double> cell_params = {10.0, 10.0, 10.0, 90.0, 90.0, 90.0};
    clipper::Cell cell = make_cell(cell_params);
    
    std::vector<gemmi::Atom*> atom_vector;
    
    // case: valid atom with isotropic U (positive value)
    {
        coot::shelx_card_info_t card = make_card("C1 1 0.5 0.5 0.5 1.0 0.02",
            {"C1", "1", "0.5", "0.5", "0.5", "1.0", "0.02"});
        
        std::string altconf = "";
        int udd_afix_handle = -1;
        int udd_non_riding_atom_flag_handle = -1;
        int udd_riding_atom_negative_u_value_handle = -1;
        bool have_udd_atoms = false;
        int current_afix = -1;
        
        gemmi::Atom *at = shelx_ins.make_atom_gemmi(card, altconf, udd_afix_handle,
                                                    udd_non_riding_atom_flag_handle,
                                                    udd_riding_atom_negative_u_value_handle,
                                                    have_udd_atoms, current_afix,
                                                    cell, atom_vector);
        
        EXPECT_TRUE(at != nullptr);
        
        if (at) {
            EXPECT_EQ(std::string(at->name), " C1 ");
            EXPECT_EQ(get_padded_element_name(*at), " C");
            EXPECT_NEAR(at->pos.x, 0.5, 1e-4);
            EXPECT_NEAR(at->pos.y, 0.5, 1e-4);
            EXPECT_NEAR(at->pos.z, 0.5, 1e-4);
            EXPECT_NEAR(at->occ, 1.0, 1e-4);
            EXPECT_NEAR(at->b_iso, 1.57914, 1e-4);
            EXPECT_EQ(at->flag, 1);  // Basic flag for coordinates
            atom_vector.push_back(at);
        }
    }
    
    // case: atom with negative U (riding atom case)
    {
        // First atom - non-riding
        coot::shelx_card_info_t card1 = make_card("C1 1 0.5 0.5 0.5 1.0 0.02",
            {"C1", "1", "0.5", "0.5", "0.5", "1.0", "0.02"});
        
        std::string altconf = "";
        int udd_afix_handle = -1;
        int udd_non_riding_atom_flag_handle = -1;
        int udd_riding_atom_negative_u_value_handle = -1;
        bool have_udd_atoms = false;
        int current_afix = -1;
        
        gemmi::Atom *at1 = shelx_ins.make_atom_gemmi(card1, altconf, udd_afix_handle,
                                                     udd_non_riding_atom_flag_handle,
                                                     udd_riding_atom_negative_u_value_handle,
                                                     have_udd_atoms, current_afix,
                                                     cell, atom_vector);
        if (at1) atom_vector.push_back(at1);
        
        // Second atom - riding with negative U
        coot::shelx_card_info_t card2 = make_card("C2 1 0.6 0.6 0.6 1.0 -0.5",
            {"C2", "1", "0.6", "0.6", "0.6", "1.0", "-0.5"});
        
        gemmi::Atom *at2 = shelx_ins.make_atom_gemmi(card2, altconf, udd_afix_handle,
                                                     udd_non_riding_atom_flag_handle,
                                                     udd_riding_atom_negative_u_value_handle,
                                                     have_udd_atoms, current_afix,
                                                     cell, atom_vector);
        
        EXPECT_TRUE(at2 != nullptr);
        
        if (at2) {
            EXPECT_EQ(std::string(at2->name), " C2 ");
            // tempFactor should be positive after taking absolute value for riding atoms
            EXPECT_NEAR(at2->b_iso, 0.789568, 1e-4);
            EXPECT_NEAR(at2->pos.x, 0.6, 1e-4);
            EXPECT_NEAR(at2->pos.y, 0.6, 1e-4);
            EXPECT_NEAR(at2->pos.z, 0.6, 1e-4);
            atom_vector.push_back(at2);
        }
    }
    
    // case: atom with anisotropic U
    {
        coot::shelx_card_info_t card_aniso = make_card("N1 2 0.3 0.3 0.3 1.0 0.02 0.01 0.015 0.005 0.003 0.004",
            {"N1", "2", "0.3", "0.3", "0.3", "1.0", "0.02", "0.01", "0.015", "0.005", "0.003", "0.004"});
        
        std::string altconf = "";
        int udd_afix_handle = -1;
        int udd_non_riding_atom_flag_handle = -1;
        int udd_riding_atom_negative_u_value_handle = -1;
        bool have_udd_atoms = false;
        int current_afix = -1;
        
        gemmi::Atom *at = shelx_ins.make_atom_gemmi(card_aniso, altconf, udd_afix_handle,
                                                    udd_non_riding_atom_flag_handle,
                                                    udd_riding_atom_negative_u_value_handle,
                                                    have_udd_atoms, current_afix,
                                                    cell, atom_vector);
        
        EXPECT_TRUE(at != nullptr);
        
        if (at) {
            EXPECT_EQ(std::string(at->name), " N1 ");
            EXPECT_EQ(get_padded_element_name(*at), " N");
            EXPECT_NEAR(at->aniso.u11, 0.02, 1e-4);
            EXPECT_NEAR(at->aniso.u22, 0.01, 1e-4);
            EXPECT_NEAR(at->aniso.u33, 0.015, 1e-4);
            EXPECT_NEAR(at->aniso.u23, 0.005, 1e-4);
            EXPECT_NEAR(at->aniso.u13, 0.003, 1e-4);
            EXPECT_NEAR(at->aniso.u12, 0.004, 1e-4);
            // Original WhatIsSet = 263 (0x00000100 | 0x00000002 | 0x00000001 | 0x00000040 | 0x00000080 | 0x00000200)
            // In gemmi, flag is a char - we'll just verify it's set for anisotropic
            EXPECT_EQ(at->flag, 1);  // Basic flag
            atom_vector.push_back(at);
        }
    }
    
    // case: invalid element (sfac_index out of range)
    {
        coot::shelx_card_info_t card_invalid = make_card("XXX 99 0.5 0.5 0.5 1.0 0.02",
            {"XXX", "99", "0.5", "0.5", "0.5", "1.0", "0.02"});
        
        std::string altconf = "";
        int udd_afix_handle = -1;
        int udd_non_riding_atom_flag_handle = -1;
        int udd_riding_atom_negative_u_value_handle = -1;
        bool have_udd_atoms = false;
        int current_afix = -1;
        
        gemmi::Atom *at = shelx_ins.make_atom_gemmi(card_invalid, altconf, udd_afix_handle,
                                                    udd_non_riding_atom_flag_handle,
                                                    udd_riding_atom_negative_u_value_handle,
                                                    have_udd_atoms, current_afix,
                                                    cell, atom_vector);
        
        EXPECT_TRUE(at == nullptr);
    }
    
    // Cleanup
    for (auto *at : atom_vector) {
        delete at;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}