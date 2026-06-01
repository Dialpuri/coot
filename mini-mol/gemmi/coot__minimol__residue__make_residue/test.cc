#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

// Helper to normalize gemmi insertion code ' ' back to MMDB-style ""
static std::string gemmi_icode_to_mmdb(gemmi::Residue const& r) {
    if (r.seqid.icode == ' ') return "";
    return std::string(1, r.seqid.icode);
}

// Helper to pad single-char element names to PDB-style 2-char format
static std::string gemmi_elem_to_mmdb(gemmi::Element e) {
    std::string s = e.name();
    if (s.size() < 2) s = " " + s;
    return s;
}

TEST(OracleTest, make_residue) {
    // Case 1: Typical case — residue with atoms (ALA with backbone atoms)
    {
        coot::minimol::residue res(42, "ALA");
        res.ins_code = "";

        // Add backbone atoms: N, CA, C
        coot::minimol::atom n_atom(" N  ", " N", clipper::Coord_orth(10.0, 20.0, 30.0), "", 0.5f, 20.0f);
        coot::minimol::atom ca_atom(" CA ", " C", clipper::Coord_orth(11.0, 21.0, 31.0), "", 0.5f, 20.0f);
        coot::minimol::atom c_atom(" C  ", " C", clipper::Coord_orth(12.0, 22.0, 32.0), "", 0.5f, 20.0f);

        res.atoms.push_back(n_atom);
        res.atoms.push_back(ca_atom);
        res.atoms.push_back(c_atom);

        auto gemmi_res = coot::minimol::make_residue_gemmi(res);

        EXPECT_NE(gemmi_res, nullptr);
        if (gemmi_res) {
            EXPECT_EQ(gemmi_res->name, "ALA");
            EXPECT_EQ(gemmi_res->seqid.num.value, 42);
            EXPECT_EQ(gemmi_icode_to_mmdb(*gemmi_res), "");
            EXPECT_EQ(gemmi_res->atoms.size(), 3u);

            // Print first atom details
            if (gemmi_res->atoms.size() > 0) {
                gemmi::Atom const* first_atom = &gemmi_res->atoms[0];
                EXPECT_NE(first_atom, nullptr);
                EXPECT_EQ(first_atom->name, " N  ");
                EXPECT_EQ(gemmi_elem_to_mmdb(first_atom->element), " N");
                EXPECT_NEAR(first_atom->pos.x, 10.0, 1e-4);
                EXPECT_NEAR(first_atom->pos.y, 20.0, 1e-4);
                EXPECT_NEAR(first_atom->pos.z, 30.0, 1e-4);
            }
        }
    }

    // Case 2: Edge case — empty residue (no atoms), should return nullptr
    {
        coot::minimol::residue res(99, "GLY");
        // No atoms added — atoms vector is empty

        auto gemmi_res = coot::minimol::make_residue_gemmi(res);

        EXPECT_EQ(gemmi_res, nullptr);
    }

    // Case 3: Residue with insertion code
    {
        coot::minimol::residue res(5, "SER");
        res.ins_code = "A";

        coot::minimol::atom n_atom(" N  ", " N", clipper::Coord_orth(5.0, 6.0, 7.0), "", 1.0f, 15.0f);
        coot::minimol::atom ca_atom(" CA ", " C", clipper::Coord_orth(5.5, 6.5, 7.5), "", 1.0f, 15.0f);

        res.atoms.push_back(n_atom);
        res.atoms.push_back(ca_atom);

        auto gemmi_res = coot::minimol::make_residue_gemmi(res);

        EXPECT_NE(gemmi_res, nullptr);
        if (gemmi_res) {
            EXPECT_EQ(gemmi_res->name, "SER");
            EXPECT_EQ(gemmi_res->seqid.num.value, 5);
            EXPECT_EQ(gemmi_icode_to_mmdb(*gemmi_res), "A");
            EXPECT_EQ(gemmi_res->atoms.size(), 2u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}