#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_methyl_Hs) {
    // Case 1: ALA residue with methyl group (CB)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Residue* ala_res = nullptr;

        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "ALA") {
                    ala_res = &res;
                    break;
                }
            }
            if (ala_res) break;
        }

        ASSERT_NE(ala_res, nullptr);

        // Verify the residue has the atoms needed for torsion_info_t (C, CA, CB)
        const gemmi::Atom* at_C  = nullptr;
        const gemmi::Atom* at_CA = nullptr;
        const gemmi::Atom* at_CB = nullptr;

        for (const auto& atom : ala_res->atoms) {
            std::string trimmed = atom.name;
            size_t start = trimmed.find_first_not_of(' ');
            if (start != std::string::npos) {
                size_t end = trimmed.find_last_not_of(' ');
                trimmed = trimmed.substr(start, end - start + 1);
            }
            if (trimmed == "C")  at_C  = &atom;
            if (trimmed == "CA") at_CA = &atom;
            if (trimmed == "CB") at_CB = &atom;
        }

        ASSERT_NE(at_C, nullptr);
        ASSERT_NE(at_CA, nullptr);
        ASSERT_NE(at_CB, nullptr);

        int seqnum = ala_res->seqid.num.value;

        // Build cid using the model chain info.
        gemmi::Chain* ala_chain = nullptr;
        for (auto& ch : model.chains) {
            for (auto& r : ch.residues) {
                if (&r == ala_res) {
                    ala_chain = &ch;
                    break;
                }
            }
            if (ala_chain) break;
        }
        ASSERT_NE(ala_chain, nullptr);

        std::string cid2 = std::string("//") + ala_chain->name + "/" + std::to_string(seqnum);
        EXPECT_EQ(cid2, "//A/21");

        int atoms_before = static_cast<int>(ala_res->atoms.size());
        EXPECT_EQ(atoms_before, 10);

        // torsion_info_t: reference atoms C, CA, CB to define the methyl orientation
        coot::reduce::torsion_info_t ti(" C  ", " CA ", " CB ", 1.09, 109.5, 180.0);

        EXPECT_EQ(ti.at_name_1, " C  ");
        EXPECT_EQ(ti.at_name_2, " CA ");
        EXPECT_EQ(ti.at_name_3, " CB ");
        EXPECT_NEAR(ti.bond_length, 1.09, 1e-4);
        EXPECT_NEAR(ti.angle_deg, 109.5, 1e-4);
        EXPECT_NEAR(ti.torsion_deg, 180.0, 1e-4);

        // Call add_methyl_Hs
        coot::reduce::atoms_with_spinnable_Hs spinables;
        coot::reduce::add_methyl_Hs_gemmi("HB1", "HB2", "HB3", ti, *ala_res, spinables);

        int atoms_after = static_cast<int>(ala_res->atoms.size());
        EXPECT_EQ(atoms_after, 13);
        EXPECT_EQ((atoms_after - atoms_before), 3);

        // Verify hydrogens were added
        const gemmi::Atom* hb1 = nullptr;
        const gemmi::Atom* hb2 = nullptr;
        const gemmi::Atom* hb3 = nullptr;

        for (const auto& atom : ala_res->atoms) {
            std::string trimmed = atom.name;
            size_t start = trimmed.find_first_not_of(' ');
            if (start != std::string::npos) {
                size_t end = trimmed.find_last_not_of(' ');
                trimmed = trimmed.substr(start, end - start + 1);
            }
            if (trimmed == "HB1") hb1 = &atom;
            if (trimmed == "HB2") hb2 = &atom;
            if (trimmed == "HB3") hb3 = &atom;
        }
        EXPECT_NE(hb1, nullptr);
        EXPECT_NE(hb2, nullptr);
        EXPECT_NE(hb3, nullptr);
    }

    // Case 2: GLY residue without CB — should do nothing
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model& model = st.models[0];
        gemmi::Residue* gly_res = nullptr;

        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.name == "GLY") {
                    gly_res = &res;
                    break;
                }
            }
            if (gly_res) break;
        }

        if (gly_res) {
            int gly_seqnum = gly_res->seqid.num.value;

            gemmi::Chain* gly_chain = nullptr;
            for (auto& ch : model.chains) {
                for (auto& r : ch.residues) {
                    if (&r == gly_res) {
                        gly_chain = &ch;
                        break;
                    }
                }
                if (gly_chain) break;
            }
            ASSERT_NE(gly_chain, nullptr);

            std::string gly_cid = std::string("//") + gly_chain->name + "/" + std::to_string(gly_seqnum);
            EXPECT_EQ(gly_cid, "//A/11");

            int gly_before = static_cast<int>(gly_res->atoms.size());
            EXPECT_EQ(gly_before, 7);

            // GLY doesn't have CB, so guard should fail and no atoms added
            coot::reduce::torsion_info_t ti_gly(" C  ", " CA ", " CB ", 1.09, 109.5, 180.0);
            coot::reduce::atoms_with_spinnable_Hs spinables;
            coot::reduce::add_methyl_Hs_gemmi("HB1", "HB2", "HB3", ti_gly, *gly_res, spinables);

            int gly_after = static_cast<int>(gly_res->atoms.size());
            EXPECT_EQ(gly_after, 7);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}