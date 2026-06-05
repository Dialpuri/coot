#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <sstream>

TEST(OracleTest, dipole) {
    // Load the structure
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // --- Case 1: LZA ligand (no partial charges -> exception) ---
    {
        // Find residue A/1299/LZA
        gemmi::CRA target_cra{nullptr, nullptr, nullptr};
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 1299 && res.name == "LZA") {
                        target_cra = gemmi::CRA{&chain, &res, nullptr};
                        break;
                    }
                }
            }
        }
        ASSERT_NE(target_cra.residue, nullptr);

        EXPECT_EQ(target_cra.chain->name, "A");
        EXPECT_EQ(target_cra.residue->seqid.num.value, 1299);
        EXPECT_EQ(target_cra.residue->name, "LZA");

        // Count atoms
        EXPECT_EQ(target_cra.residue->atoms.size(), 25u);

        // Suppress "no partial charge for" messages
        std::ostringstream nullstream;
        std::streambuf *old_cout = std::cout.rdbuf(nullstream.rdbuf());

        coot::dictionary_residue_restraints_t rest("LZA", 1299);
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> dict_res_pairs;
        dict_res_pairs.push_back(std::make_pair(rest, target_cra));

        coot::dipole::dipole_state_t dip;
        bool threw = false;
        std::string errmsg;
        try {
            coot::dipole::init_gemmi(dip, dict_res_pairs);
        } catch (const std::exception &e) {
            threw = true;
            errmsg = e.what();
        } catch (...) {
            threw = true;
            errmsg = "unknown exception";
        }

        // Restore stdout
        std::cout.rdbuf(old_cout);

        EXPECT_EQ(dip.dipole_is_good_flag, false);
        EXPECT_TRUE(threw);
        EXPECT_EQ(errmsg, "Dipole is not good for 1 residue A 1299 LZA, ");
    }

    // --- Case 2: Empty vector (triggers "No atoms" exception) ---
    {
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> empty_pairs;
        EXPECT_EQ(empty_pairs.size(), 0u);

        coot::dipole::dipole_state_t dip2;
        bool threw2 = false;
        std::string errmsg2;
        try {
            coot::dipole::init_gemmi(dip2, empty_pairs);
        } catch (const std::exception &e) {
            threw2 = true;
            errmsg2 = e.what();
        } catch (...) {
            threw2 = true;
            errmsg2 = "unknown exception";
        }

        EXPECT_TRUE(threw2);
        EXPECT_EQ(errmsg2, "No atoms in 0 residues ");
    }

    // --- Case 3: Multiple residues (LYS at 9 + ILE at 10) ---
    {
        gemmi::CRA cra1{nullptr, nullptr, nullptr};
        gemmi::CRA cra2{nullptr, nullptr, nullptr};

        // Find residues A/9/LYS and A/10/ILE
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == 9 && res.name == "LYS") {
                        cra1 = gemmi::CRA{&chain, &res, nullptr};
                    }
                    if (res.seqid.num.value == 10 && res.name == "ILE") {
                        cra2 = gemmi::CRA{&chain, &res, nullptr};
                    }
                }
            }
        }

        ASSERT_NE(cra1.residue, nullptr);
        ASSERT_NE(cra2.residue, nullptr);

        EXPECT_EQ(cra1.chain->name, "A");
        EXPECT_EQ(cra1.residue->seqid.num.value, 9);
        EXPECT_EQ(cra1.residue->name, "LYS");

        EXPECT_EQ(cra2.chain->name, "A");
        EXPECT_EQ(cra2.residue->seqid.num.value, 10);
        EXPECT_EQ(cra2.residue->name, "ILE");

        // Suppress noisy messages
        std::ostringstream nullstream3;
        std::streambuf *old_cout3 = std::cout.rdbuf(nullstream3.rdbuf());

        coot::dictionary_residue_restraints_t rest1("LYS", 9);
        coot::dictionary_residue_restraints_t rest2("ILE", 10);

        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> pairs3;
        pairs3.push_back(std::make_pair(rest1, cra1));
        pairs3.push_back(std::make_pair(rest2, cra2));

        coot::dipole::dipole_state_t dip3;
        bool threw3 = false;
        std::string errmsg3;
        try {
            coot::dipole::init_gemmi(dip3, pairs3);
        } catch (const std::exception &e) {
            threw3 = true;
            errmsg3 = e.what();
        } catch (...) {
            threw3 = true;
            errmsg3 = "unknown exception";
        }

        std::cout.rdbuf(old_cout3);

        EXPECT_EQ(dip3.dipole_is_good_flag, false);
        EXPECT_TRUE(threw3);
        EXPECT_EQ(errmsg3, "Dipole is not good for 2 residues A 9 LYS, A 10 ILE, ");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}