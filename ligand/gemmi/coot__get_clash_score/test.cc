#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <streambuf>
#include "function.hh"

// Redirect cout to suppress debug prints from get_clash_score
static std::ofstream null_stream("/dev/null");
static std::streambuf *old_cout_buf = nullptr;

void suppress_cout() { old_cout_buf = std::cout.rdbuf(null_stream.rdbuf()); }
void restore_cout() { std::cout.rdbuf(old_cout_buf); }

TEST(OracleTest, get_clash_score_gemmi) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Build rotamer_atoms from ALL atoms (matching minimol built from manager with all models)
    // Build atom_selection with only HOH atoms (matching original model 1 which had just the water)
    std::vector<gemmi::CRA> rotamer_atoms;
    std::vector<gemmi::CRA> atom_selection;

    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                gemmi::CRA cra{&chain, &res, &atom};
                rotamer_atoms.push_back(cra);
                if (res.name == "HOH") {
                    atom_selection.push_back(cra);
                }
            }
        }
    }

    // Verify inputs match oracle: atom_selection should have 1 HOH atom
    EXPECT_EQ(atom_selection.size(), 1u);

    // Case 1: water_interaction_mode=0 (waters excluded from score, collected instead)
    {
        suppress_cout();
        auto p = coot::get_clash_score_gemmi(rotamer_atoms, atom_selection, 0);
        restore_cout();

        EXPECT_NEAR(p.first, 0.0f, 1e-4);
        EXPECT_EQ(p.second.size(), 1u);
    }

    // Case 2: water_interaction_mode=1 (waters included in score)
    {
        suppress_cout();
        auto p2 = coot::get_clash_score_gemmi(rotamer_atoms, atom_selection, 1);
        restore_cout();

        EXPECT_NEAR(p2.first, 7.56915f, 1e-4);
        EXPECT_EQ(p2.second.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}