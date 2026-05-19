#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, TorsionTests) {
    // Load structure using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Find residue A/10
    const gemmi::Model& model = st.models[0];
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& c : model.chains) {
        if (c.name == "A") {
            chain_a = &c;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr) << "could not find chain A";
    
    const gemmi::Residue* res = nullptr;
    for (const auto& r : chain_a->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr) << "could not find residue 10";

    // Case 1: valid torsion with all 4 atoms present (BACKBONE torsion: N-CA-C-O)
    {
        std::vector<std::string> atom_names = {"N", "CA", "C", "O"};
        coot::atom_name_quad quad(atom_names[0], atom_names[1], atom_names[2], atom_names[3]);
        double torsion_val = quad.torsion_gemmi(*res);
        // The oracle printed 174.318, actual value is 174.31781459093162
        // Using a tolerance of 0.001 to accommodate the printed precision
        EXPECT_NEAR(torsion_val, 174.318, 0.001);
    }

    // Case 2: invalid torsion — one atom missing (missing O)
    {
        std::vector<std::string> atom_names = {"N", "CA", "C", "XX"};
        coot::atom_name_quad quad(atom_names[0], atom_names[1], atom_names[2], atom_names[3]);
        double torsion_val = quad.torsion_gemmi(*res);
        EXPECT_NEAR(torsion_val, -999.9, 1e-4);
    }

    // Case 3: another valid residue with backbone atoms (residue 20)
    const gemmi::Residue* res2 = nullptr;
    for (const auto& r : chain_a->residues) {
        if (r.seqid.num.value == 20) {
            res2 = &r;
            break;
        }
    }
    if (res2) {
        std::vector<std::string> atom_names = {"N", "CA", "C", "O"};
        coot::atom_name_quad quad(atom_names[0], atom_names[1], atom_names[2], atom_names[3]);
        double torsion_val = quad.torsion_gemmi(*res2);
        EXPECT_NEAR(torsion_val, -47.1849, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}