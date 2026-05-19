#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ShelxInsMessageForAtom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Check that loading succeeded
    EXPECT_GE(st.models.size(), 1) << "Failed to load PDB file";
    
    // Create a ShelxIns object
    coot::ShelxIns shelx_ins;
    
    // Get model
    gemmi::Model& model = st.models[0];
    
    // Find residue 10 (ILE) in chain A
    gemmi::Residue* res = nullptr;
    for (auto& chain : model.chains) {
        for (auto& r : chain.residues) {
            if (r.seqid.num.value == 10 && chain.name == "A") {
                res = &r;
                break;
            }
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Failed to get residue 10 in chain A";
    
    // Get first atom (should be N)
    ASSERT_GE(res->atoms.size(), 1) << "Residue has no atoms";
    gemmi::CRA cra{&model.chains[0], res, &res->atoms[0]};
    
    // Case 1: Call with a prefix string
    {
        std::string prefix = "WARNING:: test message for atom ";
        std::string result = shelx_ins.message_for_atom_gemmi(prefix, cra);
        EXPECT_EQ(prefix, "WARNING:: test message for atom ");
        // Check that result contains expected format
        EXPECT_NE(result.find("WARNING:: test message for atom"), std::string::npos);
        EXPECT_NE(result.find("\"A\""), std::string::npos);
        EXPECT_NE(result.find("10"), std::string::npos);
        EXPECT_NE(result.find("\"ILE\""), std::string::npos);
        // Atom name in original format with quotes - need to match what gemmi produces
        // The original MMDB version produced quotes around the atom name
        EXPECT_NE(result.find("\"N\""), std::string::npos);
    }
    
    // Case 2: Call with empty prefix
    {
        std::string prefix = "";
        std::string result = shelx_ins.message_for_atom_gemmi(prefix, cra);
        EXPECT_EQ(prefix, "");
        // Check that result contains expected format without prefix
        EXPECT_NE(result.find("\"A\""), std::string::npos);
        EXPECT_NE(result.find("10"), std::string::npos);
        EXPECT_NE(result.find("\"ILE\""), std::string::npos);
        EXPECT_NE(result.find("\"N\""), std::string::npos);
    }
    
    // Case 3: Find residue 11 (GLY) and look for an atom with altLoc
    gemmi::Residue* res2 = nullptr;
    for (auto& chain : model.chains) {
        for (auto& r : chain.residues) {
            if (r.seqid.num.value == 11 && chain.name == "A") {
                res2 = &r;
                break;
            }
        }
        if (res2) break;
    }
    
    if (res2 != nullptr) {
        bool found_altloc = false;
        for (auto& atom : res2->atoms) {
            if (atom.altloc != '\0' && atom.altloc != ' ') {
                gemmi::CRA cra2{&model.chains[0], res2, &atom};
                std::string prefix = "ALTLOC atom: ";
                std::string result = shelx_ins.message_for_atom_gemmi(prefix, cra2);
                EXPECT_NE(result.find("ALTLOC atom:"), std::string::npos);
                found_altloc = true;
                break;
            }
        }
        // Note: If no altloc is found, that's okay - the test just skips this case
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}