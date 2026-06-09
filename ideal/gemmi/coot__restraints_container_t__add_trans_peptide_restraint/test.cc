#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

static std::string residue_label(gemmi::Chain* chain, gemmi::Residue* res) {
    if (!res) return "nullptr";
    std::string s;
    if (chain) {
        s += chain->name;
    }
    s += "/";
    s += std::to_string(res->seqid.num.value);
    s += "/";
    s += res->name;
    return s;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTestGemmi, add_trans_peptide_restraint) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Get chain A
    gemmi::Chain* chain_A = nullptr;
    for (auto& ch : model.chains) {
        if (ch.name == "A") {
            chain_A = &ch;
            break;
        }
    }
    ASSERT_NE(chain_A, nullptr);

    // Find residues with seqnum 1 (MET) and 2 (GLU) in chain A
    gemmi::Residue* resA = nullptr;
    gemmi::Residue* resB = nullptr;
    for (auto& res : chain_A->residues) {
        if (res.seqid.num.value == 1) {
            resA = &res;
        }
        if (res.seqid.num.value == 2) {
            resB = &res;
        }
    }
    ASSERT_NE(resA, nullptr);
    ASSERT_NE(resB, nullptr);

    std::string labelA = residue_label(chain_A, resA);
    std::string labelB = residue_label(chain_A, resB);

    // Build a simple restraints container (default constructor)
    coot::restraints_container_t restraints;

    // === Case 1: add trans-peptide restraint between two consecutive residues ===
    unsigned int before1 = restraints.size();

    restraints.add_trans_peptide_restraint_gemmi(chain_A, resA, chain_A, resB);

    unsigned int after1 = restraints.size();

    // === Case 2: same pair again — should find existing restraint, not add new one ===
    unsigned int before2 = restraints.size();

    restraints.add_trans_peptide_restraint_gemmi(chain_A, resA, chain_A, resB);

    unsigned int after2 = restraints.size();

    // Assertions - labels
    EXPECT_EQ(labelA, "A/1/MET");
    EXPECT_EQ(labelB, "A/2/GLU");

    // Assertions for case 1: one new trans-peptide restraint added
    EXPECT_EQ(after1, before1 + 1);

    // Assertions for case 2: size should not change (existing restraint found, not duplicated)
    EXPECT_EQ(after2, before2);
}
