#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, IsNucleotideFunction) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: a protein residue (ile) — should return 0
    {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.name == "ILE" && res.seqid.num.value == 10) {
                        for (auto& atom : res.atoms) {
                            if (atom.name == "CA") {
                                short int nuc = coot::util::is_nucleotide_gemmi(res);
                                EXPECT_EQ(atom.name, "CA");
                                EXPECT_EQ(res.name, "ILE");
                                EXPECT_EQ(res.seqid.num.value, 10);
                                EXPECT_EQ(nuc, 0);
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        EXPECT_TRUE(found) << "Residue 10 not found";
    }

    // case 2: another protein residue (tyr) — should also return 0
    {
        bool found = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.name == "TYR" && res.seqid.num.value == 15) {
                        for (auto& atom : res.atoms) {
                            if (atom.name == "N") {
                                short int nuc = coot::util::is_nucleotide_gemmi(res);
                                EXPECT_EQ(atom.name, "N");
                                EXPECT_EQ(res.name, "TYR");
                                EXPECT_EQ(res.seqid.num.value, 15);
                                EXPECT_EQ(nuc, 0);
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        EXPECT_TRUE(found) << "Residue 15 not found";
    }

    // case 3: null residue (not applicable in gemmi, skip)
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}