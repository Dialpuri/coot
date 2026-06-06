#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cis_trans_conversion) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Load standard residues for CIS/TRANS templates
        gemmi::Structure st_standard = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/standard-residues.pdb");

        // Get the first model from each structure
        ASSERT_FALSE(st.models.empty());
        ASSERT_FALSE(st_standard.models.empty());

        gemmi::Model& mol = st.models[0];
        gemmi::Model& standard_residues = st_standard.models[0];

        // Case: valid N atom — e.g. //A/10/N
        {
            // Find the atom in chain A, residue 10, named N
            gemmi::Atom* at = nullptr;
            gemmi::CRA cra{};
            for (auto& chain : mol.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        if (res.seqid.num.value == 10) {
                            for (auto& atom : res.atoms) {
                                if (atom.name == "N") {
                                    at = &atom;
                                    cra = {&chain, &res, &atom};
                                    break;
                                }
                            }
                        }
                    }
                }
                if (at) break;
            }

            EXPECT_NE(at, nullptr);

            if (at) {
                bool is_N = (at->name == "N");
                EXPECT_TRUE(is_N);

                int result = coot::util::cis_trans_conversion_gemmi(cra, is_N, &mol, &standard_residues);
                EXPECT_EQ(result, 1);
            }
        }

        // Case: invalid atom
        {
            gemmi::Atom* at = nullptr;
            for (auto& chain : mol.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        if (res.seqid.num.value == 9999) {
                            for (auto& atom : res.atoms) {
                                if (atom.name == "N") {
                                    at = &atom;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (at) break;
            }

            EXPECT_EQ(at, nullptr);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}