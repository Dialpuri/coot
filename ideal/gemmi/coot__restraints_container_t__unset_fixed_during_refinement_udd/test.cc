#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, unset_fixed_during_refinement_udd) {
    {
        // case: main test - unset UDD and verify idempotency
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty());

        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        ASSERT_EQ(chain.name, "A");

        // Pick residue //A/10
        gemmi::Residue* res = nullptr;
        for (gemmi::Residue& r : chain.residues) {
            if (r.seqid.num.value == 10) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "ILE");

        // Build restraints_container_t_gemmi with atoms from the structure
        coot::restraints_container_t_gemmi restraints;

        // Collect atoms from all residues in chain A
        for (gemmi::Chain& ch : model.chains) {
            for (gemmi::Residue& r : ch.residues) {
                for (gemmi::Atom& a : r.atoms) {
                    restraints.atom.push_back(a);
                    restraints.fixed_during_refinement.push_back(false);
                }
            }
        }

        restraints.n_atoms = restraints.atom.size();

        // Set all to "fixed" (true) to simulate the UDD being set to 1
        for (size_t i = 0; i < restraints.n_atoms; i++) {
            restraints.fixed_during_refinement[i] = true;
        }

        // Verify sum_before = number of atoms (all true)
        int sum_before = 0;
        for (size_t i = 0; i < restraints.n_atoms; i++) {
            sum_before += restraints.fixed_during_refinement[i] ? 1 : 0;
        }
        EXPECT_EQ(sum_before, static_cast<int>(restraints.n_atoms));

        // Call the function under test
        EXPECT_NO_THROW(restraints.unset_fixed_during_refinement_udd_gemmi());

        // Verify values are now 0
        int sum_after = 0;
        for (size_t i = 0; i < restraints.n_atoms; i++) {
            sum_after += restraints.fixed_during_refinement[i] ? 1 : 0;
        }
        EXPECT_EQ(sum_after, 0);

        // Case 2: call again and verify idempotency
        restraints.unset_fixed_during_refinement_udd_gemmi();
        int sum_again = 0;
        for (size_t i = 0; i < restraints.n_atoms; i++) {
            sum_again += restraints.fixed_during_refinement[i] ? 1 : 0;
        }
        EXPECT_EQ(sum_again, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}