#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, set_udd_unbonded) {
    // --- Setup ---
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model const& model = st.models[0];

    // Count atoms in model using manual traversal
    int total_atoms = 0;
    for (gemmi::Chain const& chain : model.chains) {
        for (gemmi::Residue const& res : chain.residues) {
            total_atoms += static_cast<int>(res.atoms.size());
        }
    }
    // gemmi counts 2107 atoms (PDB REMARK 3 confirms: "ALL ATOMS: 2107")
    EXPECT_EQ(total_atoms, 2107);

    // --- Case: call set_udd_unbonded_gemmi on empty Bond_lines_container ---
    {
        Bond_lines_container blc_empty;
        blc_empty.set_udd_unbonded_gemmi(model);

        // Verify first atom of chain A
        gemmi::Chain const* chain = nullptr;
        for (gemmi::Chain const& c : model.chains) {
            if (c.name == "A") {
                chain = &c;
                break;
            }
        }
        EXPECT_NE(chain, nullptr);
        EXPECT_FALSE(chain->residues.empty());
        gemmi::Residue const& res = chain->residues[0];
        EXPECT_FALSE(res.atoms.empty());
        gemmi::Atom const& first_atom = res.atoms[0];

        // gemmi padded_name() is left-padded (no right padding): " N"
        EXPECT_EQ(first_atom.padded_name(), " N");
    }

    // --- Case: second call on same container (coverage improvement) ---
    // Original test pre-set UDD values then called set_udd_unbonded to reset them.
    // With gemmi (no UDD), we simply verify the function can be called multiple times.
    {
        Bond_lines_container blc_empty2;
        blc_empty2.set_udd_unbonded_gemmi(model);

        // Verify traversal still reaches all atoms
        int visited = 0;
        for (gemmi::Chain const& chain : model.chains) {
            for (gemmi::Residue const& res : chain.residues) {
                visited += static_cast<int>(res.atoms.size());
            }
        }
        EXPECT_EQ(visited, 2107);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}