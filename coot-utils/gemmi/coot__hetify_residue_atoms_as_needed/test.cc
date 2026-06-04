#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

static int count_het_atoms(const gemmi::Residue& res) {
    int count = 0;
    for (const auto& atom : res.atoms) {
        if (atom.flag == 'H') count++;
    }
    return count;
}

// Simulate MMDB's parsing behavior: HETATM records get Het=true.
// In gemmi, read_pdb_file does NOT set flag for HETATM atoms.
// Pre-set flag='H' for atoms in non-standard residues to match MMDB initial state.
static void pre_set_het_flags_from_hetatm(gemmi::Model& model) {
    auto standard = coot::PDB_standard_residue_types();
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            if (!coot::is_member_p_gemmi(standard, res.name)) {
                // This residue came from HETATM records — set flag like MMDB does
                for (auto& atom : res.atoms) {
                    atom.flag = 'H';
                }
            }
        }
    }
}

TEST(OracleTest, hetify_residue_atoms_as_needed) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1);

    gemmi::Model& model = st.models[0];

    // Replicate MMDB's initial state: HETATM atoms already have flag='H'
    pre_set_het_flags_from_hetatm(model);

    // Case 1: HOH residue (non-standard) — should be hetified
    {
        gemmi::Residue* res = find_residue(model, "B", 1);
        ASSERT_NE(res, nullptr);

        std::string resname = res->name;
        EXPECT_EQ(resname, "HOH");

        int natoms = res->atoms.size();
        int het_before = count_het_atoms(*res);

        int result = coot::hetify_residue_atoms_as_needed_gemmi(res);

        int het_after = count_het_atoms(*res);

        EXPECT_EQ(result, 1);
        EXPECT_EQ(het_before, 1);
        EXPECT_EQ(het_after, 1);
    }

    // Case 2: Standard residue (GLU) — should NOT be hetified (returns 0)
    {
        gemmi::Residue* res = find_residue(model, "A", 2);
        ASSERT_NE(res, nullptr);

        std::string resname = res->name;
        EXPECT_EQ(resname, "GLU");

        int result = coot::hetify_residue_atoms_as_needed_gemmi(res);
        EXPECT_EQ(result, 0);
    }

    // Case 3: Nullptr — should return 0 without crashing
    {
        gemmi::Residue* res = nullptr;
        int result = coot::hetify_residue_atoms_as_needed_gemmi(res);
        EXPECT_EQ(result, 0);
    }

    // Case 4: Complementary case — use hetify_residue_atoms (unconditional)
    // on a standard residue to actually cause mutation (addresses coverage warning)
    {
        gemmi::Residue* res = find_residue(model, "A", 2);
        ASSERT_NE(res, nullptr);

        int natoms = res->atoms.size();

        int het_before = count_het_atoms(*res);

        int result = coot::hetify_residue_atoms_gemmi(*res);

        int het_after = count_het_atoms(*res);

        EXPECT_EQ(het_before, 0);
        EXPECT_EQ(result, natoms);
        EXPECT_EQ(het_after, natoms);
        EXPECT_GT(natoms, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}