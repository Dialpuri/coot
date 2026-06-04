#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a gemmi residue by chain ID and sequence number
static gemmi::Residue* find_residue(
    gemmi::Model& model,
    const std::string& chain_id,
    int seqnum)
{
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == seqnum)
                return &res;
        }
    }
    return nullptr;
}

TEST(OracleTest, is_nucleotide_by_dict_dynamic_add) {
    // Setup: load PDB and initialise geometry dictionary
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Standard protein residue (should return false — not a nucleotide)
    {
        gemmi::Residue* res = find_residue(model, "A", 10);
        ASSERT_NE(res, nullptr);
        bool is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add_gemmi(res, &geom);
        EXPECT_FALSE(is_nuc);
    }

    // Case 2: Another protein residue (HOH in chain B) — should return false
    {
        gemmi::Residue* res = find_residue(model, "B", 1);
        ASSERT_NE(res, nullptr);
        bool is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add_gemmi(res, &geom);
        EXPECT_FALSE(is_nuc);
    }

    // Case 3: Temporarily rename residue to RNA nucleotide "A" to test true path
    {
        gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);
        std::string original_name = res->name;
        res->name = "A";  // RNA adenosine
        bool is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add_gemmi(res, &geom);
        res->name = original_name;  // restore
        EXPECT_TRUE(is_nuc);
    }

    // Case 4: Temporarily rename residue to DNA nucleotide "DA" to test true path
    {
        gemmi::Residue* res = find_residue(model, "A", 2);
        ASSERT_NE(res, nullptr);
        std::string original_name = res->name;
        res->name = "DA";  // DNA deoxyadenosine
        bool is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add_gemmi(res, &geom);
        res->name = original_name;  // restore
        EXPECT_TRUE(is_nuc);
    }

    // Case 5: Unknown residue name to test dynamic add path (should return false)
    {
        gemmi::Residue* res = find_residue(model, "A", 5);
        ASSERT_NE(res, nullptr);
        std::string original_name = res->name;
        res->name = "XXXX";  // unknown residue, not in dictionary
        bool is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add_gemmi(res, &geom);
        res->name = original_name;  // restore
        EXPECT_FALSE(is_nuc);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}