#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_fragment_from_atom_spec) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    // Find a CA atom to get exact name format and chain info
    std::string ca_atom_name;
    int ca_res_no = -1;
    std::string ca_chain_id;

    bool found_ca = false;
    for (gemmi::Chain& ch : model.chains) {
        if (found_ca) break;
        for (gemmi::Residue& res : ch.residues) {
            if (found_ca) break;
            for (gemmi::Atom& at : res.atoms) {
                if (at.name == "CA") {
                    ca_atom_name = at.name;
                    ca_res_no = res.seqid.num.value;
                    ca_chain_id = ch.name;
                    found_ca = true;
                    break;
                }
            }
        }
    }
    ASSERT_TRUE(found_ca) << "Could not find a CA atom in the PDB";

    // Case 1: Valid atom in a continuous stretch — expects fragment
    {
        coot::atom_spec_t atom_spec(ca_chain_id, ca_res_no, "", ca_atom_name, "");
        std::pair<std::optional<gemmi::Structure>, std::vector<coot::residue_spec_t>> result =
            coot::util::get_fragment_from_atom_spec_gemmi(atom_spec, st);

        EXPECT_TRUE(result.first.has_value());  // mol_found: true
        EXPECT_EQ(result.second.size(), 36u);  // fragment_size: 36

        if (!result.second.empty()) {
            coot::residue_spec_t first = result.second.front();
            EXPECT_EQ(first.chain_id, "A");
            EXPECT_EQ(first.res_no, 1);

            coot::residue_spec_t last = result.second.back();
            EXPECT_EQ(last.chain_id, "A");
            EXPECT_EQ(last.res_no, 36);
        }
    }

    // Case 2: Atom with invalid alt_conf — should not find the atom
    {
        coot::atom_spec_t atom_spec(ca_chain_id, ca_res_no, "", ca_atom_name, "A");
        std::pair<std::optional<gemmi::Structure>, std::vector<coot::residue_spec_t>> result =
            coot::util::get_fragment_from_atom_spec_gemmi(atom_spec, st);

        EXPECT_FALSE(result.first.has_value());  // mol_found: false
        EXPECT_EQ(result.second.size(), 0u);  // fragment_size: 0
    }

    // Case 3: Non-existent residue number — should return empty
    {
        coot::atom_spec_t atom_spec(ca_chain_id, 9999, "", ca_atom_name, "");
        std::pair<std::optional<gemmi::Structure>, std::vector<coot::residue_spec_t>> result =
            coot::util::get_fragment_from_atom_spec_gemmi(atom_spec, st);

        EXPECT_FALSE(result.first.has_value());  // mol_found: false
        EXPECT_EQ(result.second.size(), 0u);  // fragment_size: 0
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}