#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::CRA find_atom_by_name(gemmi::Model& model, std::string chain_id, int res_no, std::string atom_name_unpadded) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value != res_no) continue;
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.name == atom_name_unpadded) {
                    return {&chain, &res, &atom};
                }
            }
        }
    }
    return {nullptr, nullptr, nullptr};
}

TEST(OracleTest, is_angle_related_via_link) {
    {
        // Case 1: Two atoms from same residue (no links in standard PDB)
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];

        gemmi::CRA cra1 = find_atom_by_name(model, "A", 1, "CA");
        gemmi::CRA cra2 = find_atom_by_name(model, "A", 1, "CB");

        EXPECT_NE(cra1.atom, nullptr);
        EXPECT_NE(cra2.atom, nullptr);
        EXPECT_EQ(cra1.atom->name, "CA");
        EXPECT_EQ(cra2.atom->name, "CB");

        // No links in standard PDB
        std::vector<gemmi::Connection> connections;  // empty

        std::vector<std::pair<std::string, std::string>> bonds_for_at_1;
        std::vector<std::pair<std::string, std::string>> bonds_for_at_2;

        bool result = coot::is_angle_related_via_link_gemmi(cra1, cra2, &model, connections, bonds_for_at_1, bonds_for_at_2);
        EXPECT_FALSE(result);
    }

    {
        // Case 2: null first atom (guard clause)
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];

        gemmi::CRA cra1{nullptr, nullptr, nullptr};  // null atom
        gemmi::CRA cra2 = find_atom_by_name(model, "A", 1, "CB");

        EXPECT_NE(cra2.atom, nullptr);

        std::vector<gemmi::Connection> connections;
        std::vector<std::pair<std::string, std::string>> bonds_for_at_1;
        std::vector<std::pair<std::string, std::string>> bonds_for_at_2;

        bool result_null = coot::is_angle_related_via_link_gemmi(cra1, cra2, &model, connections, bonds_for_at_1, bonds_for_at_2);
        EXPECT_FALSE(result_null);
    }

    {
        // Case 3: atoms from different chains (second chain B doesn't exist)
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];

        gemmi::CRA cra1 = find_atom_by_name(model, "A", 1, "CA");
        gemmi::CRA cra3 = find_atom_by_name(model, "B", 1, "CA");  // Chain B doesn't exist

        EXPECT_NE(cra1.atom, nullptr);
        EXPECT_EQ(cra3.atom, nullptr);

        std::vector<gemmi::Connection> connections;
        std::vector<std::pair<std::string, std::string>> bonds_for_at_1;
        std::vector<std::pair<std::string, std::string>> bonds_for_at_3;

        bool result_diff = coot::is_angle_related_via_link_gemmi(cra1, cra3, &model, connections, bonds_for_at_1, bonds_for_at_3);
        EXPECT_FALSE(result_diff);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}