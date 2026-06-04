#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, deep_copy_this_residue_and_make_asc) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    auto& model = st.models[0];

    // Case 1: first residue of chain A (seqnum 1, MET), whole_residue_flag=1
    {
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") { chain = &c; break; }
        }
        ASSERT_NE(chain, nullptr);
        ASSERT_TRUE(chain->residues.size() > 0);

        gemmi::CRA cra{chain, &chain->residues[0], nullptr};

        std::string altconf = "";
        short int whole_residue_flag = 1;
        int atom_index_handle = -1;
        int udd_afix_handle = -1;

        auto result = coot::deep_copy_this_residue_and_make_asc_gemmi(st, cra, altconf,
                                                                      whole_residue_flag,
                                                                      atom_index_handle,
                                                                      udd_afix_handle);

        auto* new_res = result.first.get();
        atom_selection_container_t asc = result.second;

        EXPECT_NE(new_res, nullptr);
        if (new_res) {
            EXPECT_EQ(new_res->seqid.num.value, 1);
            EXPECT_EQ(new_res->name, "MET");
        }
        EXPECT_EQ(asc.n_selected_atoms, 8);
    }

    // Case 2: second residue of chain A (seqnum 2), whole_residue_flag=0, altconf="A"
    {
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") { chain = &c; break; }
        }
        ASSERT_NE(chain, nullptr);
        ASSERT_TRUE(chain->residues.size() > 1);

        gemmi::CRA cra{chain, &chain->residues[1], nullptr};

        std::string altconf = "A";
        short int whole_residue_flag = 0;
        int atom_index_handle = -1;
        int udd_afix_handle = -1;

        auto result = coot::deep_copy_this_residue_and_make_asc_gemmi(st, cra, altconf,
                                                                      whole_residue_flag,
                                                                      atom_index_handle,
                                                                      udd_afix_handle);

        auto* new_res = result.first.get();
        atom_selection_container_t asc = result.second;

        EXPECT_NE(new_res, nullptr);
        if (new_res) {
            EXPECT_EQ(new_res->seqid.num.value, 2);
        }
        EXPECT_EQ(asc.n_selected_atoms, 9);
    }

    // Case 3: first residue of chain B (seqnum 1, HOH), whole_residue_flag=1
    {
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "B") { chain = &c; break; }
        }
        ASSERT_NE(chain, nullptr);
        ASSERT_TRUE(chain->residues.size() > 0);

        gemmi::CRA cra{chain, &chain->residues[0], nullptr};

        std::string altconf = "";
        short int whole_residue_flag = 1;
        int atom_index_handle = -1;
        int udd_afix_handle = -1;

        auto result = coot::deep_copy_this_residue_and_make_asc_gemmi(st, cra, altconf,
                                                                      whole_residue_flag,
                                                                      atom_index_handle,
                                                                      udd_afix_handle);

        auto* new_res = result.first.get();
        atom_selection_container_t asc = result.second;

        EXPECT_NE(new_res, nullptr);
        if (new_res) {
            EXPECT_EQ(new_res->seqid.num.value, 1);
            EXPECT_EQ(new_res->name, "HOH");
        }
        EXPECT_EQ(asc.n_selected_atoms, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}