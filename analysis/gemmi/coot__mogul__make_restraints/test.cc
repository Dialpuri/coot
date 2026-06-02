#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mogul_make_restraints) {
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;

    // Find the GLN residue (residue 5 in chain A)
    const gemmi::Residue* res = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (const gemmi::Residue& residue : chain.residues) {
                    if (residue.seqid.num.value == 5) {
                        res = &residue;
                        break;
                    }
                }
            }
        }
    }
    ASSERT_NE(res, nullptr);

    std::string comp_id(res->name);
    EXPECT_EQ(comp_id, "GLN");
    EXPECT_EQ(res->seqid.num.value, 5);

    // Populate mogul.items with synthetic test data for GLN residue
    // GLN has atoms: N(1), CA(2), C(3), O(4), ... (1-indexed in mogul)
    std::vector<coot::mogul_item> items;
    {
        // Add bond items: N-CA (1-2), CA-C (2-3), C-O (3-4)
        {
            coot::mogul_item item;
            item.type = coot::mogul_item::BOND;
            item.idx_1 = 1; item.idx_2 = 2;
            item.median = 1.46f; item.std_dev = 0.02f;
            item.counts = 500; item.mean = 1.46f; item.z = 0.0f;
            items.push_back(item);
        }
        {
            coot::mogul_item item;
            item.type = coot::mogul_item::BOND;
            item.idx_1 = 2; item.idx_2 = 3;
            item.median = 1.53f; item.std_dev = 0.02f;
            item.counts = 500; item.mean = 1.53f; item.z = 0.0f;
            items.push_back(item);
        }
        {
            coot::mogul_item item;
            item.type = coot::mogul_item::BOND;
            item.idx_1 = 3; item.idx_2 = 4;
            item.median = 1.23f; item.std_dev = 0.02f;
            item.counts = 500; item.mean = 1.23f; item.z = 0.0f;
            items.push_back(item);
        }
        // Add angle items: N-CA-C (1-2-3), CA-C-O (2-3-4)
        {
            coot::mogul_item item;
            item.type = coot::mogul_item::ANGLE;
            item.idx_1 = 1; item.idx_2 = 2; item.idx_3 = 3;
            item.median = 110.5f; item.std_dev = 2.5f;
            item.counts = 500; item.mean = 110.5f; item.z = 0.0f;
            items.push_back(item);
        }
        {
            coot::mogul_item item;
            item.type = coot::mogul_item::ANGLE;
            item.idx_1 = 2; item.idx_2 = 3; item.idx_3 = 4;
            item.median = 120.0f; item.std_dev = 3.0f;
            item.counts = 500; item.mean = 120.0f; item.z = 0.0f;
            items.push_back(item);
        }
    }

    // Case 1: valid residue — should produce bonds and angles
    {
        coot::dictionary_residue_restraints_t r = coot::make_restraints_gemmi(res, comp_id, items, imol, geom);
        EXPECT_EQ(r.bond_restraint.size(), 3u);
        EXPECT_EQ(r.angle_restraint.size(), 2u);
    }

    // Case 2: nullptr residue — guarded path returns empty restraints
    {
        coot::dictionary_residue_restraints_t r = coot::make_restraints_gemmi(nullptr, comp_id, items, imol, geom);
        EXPECT_EQ(r.bond_restraint.size(), 0u);
        EXPECT_EQ(r.angle_restraint.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}