#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

static const char* example_pdb = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(GemmiTest, general_link_find_close_link) {
    // Initialize protein geometry
    coot::protein_geometry pg;
    pg.init_standard();

    // Load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file(example_pdb);
    EXPECT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Find chain A
    const gemmi::Chain* chain_a = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            chain_a = &chain;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Find residues 1 and 2 in chain A
    const gemmi::Residue* r1 = nullptr;
    const gemmi::Residue* r2 = nullptr;
    for (const gemmi::Residue& res : chain_a->residues) {
        if (res.seqid.num.value == 1) {
            r1 = &res;
        }
        if (res.seqid.num.value == 2) {
            r2 = &res;
        }
    }
    ASSERT_NE(r1, nullptr);
    ASSERT_NE(r2, nullptr);

    // ============================================================
    // Oracle Case 1: consecutive residues with empty link_infos
    // ============================================================
    {
        std::string comp_id_1 = r1->name;
        std::string comp_id_2 = r2->name;
        std::string group_1 = std::to_string(r1->seqid.num.value);
        std::string group_2 = std::to_string(r2->seqid.num.value);

        std::vector<coot::chem_link> link_infos = pg.matching_chem_links(comp_id_1, group_1,
                                                                          comp_id_2, group_2);
        EXPECT_EQ(link_infos.size(), 0u);

        bool order_switch = false;
        std::pair<std::string, bool> result = coot::general_link_find_close_link_gemmi(
            link_infos, r1, r2, order_switch, pg);

        EXPECT_EQ(result.first, "");
        EXPECT_FALSE(result.second);
    }

    // ============================================================
    // Oracle Case 2: reversed order with empty link_infos
    // ============================================================
    {
        std::string comp_id_1 = r1->name;
        std::string comp_id_2 = r2->name;
        std::string group_1 = std::to_string(r1->seqid.num.value);
        std::string group_2 = std::to_string(r2->seqid.num.value);

        std::vector<coot::chem_link> link_infos = pg.matching_chem_links(comp_id_1, group_1,
                                                                          comp_id_2, group_2);
        EXPECT_EQ(link_infos.size(), 0u);

        bool order_switch = false;
        std::pair<std::string, bool> result = coot::general_link_find_close_link_gemmi(
            link_infos, r2, r1, order_switch, pg);

        EXPECT_EQ(result.first, "");
        EXPECT_FALSE(result.second);
    }

    // ============================================================
    // Oracle Case 3: non-consecutive residues (A/1 and A/100)
    // ============================================================
    {
        const gemmi::Residue* r_far = nullptr;
        for (const gemmi::Residue& res : chain_a->residues) {
            if (res.seqid.num.value == 100) {
                r_far = &res;
                break;
            }
        }

        if (r_far) {
            std::string comp_id_1 = r1->name;
            std::string group_1 = std::to_string(r1->seqid.num.value);
            std::string comp_id_far = r_far->name;
            std::string group_far = std::to_string(r_far->seqid.num.value);

            std::vector<coot::chem_link> link_infos_far = pg.matching_chem_links(
                comp_id_1, group_1, comp_id_far, group_far);
            EXPECT_EQ(link_infos_far.size(), 0u);

            bool order_switch = false;
            std::pair<std::string, bool> result = coot::general_link_find_close_link_gemmi(
                link_infos_far, r1, r_far, order_switch, pg);

            EXPECT_EQ(result.first, "");
            EXPECT_FALSE(result.second);
        }
    }

    // ============================================================
    // Complementary Case: exercise core logic with non-empty link_infos
    // ============================================================
    {
        std::string comp_id_1 = r1->name;
        std::string comp_id_2 = r2->name;
        std::string group_1 = std::to_string(r1->seqid.num.value);
        std::string group_2 = std::to_string(r2->seqid.num.value);

        std::vector<coot::chem_link> link_infos_pep = pg.matching_chem_links(
            comp_id_1, group_1, comp_id_2, group_2);

        if (!link_infos_pep.empty()) {
            // Case: forward order (r1=MET, r2=GLU)
            {
                bool order_switch = false;
                std::pair<std::string, bool> result = coot::general_link_find_close_link_gemmi(
                    link_infos_pep, r1, r2, order_switch, pg);

                EXPECT_FALSE(result.first.empty());
                EXPECT_FALSE(result.second);
            }

            // Case: reversed order (r2=GLU, r1=MET) — should switch
            {
                bool order_switch = false;
                std::pair<std::string, bool> result = coot::general_link_find_close_link_gemmi(
                    link_infos_pep, r2, r1, order_switch, pg);

                EXPECT_FALSE(result.first.empty());
                EXPECT_TRUE(result.second);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}