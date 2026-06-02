#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_OH_H) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // --- Case 1: SER residue with valid OH hydrogen placement (but lookup fails) ---
    {
        gemmi::Residue* res_p = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 116) {
                    res_p = &res;
                    break;
                }
            }
            if (res_p) break;
        }
        ASSERT_NE(res_p, nullptr);
        EXPECT_EQ(res_p->name, "ALA");
        EXPECT_EQ(res_p->seqid.num.value, 116);

        std::string H_at_name = " HG ";
        std::string first_neighb = " OG ";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back(" CB ");
        std::map<std::string, std::vector<std::string>> third_neighb_map;
        std::vector<std::string> ca_vec;
        ca_vec.push_back(" CA ");
        third_neighb_map[" CB "] = ca_vec;

        double bond_length = 0.96;
        double ang_deg = 109.5;
        double torsion_deg = 180.0;

        EXPECT_EQ(second_neighb_vec.size(), 1u);
        EXPECT_EQ(third_neighb_map.size(), 1u);

        int atoms_before = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_before, 5);

        coot::reduce::add_OH_H_gemmi(H_at_name, first_neighb, second_neighb_vec, third_neighb_map,
                                     bond_length, ang_deg, torsion_deg, res_p);

        int atoms_after = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_after, 5);
        EXPECT_EQ((atoms_after - atoms_before), 0);

        // Check that HG atom does not exist
        bool h_found = false;
        for (const auto& atom : res_p->atoms) {
            if (coot::reduce::trim(atom.name) == "HG") {
                h_found = true;
                break;
            }
        }
        EXPECT_EQ(h_found, false);
    }

    // --- Case 2: Empty second_neighb_vec — should do nothing (guarded path) ---
    {
        gemmi::Residue* res_p = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 120) {
                    res_p = &res;
                    break;
                }
            }
            if (res_p) break;
        }
        ASSERT_NE(res_p, nullptr);
        EXPECT_EQ(res_p->name, "SER");

        std::string H_at_name = " HG ";
        std::string first_neighb = " OG ";
        std::vector<std::string> second_neighb_vec;  // empty!
        std::map<std::string, std::vector<std::string>> third_neighb_map;

        EXPECT_TRUE(second_neighb_vec.empty());

        int atoms_before = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_before, 6);

        coot::reduce::add_OH_H_gemmi(H_at_name, first_neighb, second_neighb_vec, third_neighb_map,
                                     0.96, 109.5, 180.0, res_p);

        int atoms_after = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_after, 6);
        EXPECT_EQ((atoms_after - atoms_before), 0);
    }

    // --- Case 3: second_neighb present but key missing from third_neighb_map ---
    {
        gemmi::Residue* res_p = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 130) {
                    res_p = &res;
                    break;
                }
            }
            if (res_p) break;
        }
        ASSERT_NE(res_p, nullptr);
        EXPECT_EQ(res_p->name, "PRO");

        std::string H_at_name = " HG ";
        std::string first_neighb = " OG ";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back(" CB ");
        std::map<std::string, std::vector<std::string>> third_neighb_map;  // empty — won't find " CB "

        EXPECT_EQ(second_neighb_vec.size(), 1u);
        EXPECT_TRUE(third_neighb_map.empty());

        int atoms_before = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_before, 7);

        coot::reduce::add_OH_H_gemmi(H_at_name, first_neighb, second_neighb_vec, third_neighb_map,
                                     0.96, 109.5, 180.0, res_p);

        int atoms_after = static_cast<int>(res_p->atoms.size());
        EXPECT_EQ(atoms_after, 7);
        EXPECT_EQ((atoms_after - atoms_before), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}