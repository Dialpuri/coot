#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps_container_t_kludge_filter) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr) << "No model";

    coot::atom_overlaps_container_t container;

    gemmi::Chain* chain = &model->chains[0];
    ASSERT_NE(chain, nullptr) << "No chain";

    // Get ASN residue (index 2)
    gemmi::Residue* asn_res = &chain->residues[2];
    ASSERT_NE(asn_res, nullptr) << "No ASN residue";
    std::string orig_asn_resname(asn_res->name);

    // Find ND2 atom of ASN (gemmi stores unpadded names)
    gemmi::Atom* at1 = nullptr;
    for (auto& a : asn_res->atoms) {
        if (a.name == "ND2") { at1 = &a; break; }
    }
    if (!at1) at1 = &asn_res->atoms[0];
    ASSERT_NE(at1, nullptr);
    std::string orig_atom_name_1(at1->name);

    // Another atom in same chain (first residue)
    gemmi::Residue* other_res = &chain->residues[0];
    ASSERT_NE(other_res, nullptr) << "No other residue";
    gemmi::Atom* at2 = &other_res->atoms[0];
    ASSERT_NE(at2, nullptr);
    std::string orig_resname_2(other_res->name);
    std::string orig_atom_name_2(at2->name);

    // Case 1: Same chain, ASN/NE2 vs NAG/C2 → reject = true
    {
        at1->name = "NE2";
        at2->name = "C2";
        other_res->name = "NAG";

        gemmi::CRA cra1{chain, asn_res, at1};
        gemmi::CRA cra2{chain, other_res, at2};

        bool result = container.kludge_filter_gemmi(cra1, cra2);
        EXPECT_EQ(asn_res->name, "ASN");
        EXPECT_EQ(at1->name, "NE2");
        EXPECT_EQ(other_res->name, "NAG");
        EXPECT_EQ(at2->name, "C2");
        EXPECT_EQ(cra1.chain->name, cra2.chain->name);
        EXPECT_TRUE(result);

        // Restore names for remaining cases
        asn_res->name = orig_asn_resname;
        at1->name = orig_atom_name_1;
        other_res->name = orig_resname_2;
        at2->name = orig_atom_name_2;
    }

    // Case 2: Different chains → reject = false
    {
        if (model->chains.size() >= 2) {
            gemmi::Chain* chain2 = &model->chains[1];
            if (!chain2->residues.empty()) {
                gemmi::Residue* res3 = &chain2->residues[0];
                if (!res3->atoms.empty()) {
                    gemmi::Atom* at3 = &res3->atoms[0];
                    gemmi::CRA cra1{chain, asn_res, at1};
                    gemmi::CRA cra3{chain2, res3, at3};
                    EXPECT_EQ(chain->name, "A");
                    EXPECT_EQ(chain2->name, "B");
                    bool result2 = container.kludge_filter_gemmi(cra1, cra3);
                    EXPECT_FALSE(result2);
                }
            }
        }
    }

    // Case 3: Same chain, wrong atom/residue names → reject = false
    {
        gemmi::CRA cra1{chain, asn_res, at1};
        gemmi::CRA cra2{chain, other_res, at2};
        bool result3 = container.kludge_filter_gemmi(cra1, cra2);
        EXPECT_EQ(asn_res->name, "ASN");
        EXPECT_EQ(at1->name, "ND2");
        EXPECT_EQ(other_res->name, "MET");
        EXPECT_EQ(at2->name, "N");
        EXPECT_FALSE(result3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}