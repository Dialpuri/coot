#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to find an atom by name in a gemmi residue
static const gemmi::Atom* find_atom(const gemmi::Residue& res, const std::string& name) {
    for (const auto& a : res.atoms) {
        if (a.name == name) {
            return &a;
        }
    }
    return nullptr;
}

TEST(OracleTest, add_methyl_Hs) {
    // case: VAL CG1 methyl (XHG1, XHG2, XHG3)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model *model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain *chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue *val_res = nullptr;
        for (auto& r : chain->residues) {
            if (r.name == "VAL") {
                if (find_atom(r, "CA") && find_atom(r, "CB") && find_atom(r, "CG1")) {
                    val_res = &r;
                    break;
                }
            }
        }
        ASSERT_NE(val_res, nullptr);

        EXPECT_EQ(std::string(val_res->name), "VAL");
        EXPECT_EQ(std::string(chain->name), "A");
        EXPECT_EQ(val_res->seqid.num.value, 7);

        coot::reduce::torsion_info_t ti("CA", "CB", "CG1", 1.09, 109.0, 180.0);

        // Count existing H atoms before call
        int h_found_before = 0;
        if (find_atom(*val_res, "XHG1")) h_found_before++;
        if (find_atom(*val_res, "XHG2")) h_found_before++;
        if (find_atom(*val_res, "XHG3")) h_found_before++;

        EXPECT_EQ(h_found_before, 0);

        coot::reduce::atoms_with_spinnable_Hs spinables;
        coot::reduce::add_methyl_Hs_gemmi("XHG1", "XHG2", "XHG3", ti, *val_res, spinables);

        int h_found_after = 0;
        const gemmi::Atom *h1 = find_atom(*val_res, "XHG1"); if (h1) h_found_after++;
        const gemmi::Atom *h2 = find_atom(*val_res, "XHG2"); if (h2) h_found_after++;
        const gemmi::Atom *h3 = find_atom(*val_res, "XHG3"); if (h3) h_found_after++;

        EXPECT_EQ(h_found_after, 3);
        EXPECT_EQ(h_found_after - h_found_before, 3);

        ASSERT_NE(h1, nullptr);
        ASSERT_NE(h2, nullptr);
        ASSERT_NE(h3, nullptr);

        EXPECT_NEAR(h1->pos.x, 28.9025, 1e-4);
        EXPECT_NEAR(h1->pos.y, 8.53075, 1e-4);
        EXPECT_NEAR(h1->pos.z, 76.378, 1e-4);

        EXPECT_NEAR(h2->pos.x, 28.4305, 1e-4);
        EXPECT_NEAR(h2->pos.y, 10.2192, 1e-4);
        EXPECT_NEAR(h2->pos.z, 76.042, 1e-4);

        EXPECT_NEAR(h3->pos.x, 28.1885, 1e-4);
        EXPECT_NEAR(h3->pos.y, 9.53219, 1e-4);
        EXPECT_NEAR(h3->pos.z, 77.6717, 1e-4);
    }

    // case: VAL CG2 methyl (YHG1, YHG2, YHG3)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model *model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain *chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue *val2 = nullptr;
        for (auto& r : chain->residues) {
            if (r.name == "VAL" && r.seqid.num.value != 7) {
                if (find_atom(r, "CA") && find_atom(r, "CB") && find_atom(r, "CG2")) {
                    val2 = &r;
                    break;
                }
            }
        }

        if (val2) {
            EXPECT_EQ(std::string(val2->name), "VAL");
            EXPECT_EQ(std::string(chain->name), "A");
            EXPECT_EQ(val2->seqid.num.value, 17);

            coot::reduce::torsion_info_t ti2("CA", "CB", "CG2", 1.09, 109.0, 180.0);

            int h_before2 = 0;
            if (find_atom(*val2, "YHG1")) h_before2++;
            if (find_atom(*val2, "YHG2")) h_before2++;
            if (find_atom(*val2, "YHG3")) h_before2++;
            EXPECT_EQ(h_before2, 0);

            coot::reduce::atoms_with_spinnable_Hs spinables;
            coot::reduce::add_methyl_Hs_gemmi("YHG1", "YHG2", "YHG3", ti2, *val2, spinables);

            int h_after2 = 0;
            const gemmi::Atom *h21 = find_atom(*val2, "YHG1"); if (h21) h_after2++;
            const gemmi::Atom *h22 = find_atom(*val2, "YHG2"); if (h22) h_after2++;
            const gemmi::Atom *h23 = find_atom(*val2, "YHG3"); if (h23) h_after2++;

            EXPECT_EQ(h_after2, 3);
            EXPECT_EQ(h_after2 - h_before2, 3);

            ASSERT_NE(h21, nullptr);
            ASSERT_NE(h22, nullptr);
            ASSERT_NE(h23, nullptr);

            EXPECT_NEAR(h21->pos.x, 17.5936, 1e-4);
            EXPECT_NEAR(h21->pos.y, 14.9264, 1e-4);
            EXPECT_NEAR(h21->pos.z, 66.9454, 1e-4);

            EXPECT_NEAR(h22->pos.x, 18.1564, 1e-4);
            EXPECT_NEAR(h22->pos.y, 15.0459, 1e-4);
            EXPECT_NEAR(h22->pos.z, 65.2555, 1e-4);

            EXPECT_NEAR(h23->pos.x, 17.047, 1e-4);
            EXPECT_NEAR(h23->pos.y, 13.7326, 1e-4);
            EXPECT_NEAR(h23->pos.z, 65.7361, 1e-4);
        }
    }

    // case: GLY (no CB/CG1) - should add 0 atoms
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        gemmi::Model *model = &st.models[0];
        ASSERT_NE(model, nullptr);

        gemmi::Chain *chain = &model->chains[0];
        ASSERT_NE(chain, nullptr);

        gemmi::Residue *gly_res = nullptr;
        for (auto& r : chain->residues) {
            if (r.name == "GLY") {
                gly_res = &r;
                break;
            }
        }

        if (gly_res) {
            coot::reduce::torsion_info_t ti("CA", "CB", "CG1", 1.09, 109.0, 180.0);

            int h_before3 = 0;
            if (find_atom(*gly_res, "ZHG1")) h_before3++;
            if (find_atom(*gly_res, "ZHG2")) h_before3++;
            if (find_atom(*gly_res, "ZHG3")) h_before3++;
            EXPECT_EQ(h_before3, 0);

            coot::reduce::atoms_with_spinnable_Hs spinables;
            coot::reduce::add_methyl_Hs_gemmi("ZHG1", "ZHG2", "ZHG3", ti, *gly_res, spinables);

            int h_after3 = 0;
            if (find_atom(*gly_res, "ZHG1")) h_after3++;
            if (find_atom(*gly_res, "ZHG2")) h_after3++;
            if (find_atom(*gly_res, "ZHG3")) h_after3++;
            EXPECT_EQ(h_after3, 0);
            EXPECT_EQ(h_after3 - h_before3, 0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}