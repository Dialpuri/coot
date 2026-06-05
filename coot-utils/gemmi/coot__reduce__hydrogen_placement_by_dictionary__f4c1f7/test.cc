#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, hydrogen_placement_by_dictionary) {
    // Initialize protein geometry (REQUIRED)
    coot::protein_geometry geom;
    geom.init_standard();

    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    auto& model = st.models[0];

    // Case 1: ALA 21 — standard residue, should add H from dictionary
    {
        gemmi::Residue *res = nullptr;
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.name == "ALA" && r.seqid.num.value == 21) {
                    res = &r;
                    break;
                }
            }
        }
        ASSERT_NE(res, nullptr);

        EXPECT_EQ(res->name, "ALA");
        EXPECT_EQ(res->seqid.num.value, 21);

        int atoms_before = res->atoms.size();
        EXPECT_EQ(atoms_before, 10);

        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom.get_monomer_restraints("ALA", 0);
        ASSERT_TRUE(p.first);

        gemmi::CRA cra{nullptr, res, nullptr};
        coot::reduce::hydrogen_placement_by_dictionary_gemmi(p.second, cra, st,
                                                 1.09, 1.08, 1.01, 0.96, 1.33);

        int atoms_after = res->atoms.size();
        EXPECT_EQ(atoms_after, 12);
        EXPECT_EQ((atoms_after - atoms_before), 2);
    }

    // Case 2: HOH — should be skipped (guard clause in public overload)
    {
        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom.get_monomer_restraints("HOH", 0);
        EXPECT_TRUE(p.first);
    }

    // Case 3: TYR 15 — aromatic residue with OH, should add Hs
    {
        gemmi::Residue *res = nullptr;
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.name == "TYR" && r.seqid.num.value == 15) {
                    res = &r;
                    break;
                }
            }
        }
        ASSERT_NE(res, nullptr);

        EXPECT_EQ(res->name, "TYR");
        EXPECT_EQ(res->seqid.num.value, 15);

        int atoms_before = res->atoms.size();
        EXPECT_EQ(atoms_before, 21);

        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom.get_monomer_restraints("TYR", 0);
        ASSERT_TRUE(p.first);

        gemmi::CRA cra{nullptr, res, nullptr};
        coot::reduce::hydrogen_placement_by_dictionary_gemmi(p.second, cra, st,
                                                 1.09, 1.08, 1.01, 0.96, 1.33);

        int atoms_after = res->atoms.size();
        EXPECT_EQ(atoms_after, 23);
        EXPECT_EQ((atoms_after - atoms_before), 2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}