#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "function.hh"

TEST(OracleTest, atom_tree_t) {
    // case: valid standard residue (MET, chain A, residue 1)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model &mol = st.models[0];

        gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : mol.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 1) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "MET");

        std::string alt_conf = "";

        coot::protein_geometry geom;
        geom.init_standard();
        auto restraints_info = geom.get_monomer_restraints(res->name, 0);
        EXPECT_TRUE(restraints_info.first);

        EXPECT_NO_THROW({
            coot::atom_tree_t tree;
            tree.atom_tree_t_gemmi(restraints_info.second, res, alt_conf);
            EXPECT_EQ(tree.name_to_index.size(), 22);
            EXPECT_EQ(tree.n_selected_atoms, 0);
            EXPECT_EQ(tree.atom_vertex_vec.size(), 8);
        });
    }

    // case: non-existent residue (chain Z does not exist)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model &mol = st.models[0];

        gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : mol.chains) {
            if (chain.name == "Z") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 9999) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        EXPECT_EQ(res, nullptr);
    }

    // case: another valid residue (chain A, residue 10)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model &mol = st.models[0];

        gemmi::Residue *res = nullptr;
        for (gemmi::Chain &chain : mol.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &r : chain.residues) {
                    if (r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "ILE");

        std::string alt_conf = "";

        coot::protein_geometry geom;
        geom.init_standard();
        auto restraints_info = geom.get_monomer_restraints(res->name, 0);
        EXPECT_TRUE(restraints_info.first);

        EXPECT_NO_THROW({
            coot::atom_tree_t tree;
            tree.atom_tree_t_gemmi(restraints_info.second, res, alt_conf);
            EXPECT_EQ(tree.name_to_index.size(), 24);
            EXPECT_EQ(tree.n_selected_atoms, 0);
            EXPECT_EQ(tree.atom_vertex_vec.size(), 8);
        });
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
