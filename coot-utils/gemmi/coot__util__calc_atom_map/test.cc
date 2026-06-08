#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

TEST(OracleTest, calc_atom_map) {
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    gemmi::Model& model = st.models[0];

    clipper::Cell cell(clipper::Cell_descr(100, 100, 100, 90, 90, 90));
    clipper::Spacegroup sg(clipper::Spacegroup::P1);
    clipper::Grid_sampling gs(50, 50, 50);

    // Case 1: Empty selection (original selected model 0 which was empty in MMDB)
    {
        std::vector<gemmi::Atom const*> atoms;
        EXPECT_EQ(atoms.size(), 0u);

        clipper::Xmap<float> xmap =
            coot::util::calc_atom_map_gemmi(atoms, cell, sg, gs);

        EXPECT_FALSE(xmap.is_null());

        EXPECT_EQ(xmap.grid_sampling().nu(), 50);
        EXPECT_EQ(xmap.grid_sampling().nv(), 50);
        EXPECT_EQ(xmap.grid_sampling().nw(), 50);

        float v0 = xmap.get_data(clipper::Coord_grid(0, 0, 0));
        float v1 = xmap.get_data(clipper::Coord_grid(25, 25, 25));
        EXPECT_NEAR(v0, 0.0f, 1e-4);
        EXPECT_NEAR(v1, 0.0f, 1e-4);
    }

    // Case 2: Select just LZA ligand (chain A, residue 1299)
    {
        std::vector<gemmi::Atom const*> atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 1299) continue;
                if (res.name != "LZA") continue;
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back(&atom);
                }
            }
        }
        EXPECT_EQ(atoms.size(), 25u);

        clipper::Xmap<float> xmap =
            coot::util::calc_atom_map_gemmi(atoms, cell, sg, gs);

        EXPECT_FALSE(xmap.is_null());

        EXPECT_EQ(xmap.grid_sampling().nu(), 50);
        EXPECT_EQ(xmap.grid_sampling().nv(), 50);
        EXPECT_EQ(xmap.grid_sampling().nw(), 50);

        float v0 = xmap.get_data(clipper::Coord_grid(0, 0, 0));
        float v1 = xmap.get_data(clipper::Coord_grid(25, 25, 25));
        EXPECT_NEAR(v0, 0.0f, 1e-4);
        EXPECT_NEAR(v1, 0.0f, 1e-4);

    }

    // Case 3: Empty selection (nonexistent residue)
    {
        std::vector<gemmi::Atom const*> atoms;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != "A") continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != 9999) continue;
                for (gemmi::Atom& atom : res.atoms) {
                    atoms.push_back(&atom);
                }
            }
        }
        EXPECT_EQ(atoms.size(), 0u);

        clipper::Xmap<float> xmap =
            coot::util::calc_atom_map_gemmi(atoms, cell, sg, gs);

        EXPECT_FALSE(xmap.is_null());

        EXPECT_EQ(xmap.grid_sampling().nu(), 50);
        EXPECT_EQ(xmap.grid_sampling().nv(), 50);
        EXPECT_EQ(xmap.grid_sampling().nw(), 50);

        float v0 = xmap.get_data(clipper::Coord_grid(0, 0, 0));
        EXPECT_NEAR(v0, 0.0f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}