#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <vector>
#include <utility>
#include <string>
#include <sstream>
#include <iomanip>

#include "api/molecules-container.hh"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, z_weighted_density_score_new) {
    molecules_container_t mc;
    mc.geometry_init_standard();

    // Load PDB using gemmi (instead of MMDB)
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Load MTZ map
    int imol_map = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);
    EXPECT_GE(imol_map, 0);
    EXPECT_TRUE(mc.is_valid_map_molecule(imol_map));

    clipper::Xmap<float> xmap = mc.get_xmap(imol_map);

    // Collect atoms from the gemmi structure
    std::vector<std::pair<gemmi::Atom*, float>> atom_weight_pairs;

    // Gather all atoms from all chains in the first model
    for (gemmi::Chain& chain : st.models[0].chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (gemmi::Atom& atom : res.atoms) {
                int z = atom.element.atomic_number();
                if (z == 0) {
                    z = 6;  // default to carbon, same as atomic_number returning -1
                }
                atom_weight_pairs.push_back(std::make_pair(&atom, static_cast<float>(z)));
            }
        }
    }

    // Case 1: All atoms from the structure
    {
        EXPECT_EQ(static_cast<int>(atom_weight_pairs.size()), 2107);
        float score = coot::util::z_weighted_density_score_new_gemmi(atom_weight_pairs, xmap);
        EXPECT_NEAR(score, 12143.9, 0.1);
    }

    // Case 2: Single atom (first atom)
    {
        std::vector<std::pair<gemmi::Atom*, float>> single_atom;
        single_atom.push_back(atom_weight_pairs[0]);
        float score = coot::util::z_weighted_density_score_new_gemmi(single_atom, xmap);
        gemmi::Atom* a = single_atom[0].first;
        std::string mmdb_name = a->padded_name();
        while (mmdb_name.size() < 4) mmdb_name += ' ';
        EXPECT_EQ(mmdb_name, " N  ");
        EXPECT_NEAR(a->pos.x, 12.334, 1e-4);
        EXPECT_NEAR(a->pos.y, 2.772, 1e-4);
        EXPECT_NEAR(a->pos.z, 74.951, 1e-4);
        EXPECT_NEAR(score, 3.09764, 1e-4);
    }

    // Case 3: Empty vector (edge case)
    {
        std::vector<std::pair<gemmi::Atom*, float>> empty_vec;
        EXPECT_EQ(empty_vec.size(), 0);
        float score = coot::util::z_weighted_density_score_new_gemmi(empty_vec, xmap);
        EXPECT_NEAR(score, 0.0, 1e-4);
    }
}