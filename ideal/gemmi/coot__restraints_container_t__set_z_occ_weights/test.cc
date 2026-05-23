#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, RestraintZOccWeights) {
    coot::restraints_container_t restraints;
    
    restraints.n_atoms = 5;
    restraints.atom.resize(5);
    restraints.do_neutron_refinement = false;
    restraints.cryo_em_mode = false;
    
    restraints.atom[0].name = " N  ";
    restraints.atom[0].element = gemmi::Element("N");
    restraints.atom[0].occ = 1.0;
    
    restraints.atom[1].name = "CA";
    restraints.atom[1].element = gemmi::Element("C");
    restraints.atom[1].occ = 1.0;
    
    restraints.atom[2].name = "C";
    restraints.atom[2].element = gemmi::Element("C");
    restraints.atom[2].occ = 1.0;
    
    restraints.atom[3].name = "BE  ";
    restraints.atom[3].element = gemmi::Element("Be");
    restraints.atom[3].occ = 0.8;
    
    restraints.atom[4].name = "CB";
    restraints.atom[4].element = gemmi::Element("C");
    restraints.atom[4].occ = 1.0;
    
    restraints.set_z_occ_weights_gemmi();
    
    EXPECT_EQ(restraints.atom_z_occ_weight.size(), 5);
    if (restraints.atom_z_occ_weight.size() > 0) {
        EXPECT_NEAR(restraints.atom_z_occ_weight[0], 7, 1e-4);
        EXPECT_NEAR(restraints.atom_z_occ_weight[1], 6, 1e-4);
        EXPECT_NEAR(restraints.atom_z_occ_weight[2], 6, 1e-4);
        EXPECT_NEAR(restraints.atom_z_occ_weight[3], 3.2, 1e-4);
        EXPECT_NEAR(restraints.atom_z_occ_weight[4], 6, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}