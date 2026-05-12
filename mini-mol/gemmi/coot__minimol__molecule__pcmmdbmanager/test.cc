#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, PdbConversionTest) {
    coot::minimol::molecule test_mol;
    
    coot::minimol::fragment frag;
    frag.fragment_id = "A";
    
    coot::minimol::residue res;
    res.name = "MET";
    res.seqnum = 1;
    res.ins_code = ' ';
    
    coot::minimol::atom atom1;
    atom1.pos = gemmi::Vec3(10.0, 20.0, 30.0);
    atom1.occupancy = 1.0;
    atom1.temperature_factor = 20.0;
    atom1.name = "CA";
    atom1.element = "C";
    atom1.altloc = ' ';
    atom1.int_user_data = -1;
    
    res.atoms.push_back(atom1);
    frag.residues.push_back(res);
    test_mol.fragments.push_back(frag);
    
    test_mol.have_cell = true;
    test_mol.mmdb_cell[0] = 53.512;
    test_mol.mmdb_cell[1] = 71.311;
    test_mol.mmdb_cell[2] = 72.101;
    test_mol.mmdb_cell[3] = 90.0;
    test_mol.mmdb_cell[4] = 90.0;
    test_mol.mmdb_cell[5] = 90.0;
    
    test_mol.have_spacegroup = true;
    test_mol.mmdb_spacegroup = "P 21 21 21";
    
    gemmi::Structure mmdb_mol = test_mol.pcmmdbmanager_gemmi();
    
    EXPECT_EQ(mmdb_mol.models.size(), 1);
    
    if (mmdb_mol.models.size() > 0) {
        gemmi::Model& model = mmdb_mol.models[0];
        EXPECT_EQ(model.chains.size(), 1);
        
        if (model.chains.size() > 0) {
            gemmi::Chain& chain = model.chains[0];
            EXPECT_EQ(chain.name, "A");
            
            EXPECT_EQ(chain.residues.size(), 1);
            
            if (chain.residues.size() > 0) {
                gemmi::Residue& res = chain.residues[0];
                EXPECT_EQ(res.name, "MET");
                EXPECT_EQ(res.seqid.num.value, 1);
                
                EXPECT_EQ(res.atoms.size(), 1);
                if (res.atoms.size() > 0) {
                    gemmi::Atom& atom = res.atoms[0];
                    EXPECT_EQ(atom.name, "CA");
                    EXPECT_EQ(atom.pos.x, 10.0);
                    EXPECT_EQ(atom.pos.y, 20.0);
                    EXPECT_EQ(atom.pos.z, 30.0);
                }
            }
        }
    }
    
    if (test_mol.have_cell) {
        EXPECT_NEAR(mmdb_mol.cell.a, 53.512, 1e-3);
        EXPECT_NEAR(mmdb_mol.cell.b, 71.311, 1e-3);
        EXPECT_NEAR(mmdb_mol.cell.c, 72.101, 1e-3);
        EXPECT_NEAR(mmdb_mol.cell.alpha, 90.0, 1e-3);
        EXPECT_NEAR(mmdb_mol.cell.beta, 90.0, 1e-3);
        EXPECT_NEAR(mmdb_mol.cell.gamma, 90.0, 1e-3);
    }
    
    if (test_mol.have_spacegroup) {
        EXPECT_EQ(mmdb_mol.spacegroup_hm, "P 21 21 21");
    }
    
    coot::minimol::molecule empty_mol;
    
    EXPECT_EQ(empty_mol.fragments.size(), 0);
    
    gemmi::Structure empty_mmdb_mol = empty_mol.pcmmdbmanager_gemmi();
    
    EXPECT_EQ(empty_mmdb_mol.models.size(), 1);
    
    if (empty_mmdb_mol.models.size() > 0) {
        gemmi::Model& empty_model = empty_mmdb_mol.models[0];
        EXPECT_EQ(empty_model.chains.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}