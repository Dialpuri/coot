#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <fstream>
#include <string>
#include <vector>
#include "function.hh"

TEST(OracleTest, write_ins_file_internal) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    int n_models = static_cast<int>(st.models.size());
    EXPECT_EQ(n_models, 1);

    int n_chains = 0;
    if (n_models > 0) {
        n_chains = static_cast<int>(st.models[0].chains.size());
    }
    EXPECT_EQ(n_chains, 2);

    int total_atoms = 0;
    if (n_models > 0) {
        for (const gemmi::Model &md : st.models) {
            for (const gemmi::Chain &ch : md.chains) {
                for (const gemmi::Residue &res : ch.residues) {
                    total_atoms += static_cast<int>(res.atoms.size());
                }
            }
        }
    }
    EXPECT_EQ(total_atoms, 4280);

    // Case 1: Write with a proper cell set (mol_is_from_shelx_ins = false)
    {
        coot::ShelxIns shelx;
        shelx.have_cell_flag = true;
        shelx.cell = clipper::Cell(clipper::Cell_descr(100.0, 100.0, 100.0, 90.0, 90.0, 90.0));
        shelx.title = "Test shelx ins file";
        shelx.sfac.push_back("C");
        shelx.sfac.push_back("N");
        shelx.sfac.push_back("O");
        shelx.sfac.push_back("S");

        std::string filename = "/tmp/oracle_test_ins1.ins";
        bool mol_is_from_shelx_ins = false;

        std::pair<int, std::string> result = shelx.write_ins_file_internal_gemmi(st.models[0], filename, mol_is_from_shelx_ins);

        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, "INFO:: SHELXL file /tmp/oracle_test_ins1.ins written.");

        // Check file was written and read first few lines
        std::ifstream fin(filename);
        ASSERT_TRUE(fin.is_open());
        {
            std::string line;
            std::getline(fin, line);
            EXPECT_EQ(line, "TITL PDB->ins");
            std::getline(fin, line);
            EXPECT_EQ(line, "CELL 1.54178  100.0000 100.0000 100.0000 90.0000 90.0000 90.0000");
            std::getline(fin, line);
            EXPECT_EQ(line, "ZERR 4         0.1000  0.1000  0.1000  0.0000  0.0000  0.0000");
            std::getline(fin, line);
            EXPECT_EQ(line, "LATT -1");
            std::getline(fin, line);
            EXPECT_EQ(line, "SYMM X+1/2, -Y+1/2, -Z");
        }
        fin.close();
    }

    // Case 2: No cell set — should fail
    {
        coot::ShelxIns shelx;
        shelx.have_cell_flag = false;  // no cell

        std::string filename = "/tmp/oracle_test_ins2.ins";
        bool mol_is_from_shelx_ins = false;

        std::pair<int, std::string> result = shelx.write_ins_file_internal_gemmi(st.models[0], filename, mol_is_from_shelx_ins);

        EXPECT_EQ(result.first, 0);
        EXPECT_EQ(result.second, "WARNING:: no cell available... failure to write ins file.");
    }

    // Case 3: With cell, mol_is_from_shelx_ins = true (occupancy not forced to 11)
    {
        coot::ShelxIns shelx;
        shelx.have_cell_flag = true;
        shelx.cell = clipper::Cell(clipper::Cell_descr(50.0, 60.0, 70.0, 90.0, 90.0, 90.0));
        shelx.title = "Test shelx ins file from shelx";

        std::string filename = "/tmp/oracle_test_ins3.ins";
        bool mol_is_from_shelx_ins = true;

        std::pair<int, std::string> result = shelx.write_ins_file_internal_gemmi(st.models[0], filename, mol_is_from_shelx_ins);

        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, "INFO:: SHELXL file /tmp/oracle_test_ins3.ins written.");

        // Check file was written
        std::ifstream fin(filename);
        ASSERT_TRUE(fin.is_open());
        {
            std::string line;
            std::getline(fin, line);
            EXPECT_EQ(line, "TITL PDB->ins");
            std::getline(fin, line);
            EXPECT_EQ(line, "CELL 1.54178  50.0000 60.0000 70.0000 90.0000 90.0000 90.0000");
            std::getline(fin, line);
            EXPECT_EQ(line, "ZERR 4         0.0500  0.0600  0.0700  0.0000  0.0000  0.0000");
            std::getline(fin, line);
            EXPECT_EQ(line, "LATT -1");
            std::getline(fin, line);
            EXPECT_EQ(line, "SYMM X+1/2, -Y+1/2, -Z");
        }
        fin.close();
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}