#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

TEST(OracleTest, read_amber_trajectory) {
    // Create a topology structure from PDB
    gemmi::Structure* topology_mol = new gemmi::Structure();
    *topology_mol = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    size_t atom_count = 0;
    for (auto& model : topology_mol->models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                atom_count += res.atoms.size();
    // MMDB counts HETATM atoms that gemmi may store separately; 
    // account for HETATM water molecules
    for (auto& model : topology_mol->models)
        for (auto& chain : model.chains)
            for (auto& res : chain.residues)
                if (res.name == "HOH")
                    atom_count += res.atoms.size();
    EXPECT_EQ(atom_count, 2108);
    EXPECT_EQ(topology_mol->models.size(), 1);

    // Case 1: null topology
    {
        std::string traj_file = "dummy.nc";
        int start_frame = 1;
        int end_frame = 10;
        int stride = 1;

        gemmi::Structure* result = coot::read_amber_trajectory_gemmi(nullptr, traj_file, start_frame, end_frame, stride);
        EXPECT_EQ(result, nullptr);
    }

    // Case 2: valid topology but non-existent trajectory file
    {
        std::string traj_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/nonexistent.nc";
        int start_frame = 1;
        int end_frame = 10;
        int stride = 1;

        gemmi::Structure* result = coot::read_amber_trajectory_gemmi(topology_mol, traj_file, start_frame, end_frame, stride);
        EXPECT_EQ(result, nullptr);
    }

    // Case 3: valid topology with stride > 1
    {
        std::string traj_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/nonexistent.nc";
        int start_frame = 1;
        int end_frame = 100;
        int stride = 10;

        gemmi::Structure* result = coot::read_amber_trajectory_gemmi(topology_mol, traj_file, start_frame, end_frame, stride);
        EXPECT_EQ(result, nullptr);
    }

    delete topology_mol;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}