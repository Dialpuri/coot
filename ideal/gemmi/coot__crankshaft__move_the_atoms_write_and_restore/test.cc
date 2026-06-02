#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "function.hh"

TEST(OracleTest, move_the_atoms_write_and_restore_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Get chain A (first chain)
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];

    size_t nResidues = chain.residues.size();
    EXPECT_EQ(nResidues, 267u);

    // Collect 9 consecutive residues starting from index 10
    std::vector<gemmi::Residue*> resids;
    for (int i = 10; i < 10 + 9 && i < static_cast<int>(nResidues); i++) {
        resids.push_back(&chain.residues[i]);
    }

    ASSERT_EQ(resids.size(), 9u);

    // Build atom_sets mirroring crankshaft_set::v
    // Original crankshaft_set v contains ALL atoms from 4 residues concatenated.
    // Each crankshaft has 4 residues, sas has 3 crankshafts.
    // We need atom_sets[i] to have at least 6 atoms so indices {2,3,4,5} are valid.
    std::vector<std::vector<gemmi::Atom*>> atom_sets(3);
    for (int cs_idx = 0; cs_idx < 3; cs_idx++) {
        // Each crankshaft covers 4 residues: resids[cs_idx*4 .. cs_idx*4+3]
        for (int r_idx = 0; r_idx < 4; r_idx++) {
            int ri = cs_idx * 4 + r_idx;
            if (ri < static_cast<int>(resids.size())) {
                for (auto& atom : resids[ri]->atoms) {
                    atom_sets[cs_idx].push_back(&atom);
                }
            }
        }
    }

    // Verify we have enough atoms per set
    for (std::size_t i = 0; i < 3; i++) {
        ASSERT_GT(atom_sets[i].size(), 5u);
    }

    // Angles matching original test
    std::vector<float> angles(3);
    angles[0] = 0.1f;
    angles[1] = 0.2f;
    angles[2] = 0.3f;

    EXPECT_EQ(angles.size(), 3u);
    EXPECT_NEAR(angles[0], 0.1f, 1e-4);
    EXPECT_NEAR(angles[1], 0.2f, 1e-4);
    EXPECT_NEAR(angles[2], 0.3f, 1e-4);

    // Record original positions of atoms that will be moved
    int indices[] = {2, 3, 4, 5};
    std::map<void*, std::string> original_pos_str;
    for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t iat = 0; iat < 4; iat++) {
            gemmi::Atom* at = atom_sets[i][indices[iat]];
            if (at) {
                clipper::Coord_orth pos = coot::co_gemmi(at);
                std::ostringstream oss;
                oss << pos[0] << " " << pos[1] << " " << pos[2];
                original_pos_str[static_cast<void*>(at)] = oss.str();
            }
        }
    }

    ASSERT_FALSE(original_pos_str.empty());

    std::string out_pdb = "/tmp/crankshaft_test_output.pdb";

    EXPECT_NO_THROW(coot::crankshaft::move_the_atoms_write_and_restore_gemmi(
        st, atom_sets, angles, out_pdb));

    // Check output PDB written
    std::ifstream outfile(out_pdb);
    bool file_exists = outfile.good();
    EXPECT_TRUE(file_exists);

    // Verify positions restored
    bool all_restored = true;
    std::map<void*, std::string> final_pos_str;
    for (auto const &kv : original_pos_str) {
        gemmi::Atom* at = static_cast<gemmi::Atom*>(kv.first);
        clipper::Coord_orth pos = coot::co_gemmi(at);
        std::ostringstream oss;
        oss << pos[0] << " " << pos[1] << " " << pos[2];
        final_pos_str[at] = oss.str();
        if (kv.second != oss.str()) all_restored = false;
    }

    EXPECT_TRUE(all_restored);

    // Sample atom info
    if (!original_pos_str.empty()) {
        auto it = original_pos_str.begin();
        gemmi::Atom* sa = static_cast<gemmi::Atom*>(it->first);
        std::string sample_name = sa->name;
        std::string sample_orig = original_pos_str[sa];
        std::string sample_final = final_pos_str[sa];

        EXPECT_FALSE(sample_name.empty());
        EXPECT_FALSE(sample_orig.empty());
        EXPECT_EQ(sample_orig, sample_final);

        std::cerr << "INFO sample_atom_name: " << sample_name << std::endl;
        std::cerr << "INFO sample_original_pos: " << sample_orig << std::endl;
        std::cerr << "INFO sample_final_pos: " << sample_final << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}