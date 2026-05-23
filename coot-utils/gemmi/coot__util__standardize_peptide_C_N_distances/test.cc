#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <fstream>
#include <iomanip>
#include <cstring>

// Helper to create a temporary PDB file with specific C-N distances
static std::string create_test_pdb(const std::string& path, const std::vector<double>& initial_distances) {
    std::ofstream ofs(path);
    
    for (size_t i = 0; i < initial_distances.size(); i++) {
        double dist = initial_distances[i];
        
        // Residue 1
        ofs << "ATOM  " << std::setw(5) << (i*10+1) << "  CA  ALA A" << std::setw(4) << (i*2+1) 
            << "    " << std::fixed << std::setprecision(3)
            << std::setw(8) << 0.000
            << std::setw(8) << 0.000
            << std::setw(8) << 0.000 << "  1.00  0.00           C" << std::endl;
        
        // C atom of residue 1 (at standard position relative to CA)
        ofs << "ATOM  " << std::setw(5) << (i*10+2) << "  C   ALA A" << std::setw(4) << (i*2+1) 
            << "    " << std::fixed << std::setprecision(3)
            << std::setw(8) << 1.525
            << std::setw(8) << 0.000
            << std::setw(8) << 0.000 << "  1.00  0.00           C" << std::endl;
        
        // Residue 2 (N atom at specified distance from C atom of residue 1)
        ofs << "ATOM  " << std::setw(5) << (i*10+3) << "  N   ALA A" << std::setw(4) << (i*2+2)
            << "    " << std::fixed << std::setprecision(3)
            << std::setw(8) << (1.525 + dist)  // N position = C position + distance
            << std::setw(8) << 0.000
            << std::setw(8) << 0.000 << "  1.00  0.00           N" << std::endl;
        
        // CA of residue 2
        ofs << "ATOM  " << std::setw(5) << (i*10+4) << "  CA  ALA A" << std::setw(4) << (i*2+2)
            << "    " << std::fixed << std::setprecision(3)
            << std::setw(8) << (1.525 + dist + 1.525)
            << std::setw(8) << 0.000
            << std::setw(8) << 0.000 << "  1.00  0.00           C" << std::endl;
    }
    ofs.close();
    return path;
}

static std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> build_c_n_pairs(gemmi::Model& model) {
    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> c_n_pairs;
    
    for (gemmi::Chain& chain : model.chains) {
        for (size_t ires = 0; ires < chain.residues.size() - 1; ires++) {
            gemmi::Residue& res1 = chain.residues[ires];
            gemmi::Residue& res2 = chain.residues[ires + 1];
            
            // Skip water molecules
            if (res1.name == "HOH" || res2.name == "HOH")
                continue;
            
            // Get C atom from first residue (exactly "C")
            gemmi::Atom* c_atom = nullptr;
            for (gemmi::Atom& atom : res1.atoms) {
                if (atom.name == "C") {
                    c_atom = &atom;
                    break;
                }
            }
            
            // Get N atom from second residue (exactly "N")
            gemmi::Atom* n_atom = nullptr;
            for (gemmi::Atom& atom : res2.atoms) {
                if (atom.name == "N") {
                    n_atom = &atom;
                    break;
                }
            }
            
            if (c_atom && n_atom) {
                c_n_pairs.push_back(std::make_pair(c_atom, n_atom));
            }
        }
    }
    return c_n_pairs;
}

static double get_c_n_distance(gemmi::Atom* c_at, gemmi::Atom* n_at) {
    clipper::Coord_orth pt_1 = coot::co_gemmi(c_at);
    clipper::Coord_orth pt_2 = coot::co_gemmi(n_at);
    return clipper::Coord_orth::length(pt_1, pt_2);
}

TEST(OracleTest, StandardizePeptideCNDistances) {
    std::string pdb_path = "/tmp/test_peptide.pdb";
    
    // Case 1: Long C-N distances (should be standardized to ~1.320)
    std::vector<double> long_distances = {2.442, 2.434, 2.463, 2.412, 2.435};
    create_test_pdb(pdb_path, long_distances);
    
    gemmi::Structure st1 = gemmi::read_pdb_file(pdb_path);
    ASSERT_GE(st1.models.size(), 1) << "Failed to load PDB file";
    
    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> C_N_pairs = build_c_n_pairs(st1.models[0]);
    
    ASSERT_EQ(C_N_pairs.size(), 5) << "Expected 5 C-N pairs";
    
    // Get distances before
    std::vector<double> distances_before;
    for (size_t i = 0; i < C_N_pairs.size(); i++) {
        distances_before.push_back(get_c_n_distance(C_N_pairs[i].first, C_N_pairs[i].second));
    }
    
    // Call the function
    coot::util::standardize_peptide_C_N_distances_gemmi(C_N_pairs);
    
    // Get distances after
    std::vector<double> distances_after;
    for (size_t i = 0; i < C_N_pairs.size(); i++) {
        distances_after.push_back(get_c_n_distance(C_N_pairs[i].first, C_N_pairs[i].second));
    }
    
    // Check that all long distances were standardized to ~1.320
    for (size_t i = 0; i < distances_after.size(); i++) {
        EXPECT_NEAR(distances_after[i], 1.320, 0.01) 
            << "Pair " << i << " should be standardized to 1.320, got " << distances_after[i];
    }
    
    // Case 2: Already short C-N distances (should remain unchanged at ~1.320)
    std::vector<double> short_distances = {1.320, 1.320, 1.320};
    create_test_pdb(pdb_path, short_distances);
    
    gemmi::Structure st2 = gemmi::read_pdb_file(pdb_path);
    ASSERT_GE(st2.models.size(), 1) << "Failed to load second PDB file";
    
    std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> C_N_pairs2 = build_c_n_pairs(st2.models[0]);
    
    // Get distances before (should already be ~1.320)
    std::vector<double> distances_before2;
    for (size_t i = 0; i < C_N_pairs2.size(); i++) {
        distances_before2.push_back(get_c_n_distance(C_N_pairs2[i].first, C_N_pairs2[i].second));
    }
    
    // Call the function
    coot::util::standardize_peptide_C_N_distances_gemmi(C_N_pairs2);
    
    // Get distances after
    std::vector<double> distances_after2;
    for (size_t i = 0; i < C_N_pairs2.size(); i++) {
        distances_after2.push_back(get_c_n_distance(C_N_pairs2[i].first, C_N_pairs2[i].second));
    }
    
    // Check that short distances remain at ~1.320
    for (size_t i = 0; i < distances_after2.size(); i++) {
        EXPECT_NEAR(distances_after2[i], 1.320, 0.01)
            << "Pair " << i << " should remain at 1.320, got " << distances_after2[i];
    }
    
    // Clean up
    std::remove(pdb_path.c_str());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}