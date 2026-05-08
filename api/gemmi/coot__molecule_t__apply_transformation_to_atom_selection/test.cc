#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>
#include "function.hh"

TEST(OracleTest, ApplyTransformationToAtomSelection) {
    // ---- Test Case 1: Identity transformation on chain A atoms ----
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Select all atoms in chain A
        auto atoms_a = coot::molecule_t::parse_selection("//A", st);
        int n_atoms_chain_A = atoms_a.size();
        EXPECT_GT(n_atoms_chain_A, 0);
        
        // Identity transformation
        gemmi::Vec3 rotation_centre(0.0, 0.0, 0.0);
        gemmi::Transform rtop;
        rtop.mat = gemmi::Mat33(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        rtop.vec = gemmi::Vec3(0.0, 0.0, 0.0);
        
        // Store original coords for first atom in residue A/10
        gemmi::Residue* res10 = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10) {
                            res10 = &res;
                            break;
                        }
        
        gemmi::Vec3 original_coord(0.0, 0.0, 0.0);
        if (res10 && !res10->atoms.empty())
            original_coord = res10->atoms[0].pos;
        
        int n_moved = coot::molecule_t::apply_transformation_to_atom_selection_gemmi(
            "//A", n_atoms_chain_A, rotation_centre, rtop, st);
        EXPECT_EQ(n_moved, n_atoms_chain_A);
        
        // Verify coordinates unchanged after identity transform
        gemmi::Residue* res10_after = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10) {
                            res10_after = &res;
                            break;
                        }
        
        if (res10_after && !res10_after->atoms.empty()) {
            gemmi::Vec3 after_coord = res10_after->atoms[0].pos;
            EXPECT_FLOAT_EQ(after_coord.x, original_coord.x);
            EXPECT_FLOAT_EQ(after_coord.y, original_coord.y);
            EXPECT_FLOAT_EQ(after_coord.z, original_coord.z);
        }
    }
    
    // ---- Test Case 2: Apply translation to chain A ----
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Select all atoms in chain A
        auto atoms_a = coot::molecule_t::parse_selection("//A", st);
        int n_atoms_chain_A = atoms_a.size();
        EXPECT_GT(n_atoms_chain_A, 0);
        
        // Apply translation of (10, 0, 0)
        gemmi::Vec3 translation(10.0, 0.0, 0.0);
        gemmi::Transform rtop;
        rtop.mat = gemmi::Mat33(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        rtop.vec = translation;
        gemmi::Vec3 rotation_centre(0.0, 0.0, 0.0);
        
        int n_moved = coot::molecule_t::apply_transformation_to_atom_selection_gemmi(
            "//A", n_atoms_chain_A, rotation_centre, rtop, st);
        EXPECT_EQ(n_moved, n_atoms_chain_A);
        
        // Verify that atom[0] in residue A/10 has shifted by +10 in x
        gemmi::Residue* res10 = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10) {
                            res10 = &res;
                            break;
                        }
        
        if (res10 && !res10->atoms.empty()) {
            gemmi::Vec3 after_coord = res10->atoms[0].pos;
            EXPECT_NEAR(after_coord.x, 36.018, 0.001);
            EXPECT_NEAR(after_coord.y, 12.491, 0.001);
            EXPECT_NEAR(after_coord.z, 68.556, 0.001);
        }
    }
    
    // ---- Test Case 3: Apply transformation to single residue A/10 with identity about custom centre ----
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Select atoms for residue A/10
        auto atoms_res = coot::molecule_t::parse_selection("//A/10", st);
        int n_sel = atoms_res.size();
        EXPECT_EQ(n_sel, 8); // Oracle shows 8 atoms in residue A/10
        
        // Identity transformation about a custom centre
        // For identity rotation and zero translation:
        // new_coord = coord - 2*rotation_centre
        gemmi::Vec3 rotation_centre(10.0, 20.0, 30.0);
        gemmi::Transform rtop;
        rtop.mat = gemmi::Mat33(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        rtop.vec = gemmi::Vec3(0.0, 0.0, 0.0);
        
        int n_moved = coot::molecule_t::apply_transformation_to_atom_selection_gemmi(
            "//A/10", n_sel, rotation_centre, rtop, st);
        EXPECT_EQ(n_moved, n_sel);
        
        // With identity rotation and zero translation, coordinates are transformed as:
        // new_coord = coord - 2*rotation_centre
        gemmi::Residue* res10 = nullptr;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                if (chain.name == "A")
                    for (auto& res : chain.residues)
                        if (res.seqid.num.value == 10) {
                            res10 = &res;
                            break;
                        }
        
        if (res10 && !res10->atoms.empty()) {
            gemmi::Vec3 coord_before = res10->atoms[0].pos;
            // Original coord (26.018, 12.491, 68.556) - 2*(10, 20, 30) = (6.018, -27.509, 8.556)
            EXPECT_NEAR(res10->atoms[0].pos.x, 6.018, 0.01);
            EXPECT_NEAR(res10->atoms[0].pos.y, -27.509, 0.01);
            EXPECT_NEAR(res10->atoms[0].pos.z, 8.556, 0.01);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}