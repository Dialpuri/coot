#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// case: MET residue 10 in chain A
TEST(OracleTest, PrimitiveChiAnglesMET) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get model and chain
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    
    // Find residue with seqid 10
    gemmi::Residue* res = nullptr;
    for (auto& r : chain.residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    
    ASSERT_NE(res, nullptr) << "Residue 10 not found";
    
    // Verify residue info (gemmi uses unpadded atom names)
    EXPECT_EQ(chain.name, "A");
    EXPECT_EQ(res->seqid.num.value, 10);
    
    // Verify atom count and names (gemmi uses unpadded names)
    size_t n_atoms = res->atoms.size();
    EXPECT_EQ(n_atoms, 8);
    
    // Verify atom names (gemmi uses unpadded names like "N", "CA", etc.)
    std::vector<std::string> expected_names = {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"};
    for (size_t i = 0; i < n_atoms; ++i) {
        EXPECT_EQ(res->atoms[i].name, expected_names[i]);
    }
    
    // Construct atom_name_quads (using unpadded names for gemmi)
    std::vector<coot::primitive_chi_angles::atom_name_quad> atom_name_quads;
    atom_name_quads.push_back(coot::primitive_chi_angles::atom_name_quad("N", "CA", "CB", "CG1"));
    atom_name_quads.push_back(coot::primitive_chi_angles::atom_name_quad("N", "CA", "CB", "CG2"));
    
    // Call get_quads_using_altconfs
    std::vector<coot::primitive_chi_angles::alt_confed_atom_index_quad> result =
        coot::primitive_chi_angles::get_quads_using_altconfs(atom_name_quads, *res);
    
    // MET residue 10 has no alt confs, so result should be empty
    EXPECT_EQ(result.size(), 0u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}