#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

using namespace std;

TEST(OracleTest, get_atom_index_pairs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    const gemmi::Model& model = st.models[0];
    const gemmi::Chain& chain_a = model.chains[0];
    const gemmi::Residue& res = chain_a.residues[9];
    
    EXPECT_EQ(chain_a.name, "A");
    EXPECT_EQ(res.seqid.num.value, 10);
    EXPECT_EQ(res.name, "ILE");
    
    const auto& atoms = res.atoms;
    EXPECT_EQ(atoms.size(), 8);
    
    {
        std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs;
        atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("N", "CA"));
        atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CA", "CB"));
        atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CB", "CG1"));

        std::vector<coot::chi_angles::atom_index_pair> index_pairs = 
            coot::chi_angles::get_atom_index_pairs_gemmi(atom_name_pairs, res);

        EXPECT_EQ(index_pairs.size(), 3);
        EXPECT_EQ(index_pairs[0].index1, 0);
        EXPECT_EQ(index_pairs[0].index2, 1);
        EXPECT_EQ(index_pairs[1].index1, 1);
        EXPECT_EQ(index_pairs[1].index2, 4);
        EXPECT_EQ(index_pairs[2].index1, 4);
        EXPECT_EQ(index_pairs[2].index2, 5);
    }

    {
        std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs2;
        atom_name_pairs2.push_back(coot::chi_angles::atom_name_pair("ZZZ", "CA"));

        std::vector<coot::chi_angles::atom_index_pair> index_pairs2 = 
            coot::chi_angles::get_atom_index_pairs_gemmi(atom_name_pairs2, res);

        EXPECT_EQ(index_pairs2.size(), 0);
    }

    {
        std::vector<coot::chi_angles::atom_name_pair> empty_pairs;
        std::vector<coot::chi_angles::atom_index_pair> empty_result = 
            coot::chi_angles::get_atom_index_pairs_gemmi(empty_pairs, res);

        EXPECT_EQ(empty_result.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}