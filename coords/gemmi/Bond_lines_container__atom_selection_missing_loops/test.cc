#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MissingLoopsBondLines) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    EXPECT_GT(st.models.size(), 0);
    
    gemmi::Model& model = st.models[0];
    EXPECT_EQ(model.chains.size(), 2);
    
    int max_gap = 0;
    for (const auto& chain : model.chains) {
        if (chain.residues.size() < 2) continue;
        for (size_t i = 1; i < chain.residues.size(); ++i) {
            int res_no_prev = chain.residues[i-1].seqid.num.value;
            int res_no_this = chain.residues[i].seqid.num.value;
            int delta = res_no_this - res_no_prev;
            if (delta > 1) {
                if (delta > max_gap) max_gap = delta;
            }
        }
    }
    EXPECT_EQ(max_gap, 20);
    
    size_t n_atoms = 0;
    for (const auto& chain : model.chains) {
        for (const auto& res : chain.residues) {
            n_atoms += res.atoms.size();
        }
    }
    EXPECT_GT(n_atoms, 0);
    
    coot::atom_selection_container_t asc;
    asc.structure = &st;
    asc.model = &model;
    asc.n_selected_atoms = static_cast<int>(n_atoms);
    asc.UDDAtomIndexHandle = 0;
    asc.mol = 0;
    
    int udd_atom_index_handle = 0;
    int udd_fixed_during_refinement_handle = 0;
    
    coot::Bond_lines_container blc;
    blc.atom_selection_missing_loops_gemmi(asc, udd_atom_index_handle, udd_fixed_during_refinement_handle);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}