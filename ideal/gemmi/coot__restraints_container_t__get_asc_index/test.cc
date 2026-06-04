#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, restraints_container_get_asc_index) {
    // Setup: load PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Collect CRA from chain A, residues 1, 2, 3 (matching original test setup)
    std::vector<gemmi::CRA> cras;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue& residue : chain.residues) {
            int seq_num = residue.seqid.num.value;
            if (seq_num < 1 || seq_num > 3) continue;
            for (gemmi::Atom& atom : residue.atoms) {
                cras.push_back({&chain, &residue, &atom});
            }
        }
    }

    // Case 1: valid atom spec - CA of residue 1, chain A
    {
        coot::atom_spec_t spec;
        spec.atom_name = "CA";
        spec.alt_conf = "";
        spec.res_no = 1;
        spec.ins_code = "";
        spec.chain_id = "A";

        int idx = coot::get_asc_index_gemmi(cras, spec);
        EXPECT_EQ(idx, 1);
    }

    // Case 2: valid atom spec - N of residue 2, chain A
    {
        coot::atom_spec_t spec;
        spec.atom_name = "N";
        spec.alt_conf = "";
        spec.res_no = 2;
        spec.ins_code = "";
        spec.chain_id = "A";

        int idx = coot::get_asc_index_gemmi(cras, spec);
        EXPECT_EQ(idx, 8);
    }

    // Case 3: non-existent atom - invalid residue number
    {
        coot::atom_spec_t spec;
        spec.atom_name = "CA";
        spec.alt_conf = "";
        spec.res_no = 9999;
        spec.ins_code = "";
        spec.chain_id = "A";

        int idx = coot::get_asc_index_gemmi(cras, spec);
        EXPECT_EQ(idx, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}