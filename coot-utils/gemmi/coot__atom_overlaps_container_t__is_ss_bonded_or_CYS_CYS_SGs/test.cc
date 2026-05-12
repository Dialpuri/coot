#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_ss_bonded_or_CYS_CYS_SGs) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    gemmi::Atom* sg_118 = nullptr;
    gemmi::Residue* res_118 = nullptr;
    for (gemmi::Model& model : st.models)
    for (gemmi::Chain& chain : model.chains)
    for (gemmi::Residue& res : chain.residues) {
        if (res.name == "CYS" && res.seqid.num.value == 118) {
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.name == "SG") {
                    sg_118 = &atom;
                    res_118 = &res;
                    break;
                }
            }
        }
    }
    
    gemmi::Atom* sg_177 = nullptr;
    gemmi::Residue* res_177 = nullptr;
    for (gemmi::Model& model : st.models)
    for (gemmi::Chain& chain : model.chains)
    for (gemmi::Residue& res : chain.residues) {
        if (res.name == "CYS" && res.seqid.num.value == 177) {
            for (gemmi::Atom& atom : res.atoms) {
                if (atom.name == "SG") {
                    sg_177 = &atom;
                    res_177 = &res;
                    break;
                }
            }
        }
    }
    
    coot::atom_overlaps_container_t container;
    
    if (sg_118 && sg_177 && res_118 && res_177) {
        bool result = container.is_ss_bonded_or_CYS_CYS_SGs_gemmi(*sg_118, *res_118, *sg_177, *res_177);
        EXPECT_TRUE(result);
    } else {
        fprintf(stderr, "SG atoms not found for test\n");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}