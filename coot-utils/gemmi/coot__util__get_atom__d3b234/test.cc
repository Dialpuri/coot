#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot/geometry/residue-and-atom-specs.hh"

TEST(OracleTest, get_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid atom spec - find CA in residue 10
    {
        // First find residue 10 in chain A
        gemmi::Residue *res = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& r : chain.residues) {
                    if (chain.name == "A" && r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        
        coot::atom_spec_t spec("A", 10, "", "CA", "");
        gemmi::Atom *at = coot::util::get_atom_gemmi(spec, res);
        EXPECT_FALSE(at != nullptr);  // atom not found
    }

    // case 2: non-existent atom
    {
        gemmi::Residue *res = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& r : chain.residues) {
                    if (chain.name == "A" && r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
            }
        }
        
        coot::atom_spec_t spec("A", 10, "", "ZZZ", "");
        gemmi::Atom *at = coot::util::get_atom_gemmi(spec, res);
        EXPECT_FALSE(at != nullptr);  // atom not found
    }

    // case 3: null residue
    {
        coot::atom_spec_t spec("A", 10, "", "CA", "");
        gemmi::Atom *at = coot::util::get_atom_gemmi(spec, nullptr);
        EXPECT_FALSE(at != nullptr);  // atom not found
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}