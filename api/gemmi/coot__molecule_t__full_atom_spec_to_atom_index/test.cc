#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/residue-and-atom-specs.hh"

TEST(OracleTest, full_atom_spec_to_atom_index) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid atom — CA of residue 1 in chain A
    {
        std::string chain = "A";
        int resno = 1;
        std::string ins_code = "";
        std::string atom_name = "CA";
        std::string alt_conf = "";

        coot::atom_spec_t spec;
        spec.chain_id = chain;
        spec.res_no = resno;
        spec.ins_code = ins_code;
        spec.atom_name = atom_name;
        spec.alt_conf = alt_conf;

        int idx = coot::molecule_t::full_atom_spec_to_atom_index_gemmi(spec, st);

        EXPECT_EQ(idx, 1);
    }

    // Case 2: non-existent atom — atom "XYZ" of residue 1 in chain A
    {
        std::string chain = "A";
        int resno = 1;
        std::string ins_code = "";
        std::string atom_name = "XYZ";
        std::string alt_conf = "";

        coot::atom_spec_t spec;
        spec.chain_id = chain;
        spec.res_no = resno;
        spec.ins_code = ins_code;
        spec.atom_name = atom_name;
        spec.alt_conf = alt_conf;

        int idx = coot::molecule_t::full_atom_spec_to_atom_index_gemmi(spec, st);

        EXPECT_EQ(idx, -1);
    }

    // Case 3: via atom_spec_t overload — N of residue 50 (ARG) in chain A
    {
        coot::atom_spec_t spec;
        spec.chain_id = "A";
        spec.res_no = 50;
        spec.ins_code = "";
        spec.atom_name = "N";
        spec.alt_conf = "";

        int idx = coot::molecule_t::full_atom_spec_to_atom_index_gemmi(spec, st);

        EXPECT_EQ(idx, 309);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}