#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

struct dummy_geometry {
  struct residue_info_t { std::string group; };
  struct restraints_t { residue_info_t residue_info; };
  std::pair<int, restraints_t> get_monomer_restraints(const std::string& name, int) const {
    restraints_t r;
    if (name == "A" || name == "DA" || name == "G" || name == "DG" ||
        name == "C" || name == "DC" || name == "U" || name == "DT") {
      r.residue_info.group = (name[0] == 'D') ? "DNA" : "RNA";
    } else {
      r.residue_info.group = "unknown";
    }
    return std::make_pair(0, r);
  }
};

TEST(OracleTest, is_nucleotide_by_dict) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    dummy_geometry geom;

    // case: protein residue (ILE) should return false
    {
        bool found = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == 10 && chain.name == "A") {
                        bool result = coot::util::is_nucleotide_by_dict_gemmi(res, "ILE");
                        EXPECT_FALSE(result);
                        found = true;
                    }
        EXPECT_TRUE(found) << "Residue ILE at chain A, seq 10 not found";
    }

    // case: water residue (HOH) should return false
    {
        bool found = false;
        for (auto& model : st.models)
            for (auto& chain : model.chains)
                for (auto& res : chain.residues)
                    if (res.seqid.num.value == 1 && chain.name == "B") {
                        bool result = coot::util::is_nucleotide_by_dict_gemmi(res, "HOH");
                        EXPECT_FALSE(result);
                        found = true;
                    }
        EXPECT_TRUE(found) << "Residue HOH at chain B, seq 1 not found";
    }

    // case: RNA nucleotide (A) should return true
    {
        gemmi::Residue res;
        res.name = "A";
        bool result = coot::util::is_nucleotide_by_dict_gemmi(res, "A");
        EXPECT_TRUE(result);
    }

    // case: DNA nucleotide (DA) should return true
    {
        gemmi::Residue res;
        res.name = "DA";
        bool result = coot::util::is_nucleotide_by_dict_gemmi(res, "DA");
        EXPECT_TRUE(result);
    }

    // case: protein residue (ALA) should return false
    {
        gemmi::Residue res;
        res.name = "ALA";
        bool result = coot::util::is_nucleotide_by_dict_gemmi(res, "ALA");
        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}