#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue in the first model by chain name and seqid
static gemmi::Residue* find_residue(gemmi::Model& model, char chain_id, int seqnum) {
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name != std::string(1, chain_id)) continue;
        for (gemmi::Residue& res : ch.residues) {
            if (res.seqid.num.value == seqnum)
                return &res;
        }
    }
    return nullptr;
}

// Pad atom name to MMDB GetAtomName style: 1 space + name, right-padded to 4
static std::string pad4(const std::string& s) {
    std::string r = " " + s;
    while (r.size() < 4) r += ' ';
    return r;
}

TEST(OracleTest, ordered_residue_atoms) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  EXPECT_FALSE(st.models.empty());
  gemmi::Model& model = st.models[0];

  // Case 1: VAL (has CG1/CG2 ordering) - residue //A/7
  {
    gemmi::Residue* res = find_residue(model, 'A', 7);
    EXPECT_NE(res, nullptr);
    EXPECT_EQ(res->name, "VAL");

    std::vector<gemmi::Atom*> ordered = coot::rotamer::ordered_residue_atoms_gemmi(*res);

    EXPECT_EQ(ordered.size(), 7u);

    std::string atom_names;
    for (unsigned int i = 0; i < ordered.size(); i++) {
      if (i > 0) atom_names += ", ";
      atom_names += pad4(ordered[i]->name);
    }
    EXPECT_EQ(atom_names, " N  ,  CA ,  C  ,  O  ,  CB ,  CG1,  CG2");
  }

  // Case 2: SER (no special ordering) - residue //A/46
  {
    gemmi::Residue* res = find_residue(model, 'A', 46);
    EXPECT_NE(res, nullptr);
    EXPECT_EQ(res->name, "SER");

    std::vector<gemmi::Atom*> ordered = coot::rotamer::ordered_residue_atoms_gemmi(*res);

    EXPECT_EQ(ordered.size(), 6u);

    std::string atom_names;
    for (unsigned int i = 0; i < ordered.size(); i++) {
      if (i > 0) atom_names += ", ";
      atom_names += pad4(ordered[i]->name);
    }
    EXPECT_EQ(atom_names, " N  ,  CA ,  C  ,  O  ,  CB ,  OG ");
  }

  // Case 3: PRO (more than 2 atoms to order) - residue //A/61
  {
    gemmi::Residue* res = find_residue(model, 'A', 61);
    EXPECT_NE(res, nullptr);
    EXPECT_EQ(res->name, "PRO");

    std::vector<gemmi::Atom*> ordered = coot::rotamer::ordered_residue_atoms_gemmi(*res);

    EXPECT_EQ(ordered.size(), 7u);

    std::string atom_names;
    for (unsigned int i = 0; i < ordered.size(); i++) {
      if (i > 0) atom_names += ", ";
      atom_names += pad4(ordered[i]->name);
    }
    EXPECT_EQ(atom_names, " N  ,  CA ,  CB ,  CG ,  CD ,  C  ,  O  ");
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}