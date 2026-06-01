#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_atom/gemmi/function.cc"

// Helper: find the Chain* that contains the given Residue*
static gemmi::Chain* find_chain_of_residue(gemmi::Model &model, gemmi::Residue *res) {
  for (auto &chain : model.chains) {
    for (auto &r : chain.residues) {
      if (&r == res) return &chain;
    }
  }
  return nullptr;
}

// Helper: find the Chain* that contains the given Atom*
static gemmi::Chain* find_chain_of_atom(gemmi::Model &model, gemmi::Atom *atom) {
  for (auto &chain : model.chains) {
    for (auto &res : chain.residues) {
      for (auto &a : res.atoms) {
        if (&a == atom) return &chain;
      }
    }
  }
  return nullptr;
}

// Helper: find the Residue* that contains the given Atom*
static gemmi::Residue* find_residue_of_atom(gemmi::Model &model, gemmi::Atom *atom) {
  for (auto &chain : model.chains) {
    for (auto &res : chain.residues) {
      for (auto &a : res.atoms) {
        if (&a == atom) return &res;
      }
    }
  }
  return nullptr;
}

TEST(OracleTest, atom_is_neighbour_mainchain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Get a reference residue: //A/10 (ILE)
    gemmi::Residue *ref_res = coot::cid_to_residue_gemmi("//A/10", model);
    EXPECT_NE(ref_res, nullptr);
    gemmi::Chain *ref_chain = find_chain_of_residue(model, ref_res);
    EXPECT_NE(ref_chain, nullptr);

    gemmi::CRA ref_cra{ref_chain, ref_res, nullptr};

    // Case 1: same residue, mainchain atom CA -> result: true
    {
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi("//A/10/CA", st);
        EXPECT_NE(at, nullptr);
        gemmi::Chain *at_chain = find_chain_of_atom(model, at);
        gemmi::Residue *at_res = find_residue_of_atom(model, at);
        gemmi::CRA atom_cra{at_chain, at_res, at};
        bool result = coot::daca::atom_is_neighbour_mainchain_gemmi(atom_cra, ref_cra);
        EXPECT_EQ(result, true);
    }

    // Case 2: adjacent residue (res+1), mainchain atom N -> result: true
    {
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi("//A/11/N", st);
        EXPECT_NE(at, nullptr);
        gemmi::Chain *at_chain = find_chain_of_atom(model, at);
        gemmi::Residue *at_res = find_residue_of_atom(model, at);
        gemmi::CRA atom_cra{at_chain, at_res, at};
        bool result = coot::daca::atom_is_neighbour_mainchain_gemmi(atom_cra, ref_cra);
        EXPECT_EQ(result, true);
    }

    // Case 3: same residue, sidechain atom CB -> result: false
    {
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi("//A/10/CB", st);
        EXPECT_NE(at, nullptr);
        gemmi::Chain *at_chain = find_chain_of_atom(model, at);
        gemmi::Residue *at_res = find_residue_of_atom(model, at);
        gemmi::CRA atom_cra{at_chain, at_res, at};
        bool result = coot::daca::atom_is_neighbour_mainchain_gemmi(atom_cra, ref_cra);
        EXPECT_EQ(result, false);
    }

    // Case 4: far residue (res+10), mainchain atom CA -> result: false
    {
        gemmi::Atom *at = coot::molecule_t::cid_to_atom_gemmi("//A/20/CA", st);
        EXPECT_NE(at, nullptr);
        gemmi::Chain *at_chain = find_chain_of_atom(model, at);
        gemmi::Residue *at_res = find_residue_of_atom(model, at);
        gemmi::CRA atom_cra{at_chain, at_res, at};
        bool result = coot::daca::atom_is_neighbour_mainchain_gemmi(atom_cra, ref_cra);
        EXPECT_EQ(result, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}