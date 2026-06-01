#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

// Helper: find residue by chain name and sequence number
static gemmi::Residue* find_residue(gemmi::Model& model,
                                     const std::string& chain_name,
                                     int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num)
                    return &res;
            }
        }
    }
    return nullptr;
}

TEST(ModBondAddGemmiTest, valid_atoms_adds_restraint) {
    // Residue A/10 is ILE with atoms: N, CA, C, O, CB, CG1, CG2, CD1
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Residue* res_p = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res_p, nullptr);
    EXPECT_EQ(res_p->name, "ILE");

    restraints_container_t rc;

    // Case: add a bond between N and CA (both exist in ILE A/10)
    chem_mod_bond mod_bond("add", "N", "CA", "", 3.00, 0.05);
    int before = static_cast<int>(rc.restraints_vec.size());
    EXPECT_EQ(before, 0);

    rc.mod_bond_add_gemmi(mod_bond, *res_p);

    int after = static_cast<int>(rc.restraints_vec.size());
    EXPECT_EQ(after, 1);

    // Verify restraint type
    EXPECT_EQ(rc.restraints_vec[0][0], BOND_RESTRAINT);
}

TEST(ModBondAddGemmiTest, invalid_atom_no_restraint) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Residue* res_p = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res_p, nullptr);

    restraints_container_t rc;

    // Case: invalid atom name ZZZ — should not add a restraint
    chem_mod_bond mod_bond("add", "N", "ZZZ", "", 3.00, 0.05);
    int before = static_cast<int>(rc.restraints_vec.size());
    EXPECT_EQ(before, 0);

    rc.mod_bond_add_gemmi(mod_bond, *res_p);

    int after = static_cast<int>(rc.restraints_vec.size());
    EXPECT_EQ(after, 0);
}

TEST(ModBondAddGemmiTest, og1_not_in_ile_no_restraint) {
    // Original test case 1: N and OG1 in residue A/10 (ILE, not THR)
    // OG1 does not exist in ILE, so no restraint should be added
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Residue* res_p = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res_p, nullptr);

    restraints_container_t rc;
    chem_mod_bond mod_bond("add", "N", "OG1", "", 3.00, 0.05);
    EXPECT_EQ(rc.restraints_vec.size(), 0u);

    rc.mod_bond_add_gemmi(mod_bond, *res_p);

    // OG1 doesn't exist in ILE — no restraint added
    EXPECT_EQ(rc.restraints_vec.size(), 0u);
}

TEST(ModBondAddGemmiTest, bonded_atom_indices_populated) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Residue* res_p = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res_p, nullptr);

    restraints_container_t rc;
    chem_mod_bond mod_bond("add", "N", "CA", "", 3.00, 0.05);
    rc.mod_bond_add_gemmi(mod_bond, *res_p);

    // Find indices of N and CA atoms in the residue
    int n_idx = -1, ca_idx = -1;
    for (size_t i = 0; i < res_p->atoms.size(); i++) {
        std::string name = trim_atom_name(res_p->atoms[i].name);
        if (name == "N")  n_idx  = static_cast<int>(i);
        if (name == "CA") ca_idx = static_cast<int>(i);
    }
    ASSERT_NE(n_idx, -1);
    ASSERT_NE(ca_idx, -1);

    // bonded_atom_indices should have each other
    EXPECT_TRUE(rc.bonded_atom_indices[n_idx].count(ca_idx) > 0);
    EXPECT_TRUE(rc.bonded_atom_indices[ca_idx].count(n_idx) > 0);
}

TEST(ModBondAddGemmiTest, different_altlocs_rejected) {
    // Build a residue with atoms that have different non-space altlocs
    gemmi::Residue res;
    res.seqid = gemmi::SeqId(10, ' ');
    res.name = "UNK";

    gemmi::Atom a1;
    a1.name = " N  ";
    a1.altloc = 'A';
    a1.pos = gemmi::Position(0, 0, 0);
    a1.occ = 0.5;
    a1.b_iso = 0;
    res.atoms.push_back(a1);

    gemmi::Atom a2;
    a2.name = "CA ";
    a2.altloc = 'B';
    a2.pos = gemmi::Position(1, 0, 0);
    a2.occ = 0.5;
    a2.b_iso = 0;
    res.atoms.push_back(a2);

    restraints_container_t rc;
    chem_mod_bond mod_bond("add", "N", "CA", "", 3.00, 0.05);
    rc.mod_bond_add_gemmi(mod_bond, res);

    // Different altlocs (A vs B) should NOT create a bond
    EXPECT_EQ(rc.restraints_vec.size(), 0u);
}

TEST(ModBondAddGemmiTest, same_altlocs_accepted) {
    // Atoms with the same non-space altloc SHOULD be bonded
    gemmi::Residue res;
    res.seqid = gemmi::SeqId(10, ' ');
    res.name = "UNK";

    gemmi::Atom a1;
    a1.name = " N  ";
    a1.altloc = 'A';
    a1.pos = gemmi::Position(0, 0, 0);
    a1.occ = 0.5;
    a1.b_iso = 0;
    res.atoms.push_back(a1);

    gemmi::Atom a2;
    a2.name = "CA ";
    a2.altloc = 'A';
    a2.pos = gemmi::Position(1, 0, 0);
    a2.occ = 0.5;
    a2.b_iso = 0;
    res.atoms.push_back(a2);

    restraints_container_t rc;
    chem_mod_bond mod_bond("add", "N", "CA", "", 3.00, 0.05);
    rc.mod_bond_add_gemmi(mod_bond, res);

    // Same altloc (A) should create a bond
    EXPECT_EQ(rc.restraints_vec.size(), 1u);
}

TEST(ModBondAddGemmiTest, one_no_altloc_accepted) {
    // One atom with no altloc (' ') and one with altloc — should still bond
    gemmi::Residue res;
    res.seqid = gemmi::SeqId(10, ' ');
    res.name = "UNK";

    gemmi::Atom a1;
    a1.name = " N  ";
    a1.altloc = ' ';
    a1.pos = gemmi::Position(0, 0, 0);
    a1.occ = 1.0;
    a1.b_iso = 0;
    res.atoms.push_back(a1);

    gemmi::Atom a2;
    a2.name = "CA ";
    a2.altloc = 'A';
    a2.pos = gemmi::Position(1, 0, 0);
    a2.occ = 0.5;
    a2.b_iso = 0;
    res.atoms.push_back(a2);

    restraints_container_t rc;
    chem_mod_bond mod_bond("add", "N", "CA", "", 3.00, 0.05);
    rc.mod_bond_add_gemmi(mod_bond, res);

    EXPECT_EQ(rc.restraints_vec.size(), 1u);
}

} // namespace coot

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}