#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <map>
#include "function.hh"

TEST(OracleTest, invert_deloc_bond_displacement_vector) {
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Get model 0 (gemmi is 0-indexed; MMDB's model 1)
    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    // Find chain A
    const gemmi::Chain *chain = nullptr;
    for (auto &ch : model->chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Find residue with seqid 4 (PHE 4)
    const gemmi::Residue *res = nullptr;
    for (auto &r : chain->residues) {
        if (r.seqid.num.value == 4) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr);

    // Build atom table as vector of pointers (gemmi equivalent of mmdb::PPAtom)
    std::vector<const gemmi::Atom*> atom_table;
    for (auto &a : res->atoms) {
        atom_table.push_back(&a);
    }
    int n_atoms = static_cast<int>(atom_table.size());

    ASSERT_GE(n_atoms, 3);

    // Build atom name to index map using atom_id_mmdb_expand to match _4c() format
    std::map<std::string, int> atom_name_map;
    for (int i = 0; i < n_atoms; i++) {
        std::string raw = atom_table[i]->name;
        size_t start = raw.find_first_not_of(" ");
        size_t end = raw.find_last_not_of(" ");
        std::string trimmed;
        if (start != std::string::npos)
            trimmed = raw.substr(start, end - start + 1);
        atom_name_map[coot::atom_id_mmdb_expand(trimmed)] = i;
    }

    // Pick iat_1 = index of "CA", iat_2 = index of "C"
    int iat_1 = -1, iat_2 = -1;
    auto it_ca = atom_name_map.find(coot::atom_id_mmdb_expand("CA"));
    auto it_c = atom_name_map.find(coot::atom_id_mmdb_expand("C"));
    if (it_ca != atom_name_map.end()) iat_1 = it_ca->second;
    if (it_c != atom_name_map.end()) iat_2 = it_c->second;

    ASSERT_NE(iat_1, -1);
    ASSERT_NE(iat_2, -1);

    // Create a default Bond_lines_container (the method doesn't use member state)
    Bond_lines_container bonds;

    // ---------------------------------------------------------------
    // Case 1: deloc bond restraint, vect=(1,0,0)
    //   dot(vect, diff) = -1.045 < 0 → result=true
    // ---------------------------------------------------------------
    {
        clipper::Coord_orth vect(1, 0, 0);

        coot::dict_bond_restraint_t br("CA", "N", "deloc");

        std::vector<coot::dict_bond_restraint_t> restraints;
        restraints.push_back(br);

        bool result = bonds.invert_deloc_bond_displacement_vector_gemmi(
            vect, iat_1, iat_2, atom_table, n_atoms, restraints);

        EXPECT_TRUE(result);
    }

    // ---------------------------------------------------------------
    // Case 2: deloc bond restraint, vect=(-1,0,0)
    //   dot(vect, diff) = 1.045 >= 0 → result=false
    // ---------------------------------------------------------------
    {
        clipper::Coord_orth vect(-1, 0, 0);

        coot::dict_bond_restraint_t br("CA", "N", "deloc");

        std::vector<coot::dict_bond_restraint_t> restraints;
        restraints.push_back(br);

        bool result = bonds.invert_deloc_bond_displacement_vector_gemmi(
            vect, iat_1, iat_2, atom_table, n_atoms, restraints);

        EXPECT_FALSE(result);
    }

    // ---------------------------------------------------------------
    // Case 3: no "deloc" restraints → false (default)
    // ---------------------------------------------------------------
    {
        clipper::Coord_orth vect(1, 0, 0);

        coot::dict_bond_restraint_t br("CA", "N", "single");

        std::vector<coot::dict_bond_restraint_t> restraints;
        restraints.push_back(br);

        bool result = bonds.invert_deloc_bond_displacement_vector_gemmi(
            vect, iat_1, iat_2, atom_table, n_atoms, restraints);

        EXPECT_FALSE(result);
    }

    // ---------------------------------------------------------------
    // Case 4: empty restraints → false (default)
    // ---------------------------------------------------------------
    {
        clipper::Coord_orth vect(1, 0, 0);
        std::vector<coot::dict_bond_restraint_t> restraints; // empty

        bool result = bonds.invert_deloc_bond_displacement_vector_gemmi(
            vect, iat_1, iat_2, atom_table, n_atoms, restraints);

        EXPECT_FALSE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}