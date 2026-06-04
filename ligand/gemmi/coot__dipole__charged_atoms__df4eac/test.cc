#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string pad_4c(const std::string &n) {
    if (n.size() == 1) return " " + n + "  ";
    if (n.size() == 2) return " " + n + " ";
    if (n.size() == 3) return n + " ";
    if (n.size() >= 4) return n.substr(0, 4);
    return "    ";
}

TEST(OracleTest, charged_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get first model and first residue
    ASSERT_GE(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];
    ASSERT_GE(model.chains.size(), 1u);
    gemmi::Chain &chain = model.chains[0];
    ASSERT_GE(chain.residues.size(), 1u);
    gemmi::Residue &res = chain.residues[0];

    ASSERT_EQ(chain.name, "A");
    ASSERT_EQ(res.seqid.num.value, 1);
    ASSERT_EQ(res.name, "MET");

    // Get atoms from residue
    auto &atoms = res.atoms;

    // Case 1: normal — all atoms have partial charges
    {
        coot::dictionary_residue_restraints_t rest;

        for (size_t i = 0; i < atoms.size(); i++) {
            coot::dict_atom da;
            da.atom_id = atoms[i].name;
            da.atom_id_4c = atoms[i].name;
            da.partial_charge = std::make_pair(true, 0.1f * static_cast<float>(i + 1));
            rest.atom_info.push_back(da);
        }

        ASSERT_EQ(static_cast<int>(rest.atom_info.size()), 8);

        std::vector<std::pair<gemmi::Atom *, float>> result = coot::dipole::charged_atoms_gemmi(&res, rest);

        EXPECT_EQ(result.size(), 8u);

        EXPECT_EQ(pad_4c(result[0].first->name), std::string(" N  "));
        EXPECT_NEAR(result[0].second, 0.1, 1e-4);

        EXPECT_EQ(pad_4c(result[1].first->name), std::string(" CA "));
        EXPECT_NEAR(result[1].second, 0.2, 1e-4);

        EXPECT_EQ(pad_4c(result[2].first->name), std::string(" C  "));
        EXPECT_NEAR(result[2].second, 0.3, 1e-4);

        EXPECT_EQ(pad_4c(result[3].first->name), std::string(" O  "));
        EXPECT_NEAR(result[3].second, 0.4, 1e-4);

        EXPECT_EQ(pad_4c(result[4].first->name), std::string(" CB "));
        EXPECT_NEAR(result[4].second, 0.5, 1e-4);

        EXPECT_EQ(pad_4c(result[5].first->name), std::string(" CG "));
        EXPECT_NEAR(result[5].second, 0.6, 1e-4);

        EXPECT_EQ(pad_4c(result[6].first->name), std::string(" SD "));
        EXPECT_NEAR(result[6].second, 0.7, 1e-4);

        EXPECT_EQ(pad_4c(result[7].first->name), std::string(" CE "));
        EXPECT_NEAR(result[7].second, 0.8, 1e-4);
    }

    // Case 2: empty restraints — should return empty
    {
        coot::dictionary_residue_restraints_t rest_empty;
        std::vector<std::pair<gemmi::Atom *, float>> result2 = coot::dipole::charged_atoms_gemmi(&res, rest_empty);
        EXPECT_EQ(result2.size(), 0u);
    }

    // Case 3: partial charges disabled (partial_charge.first = false)
    {
        coot::dictionary_residue_restraints_t rest_no_charge;
        for (size_t i = 0; i < atoms.size(); i++) {
            coot::dict_atom da;
            da.atom_id = atoms[i].name;
            da.atom_id_4c = atoms[i].name;
            da.partial_charge = std::make_pair(false, 0.0f);
            rest_no_charge.atom_info.push_back(da);
        }
        std::vector<std::pair<gemmi::Atom *, float>> result3 = coot::dipole::charged_atoms_gemmi(&res, rest_no_charge);
        EXPECT_EQ(result3.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}