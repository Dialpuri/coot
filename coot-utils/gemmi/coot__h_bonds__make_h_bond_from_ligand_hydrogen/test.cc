#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find atom by chain name, residue seqnum, and atom name (trimmed)
static const gemmi::Atom*
find_atom(const gemmi::Model &model, char chain_name, int seqnum, const std::string &atom_name) {
    std::string trimmed = atom_name;
    // trim whitespace
    size_t start = trimmed.find_first_not_of(" \t");
    size_t end   = trimmed.find_last_not_of(" \t");
    if (start == std::string::npos) return nullptr;
    trimmed = trimmed.substr(start, end - start + 1);

    for (const auto &chain : model.chains) {
        if (chain.name != std::string(1, chain_name)) continue;
        for (const auto &res : chain.residues) {
            if (res.seqid.num.value != seqnum) continue;
            for (const auto &atom : res.atoms) {
                std::string aname = atom.name;
                size_t astart = aname.find_first_not_of(" \t");
                size_t aend   = aname.find_last_not_of(" \t");
                if (astart == std::string::npos) continue;
                aname = aname.substr(astart, aend - astart + 1);
                if (aname == trimmed)
                    return &atom;
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, make_h_bond_from_ligand_hydrogen) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model &model = st.models[0];

    // Get atoms we'll reuse
    // //A/2/H   -> chain A, seqnum 2, atom "H"
    const gemmi::Atom *at_1 = find_atom(model, 'A', 2, "H");
    // //A/3/OD1 -> chain A, seqnum 3, atom "OD1"
    const gemmi::Atom *at_2 = find_atom(model, 'A', 3, "OD1");
    // //A/2/N   -> chain A, seqnum 2, atom "N"
    const gemmi::Atom *donor = find_atom(model, 'A', 2, "N");
    // //A/3/CG  -> chain A, seqnum 3, atom "CG"
    const gemmi::Atom *aa = find_atom(model, 'A', 3, "CG");

    ASSERT_NE(at_1, nullptr);
    ASSERT_NE(at_2, nullptr);

    // Case 1: with neighbours — attempt a valid H-bond (angles >= 90, dist < 3.9)
    {
        std::vector<std::pair<const gemmi::Atom *, float>> nb_1;
        std::vector<std::pair<const gemmi::Atom *, float>> nb_2;
        if (donor) nb_1.push_back({donor, 1.0f});
        if (aa)    nb_2.push_back({aa, 1.0f});

        auto result = coot::h_bonds::make_h_bond_from_ligand_hydrogen_gemmi(at_1, at_2, nb_1, nb_2);

        // LHS accessor changed: MMDB GetAtomName() -> gemmi atom.name (trimmed)
        EXPECT_STREQ(at_1->name.c_str(), "H");
        EXPECT_STREQ(at_2->name.c_str(), "OD1");
        EXPECT_EQ(nb_1.size(), 1u);
        EXPECT_EQ(nb_2.size(), 1u);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second.dist, 4.59178, 1e-4);
        EXPECT_NEAR(result.second.angle_1, 88.664, 1e-3);
        EXPECT_NEAR(result.second.angle_2, -1.0, 1e-4);
        EXPECT_NEAR(result.second.angle_3, 99.0782, 1e-4);
        if (result.second.donor)
            EXPECT_STREQ(result.second.donor->name.c_str(), "N");
        if (result.second.acceptor_neigh)
            EXPECT_STREQ(result.second.acceptor_neigh->name.c_str(), "CG");
    }

    // Case 2: empty neighbour vectors → no donor-acceptor distance check → false
    {
        std::vector<std::pair<const gemmi::Atom *, float>> nb_1;
        std::vector<std::pair<const gemmi::Atom *, float>> nb_2;

        auto result = coot::h_bonds::make_h_bond_from_ligand_hydrogen_gemmi(at_1, at_2, nb_1, nb_2);

        // LHS accessor changed: MMDB GetAtomName() -> gemmi atom.name (trimmed)
        EXPECT_STREQ(at_1->name.c_str(), "H");
        EXPECT_STREQ(at_2->name.c_str(), "OD1");
        EXPECT_EQ(nb_1.size(), 0u);
        EXPECT_EQ(nb_2.size(), 0u);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second.dist, 4.59178, 1e-4);
    }

    // Case 3: atoms too far apart (no good donor-acceptor distance)
    {
        const gemmi::Atom *h_far = find_atom(model, 'A', 2, "HZ1");
        const gemmi::Atom *o_far = find_atom(model, 'A', 10, "O");

        if (h_far && o_far) {
            std::vector<std::pair<const gemmi::Atom *, float>> nb_1;
            std::vector<std::pair<const gemmi::Atom *, float>> nb_2;
            // NZ is the donor for HZ1
            const gemmi::Atom *nz = find_atom(model, 'A', 2, "NZ");
            if (nz) nb_1.push_back({nz, 1.0f});
            // CB of residue 10 as AA neighbour
            const gemmi::Atom *cb10 = find_atom(model, 'A', 10, "CB");
            if (cb10) nb_2.push_back({cb10, 1.0f});

            auto result = coot::h_bonds::make_h_bond_from_ligand_hydrogen_gemmi(h_far, o_far, nb_1, nb_2);

            // LHS accessor changed: MMDB GetAtomName() -> gemmi atom.name (trimmed)
            EXPECT_STREQ(h_far->name.c_str(), "HZ1");
            EXPECT_STREQ(o_far->name.c_str(), "O");
            EXPECT_FALSE(result.first);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}