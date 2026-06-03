#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::CRA find_atom_by_cid(gemmi::Model& model, const std::string& cid) {
    // Parse CID like "//A/10/N" → chain="A", seqnum=10, atom_name="N"
    std::string chain_id;
    int seqnum = 0;
    std::string atom_name;

    // Simple parser for "//A/10/N"
    size_t pos = 0;
    // Skip leading slashes
    while (pos < cid.size() && cid[pos] == '/') pos++;

    // Read chain ID
    size_t start = pos;
    while (pos < cid.size() && cid[pos] != '/') pos++;
    chain_id = cid.substr(start, pos - start);
    pos++; // skip '/'

    // Read seqnum
    start = pos;
    while (pos < cid.size() && cid[pos] != '/') pos++;
    seqnum = std::stoi(cid.substr(start, pos - start));
    pos++; // skip '/'

    // Read atom name (rest of string)
    atom_name = cid.substr(pos);

    for (auto& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (auto& res : chain.residues) {
            if (res.seqid.num.value != seqnum) continue;
            for (auto& atom : res.atoms) {
                if (atom.name == atom_name) {
                    return gemmi::CRA{&chain, &res, &atom};
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

TEST(OracleTest, draw_these_atom_contacts) {
    // Load molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Initialize protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    Bond_lines_container container;

    gemmi::Model& model = st.models[0];

    // Case 1: Same residue - atoms from residue A/10 (ILE) — should return false
    {
        gemmi::CRA cra1 = find_atom_by_cid(model, "//A/10/N");
        gemmi::CRA cra2 = find_atom_by_cid(model, "//A/10/CA");
        bool result = false;
        if (cra1.atom && cra2.atom) {
            result = container.draw_these_atom_contacts_gemmi(cra1, cra2, &geom);
        }
        EXPECT_FALSE(result);
    }

    // Case 2: Different chains - atom from A and atom from B — should return true
    {
        gemmi::CRA cra1 = find_atom_by_cid(model, "//A/10/N");
        gemmi::CRA cra2 = find_atom_by_cid(model, "//B/1/O");
        bool result = false;
        if (cra1.atom && cra2.atom) {
            result = container.draw_these_atom_contacts_gemmi(cra1, cra2, &geom);
        }
        EXPECT_TRUE(result);
    }

    // Case 3: Same chain, seq diff > 1 - A/10 and A/12 — should return true
    {
        gemmi::CRA cra1 = find_atom_by_cid(model, "//A/10/N");
        gemmi::CRA cra2 = find_atom_by_cid(model, "//A/12/N");
        bool result = false;
        if (cra1.atom && cra2.atom) {
            result = container.draw_these_atom_contacts_gemmi(cra1, cra2, &geom);
        }
        EXPECT_TRUE(result);
    }

    // Case 4: Same chain, neighboring, main chain atoms (N, N) — should return false
    {
        gemmi::CRA cra1 = find_atom_by_cid(model, "//A/10/N");
        gemmi::CRA cra2 = find_atom_by_cid(model, "//A/11/N");
        bool result = false;
        if (cra1.atom && cra2.atom) {
            result = container.draw_these_atom_contacts_gemmi(cra1, cra2, &geom);
        }
        EXPECT_FALSE(result);
    }

    // Case 5: Same chain, neighboring, side chain atoms (CB, CA) — should return true
    {
        gemmi::CRA cra1 = find_atom_by_cid(model, "//A/10/CB");
        gemmi::CRA cra2 = find_atom_by_cid(model, "//A/11/CA");
        bool result = false;
        if (cra1.atom && cra2.atom) {
            result = container.draw_these_atom_contacts_gemmi(cra1, cra2, &geom);
        }
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}