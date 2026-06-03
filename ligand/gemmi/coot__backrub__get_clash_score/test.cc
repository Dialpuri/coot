#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <cmath>
#include <vector>
#include <string>

TEST(OracleTest, get_clash_score) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];

    // Find 3 consecutive residues in chain A
    gemmi::Residue *prev_res = nullptr, *this_res = nullptr, *next_res = nullptr;
    gemmi::Chain *chain_ptr = nullptr;
    for (auto &chain : model.chains) {
        if (chain.name != "A") continue;
        if (chain.residues.size() >= 3) {
            prev_res = &chain.residues[0]; // resno 1
            this_res = &chain.residues[1]; // resno 2
            next_res = &chain.residues[2]; // resno 3
            chain_ptr = &chain;
            break;
        }
    }
    ASSERT_NE(this_res, nullptr);
    ASSERT_NE(chain_ptr, nullptr);

    std::string chain_id_str = chain_ptr->name;
    int resno_prev = prev_res->seqid.num.value;
    int resno_this = this_res->seqid.num.value;
    int resno_next = next_res->seqid.num.value;

    EXPECT_EQ(chain_id_str, "A");
    EXPECT_EQ(resno_prev, 1);
    EXPECT_EQ(resno_this, 2);
    EXPECT_EQ(resno_next, 3);

    // Center of this residue
    int natoms = static_cast<int>(this_res->atoms.size());
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < natoms; i++) {
        cx += this_res->atoms[i].pos.x;
        cy += this_res->atoms[i].pos.y;
        cz += this_res->atoms[i].pos.z;
    }
    cx /= natoms;
    cy /= natoms;
    cz /= natoms;

    EXPECT_NEAR(cx, 13.323, 1e-2);
    EXPECT_NEAR(cy, 3.69133, 1e-2);
    EXPECT_NEAR(cz, 79.8516, 1e-2);

    // Sphere selection: collect all atoms within 20.0 Å of center
    float rrr = 20.0f;
    double rrr_sq = rrr * rrr;
    std::vector<gemmi::CRA> sphere_atoms;
    for (auto &ch : model.chains) {
        for (auto &res : ch.residues) {
            for (auto &atom : res.atoms) {
                double dx = atom.pos.x - cx;
                double dy = atom.pos.y - cy;
                double dz = atom.pos.z - cz;
                double dist_sq = dx * dx + dy * dy + dz * dz;
                if (dist_sq <= rrr_sq) {
                    sphere_atoms.push_back({&ch, &res, &atom});
                }
            }
        }
    }
    int n_sphere_atoms = static_cast<int>(sphere_atoms.size());

    EXPECT_FLOAT_EQ(rrr, 20.0f);
    EXPECT_EQ(n_sphere_atoms, 347);

    // Count how many sphere atoms would be skipped vs counted
    int skipped = 0, counted = 0;
    for (int i = 0; i < n_sphere_atoms; i++) {
        gemmi::CRA cra = sphere_atoms[i];
        std::string ach = cra.chain ? cra.chain->name : "";
        int are = cra.residue ? cra.residue->seqid.num.value : 0;
        bool skip = (ach == chain_id_str && (are == resno_prev || are == resno_this || are == resno_next));
        if (skip) skipped++;
        else counted++;
    }
    EXPECT_EQ(skipped, 25);
    EXPECT_EQ(counted, 322);

    // Find first non-skipped sphere atom
    gemmi::CRA target_atom{nullptr, nullptr, nullptr};
    for (int i = 0; i < n_sphere_atoms; i++) {
        gemmi::CRA cra = sphere_atoms[i];
        std::string ach = cra.chain ? cra.chain->name : "";
        int are = cra.residue ? cra.residue->seqid.num.value : 0;
        bool skip = (ach == chain_id_str && (are == resno_prev || are == resno_this || are == resno_next));
        if (!skip) { target_atom = cra; break; }
    }
    ASSERT_NE(target_atom.atom, nullptr);
    EXPECT_EQ(target_atom.residue->seqid.num.value, 4);

    // === CASE 1: Place test atom RIGHT ON TOP of a non-skipped sphere atom ===
    {
        coot::minimol::molecule minimol_close;
        coot::minimol::fragment frag_close;
        frag_close.fragment_id = chain_id_str;
        coot::minimol::residue min_res_close(resno_this, this_res->name);

        coot::minimol::atom close_atom1;
        close_atom1.name = "CB";
        close_atom1.element = " C";
        close_atom1.pos = clipper::Coord_orth(target_atom.atom->pos.x + 0.1, target_atom.atom->pos.y, target_atom.atom->pos.z);
        min_res_close.atoms.push_back(close_atom1);

        EXPECT_NEAR(close_atom1.pos.x(), 17.525, 1e-2);
        EXPECT_NEAR(close_atom1.pos.y(), 3.398, 1e-2);
        EXPECT_NEAR(close_atom1.pos.z(), 77.6, 1e-2);

        frag_close.residues.push_back(min_res_close);
        minimol_close.fragments.push_back(frag_close);

        int water_interaction_mode = 0;
        std::pair<float, std::vector<gemmi::CRA>> result_close =
            coot::backrub::get_clash_score_gemmi(
                minimol_close, sphere_atoms, chain_id_str, resno_prev, resno_this, resno_next, water_interaction_mode);

        EXPECT_NEAR(result_close.first, 10051.4, 0.1);
        EXPECT_EQ(result_close.second.size(), 0u);
    }

    // === CASE 2: Far away (no clashes) ===
    {
        coot::minimol::molecule minimol_far;
        coot::minimol::fragment frag_far;
        frag_far.fragment_id = chain_id_str;
        coot::minimol::residue min_res_far(resno_this, this_res->name);

        coot::minimol::atom far_atom1;
        far_atom1.name = "CB";
        far_atom1.element = " C";
        far_atom1.pos = clipper::Coord_orth(cx + 200.0, cy + 200.0, cz + 200.0);
        min_res_far.atoms.push_back(far_atom1);

        EXPECT_NEAR(far_atom1.pos.x(), 213.323, 1e-2);
        EXPECT_NEAR(far_atom1.pos.y(), 203.691, 1e-2);
        EXPECT_NEAR(far_atom1.pos.z(), 279.852, 1e-2);

        frag_far.residues.push_back(min_res_far);
        minimol_far.fragments.push_back(frag_far);

        int water_interaction_mode = 0;
        std::pair<float, std::vector<gemmi::CRA>> result_far =
            coot::backrub::get_clash_score_gemmi(
                minimol_far, sphere_atoms, chain_id_str, resno_prev, resno_this, resno_next, water_interaction_mode);

        EXPECT_NEAR(result_far.first, 0.0, 1e-4);
        EXPECT_EQ(result_far.second.size(), 0u);
    }

    // === CASE 3: Empty molecule ===
    {
        coot::minimol::molecule minimol_empty;
        coot::minimol::fragment frag_empty;
        frag_empty.fragment_id = chain_id_str;
        coot::minimol::residue min_res_empty(resno_this, this_res->name);
        frag_empty.residues.push_back(min_res_empty);
        minimol_empty.fragments.push_back(frag_empty);

        int water_interaction_mode = 0;
        std::pair<float, std::vector<gemmi::CRA>> result_empty =
            coot::backrub::get_clash_score_gemmi(
                minimol_empty, sphere_atoms, chain_id_str, resno_prev, resno_this, resno_next, water_interaction_mode);

        EXPECT_NEAR(result_empty.first, 0.0, 1e-4);
        EXPECT_EQ(result_empty.second.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}