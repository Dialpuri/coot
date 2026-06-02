#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include "function.hh"

// Trim whitespace from atom name (CIF files have unpadded names)
static std::string trim_atom(const std::string &s) {
   auto a = s.find_first_not_of(" \t\r\n");
   if (a == std::string::npos) return "";
   auto b = s.find_last_not_of(" \t\r\n");
   return s.substr(a, b - a + 1);
}

// Helper: find residue by chain name and sequence number
static const gemmi::Residue* find_residue(
    const gemmi::Model &model, const std::string &chain_name, int res_seq) {
   for (const auto &chain : model.chains) {
      if (chain.name == chain_name) {
         for (const auto &res : chain.residues) {
            if (res.seqid.num.value == res_seq) {
               return &res;
            }
         }
      }
   }
   return nullptr;
}

// Helper: find atom by trimmed name
static const gemmi::Atom* find_atom(const gemmi::Residue *res, const std::string &target) {
   for (const auto &at : res->atoms) {
      if (trim_atom(at.name) == target) return &at;
   }
   return nullptr;
}

// Helper: cross product
static clipper::Coord_orth cross(const clipper::Coord_orth &a, const clipper::Coord_orth &b) {
   return clipper::Coord_orth(
       a[1]*b[2] - a[2]*b[1],
       a[2]*b[0] - a[0]*b[2],
       a[0]*b[1] - a[1]*b[0]);
}

// Helper: dot product
static double dot(const clipper::Coord_orth &a, const clipper::Coord_orth &b) {
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Helper: normalize
static clipper::Coord_orth norm_vec(const clipper::Coord_orth &v) {
   double len = std::sqrt(v.lengthsq());
   if (len < 1e-10) return clipper::Coord_orth(0,0,0);
   return clipper::Coord_orth(v[0]/len, v[1]/len, v[2]/len);
}

TEST(OracleTest, side_chain_densities_sample_map) {
    // Load structure using gemmi (CIF file)
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Get a residue with CB — use chain A, residue 10 (standard amino acid)
    const gemmi::Residue *res = find_residue(model, "A", 10);
    ASSERT_NE(res, nullptr);

    // Create side_chain_densities-like parameters
    int n_steps = 1;
    float grid_box_radius = 2.0f;

    // Get CB position and axes
    clipper::Coord_orth cb_pt(0, 0, 0);
    std::vector<clipper::Coord_orth> axes;

    const gemmi::Atom *cb_atom = find_atom(res, "CB");
    const gemmi::Atom *n_atom  = find_atom(res, "N");
    const gemmi::Atom *ca_atom = find_atom(res, "CA");

    if (cb_atom && n_atom && ca_atom) {
        cb_pt    = coot::co_gemmi(cb_atom);
        clipper::Coord_orth n_pt = coot::co_gemmi(n_atom);
        clipper::Coord_orth ca_pt = coot::co_gemmi(ca_atom);

        // Axis 1: CA -> CB
        clipper::Coord_orth axis1 = norm_vec(cb_pt - ca_pt);

        // Axis 2: N -> CA, projected perpendicular to axis1
        clipper::Coord_orth axis2 = ca_pt - n_pt;
        double d = dot(axis2, axis1);
        axis2 = clipper::Coord_orth(
            axis2[0] - axis1[0]*d,
            axis2[1] - axis1[1]*d,
            axis2[2] - axis1[2]*d);
        axis2 = norm_vec(axis2);

        // Axis 3: cross(axis1, axis2)
        clipper::Coord_orth axis3 = norm_vec(cross(axis1, axis2));

        axes.push_back(axis1);
        axes.push_back(axis2);
        axes.push_back(axis3);
    }

    EXPECT_EQ(axes.size(), 3u);

    // Null xmap (no map loaded — density will be 0 everywhere)
    clipper::Xmap<float> xmap; // default-constructed = null

    // ---- Case 1: Valid residue with axes (SAMPLE_FOR_RESIDUE mode) ----
    {
        side_chain_mode_t mode = MODE_SAMPLE_FOR_RESIDUE;
        std::set<int> useable_grid_points;
        for (int i = 0; i < 27; ++i) useable_grid_points.insert(i);

        density_box_t_gemmi db = coot::side_chain_densities::sample_map_gemmi(
            res, nullptr, mode, cb_pt, axes, xmap, "",
            grid_box_radius, n_steps, useable_grid_points);

        EXPECT_FALSE(db.empty());
        EXPECT_NEAR(db.mean, 0, 1e-4);
        EXPECT_NEAR(db.var,  0, 1e-4);
    }

    // ---- Case 2: Null residue — guard clause ----
    {
        side_chain_mode_t mode = MODE_SAMPLE_FOR_RESIDUE;
        std::set<int> useable_grid_points;

        density_box_t_gemmi db = coot::side_chain_densities::sample_map_gemmi(
            nullptr, nullptr, mode, cb_pt, axes, xmap, "",
            grid_box_radius, n_steps, useable_grid_points);

        EXPECT_TRUE(db.empty());
    }

    // ---- Case 3: Empty axes — guard clause ----
    {
        std::vector<clipper::Coord_orth> empty_axes;
        side_chain_mode_t mode = MODE_SAMPLE_FOR_RESIDUE;
        std::set<int> useable_grid_points;

        density_box_t_gemmi db = coot::side_chain_densities::sample_map_gemmi(
            res, nullptr, mode, cb_pt, empty_axes, xmap, "",
            grid_box_radius, n_steps, useable_grid_points);

        EXPECT_TRUE(db.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}