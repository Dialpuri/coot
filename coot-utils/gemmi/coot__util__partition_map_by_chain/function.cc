#include "function.hh"
#include <gemmi/resinfo.hpp>

namespace coot {
namespace util {

std::vector<std::pair<std::string, clipper::Xmap<float>>>
partition_map_by_chain_gemmi(const clipper::Xmap<float>& xmap,
                             const gemmi::Model* model_p,
                             std::string* state_string_p) {

   std::vector<std::pair<std::string, clipper::Xmap<float>>> v;

   if (!model_p) return v;

   // Lambda: make_side_chain_centre
   auto make_side_chain_centre = [] (const gemmi::Residue& residue) {
      bool status = false;
      clipper::Coord_orth centre;
      int count = 0;
      for (const auto& atom : residue.atoms) {
         // In gemmi, TER is handled by chain boundaries - no TER atoms to skip
         std::string atom_name(atom.name);
         if (atom_name == " N  " || atom_name == " H  " || atom_name == " C  " ||
             atom_name == " CA " || atom_name == " O  " || atom_name == " HA ")
            continue;
         centre += coot::co_gemmi(&atom);
         count += 1;
      }
      if (count > 0) {
         double sf = 1.0 / static_cast<double>(count);
         centre = clipper::Coord_orth(centre.x() * sf, centre.y() * sf, centre.z() * sf);
         status = true;
      }
      return std::make_pair(status, centre);
   };

   // Lambda: make_reference_points_for_chains
   auto make_reference_points_for_chains = [&make_side_chain_centre] (
       const gemmi::Model& model,
       const clipper::Cell& cell,
       const clipper::Grid_sampling& gs) {
      std::map<std::string, std::vector<clipper::Coord_grid>> coordinates_grid_points_map;

      for (const auto& chain : model.chains) {
         std::string chain_id(chain.name);
         for (const auto& residue : chain.residues) {
            for (const auto& atom : residue.atoms) {
               std::string name(atom.name);
               if (name == " CA " || name == " P  " || name == " C4'" ||
                   name == " N9 " || name == " N1 " || name == " C4 ") {
                  clipper::Coord_orth co = coot::co_gemmi(&atom);
                  clipper::Coord_frac cf = co.coord_frac(cell);
                  clipper::Coord_grid cg = cf.coord_grid(gs);
                  coordinates_grid_points_map[chain_id].push_back(cg);
               }
            }

            // Check if amino acid using chemcomp database
            auto comp = gemmi::find_tabulated_residue(residue.name);
            if (comp.is_amino_acid()) {
               auto side_chain_centre = make_side_chain_centre(residue);
               if (side_chain_centre.first) {
                  clipper::Coord_frac cf = side_chain_centre.second.coord_frac(cell);
                  clipper::Coord_grid cg = cf.coord_grid(gs);
                  coordinates_grid_points_map[chain_id].push_back(cg);
               }
            }
         }
      }
      return coordinates_grid_points_map;
   };

   // Lambda: split_biggest_chain
   auto split_biggest_chain = [] (std::map<std::string, std::vector<clipper::Coord_grid>>* rp_p) {
      std::string biggest_chain;
      unsigned int n_biggest = 0;
      for (const auto& item : *rp_p) {
         unsigned int n = item.second.size();
         if (n > n_biggest) {
            n_biggest = n;
            biggest_chain = item.first;
         }
      }
      if (n_biggest > 0) {
         std::vector<clipper::Coord_grid>& v_vec = (*rp_p)[biggest_chain];
         size_t half_size = n_biggest / 2;
         if (n_biggest % 2 != 0) half_size += 1;
         std::vector<clipper::Coord_grid> new_v;
         for (unsigned int i = half_size; i < v_vec.size(); i++)
            new_v.push_back(v_vec[i]);
         v_vec.resize(half_size);
         std::string new_chain_id = biggest_chain + "+";
         (*rp_p)[new_chain_id] = new_v;
      }
   };

   std::cout << "Making reference points for chains" << std::endl;
   if (state_string_p) *state_string_p = "Making reference points for chains";

   clipper::Cell cell = xmap.cell();
   clipper::Spacegroup sg = xmap.spacegroup();
   clipper::Grid_sampling gs = xmap.grid_sampling();

   std::map<std::string, std::vector<clipper::Coord_grid>> rp =
       make_reference_points_for_chains(*model_p, cell, gs);

   split_biggest_chain(&rp);

   std::vector<std::string> chain_ids;
   for (const auto& item : rp)
      chain_ids.push_back(item.first);

   clipper::Xmap<std::map<std::string, int>> distance_map;
   distance_map.init(sg, cell, gs);

   if (state_string_p) *state_string_p = "Filling distance map with initial values";
   std::map<std::string, int> starting_distance_map;
   for (const auto& item : chain_ids)
      starting_distance_map[item] = 999999;

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = distance_map.first(); !ix.last(); ix.next()) {
      clipper::Coord_grid cg = ix.coord();
      float f = xmap.get_data(cg);
      if (f != 0.0)
         distance_map[ix] = starting_distance_map;
   }

   auto manhattan_check = [] (const std::string& chain_id,
                              const std::vector<clipper::Coord_grid>& reference_points,
                              clipper::Xmap<std::map<std::string, int>>* distance_map_p,
                              std::string* info_string_p) {

      if (info_string_p) *info_string_p = "Distance check for " + chain_id;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = distance_map_p->first(); !ix.last(); ix.next()) {
         if ((*distance_map_p)[ix].empty()) continue;
         auto& current_dist = (*distance_map_p)[ix][chain_id];
         clipper::Coord_grid cg = ix.coord();
         for (unsigned int i = 0; i < reference_points.size(); i++) {
            const auto& ref_grid_pt = reference_points[i];
            int d_u = ref_grid_pt.u() - cg.u();
            int d_v = ref_grid_pt.v() - cg.v();
            int d_w = ref_grid_pt.w() - cg.w();
            int manhattan_dist = abs(d_u) + abs(d_v) + abs(d_w);
            if (manhattan_dist < current_dist) {
               current_dist = manhattan_dist;
            }
         }
      }
   };

   std::vector<std::thread> threads;
   for (const auto& chain_id : chain_ids) {
      const std::vector<clipper::Coord_grid>& reference_points = rp[chain_id];
      threads.push_back(std::thread(manhattan_check, chain_id, reference_points,
                                    &distance_map, state_string_p));
   }

   if (state_string_p) *state_string_p = "Joining threads...";
   for (auto& thread : threads) thread.join();

   // now extract each of the maps for each chain
   clipper::Xmap<std::string> chain_map;
   chain_map.init(sg, cell, gs);
   for (ix = distance_map.first(); !ix.last(); ix.next()) {
      const auto& dist_map = distance_map[ix];
      std::string best_chain_id;
      int dist_best = 999999;
      for (const auto& item : dist_map) {
         if (item.second < dist_best) {
            dist_best = item.second;
            best_chain_id = item.first;
         }
      }
      // chain-ids ending in "+" really are part of another chain
      if (best_chain_id.size() == 2 && best_chain_id[1] == '+')
         best_chain_id = std::string(1, best_chain_id[0]);
      chain_map[ix] = best_chain_id;
   }

   for (const auto& chain_id : chain_ids) {
      // chain-ids ending in "+" really are part of another chain
      if (chain_id.size() == 2 && chain_id[1] == '+')
         continue;

      if (state_string_p) *state_string_p = "Constructing map for chain " + chain_id;

      clipper::Xmap<float> map_for_chain;
      map_for_chain.init(sg, cell, gs);
      for (ix = chain_map.first(); !ix.last(); ix.next()) {
         const auto& chain_for_this_grid_point = chain_map[ix];
         clipper::Coord_grid cg = ix.coord();
         if (chain_for_this_grid_point == chain_id) {
            float f = xmap.get_data(cg);
            map_for_chain.set_data(cg, f);
         } else {
            map_for_chain.set_data(cg, 0.0f);
         }
      }
      v.push_back(std::make_pair(chain_id, map_for_chain));
   }

   return v;
}

} // namespace util
} // namespace coot