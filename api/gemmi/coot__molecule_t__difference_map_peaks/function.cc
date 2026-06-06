#include "function.hh"

namespace coot {

namespace {

double distance_sq(const clipper::Coord_orth& co, const gemmi::Atom& atom) {
    double dx = atom.pos.x - co.x();
    double dy = atom.pos.y - co.y();
    double dz = atom.pos.z - co.z();
    return dx * dx + dy * dy + dz * dz;
}

coot::residue_spec_t get_closest_residue_to_point(const gemmi::Structure& st,
                                                    const clipper::Coord_orth& co) {
    coot::residue_spec_t spec;

    double d_best = 999999999999.0;

    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                for (const gemmi::Atom& atom : residue.atoms) {
                    bool is_ter = (atom.name == "TER");
                    if (!is_ter) {
                        double dd = distance_sq(co, atom);
                        if (dd < d_best) {
                            d_best = dd;
                            spec.chain_id = chain.name;
                            spec.res_no = residue.seqid.num.value;
                            spec.ins_code = std::string(1, residue.seqid.icode);
                        }
                    }
                }
            }
        }
    }

    return spec;
}

} // anonymous namespace

std::vector<coot::molecule_t::interesting_place_t>
difference_map_peaks_gemmi(coot::molecule_t& mol, const gemmi::Structure& st, float n_rmsd) {

   auto make_button_label = [] (unsigned int i, const std::pair<clipper::Coord_orth, float> &peak) {
      std::string s = "Index ";
      s += coot::util::int_to_string(i);
      s += " Position: (";
      s += coot::util::float_to_string_using_dec_pl(peak.first.x(), 2);
      s += std::string(", ");
      s += coot::util::float_to_string_using_dec_pl(peak.first.y(), 2);
      s += std::string(", ");
      s += coot::util::float_to_string_using_dec_pl(peak.first.z(), 2);
      s += std::string(") Height ");
      s += coot::util::float_to_string_using_dec_pl(peak.second, 2);
      return s;
   };

   unsigned int n_limit = 100;

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (!st.models.empty()) {
      coot::peak_search ps(mol.xmap);
      float rmsd = mol.get_map_rmsd_approx();
      bool skip_symmetry_check = false;
      std::vector<std::pair<clipper::Coord_orth, float> > peaks = ps.get_peaks(mol.xmap, n_rmsd, true, true);
      for (unsigned int i=0; i<peaks.size(); i++) {
         const auto &peak = peaks[i];
         float badness = 7.0f * std::abs(peak.second/rmsd);
         std::string button_label = make_button_label(i, peak);
         coot::molecule_t::interesting_place_t dmp("difference-map-peak", peak.first, button_label);
         dmp.residue_spec = get_closest_residue_to_point(st, peak.first);
         dmp.set_feature_value(peak.second);
         dmp.set_badness_value(badness);
         v.push_back(dmp);
      }

      if (v.size() <= n_limit) {
         auto sorter_1 = [] (const coot::molecule_t::interesting_place_t &i1, const coot::molecule_t::interesting_place_t &i2) {
            return i1.feature_value < i2.feature_value;
         };
         std::sort(v.begin(), v.end(), sorter_1);
      } else {
         auto sorter_1 = [] (const coot::molecule_t::interesting_place_t &i1, const coot::molecule_t::interesting_place_t &i2) {
            return i1.feature_value < i2.feature_value;
         };
         auto sorter_2 = [] (const coot::molecule_t::interesting_place_t &i1, const coot::molecule_t::interesting_place_t &i2) {
            return fabsf(i2.feature_value) < fabsf(i1.feature_value);
         };
         std::sort(v.begin(), v.end(), sorter_2);
         v.resize(n_limit);
         std::sort(v.begin(), v.end(), sorter_1);
      }

   } else {
      std::cout << "ERROR:: difference_map_peaks_gemmi() null mol" << std::endl;
   }
   return v;
}

} // namespace coot