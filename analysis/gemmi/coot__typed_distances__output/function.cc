#include "function.hh"

namespace coot {
namespace typed_distances {

void output_gemmi(
    const std::map<const gemmi::Residue*, std::map<int, std::vector<float>>>& residue_distances_map,
    const std::map<const gemmi::Residue*, std::vector<const gemmi::Residue*>>& residues_within_window
) {
    bool do_windowed_residue_table = true;

    if (do_windowed_residue_table) {
        std::cout << "--------------- windowed residue table -----------" << std::endl;

        std::cout << "residue_distances_map size " << residue_distances_map.size() << std::endl;

        for (auto it = residue_distances_map.begin(); it != residue_distances_map.end(); ++it) {
            const gemmi::Residue& res = *it->first;
            if (res.name != "HOH") {
                auto wit = residues_within_window.find(it->first);
                if (wit == residues_within_window.end()) {
                    continue;
                }
                const std::vector<const gemmi::Residue*>& riw = wit->second;
                std::vector<float> sums(6, 0.0);
                std::vector<unsigned int> n_distances(6, 0);
                
                for (int j = 0; j < 6; j++) {
                    for (std::size_t ir = 0; ir < riw.size(); ir++) {
                        const gemmi::Residue* res_inner = riw[ir];
                        auto it_r = residue_distances_map.find(res_inner);
                        if (it_r == residue_distances_map.end()) {
                            std::cout << "Oopps - bad residue for map " << std::endl;
                        } else {
                            // happy path
                            const std::map<int, std::vector<float>>& m = it_r->second;
                            auto it_inner = m.find(j);
                            if (it_inner == m.end()) {
                                if (false)
                                    std::cout << "Oopps failed to find bin-index " << j << " "
                                              << residue_to_string(*res_inner) << std::endl;
                            } else {
                                const std::vector<float>& v = it_inner->second;
                                n_distances[j] += v.size();
                                for (std::size_t i = 0; i < v.size(); i++) {
                                    sums[j] += v[i];
                                }
                            }
                        }
                    }
                }

                int total_n_distances = 0;
                std::vector<float> fractions(6, 0.0);
                std::vector<float> means(6, 0.0);
                for (int j = 0; j < 6; j++)
                    total_n_distances += n_distances[j];

                for (int j = 0; j < 6; j++) {
                    fractions[j] = static_cast<float>(n_distances[j]) / static_cast<float>(total_n_distances);
                }

                std::cout << "Sums: " << residue_to_string(res);
                for (int j = 0; j < 6; j++)
                    std::cout << " type " << j << ":  " << n_distances[j] << " " << fractions[j] << " ";
                std::cout << "\n";
            }
        }
    }
}

} // namespace typed_distances
} // namespace coot