#include "function.hh"
#include <fstream>
#include <iostream>
#include <cstdlib>

namespace coot::side_chain_densities {

void write_density_box_gemmi(
    float* density_box,
    int    n_steps,
    gemmi::CRA cra,
    const std::string& id,
    const std::string& rotamer_name) 
{
    if (!cra.residue) return;

    std::string res_name = cra.residue->name;
    std::string rot_name = rotamer_name;

    if (rot_name.empty()) return;

    std::string dir = "side-chain-data";
    std::string rot_dir = dir + "/" + res_name + "/" + rot_name;

    // Chain ID from cra.chain->name
    std::string chain_id = "A";
    if (cra.chain) {
        chain_id = cra.chain->name;
    }
    int seq_num = cra.residue->seqid.num.value;

    std::string file_name = rot_dir + "/" + id + "-" + chain_id + "-" + std::to_string(seq_num) + ".tab";

    // Do not create directory — matching original util::create_directory behavior
    // which does not create nested parent directories, so file open fails silently.

    std::ofstream f(file_name.c_str());
    if (f) {
        int n_per_side = 2 * n_steps + 1;
        int n_box_vol = n_per_side * n_per_side * n_per_side;
        for (int i = 0; i < n_box_vol; i++) {
            float v = density_box[i];
            if (std::isnan(v)) {
                std::cout << "ERROR:: " << file_name << " found a nan " << i << std::endl;
                exit(1);
            } else {
                f << v << " ";
                if (i % n_per_side == 0)
                    f << "\n";
            }
        }
        f << "\n";
    } else {
        std::cout << "WARNING:: cannot open file " << file_name << std::endl;
    }
}

} // namespace coot::side_chain_densities