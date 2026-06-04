#pragma once

#include <gemmi/math.hpp>
#include <gemmi/model.hpp>
#include <ostream>
#include <string>

namespace coot {

struct molecule_extents_t_gemmi {
  gemmi::Vec3 front;
  gemmi::Vec3 back;
  gemmi::Vec3 left;
  gemmi::Vec3 right;
  gemmi::Vec3 top;
  gemmi::Vec3 bottom;
};

inline molecule_extents_t_gemmi compute_extents_gemmi(const gemmi::Model& model, double padding = 0.0) {
  molecule_extents_t_gemmi e;
  double min_x = 1e9, max_x = -1e9;
  double min_y = 1e9, max_y = -1e9;
  double min_z = 1e9, max_z = -1e9;

  for (const auto& chain : model.chains) {
    for (const auto& res : chain.residues) {
      for (const auto& atom : res.atoms) {
        double x = atom.pos.x;
        double y = atom.pos.y;
        double z = atom.pos.z;
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
        if (z < min_z) min_z = z;
        if (z > max_z) max_z = z;
      }
    }
  }

  min_x -= padding; max_x += padding;
  min_y -= padding; max_y += padding;
  min_z -= padding; max_z += padding;

  double center_x = (min_x + max_x) / 2.0;
  double center_y = (min_y + max_y) / 2.0;
  double center_z = (min_z + max_z) / 2.0;

  e.front  = gemmi::Vec3(center_x, center_y, min_z);
  e.back   = gemmi::Vec3(center_x, center_y, max_z);
  e.left   = gemmi::Vec3(min_x,  center_y, center_z);
  e.right  = gemmi::Vec3(max_x,  center_y, center_z);
  e.top    = gemmi::Vec3(center_x, max_y,  center_z);
  e.bottom = gemmi::Vec3(center_x, min_y,  center_z);

  return e;
}

inline std::ostream& operator_stream_insert_gemmi(std::ostream &s, const molecule_extents_t_gemmi &e) {
  s << "front:  (" << e.front.x << "," << e.front.y << "," << e.front.z  << ")" << std::endl;
  s << "back :  (" << e.back.x  << "," << e.back.y  << "," << e.back.z   << ")" << std::endl;
  s << "left :  (" << e.left.x  << "," << e.left.y  << "," << e.left.z   << ")" << std::endl;
  s << "right:  (" << e.right.x << "," << e.right.y << "," << e.right.z  << ")" << std::endl;
  s << "top  :  (" << e.top.x   << "," << e.top.y   << "," << e.top.z    << ")" << std::endl;
  s << "bottom: (" << e.bottom.x<< "," << e.bottom.y<< "," << e.bottom.z << ")" << std::endl;
  return s;
}

} // namespace coot