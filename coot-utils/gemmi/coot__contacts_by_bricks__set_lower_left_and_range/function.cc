#include "function.hh"

void coot::contacts_by_bricks::set_lower_left_and_range_gemmi(
    gemmi::Atom const* atoms[],
    int n_atoms,
    float brick_size,
    float lower_left[3],
    int range[3]) {

   // Initialize lower_left to a large value
   for (int i = 0; i < 3; i++)
      lower_left[i] = 9999.9f;

   // Initialize range to 0
   for (int i = 0; i < 3; i++)
      range[i] = 0;

   unsigned int n = static_cast<unsigned int>(n_atoms);

   // Find minimum coordinates
   for (unsigned int i = 0; i < n; i++) {
      gemmi::Atom const* atom = atoms[i];
      float pos[3];
      pos[0] = atom->pos.x;
      pos[1] = atom->pos.y;
      pos[2] = atom->pos.z;
      for (int j = 0; j < 3; j++)
         if (pos[j] < lower_left[j])
            lower_left[j] = pos[j];
   }

   // Subtract 6.0 from each component
   for (int i = 0; i < 3; i++)
      lower_left[i] -= 6.0f;

   // Compute brick ranges
   float inv_brick_size = 1.0f / brick_size;
   for (unsigned int i = 0; i < n; i++) {
      gemmi::Atom const* atom = atoms[i];
      float pos[3];
      pos[0] = atom->pos.x;
      pos[1] = atom->pos.y;
      pos[2] = atom->pos.z;
      for (int j = 0; j < 3; j++) {
         float f = (pos[j] - lower_left[j]) * inv_brick_size;
         int brick_idx = static_cast<int>(f);
         if (brick_idx > range[j]) {
            range[j] = brick_idx;
         }
      }
   }

   // Add 1 to each range component
   for (int i = 0; i < 3; i++)
      range[i] += 1;
}