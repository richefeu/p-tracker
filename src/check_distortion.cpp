// p-tracker
// Copyright (C) 2023  Vincent Richefeu, Gael Combe
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// With this procedure (for 1G2E), the idea is to follow pairs of points whose distance is meant to remain 
// constant during the loading of the granular sample. To do this, the points are placed on the largest 
// grains so that their tracked-patterns do not overlap. We then follow the variation in average distance 
// <Delta S>, which must be as close as possible to as close as possible to 0, without drifting. 

#include "p-tracker.hpp"

// The procedure to check distortion
void ptracker::check_distortion() {
  std::cout << "*** PROCEDURE CHECK DISTORTION ***" << std::endl;

  char dic_filename[256];
  snprintf(dic_filename, 256, dic_name, iref);
  read_grains(dic_filename);

  std::ofstream file("len_variation.txt");
  double xbeg, ybeg;
  double xend, yend;
  double sx, sy;

  std::vector<double> initial_length;
  std::vector<double> initial_length_corrected;
  for (int igrain = 0; igrain < num_grains; igrain += 2) {
    xbeg = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
    ybeg = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
    xend = (double)(grain[igrain + 1].refcoord_xpix) + grain[igrain + 1].dx;
    yend = (double)(grain[igrain + 1].refcoord_ypix) + grain[igrain + 1].dy;
    sx = xend - xbeg;
    sy = yend - ybeg;
    initial_length.push_back(sqrt(sx * sx + sy * sy));

    undistor(disto_parameters, (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx,
             (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy, xbeg, ybeg);
    undistor(disto_parameters, (double)(grain[igrain + 1].refcoord_xpix) + grain[igrain + 1].dx,
             (double)(grain[igrain + 1].refcoord_ypix) + grain[igrain + 1].dy, xend, yend);
    sx = xend - xbeg;
    sy = yend - ybeg;
    initial_length_corrected.push_back(sqrt(sx * sx + sy * sy));
  }

  for (int idic = ibeg; idic <= iend; idic += iinc) {
    snprintf(dic_filename, 256, dic_name, idic);
    read_grains(dic_filename);

    // current length (raw and undistorted)
    std::vector<double> current_length;
    std::vector<double> current_length_corrected;
    for (int igrain = 0; igrain < num_grains; igrain += 2) {
      xbeg = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
      ybeg = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
      xend = (double)(grain[igrain + 1].refcoord_xpix) + grain[igrain + 1].dx;
      yend = (double)(grain[igrain + 1].refcoord_ypix) + grain[igrain + 1].dy;
      sx = xend - xbeg;
      sy = yend - ybeg;
      current_length.push_back(sqrt(sx * sx + sy * sy));

      undistor(disto_parameters, (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx,
               (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy, xbeg, ybeg);
      undistor(disto_parameters, (double)(grain[igrain + 1].refcoord_xpix) + grain[igrain + 1].dx,
               (double)(grain[igrain + 1].refcoord_ypix) + grain[igrain + 1].dy, xend, yend);
      sx = xend - xbeg;
      sy = yend - ybeg;
      current_length_corrected.push_back(sqrt(sx * sx + sy * sy));
    }

    // length variation (raw and undistorted), and averaging
    std::vector<double> DL_L;
    std::vector<double> DL_L_corrected;
    double mean_DL_L = 0.0;
    double mean_DL_L_corrected = 0.0;
    for (int iseg = 0; iseg < initial_length.size(); iseg++) {
      double variation = (current_length[iseg] - initial_length[iseg]) / initial_length[iseg];
      mean_DL_L += variation;
      DL_L.push_back(variation);
      double variation_corrected =
          (current_length_corrected[iseg] - initial_length_corrected[iseg]) / initial_length_corrected[iseg];
      mean_DL_L_corrected += variation_corrected;
      DL_L_corrected.push_back(variation_corrected);
    }
    mean_DL_L /= (double)initial_length.size();
    mean_DL_L_corrected /= (double)initial_length.size();

    file << idic << ' ' << mean_DL_L << ' ' << mean_DL_L_corrected << std::endl;
  } // end loop over dic files
}