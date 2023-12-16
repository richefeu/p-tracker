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

#include "p-tracker.hpp"

void ptracker::compute_fluctuations() {
  std::cout << "*** COMPUTE FLUCTUATIONS ***" << std::endl;
  
  int dic_num_window = 1;           // TODO ce sera un parametre entré par l'utilisateur
  bool disto_para_available = true; // TODO ce sera un parametre mis à true si un param de disto est lu

  for (int idic = ibeg; idic <= iend - dic_num_window; idic += iinc) {

    // get the origin position (and undistort it)
    char dic_filename[256];
    snprintf(dic_filename, 256, dic_name, idic);
    read_grains(dic_filename);
    std::vector<double> xfrom(num_grains), yfrom(num_grains);
    for (int igrain = 0; igrain < num_grains; igrain++) {
      xfrom[igrain] = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
      yfrom[igrain] = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
      if (disto_para_available) {
        double xu, yu;
        undistor(disto_parameters, xfrom[igrain], yfrom[igrain], xu, yu);
        xfrom[igrain] = xu;
        yfrom[igrain] = yu;
      }
    }
    
    // get the arrival position (and undistort it)
    snprintf(dic_filename, 256, dic_name, idic+dic_num_window);
    read_grains(dic_filename);
    std::vector<double> xto(num_grains), yto(num_grains);
    for (int igrain = 0; igrain < num_grains; igrain++) {
      xto[igrain] = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
      yto[igrain] = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
      if (disto_para_available) {
        double xu, yu;
        undistor(disto_parameters, xto[igrain], yto[igrain], xu, yu);
        xto[igrain] = xu;
        yto[igrain] = yu;
      }
    }

    // ...
    
    for (int igrain = 0; igrain < num_grains; igrain++) {
      // get the reduced origin position of the grain
      /*
      frame.dxBL = ...;
      frame.dyBL = ...;
      frame.dxBR = ...;
      frame.dyBR = ...;
      frame.dxTL = ...;
      frame.dyTL = ...;
      frame.dxTR = ...;
      frame.dyTR = ...;
      */
      
      double s, t, s2, t2;
      int nok = frame.inverseBilinear(xfrom[igrain], yfrom[igrain], s, t, s2, t2);
      if (nok == 1) {
        double interp_dx, interp_dy;
        frame.bilinear(s, t, interp_dx, interp_dx);
      }
    }


    /// <<<<<
    /*
    for (int i = 4; i < data.size(); i++) {
      nok = inverseBilinear(data[0].x, data[0].y, data[3].x, data[3].y, data[1].x, data[1].y, data[2].x, data[2].y,
                            data[i].x, data[i].y, &s, &t, &s2, &t2);
      bilinear(data[0].dx, data[0].dy, data[3].dx, data[3].dy, data[1].dx, data[1].dy, data[2].dx, data[2].dy, s, t,
               &(data[i].dxmean), &(data[i].dymean));
      out << data[i].x << "\t" << data[i].y << "\t" << data[i].dx << "\t" << data[i].dy << "\t" << data[i].dxmean
          << "\t" << data[i].dymean << endl;
    }
    */
    /// <<<<<<
  }
}
