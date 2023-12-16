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

#pragma once

class ptracker;
 
struct infered_grain {
  double xref;
  double yref;
  
  double x;
  double y;
  double rot;
  
  double dx;
  double dy;
  double drot;
  
  double xhom;
  double yhom;
};

struct sample {
  std::vector <infered_grain> Scheenb;
};

struct frame_1g2e {
  
  ptracker * ptrk;
  bool is_1g2e{false};
  
  size_t iBL{0};
  size_t iBR{0};
  size_t iTL{0};
  size_t iTR{0};
  
  double xBL, yBL;
  double xBR, yBR;
  double xTR, yTR;
  double xTL, yTL;
  
  double dxBL, dyBL;
  double dxBR, dyBR;
  double dxTR, dyTR;
  double dxTL, dyTL;
  
  std::vector<size_t> ID_pts_fix;
  std::vector<size_t> ID_pts_cam_movs;
  
  
  
  double tolerance{1e-10};

  void try_to_plug(ptracker * t_ptrk);
  
  // fluctuations
  int inverseBilinear(double x, double y, double &sout, double &tout, double &s2out, double &t2out);
  void bilinear(double s, double t, double &x, double &y);
  
  // rescale
  
  // correction of camera micro-movements
  
  // 
  

private:
  bool equals(double a, double b);
  double cross2(double x0, double y0, double x1, double y1);
  bool in_range(double val, double range_min, double range_max);
};
