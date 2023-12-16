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

#include "frame_1g2e.hpp"
#include "p-tracker.hpp"

// TODO rename frame_1g2e

void frame_1g2e::try_to_plug(ptracker *t_ptrk) {
  ptrk = t_ptrk;

  std::vector<size_t> cornerIds;
  for (int i = 0; i < ptrk->num_grains; i++) {
    // R = 1.0 -> corner
    // R = 2.1 -> scaling point
    // R = 2.2 -> non-moving point
    if (ptrk->grain[i].radius_pix == 1.0) {
      cornerIds.push_back((size_t)i);
    }
  }

  if (cornerIds.size() != 4) {
    is_1g2e = false;
    return;
  }

  is_1g2e = true;
  // sort from bottom to top (remember increasing y if towards the bottom)
  std::sort(cornerIds.begin(), cornerIds.end(),
            [ptrk](size_t a, size_t b) { return (ptrk->grain[a].refcoord_ypix > ptrk->grain[b].refcoord_ypix); });
  
  if (ptrk->grain[cornerIds[0]].refcoord_xpix > ptrk->grain[cornerIds[1]].refcoord_xpix) {
    size_t tmp = cornerIds[0];
    cornerIds[0] = cornerIds[1];
    cornerIds[1] = tmp;
  }
  
  if (ptrk->grain[cornerIds[2]].refcoord_xpix > ptrk->grain[cornerIds[2]].refcoord_xpix) {
    size_t tmp = cornerIds[2];
    cornerIds[2] = cornerIds[3];
    cornerIds[3] = tmp;
  }
  
  iBL = cornerIds[0];
  iBR = cornerIds[1];
  iTL = cornerIds[2];
  iTR = cornerIds[3];
  
  // here the reference image is supposed to have never changed 
  dxBL = ptrk->grain[iBL].dx; 
  dxBR = ptrk->grain[iBR].dx; 
  dxTL = ptrk->grain[iTL].dx; 
  dxTR = ptrk->grain[iTR].dx; 
  
  dyBL = ptrk->grain[iBL].dy; 
  dyBR = ptrk->grain[iBR].dy; 
  dyTL = ptrk->grain[iTL].dy; 
  dyTR = ptrk->grain[iTR].dy;
  
  xBL = (double)(ptrk->grain[iBL].refcoord_xpix) + dxBL; 
  xBR = (double)(ptrk->grain[iBR].refcoord_xpix) + dxBR; 
  xTL = (double)(ptrk->grain[iTL].refcoord_xpix) + dxTL; 
  xTR = (double)(ptrk->grain[iTR].refcoord_xpix) + dxTR; 
  
  yBL = (double)(ptrk->grain[iBL].refcoord_ypix) + dyBL; 
  yBR = (double)(ptrk->grain[iBR].refcoord_ypix) + dyBR; 
  yTL = (double)(ptrk->grain[iTL].refcoord_ypix) + dyTL; 
  yTR = (double)(ptrk->grain[iTR].refcoord_ypix) + dyTR; 
}

bool frame_1g2e::equals(double a, double b) {
  // return (a == b) || ((a <= (b + tolerance)) && (a >= (b - tolerance)));
  return (fabs(a - b) < tolerance);
}

double frame_1g2e::cross2(double xa, double ya, double xb, double yb) { return xa * yb - ya * xb; }

bool frame_1g2e::in_range(double val, double range_min, double range_max) {
  return ((val + tolerance) >= range_min) && ((val - tolerance) <= range_max);
}

// Returns number of solutions found.  If there is one valid solution, it will be put in s and t
//
//   pTL+------------+ pTR
//      |            |
//    t +     x (p)  |
//      |            |
//   pBL+-----+------+ pBR
//            s
int frame_1g2e::inverseBilinear(double x, double y, double &sout, double &tout, double &s2out, double &t2out) {
  int t_valid, t2_valid;

  double a = cross2(xBL - x, yBL - y, xBL - xTL, yBL - yTL);
  double b1 = cross2(xBL - x, yBL - y, xBR - xTR, yBR - yTR);
  double b2 = cross2(xBR - x, yBR - y, xBL - xTL, yBL - yTL);
  double c = cross2(xBR - x, yBR - y, xBR - xTR, yBR - yTR);
  double b = 0.5 * (b1 + b2);

  double s, s2, t, t2;

  double am2bpc = a - 2 * b + c;
  // this is how many valid s values we have
  int num_valid_s = 0;

  if (equals(am2bpc, 0.0)) {
    if (equals(a - c, 0.0)) {
      // Looks like the input is a line
      // You could set s=0.5 and solve for t if you wanted to
      return 0;
    }
    s = a / (a - c);
    if (in_range(s, 0.0, 1.0))
      num_valid_s = 1;
  } else {
    double sqrtbsqmac = sqrt(b * b - a * c);
    s = ((a - b) - sqrtbsqmac) / am2bpc;
    s2 = ((a - b) + sqrtbsqmac) / am2bpc;
    num_valid_s = 0;
    if (in_range(s, 0.0, 1.0)) {
      num_valid_s++;
      if (in_range(s2, 0.0, 1.0))
        num_valid_s++;
    } else {
      if (in_range(s2, 0.0, 1.0)) {
        num_valid_s++;
        s = s2;
      }
    }
  }

  if (num_valid_s == 0)
    return 0;

  t_valid = 0;
  if (num_valid_s >= 1) {
    double tdenom_x = (1 - s) * (xBL - xTL) + s * (xBR - xTR);
    double tdenom_y = (1 - s) * (yBL - yTL) + s * (yBR - yTR);
    t_valid = 1;
    if (equals(tdenom_x, 0.0) && equals(tdenom_y, 0.0)) {
      t_valid = 0;
    } else {
      // Choose the more robust denominator
      if (fabs(tdenom_x) > fabs(tdenom_y)) {
        t = ((1 - s) * (xBL - x) + s * (xBR - x)) / (tdenom_x);
      } else {
        t = ((1 - s) * (yBL - y) + s * (yBR - y)) / (tdenom_y);
      }
      if (!in_range(t, 0.0, 1.0))
        t_valid = 0;
    }
  }

  // Same thing for s2 and t2
  t2_valid = 0;
  if (num_valid_s == 2) {
    double tdenom_x = (1 - s2) * (xBL - xTL) + s2 * (xBR - xTR);
    double tdenom_y = (1 - s2) * (yBL - yTL) + s2 * (yBR - yTR);
    t2_valid = 1;
    if (equals(tdenom_x, 0.0) && equals(tdenom_y, 0.0)) {
      t2_valid = 0;
    } else {
      // Choose the more robust denominator
      if (fabs(tdenom_x) > fabs(tdenom_y)) {
        t2 = ((1 - s2) * (xBL - x) + s2 * (xBR - x)) / (tdenom_x);
      } else {
        t2 = ((1 - s2) * (yBL - y) + s2 * (yBR - y)) / (tdenom_y);
      }
      if (!in_range(t2, 0.0, 1.0))
        t2_valid = 0;
    }
  }

  // Final cleanup
  if (t2_valid && !t_valid) {
    s = s2;
    t = t2;
    t_valid = t2_valid;
    t2_valid = 0;
  }

  // Output
  if (t_valid) {
    sout = s;
    tout = t;
  }

  if (t2_valid) {
    s2out = s2;
    t2out = t2;
  }

  return t_valid + t2_valid;
}

void frame_1g2e::bilinear(double s, double t, double &x, double &y) {
  x = t * (s * dxTR + (1 - s) * dxTL) + (1 - t) * (s * dxBR + (1 - s) * dxBL);
  y = t * (s * dyTR + (1 - s) * dyTL) + (1 - t) * (s * dyBR + (1 - s) * dyBL);
}
