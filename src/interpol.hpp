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

#include <cmath>
#include <cstdint>
#include <vector>

struct image_interpolator {
  // In this pur virtual function, x and y are subpixel coordinates
  virtual double getValue(std::vector<std::vector<uint16_t>> &im, double, double) = 0;
};

struct image_interpLinear : public image_interpolator {
  double LinearInterpolate(double p[2], double x);
  double getValue(std::vector<std::vector<uint16_t>> &im, double x, double y);
};

struct image_interpCubic : public image_interpolator {
  double CubicInterpolate(double p[4], double x);
  double getValue(std::vector<std::vector<uint16_t>> &im, double x, double y);
};

struct image_interpQuintic : public image_interpolator {
  double QuinticInterpolate(double p[4], double x);
  double getValue(std::vector<std::vector<uint16_t>> &im, double x, double y);
};

