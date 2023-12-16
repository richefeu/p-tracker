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

#define P_TRACKER_VERSION "0.8.0"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#define MAX_NUMBER_THREADS 48

#if defined(_OPENMP)
#include <omp.h>
#endif

// images
#include "libraw.h"
#include <tiffio.h>

// things from toofus
#include "powell_nr3.hpp"

// local includes
#include "NCC_optimizer_functor.hpp"
#include "grain_type_2D.hpp"
#include "image_data.hpp"
#include "interpol.hpp"
#include "relative_coord_type.hpp"
#include "search_zone_type.hpp"
#include "subpix_coord_type.hpp"
#include "grain_pair.hpp"
#include "EquiProj_optimizer_functor.hpp"
//#include "frame_1g2e.hpp"

class ptracker {
public:
  // === INTERFACE (FLAGS) ===
  std::string procedure{"particle_tracking"};
  int verbose_level{0};
  int rotations{1};
  int progress{0};

  // this allows to define special pattern according to a specified radius
  // A negative radius means "all the grains"
  double targetRadiusPattern = -1.0;

  // === IMAGE ===
  std::vector<std::vector<std::vector<uint16_t>>> image; // image[0 or 1][dimx][dimy]
  std::vector<imageDATA> imageData;                      // imageData[0 or 1]
  int im_index_ref{0};
  int im_index_current{1}; // Index of image origin and target in table image
  int RawImages{0};        // 0 = TIFF, 1 = RAW
  int DemosaicModel{0};    // for libraw
  std::map<int, std::string> DemosaicModelName;
  int rescaleGrayLevels{0};
  bool require_precomputations{true}; // If true, remake the pre-computations
  int dimx{0};                        // image width
  int dimy{0};                        // image height

  // === FILES ===
  char image_name[256]; // Image name in "printf" style. Example: "./quelquepart/toto%04d.tif"
  char dic_name[256];   // DIC-file name in "printf" style. Example: "./quelquepart/dic_out_%d.txt"
  int iref{0};          // Control of file numbers
  int ibeg{0};          // Control of file numbers
  int iend{0};          // Control of file numbers
  int iinc{0};          // Control of file numbers
  int iraz{0};          // Control of file numbers
  int idelta{0};        // Control of file numbers

  // === TRACKED POINTS / PATTERNS ===
  std::vector<grain_type_2D> grain;
  int num_grains = 0;

  // === RESCUE ===
  int rescue_level{2}; // in range [0 2].
                       // A negative value means no pixel-resolution tracking (go directly to subpixel optimization)
  std::vector<int> to_be_rescued;
  int num_to_be_rescued{0};
  double NCC_min{0.7}; // The NCC-value above which a first rescue is mandatory

  // === SUPER-RESCUES ===
  std::vector<int> to_be_super_rescued;
  int num_to_be_super_rescued{0};
  double NCC_min_super{0.5}; // The NCC-value above which a super_rescue is mandatory

  // === SUB-PIXEL ===
  int subpixel{1};
  double subpix_tol{1e-8};
  double initial_direction_dx{0.01};
  double initial_direction_dy{0.01};
  double initial_direction_dr{0.01};
  double subpix_displacement_max{3.0};
  double overflow_stiffness{0.0};

  image_interpLinear IMAGE_INTERPOLATOR_LINEAR;
  image_interpCubic IMAGE_INTERPOLATOR_CUBIC;
  image_interpQuintic IMAGE_INTERPOLATOR_QUINTIC;
  image_interpolator *IMAGE_INTERPOLATOR{&IMAGE_INTERPOLATOR_CUBIC};

  // === MULTI-THREADING ===
  int wanted_num_threads{4};

  // A table that hold the grain number being processed for each thread.
  // It is required in subpixel tracking.
  // Do not hesitate to increase the maximum number of threads if necessary
  int igrain_of_thread[MAX_NUMBER_THREADS];

  // === NEIGHBORS ===
  int use_neighbour_list{1};
  int num_neighbour_max{40};
  double neighbour_dist_pix{250.0}; // max distance between points to find neighbours
  int period_rebuild_neighbour_list{50};

  // === SEARCH ZONES ===
  search_zone_type search_zone;              // For non sub-pixel tracking
  search_zone_type search_zone_rescue;       // For points with bad mathing
  search_zone_type search_zone_super_rescue; // For points with very bad mathing

  // === Correction of distortion ===
  std::vector<grain_pair> grain_pairs;
  double equiProj_dist_min{0.0};
  double equiProj_dist_max{1000.0};
  
  std::vector<int> image_numbers_corrDisto;
  std::vector<double> imposed_displ_x_pix, imposed_displ_y_pix;
  std::vector<double> orig_x, orig_y;
  std::vector<std::vector<double>> dx_corrDisto, dy_corrDisto, NCC_subpix_corrDisto;
  std::vector<double> disto_parameters;
  std::vector<double> disto_parameters_perturb;
  // 0  -> xc (x center)
  // 1  -> yc (y center)
  // 2  -> K1
  // 3  -> K2
  // 4  -> K3
  // 5  -> P1
  // 6  -> P2
  // 7  -> P3
  // 8  -> xc parallax
  // 9  -> yc parallax
  // 10  -> theta1 (with respect to horizontal)
  // 11  -> parallax1 (0.0 coresponds to no correction in direction 1)
  // 12 -> theta2 (with respect to vertical)
  // 13 -> parallax2 (0.0 coresponds to no correction in direction 2)
  
  // === Check distortion ===
  int nbseg{2};

  // === 1g2e-related things ===
  //frame_1g2e frame;
  
  // === METHODS ===
  ptracker();

  void header();
  int init(int argc, char *argv[]);

  // command file
  int read_data(const char *name);
  void run();

  // images
  void read_image(int i, int num, bool first_time = false);
  void read_tiff_image(int i, const char *name, bool first_time = false);
  void read_raw_image(int i, const char *name, bool first_time = false);

  // Procedures
  void particle_tracking();
  void correction_distortion();
  void check_distortion();
  void compute_fluctuations(); 
  
  // particle tracking
  void do_precomputations();
  void follow_pattern_pixel(int igrain);
  void follow_pattern_rescue_pixel(int igrain);
  void follow_pattern_super_rescue_pixel(int igrain);
  void follow_pattern_subpixel_xyR(int igrain);
  double NCC_to_minimize_xyR(std::vector<double> &X);
  
  // correction of distortion
  void undistor(std::vector<double> &X, const double xd, const double yd, double &xu, double &yu);
  void precompute_paires();
  double disto_to_minimize_Equiproj(std::vector<double> &X);
  
  // check distortion
  
  // compute fluctuations
  
  // other
  void make_grid(double xmin, double xmax, double ymin, double ymax, int nx, int ny, int aleaMax);
  
  int read_grains(const char *name, bool ptracker_format = true);
  int save_grains(int num);
  int save_grains(const char *name, int num, bool simpleVersion = false);
  
  void find_neighbours(int igrain);
  void make_circ_pattern(int igrain, int radius_pix);
  void make_ring_pattern(int igrain, int radius_IN_pix, int radius_OUT_pic);
  void make_rect_pattern(int igrain, int half_width_pix, int half_height_pix);
  void make_custom_pattern(const char *name);

  void mask_rect(int xmin, int xmax, int ymin, int ymax);

private:
  void rotate_pixel_pattern(int igrain, int i, double c, double s, int *xpixel, int *ypixel);
  double get_time();
  void loadbar(size_t x, size_t n, size_t w = 50);
  std::string timestamp2string(time_t rawtime);
};
