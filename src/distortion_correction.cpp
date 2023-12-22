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

// http://en.wikipedia.org/wiki/Distortion_(optics)
void ptracker::undistor(std::vector<double>& X, const double xd, const double yd, double& xu, double& yu) {
  double Xc = X[0];
  double Yc = X[1];
  double K1 = X[2];
  double K2 = X[3];
  double K3 = X[4];
  double P1 = X[5];
  double P2 = X[6];
  double P3 = X[7];

  double dx = xd - Xc;
  double dy = yd - Yc;
  double r2 = dx * dx + dy * dy;
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  // Brown's lens distortion formula
  // This is what is written in the wikipedia page (December, 2023)
  xu = xd + dx * (K1 * r2 + K2 * r4 + K3 * r6) + (P1 * (r2 + 2.0 * dx * dx) + 2.0 * P2 * dx * dy) * (1.0 + P3 * r2);
  yu = yd + dy * (K1 * r2 + K2 * r4 + K3 * r6) + (2.0 * P1 * dx * dy + P2 * (r2 + 2.0 * dy * dy)) * (1.0 + P3 * r2);

  // ======================================

  double Xcp = X[8];
  double Ycp = X[9];
  double theta1 = X[10];
  double parallax1 = X[11];
  double theta2 = X[12];
  double parallax2 = X[13];

  // Parallax correction of the undistorted positions
  double nx1 = cos(theta1);
  double ny1 = -sin(theta1);
  double nx2 = -sin(theta2);
  double ny2 = cos(theta2);
  double dX = xu - Xcp;
  double dY = yu - Ycp;
  double S1 = (1.0 + parallax1) * (dX * nx1 + dY * ny1);
  double S2 = (1.0 + parallax2) * (dX * nx2 + dY * ny2);
  xu = Xcp + S1 * nx1 + S2 * nx2;
  yu = Ycp + S1 * ny1 + S2 * ny2;
}

// This function constructs pairs of points on a grid that are more than
// equiProj_dist_min and less than equiProj_dist_max apart.
void ptracker::precompute_paires() {
  grain_pairs.clear();

  for (int igrain = 0; igrain < num_grains; igrain++) {
    double xi = grain[igrain].refcoord_xpix;
    double yi = grain[igrain].refcoord_ypix;
    for (int jgrain = igrain + 1; jgrain < num_grains; jgrain++) {
      double xj = grain[jgrain].refcoord_xpix;
      double yj = grain[jgrain].refcoord_ypix;
      double xij = (xj - xi);
      double yij = (yj - yi);
      double dij = sqrt(xij * xij + yij * yij);
      if (dij >= equiProj_dist_min && dij <= equiProj_dist_max) {
        grain_pairs.push_back(grain_pair(igrain, jgrain));
      }
    }
  }

  std::cout << "Number of points for equiprojectivity computations: " << grain_pairs.size() << std::endl;
}

/*
// The function to minimise if we want to use the equiprojectivity criterion.
double ptracker::disto_to_minimize_Equiproj(std::vector<double>& X) {
  double A1x, A1y, B1x, B1y;  // undistorted positions on image 1
  double A2x, A2y, B2x, B2y;  // undistorted positions on image 2
  double AAx, AAy, BBx, BBy, ABx, ABy, AB;

  double func = 0.0;
  for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {
    for (size_t il = 0; il < grain_pairs.size(); il++) {
      int igrain = grain_pairs[il].i;
      int jgrain = grain_pairs[il].j;
      if (NCC_subpix_corrDisto[i_image][igrain] < NCC_min || NCC_subpix_corrDisto[i_image][jgrain] < NCC_min) continue;

      // undistor(X, (double)(grain[igrain].refcoord_xpix),
      // (double)(grain[igrain].refcoord_ypix), A1x, A1y);
      A1x = orig_x[igrain];
      A1y = orig_y[igrain];
      undistor(X, (double)(grain[igrain].refcoord_xpix + dx_corrDisto[i_image][igrain]),
               (double)(grain[igrain].refcoord_ypix + dy_corrDisto[i_image][igrain]), A2x, A2y);

      // undistor(X, (double)(grain[jgrain].refcoord_xpix),
      // (double)(grain[jgrain].refcoord_ypix), B1x, B1y);
      B1x = orig_x[jgrain];
      B1y = orig_y[jgrain];
      undistor(X, (double)(grain[jgrain].refcoord_xpix + dx_corrDisto[i_image][jgrain]),
               (double)(grain[jgrain].refcoord_ypix + dy_corrDisto[i_image][jgrain]), B2x, B2y);

      ABx = B1x - A1x;
      ABy = B1y - A1y;
      AB = sqrt(ABx * ABx + ABy * ABy);
      ABx /= AB;
      ABy /= AB;
      BBx = B2x - B1x;
      BBy = B2y - B1y;
      AAx = A2x - A1x;
      AAy = A2y - A1y;

      // version sum(err^2)
      double err = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy));
      func += err * err;
    }
  }

  return (func);
}
*/

// The procedure to correct the image distortions
void ptracker::correction_distortion() {
  std::cout << '\n';
  std::cout << "*** PROCEDURE CORRECTION OF DISTORTION ***" << '\n';
  std::cout << std::endl;

  // Reserve memory
  dx_corrDisto.resize(image_numbers_corrDisto.size());
  dy_corrDisto.resize(image_numbers_corrDisto.size());
  NCC_subpix_corrDisto.resize(image_numbers_corrDisto.size());
  for (size_t i = 0; i < image_numbers_corrDisto.size(); i++) {
    dx_corrDisto[i].resize(num_grains);
    dy_corrDisto[i].resize(num_grains);
    NCC_subpix_corrDisto[i].resize(num_grains);
  }

  orig_x.resize(num_grains);
  orig_y.resize(num_grains);
  for (int igrain = 0; igrain < num_grains; igrain++) {
    orig_x[igrain] = (double)(grain[igrain].refcoord_xpix);
    orig_y[igrain] = (double)(grain[igrain].refcoord_ypix);
  }

  int igrain;
  double index_ini = 0.0, index_end = 0.0, index_gain = 0.0;

  std::ofstream logfile("corrDisto_log.txt");

  read_image(0, image_numbers_corrDisto[0], true);  // Read reference image
  iref = image_numbers_corrDisto[0];
  do_precomputations();
  precompute_paires();

  std::cout << "Tracking displacements... " << std::endl;
  double tbeg;
  for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {

    read_image(1, image_numbers_corrDisto[i_image]);

    // Accounting for the estimated displacement
    // (it allows the correlation to perform faster, with a smaller search zone)
    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].dx = imposed_displ_x_pix[i_image];
      grain[igrain].dy = imposed_displ_y_pix[i_image];
    }

    // ======== Performing a DIC
    fprintf(stdout, "> Follow %d grains... (from image number %d to number %d)\n", num_grains, iref,
            image_numbers_corrDisto[i_image]);
    fflush(stdout);
    tbeg = get_time();

    progress = 0;
#pragma omp parallel for
    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].reset();
      follow_pattern_pixel(igrain);
#pragma omp critical
      { loadbar(++progress, grain.size()); }
    }
    std::cerr << std::endl;
    fprintf(stdout, "< [DONE in %f seconds]\n", get_time() - tbeg);

    // NO RESCUE IS ATTEMPTED !!

    // ======== Sub-pixel
    fprintf(stdout, "> Sub-pixel resolution... \n");
    fflush(stdout);
    tbeg = get_time();

    progress = 0;
#pragma omp parallel for
    for (igrain = 0; igrain < num_grains; igrain++) {
      follow_pattern_subpixel_xyR(igrain);
#pragma omp critical
      { loadbar(++progress, grain.size()); }
    }
    std::cerr << std::endl;
    fprintf(stdout, "< [DONE in %f seconds]\n", get_time() - tbeg);

    int nbFailed = 0;
    for (igrain = 0; igrain < num_grains; igrain++) {
      dx_corrDisto[i_image][igrain] = grain[igrain].dx;
      dy_corrDisto[i_image][igrain] = grain[igrain].dy;
      NCC_subpix_corrDisto[i_image][igrain] = grain[igrain].NCC_subpix;
      if (grain[igrain].NCC_subpix < NCC_min) nbFailed++;
    }
    std::cout << "Failure rate = " << 100.0 * (double)nbFailed / (double)num_grains << "%" << std::endl;

    // save_grains(i_image);
    save_grains(image_numbers_corrDisto[i_image]);

  }  // for loop i_image
  std::cout << "Tracking of points done" << std::endl;

  // The parameters (initial guess) and perturbations are set by the user
  if (disto_parameters[0] < 0.0) disto_parameters[0] = 0.5 * (double)dimx;
  if (disto_parameters[1] < 0.0) disto_parameters[1] = 0.5 * (double)dimy;
  if (disto_parameters[8] < 0.0) disto_parameters[8] = 0.5 * (double)dimx;
  if (disto_parameters[9] < 0.0) disto_parameters[9] = 0.5 * (double)dimy;

  // Minimization
  int npairs = grain_pairs.size();
  int nb_image_pairs = image_numbers_corrDisto.size() - 1;
  std::cout << "Numbers of position pairs used to undistor: " << npairs << std::endl;
  std::cout << "Numbers of image pairs: " << nb_image_pairs << std::endl;

  if (optim_algo == POWELL) {
    
    // ==================================
    //.    POWELL
    // ==================================

    // With this functor, both origin and target points are undistorted
    struct opti_functor {
      ptracker* ptrk;
      double operator()(std::vector<double>& X) {
        double A1x, A1y, B1x, B1y;  // undistorted positions on image 1
        double A2x, A2y, B2x, B2y;  // undistorted positions on image 2
        double AAx, AAy, BBx, BBy, ABx, ABy, AB;

        double valret = 0.0;
        for (size_t i_image = 1; i_image < ptrk->image_numbers_corrDisto.size(); i_image++) {
          for (size_t il = 0; il < ptrk->grain_pairs.size(); il++) {
            int igrain = ptrk->grain_pairs[il].i;
            int jgrain = ptrk->grain_pairs[il].j;
            if (ptrk->NCC_subpix_corrDisto[i_image][igrain] < ptrk->NCC_min ||
                ptrk->NCC_subpix_corrDisto[i_image][jgrain] < ptrk->NCC_min)
              continue;

            ptrk->undistor(X, (double)(ptrk->grain[igrain].refcoord_xpix), (double)(ptrk->grain[igrain].refcoord_ypix),
                           A1x, A1y);
            ptrk->undistor(X, (double)(ptrk->grain[igrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][igrain]),
                           (double)(ptrk->grain[igrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][igrain]), A2x, A2y);

            ptrk->undistor(X, (double)(ptrk->grain[jgrain].refcoord_xpix), (double)(ptrk->grain[jgrain].refcoord_ypix),
                           B1x, B1y);
            ptrk->undistor(X, (double)(ptrk->grain[jgrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][jgrain]),
                           (double)(ptrk->grain[jgrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][jgrain]), B2x, B2y);

            ABx = B1x - A1x;
            ABy = B1y - A1y;
            AB = sqrt(ABx * ABx + ABy * ABy);
            ABx /= AB;
            ABy /= AB;
            BBx = B2x - B1x;
            BBy = B2y - B1y;
            AAx = A2x - A1x;
            AAy = A2y - A1y;

            double err = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy));
            valret += err * err;
          }
        }

        return valret;
      }
    };

    opti_functor func;
    func.ptrk = this;
    Powell<opti_functor> powell(func, 1e-12);
    index_ini = func(disto_parameters) / (double)(nb_image_pairs);
    std::cout << "initial index = " << index_ini << std::endl;

    std::cout << "START minimizing equiprojectivity criterion (with POWELL)" << std::endl;
    disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);
    std::cout << "END minimizing equiprojectivity criterion (with POWELL)\n" << std::endl;

    index_end = func(disto_parameters) / (double)(nb_image_pairs);
    index_gain = index_ini / index_end;

  } else if (optim_algo == POWELL_ITER) {
    
    // ==================================
    //.    POWELL_ITER
    // ==================================

    struct opti_functor {
      ptracker* ptrk;
      double operator()(std::vector<double>& X) {
        double A1x, A1y, B1x, B1y;  // undistorted positions on image 1
        double A2x, A2y, B2x, B2y;  // undistorted positions on image 2
        double AAx, AAy, BBx, BBy, ABx, ABy, AB;

        double valret = 0.0;
        for (size_t i_image = 1; i_image < ptrk->image_numbers_corrDisto.size(); i_image++) {
          for (size_t il = 0; il < ptrk->grain_pairs.size(); il++) {
            int igrain = ptrk->grain_pairs[il].i;
            int jgrain = ptrk->grain_pairs[il].j;
            if (ptrk->NCC_subpix_corrDisto[i_image][igrain] < ptrk->NCC_min ||
                ptrk->NCC_subpix_corrDisto[i_image][jgrain] < ptrk->NCC_min)
              continue;

            A1x = ptrk->orig_x[igrain];
            A1y = ptrk->orig_y[igrain];
            ptrk->undistor(X, (double)(ptrk->grain[igrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][igrain]),
                           (double)(ptrk->grain[igrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][igrain]), A2x, A2y);

            B1x = ptrk->orig_x[jgrain];
            B1y = ptrk->orig_y[jgrain];
            ptrk->undistor(X, (double)(ptrk->grain[jgrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][jgrain]),
                           (double)(ptrk->grain[jgrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][jgrain]), B2x, B2y);

            ABx = B1x - A1x;
            ABy = B1y - A1y;
            AB = sqrt(ABx * ABx + ABy * ABy);
            ABx /= AB;
            ABy /= AB;
            BBx = B2x - B1x;
            BBy = B2y - B1y;
            AAx = A2x - A1x;
            AAy = A2y - A1y;

            double err = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy));
            valret += err * err;
          }
        }

        return valret;
      }
    };

    opti_functor func;
    func.ptrk = this;
    Powell<opti_functor> powell(func, 1e-12);

    index_ini = func(disto_parameters) / (double)(nb_image_pairs);
    std::cout << "initial index = " << index_ini << std::endl;

    std::cout << "START minimizing equiprojectivity criterion (with POWELL_ITER)" << std::endl;
    std::ofstream convergeFile("convergence.txt");
    std::cout << std::scientific << std::setprecision(15);
    convergeFile << std::scientific << std::setprecision(15);
    std::vector<std::vector<double>> dx0_corrDisto = dx_corrDisto;
    std::vector<std::vector<double>> dy0_corrDisto = dy_corrDisto;
    for (int iter = 1; iter <= powell_num_iter; iter++) {
      disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);
      double idx = func(disto_parameters) / (double)(nb_image_pairs);
      std::cout << "iteration " << iter << ", index = " << idx << std::endl;
      convergeFile << iter << ' ' << idx << ' ' << disto_parameters[0] << ' ' << disto_parameters[1] << ' '
                   << disto_parameters[2] << ' ' << disto_parameters[3] << ' ' << disto_parameters[4] << ' '
                   << disto_parameters[5] << ' ' << disto_parameters[6] << ' ' << disto_parameters[7] << ' '
                   << disto_parameters[8] << ' ' << disto_parameters[9] << ' ' << disto_parameters[10] << ' '
                   << disto_parameters[11] << ' ' << disto_parameters[12] << ' ' << disto_parameters[13] << ' '
                   << std::endl;

      double xu, yu, dxcorr, dycorr;
      for (int igrain = 0; igrain < num_grains; igrain++) {
        undistor(disto_parameters, (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), xu,
                 yu);
        dxcorr = xu - (double)(grain[igrain].refcoord_xpix);
        dycorr = yu - (double)(grain[igrain].refcoord_ypix);
        for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {
          dx_corrDisto[i_image][igrain] = dx0_corrDisto[i_image][igrain] + dxcorr;
          dy_corrDisto[i_image][igrain] = dy0_corrDisto[i_image][igrain] + dycorr;
        }
        orig_x[igrain] = xu;
        orig_y[igrain] = yu;
      }
    }
    std::cout << "END minimizing equiprojectivity criterion (with POWELL_ITER)\n" << std::endl;

    index_end = func(disto_parameters) / (double)(nb_image_pairs);
    index_gain = index_ini / index_end;

  } else if (optim_algo == LEVMAR) {
    
    // ==================================
    //.    LEVENBERG-MARQUARDT
    // ==================================

    const int m = (int)((image_numbers_corrDisto.size() - 1) * grain_pairs.size());
    const int n = (int)disto_parameters.size();

    int lwa = levmar::get_lwa(m, n);
    std::vector<int> iwa(n);
    std::vector<double> fvec(m);
    std::vector<double> wa(lwa);

    auto fcn = [](void* p, int m, int n, const double* x, double* fvec, int iflag) {
      // remove warning during compilation
      (void)iflag;
      (void)m;

      ptracker* ptrk = (ptracker*)p;
      std::vector<double> X(n);
      for (int ip = 0; ip < n; ip++) {
        X[ip] = x[ip];
      }

      double A1x, A1y, B1x, B1y;  // undistorted positions on image 1
      double A2x, A2y, B2x, B2y;  // undistorted positions on image 2
      double AAx, AAy, BBx, BBy, ABx, ABy, AB;

      size_t n_images = ptrk->image_numbers_corrDisto.size();
      size_t n_pairs = ptrk->grain_pairs.size();
      double retval = 0.0;
      for (size_t i_image = 1; i_image < n_images; i_image++) {
        for (size_t il = 0; il < n_pairs; il++) {
          int igrain = ptrk->grain_pairs[il].i;
          int jgrain = ptrk->grain_pairs[il].j;
          if (ptrk->NCC_subpix_corrDisto[i_image][igrain] < ptrk->NCC_min ||
              ptrk->NCC_subpix_corrDisto[i_image][jgrain] < ptrk->NCC_min)
            continue;

          ptrk->undistor(X, (double)(ptrk->grain[igrain].refcoord_xpix), (double)(ptrk->grain[igrain].refcoord_ypix),
                         A1x, A1y);
          ptrk->undistor(X, (double)(ptrk->grain[igrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][igrain]),
                         (double)(ptrk->grain[igrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][igrain]), A2x, A2y);

          ptrk->undistor(X, (double)(ptrk->grain[jgrain].refcoord_xpix), (double)(ptrk->grain[jgrain].refcoord_ypix),
                         B1x, B1y);
          ptrk->undistor(X, (double)(ptrk->grain[jgrain].refcoord_xpix + ptrk->dx_corrDisto[i_image][jgrain]),
                         (double)(ptrk->grain[jgrain].refcoord_ypix + ptrk->dy_corrDisto[i_image][jgrain]), B2x, B2y);

          ABx = B1x - A1x;
          ABy = B1y - A1y;
          AB = sqrt(ABx * ABx + ABy * ABy);
          ABx /= AB;
          ABy /= AB;
          BBx = B2x - B1x;
          BBy = B2y - B1y;
          AAx = A2x - A1x;
          AAy = A2y - A1y;
          double err = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy));
            
          fvec[(i_image - 1) * n_pairs + il] = err;
          retval += err * err;
        }
      }

      return retval;
    };

    index_ini = fcn(this, m, n, &disto_parameters[0], &fvec[0], 0) / (double)(nb_image_pairs);
    std::cout << "initial index = " << index_ini << std::endl;

    std::cout << "START minimizing equiprojectivity criterion (with LEVMAR)" << std::endl;
    int info = levmar::lmdif1(fcn, this, m, n, &disto_parameters[0], &fvec[0], DBL_EPSILON, &iwa[0], &wa[0], lwa);
    std::cout << "END minimizing equiprojectivity criterion (with LEVMAR)\n" << std::endl;
    std::cout << "exit state: " << levmar::getInfo(info) << std::endl;

    index_end = fcn(this, m, n, &(disto_parameters[0]), &fvec[0], 0) / (double)(nb_image_pairs);
    index_gain = index_ini / index_end;
  }
  
  // ====================================================
  // OPTIMISATION DONE 
  // ====================================================

  logfile << std::scientific << std::setprecision(15);
  logfile << "Initial distorsion index " << index_ini << '\n';
  logfile << "Final distorsion index " << index_end << '\n';
  logfile << "Index gain (>1 is good) " << index_gain << '\n';
  logfile << "xc_corrDistor " << disto_parameters[0] << '\n';
  logfile << "yc_corrDistor " << disto_parameters[1] << '\n';
  logfile << "K1_corrDistor " << disto_parameters[2] << '\n';
  logfile << "K2_corrDistor " << disto_parameters[3] << '\n';
  logfile << "K3_corrDistor " << disto_parameters[4] << '\n';
  logfile << "P1_corrDistor " << disto_parameters[5] << '\n';
  logfile << "P2_corrDistor " << disto_parameters[6] << '\n';
  logfile << "P3_corrDistor " << disto_parameters[7] << '\n';
  logfile << "xc_corrParallax     " << disto_parameters[8] << '\n';
  logfile << "yc_corrParallax     " << disto_parameters[9] << '\n';
  logfile << "angle1_corrParallax " << disto_parameters[10] << '\n';
  logfile << "A1_corrParallax     " << disto_parameters[11] << '\n';
  logfile << "angle2_corrParallax " << disto_parameters[12] << '\n';
  logfile << "A2_corrParallax     " << disto_parameters[13] << std::endl;

  std::cout << std::scientific << std::setprecision(15);
  std::cout << "    Initial distorsion index " << index_ini << '\n';
  std::cout << "      Final distorsion index " << index_end << '\n';
  std::cout << "gain (index_ini / index_end) " << index_gain << " (>1 is good)" << '\n';
  std::cout << "xc_corrDistor " << disto_parameters[0] << '\n';
  std::cout << "yc_corrDistor " << disto_parameters[1] << '\n';
  std::cout << "K1_corrDistor " << disto_parameters[2] << '\n';
  std::cout << "K2_corrDistor " << disto_parameters[3] << '\n';
  std::cout << "K3_corrDistor " << disto_parameters[4] << '\n';
  std::cout << "P1_corrDistor " << disto_parameters[5] << '\n';
  std::cout << "P2_corrDistor " << disto_parameters[6] << '\n';
  std::cout << "P3_corrDistor " << disto_parameters[7] << '\n';
  std::cout << "xc_corrParallax     " << disto_parameters[8] << '\n';
  std::cout << "yc_corrParallax     " << disto_parameters[9] << '\n';
  std::cout << "angle1_corrParallax " << disto_parameters[10] << '\n';
  std::cout << "A1_corrParallax     " << disto_parameters[11] << '\n';
  std::cout << "angle2_corrParallax " << disto_parameters[12] << '\n';
  std::cout << "A2_corrParallax     " << disto_parameters[13] << std::endl;

  // === CREATE DATA FILE FOR PLOTTING THE FIELD OF CORRECTION DISPLACEMENTS ===
  //
  // in gnuplot
  //
  // ampl = 10.0
  // set size ratio -1
  // plot "errors.txt" u 1:(-$2):($3*ampl):(-$4*ampl) w vec lc 3 notitle

  std::ofstream errorFile("errors.txt");
  double xerrormax = -1e20;
  double yerrormax = -1e20;
  double errormax = -1e20;
  double xd, yd, xu, yu;
  double dxd = (double)dimx / 50.0;
  double dyd = (double)dimy / 50.0;

  if (dxd < 1.0) {
    dxd = 1.0;
  }
  if (dyd < 1.0) {
    dyd = 1.0;
  }

  for (xd = 0.0; xd < dimx; xd += dxd) {
    for (yd = 0.0; yd < dimy; yd += dyd) {
      undistor(disto_parameters, xd, yd, xu, yu);
      double xerr = xu - xd;
      double yerr = yu - yd;
      double err = sqrt(xerr * xerr + yerr * yerr);

      if (fabs(xerr) > xerrormax) {
        xerrormax = fabs(xerr);
      }
      if (fabs(yerr) > yerrormax) {
        yerrormax = fabs(yerr);
      }
      if (err > errormax) {
        errormax = err;
      }

      errorFile << xd << ' ' << yd << ' ' << xerr << ' ' << yerr << std::endl;
    }
    errorFile << std::endl;
  }
  logfile << "max |xu - xd| " << xerrormax << '\n';
  logfile << "max |yu - yd| " << yerrormax << '\n';
  logfile << "errormax " << errormax << std::endl;
  std::cout << "max |xu - xd| " << xerrormax << '\n';
  std::cout << "max |yu - yd| " << yerrormax << '\n';
  std::cout << "errormax      " << errormax << std::endl;

  // === CREATE DATA FILE FOR PLOTTING THE (UN)DISTORTED FRAME ===
  //
  // in gnuplot
  //
  // ampl = 10.0
  // set size ratio -1
  // plot "distobox.txt" u 1:(-$2) w l lw 2 lc 1 notitle, "" u
  // ($1+$3*ampl):(-$2-$4*ampl) w l lw 2 lc 1 notitle

  std::ofstream errorBoxFile("distobox.txt");
  yd = 0;
  for (xd = 0; xd < dimx; xd += dxd) {
    undistor(disto_parameters, xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = dimx - 1;
  for (yd = 0; yd < dimy; yd += dyd) {
    undistor(disto_parameters, xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  yd = dimy - 1;
  for (xd = dimx - 1; xd >= 0; xd -= dxd) {
    undistor(disto_parameters, xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = 0;
  for (yd = dimy - 1; yd >= 0; yd -= dyd) {
    undistor(disto_parameters, xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }

  xd = yd = 0;
  undistor(disto_parameters, xd, yd, xu, yu);
  errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
}
