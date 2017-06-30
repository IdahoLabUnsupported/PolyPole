// Copyright 2017 Battelle Energy Alliance, LLC
// Licensed under the Mozilla Public License version 2.0

// @author: Davide Pizzocri (davide.pizzocri@gmail.com)
// @author: Giovanni Pastore (giovanni.pastore@inl.gov)
// @author: Cristian Rabiti (cristian.rabiti@inl.gov)

#include "PolyPole.h"
#include <vector>
#include <cmath>

using namespace std;

namespace polypole
{
  double PolyPole1( double dt_sec,
                    double grn_radius,
                    double temp0,
                    double temp1,
                    double frate0,
                    double frate1,
                    double grn_gas_initial )
  {
    // -- Input variables
    
    // dt_sec     : time step (s)
    // grn_radius : grain radius (m)
    // temp0      : temperature at the previous time step (K)
    // temp1      : temperature at the current time step (K)
    // frate0     : fission rate at the beginning of the time step (fissions/m^3-s)
    // frate1     : fission rate at the end of the time step (fissions/m^3-s)
    
    // -- Order of the corrective polynomial
    
    unsigned short int pol_grad = 2;
    unsigned short int G = pol_grad;

    // -- Sampling
    //    The diffusion and the gas production rates
    //    are sampled along the time step length

    unsigned short int M = pol_grad + 1;
    unsigned short int m;

    double T = dt_sec;
    double delta = T / double(M-1);
    double time[M];
    double temperature[M];
    double fissionrate[M];
    double fiss_source[M];
    double diff_coeff[M];

    double fiss_source_av = 0.0;
    double diff_coeff_av  = 0.0;

    for ( m = 0; m < M; m++ )
    {
      time[m] = delta * double(m);

      temperature[m] = ( temp1  - temp0 )  / T * ( time[m] - time[0] ) + temp0;
      fissionrate[m] = ( frate1 - frate0 ) / T * ( time[m] - time[0] ) + frate0;

      // The two following lines may have to be modified depending on the specific code
      // Call functions that return the gas production rate (atoms/m^3-s) and the
      // effective diffusion coefficient (m^2/s)

      fiss_source[m] = polypole::production_rate( fissionrate[m] );
      diff_coeff[m]  = polypole::diffusion_coeff( temperature[m] );

      // -- Average of the sampled values

      if ( m == 0 || m == M-1 )
      {
        fiss_source_av += fiss_source[m] / 2.0;
        diff_coeff_av  += diff_coeff[m]  / 2.0;
      }
      else
      {
        fiss_source_av += fiss_source[m];
        diff_coeff_av  += diff_coeff[m];
      }
    }

    fiss_source_av /= double(M-1);
    diff_coeff_av  /= double(M-1);

    // -- Calculation of the modes

    unsigned short int n_limit = 1000;
    unsigned short int n = 0;
    unsigned short int np1 = 0;
    double remainder = 0.0;
    const double tol = 1.0e-07;

    double grn_gas = 0.0;
    double grn_gas_old = 0.0;

    // For the mode variables below, the algorithm needs both current and previous time step values.
    // The way to store this information may depend on the structure of your code.
    // The maximum number of modes considered (1000) was found to be adequate through testing.
    // Changing it is possible, although not recommended. 
    static std::vector<double> mode(1000,0.0);
    static std::vector<double> mode_old(1000,0.0);

    //    Loop over modes

    const double pi = 3.141593;

    const double proj_coeff  = 2.0 * sqrt( 2.0 * pow(grn_radius, 3) / pi );
    const double eigen_coeff = pow(pi / grn_radius, 2);

    double fiss_source_mode[M];
    double diff_coeff_mode[M];

    for ( n = 0; n < n_limit; n++ )
    {
      np1 = n + 1;

      // Storage - initial conditions
      grn_gas_old = grn_gas;
      mode_old[n] = mode[n];

      // Source projection on modes
      // Diffusion coefficient by laplacian eigenvalues

      const double spatial_projection = - pow(-1, np1) / double(np1) * proj_coeff;
      const double laplacian_eigenvalue = pow(np1, 2) * eigen_coeff;

      for ( m = 0; m < M; m++ )
      {
        fiss_source_mode[m] = spatial_projection   * fiss_source[m];
        diff_coeff_mode[m]  = laplacian_eigenvalue * diff_coeff[m];
      }

      const double source = spatial_projection   * fiss_source_av;
      const double pole   = laplacian_eigenvalue * diff_coeff_av;

      // Asymptotic behavior

      const double asym_sol  = source / pole;
      const double asym_diff = mode_old[n] - asym_sol;

      // Cast of the matrix
      // Cast of the RHS

      double coeff_matrix[G*G];
      double rhs[G];
      unsigned short int gg = 0;
      unsigned short int g  = 0;
      unsigned short int col, row;

      for ( row = 1; row <= G; row++ )
      {
        double exp_row = asym_diff * exp(- pole * time[row]);

        for ( col = 1; col <= G; col++ )
        {
          const double time_col = pow(time[row], col);
          coeff_matrix[gg] = (asym_sol + exp_row) * double(col) * pow(time[row], col-1)
                           - pole * exp_row * time_col
                           + diff_coeff_mode[row] * (asym_sol + exp_row) * time_col;
          gg++;
        }

        rhs[g] = fiss_source_mode[row]
               + pole * exp_row
               - diff_coeff_mode[row] * (asym_sol + exp_row);
        g++;
      }

      // Inversion (2x2)

      const double a = coeff_matrix[0];
      const double b = coeff_matrix[1];
      const double c = coeff_matrix[2];
      const double d = coeff_matrix[3];

      const double  det = a * d - b * c;
            double _det = 0.0;

      if ( det != 0.0 )
      {
        _det = 1.0 / det;
        coeff_matrix[0] = + _det * d;
        coeff_matrix[1] = - _det * b;
        coeff_matrix[2] = - _det * c;
        coeff_matrix[3] = + _det * a;
      }

      // System solution

      double sol[2] = {rhs[0], rhs[1]};

      sol[0] = coeff_matrix[0] * rhs[0]
             + coeff_matrix[1] * rhs[1];

      sol[1] = coeff_matrix[2] * rhs[0]
             + coeff_matrix[3] * rhs[1];

      // -- Polynomial evaluation @ T

      double poly_value = 1.0;

      for ( g = 0; g < G; g++ )
        poly_value += sol[g] * pow(T, g+1);

      // -- Solution reconstruction

      mode[n] = (asym_sol + asym_diff * exp(-pole * T)) * poly_value;

      // -- Spatial average

      const double spatial_average = - 3.0 * pow(-1, np1) / (double(np1) * pi * grn_radius * sqrt(2.0 * pi * grn_radius));

      grn_gas += spatial_average * mode[n];

      // -- D'Alembert series remainder estimation

      remainder = std::abs(grn_gas - grn_gas_old) / grn_gas;
      if ( remainder < tol )  break;
    }

    // -- Reset all non-calculated modes to zero

    for ( unsigned short int k = n_limit-1; k > n; k-- )
      mode[k] = 0.0;

    // -- Check for convergence

    // If the maximum number of modes allowed is not enough to reach convergence
    // with the used tolerance, a simplified solution is used.
    // This can happen if the diffusion coefficient is very small.

    if( n == n_limit-1 && remainder > tol )
    {
      const double dtau = 15.0 * diff_coeff_av * T / pow(grn_radius, 2);
      const double grn_gas_asym = fiss_source_av * pow(grn_radius, 2) / (15.0 * diff_coeff_av);

      double fact = 1.0;
      if ( dtau < 1.0e-07 )
        fact = dtau;
      else
        fact = 1.0 - exp(-dtau);

      grn_gas += (grn_gas_asym - grn_gas_initial) * fact;

      for ( n = 0; n < n_limit; n++ )
      {
        np1 = n + 1;
        mode[n] = (- pow( - 1.0, np1) / double(np1)) * sqrt(8.0 * pow(grn_radius, 3) / pi) * grn_gas;
      }
    }

    // -- Output variables
    return grn_gas;
  }

  // The following function is provided for example's sake.
  // Another function can be used instead that returns the diffusion coefficient.

  double diffusion_coeff( double temperature )
  {
    // It is worth noting that the diffusion coefficient is defined here as effective, i.e., a pure diffusion problem is considered
    // (concomitant reactions such as trapping and re-solution at bubbles are incorporated in the diffusion coefficient).
    // See PolyPole-2 algorithm for the solution of the general diffusion-reaction problem.

    const double D = 5.0e-08 * exp( - 40262.0 / temperature );  // Lassmann and Benk, JNM 280, 2000

    return D;
  }

  // The following function is provided for example's sake.
  // Another function can be used instead that returns the gas production rate.

  double production_rate( double fissionrate )
  {
    const double                   yield = 0.30;
    const double F = fissionrate * yield;

    return F;
  }

}
