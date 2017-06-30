// Copyright 2017 Battelle Energy Alliance, LLC
// Licensed under the Mozilla Public License version 2.0

// @author: Davide Pizzocri (davide.pizzocri@gmail.com)
// @author: Giovanni Pastore (giovanni.pastore@inl.gov)
// @author: Cristian Rabiti (cristian.rabiti@inl.gov)

#ifndef _POLYPOLE_H
#define _POLYPOLE_H

namespace polypole
{
  // Numerical algorithm for the solution of the equilibrium diffusion equation
  // Pizzocri et al., JNM 478 (2016) 333-342
  //
  // The algorithm requires auxiliary functions in order to sample
  // diffusion and gas production rates along the time step
  // These functions are provided in this namespace,
  // intending that they are extendable or replaceable

  double PolyPole1( double dt_sec,
                    double grn_radius,
                    double temp0,
                    double temp1,
                    double frate0,
                    double frate1,
                    double grn_gas_initial );

  // Functions calculating the effective diffusion coefficient and the gas production rate
  // May not be necessary as your code may already include functions that can be coupled to
  // PolyPole to provide the algorithm with these input quantities.

  double diffusion_coeff( double temperature );
  double production_rate( double fissionrate );
}

#endif
