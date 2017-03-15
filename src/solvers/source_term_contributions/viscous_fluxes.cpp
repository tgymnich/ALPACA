/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "viscous_fluxes.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "index_transformations.h"
#include "stencils/differentiation_utilities.h"

using ViscousFluxesDerivativeStencil = DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil>::type;
using ViscousFluxesReconstructionStencil = ReconstructionStencilSetup::Concretize<viscous_fluxes_reconstruction_stencil>::type;

/**
 * @brief The constructor for the ViscousFluxes class.
 * @param material_manager The material manager which contains information about the viscosities.
 */
ViscousFluxes::ViscousFluxes(MaterialManager const& material_manager) :
   material_manager_(material_manager)
{
   // Empty constructor besides initializer list
}

/**
 * @brief Computes cell face fluxes due to viscosity.
 * @param mat_block A pair containing a block and its material.
 * @param dissipative_flux_x, dissipative_flux_y, dissipative_flux_z Reference to the face fluxes (indirect return parameter).
 * @param cell_size The cell size.
 */
void ViscousFluxes::ComputeFluxes( const std::pair<const MaterialName, Block> &mat_block,
   double (&dissipative_flux_x)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
   double (&dissipative_flux_y)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
   double (&dissipative_flux_z)[FF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1], double const cell_size ) const {

   // y and z velocity buffers may not be available. As workaround use the x velocity buffer in these cases.
   //Tthis is legal since the respective gradients are not computed/used anyway
   double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityX);
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::DIM() != Dimension::One   ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityY) : u;
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::DIM() == Dimension::Three ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityZ) : u;

   /**
    * Description for the positions of the Array:
    * [CC::TCX()]    [CC::TCY()]    [CC::TCZ()]    [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Velocity gradient: du_i / dx_j
    */
   double velocity_gradient[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())][DTI(CC::DIM())];
   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  velocity_gradient[i][j][k][r][s] = 0.0;
               }
            }
         } //k
      } //j
   } //i

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Cell face x/y/z   Velocity gradient: du_i / dx_j
    */
   double velocity_gradient_at_cell_faces[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())];

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Cell face x/y/z   Velocity in x/y/z direction
    */
   double velocity_at_cell_faces[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())];

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  tau_ij
    */
   double tau[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())];

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  for( unsigned int t = 0; t < DTI(CC::DIM()); ++t ) {
                     velocity_gradient_at_cell_faces[i][j][k][r][s][t] = 0.0;
                  }
                  tau[i][j][k][r][s] = 0.0;
                  velocity_at_cell_faces[i][j][k][r][s] = 0.0;
               }
            }
         }
      }
   }

   ComputeVelocityGradient(u, v, w, cell_size, velocity_gradient);
   ComputeVelocityGradientAtCellFaces(velocity_gradient, cell_size, u, v, w, velocity_gradient_at_cell_faces);
   ComputeVelocityAtCellFaces(u, v, w, cell_size, velocity_at_cell_faces);
   ComputeTau(velocity_gradient_at_cell_faces, material_manager_.GetViscosity(mat_block.first), tau);

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            double energy_flux_x = 0.0;
            double energy_flux_y = 0.0;
            double energy_flux_z = 0.0;

            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               dissipative_flux_x[ETI(FF::AME()[r])][i][j][k] -= tau[i][j][k][0][r];
               if constexpr( CC::DIM() != Dimension::One )   dissipative_flux_y[ETI(FF::AME()[r])][i][j][k] -= tau[i][j][k][1][r];
               if constexpr( CC::DIM() == Dimension::Three ) dissipative_flux_z[ETI(FF::AME()[r])][i][j][k] -= tau[i][j][k][2][r];

               energy_flux_x += tau[i][j][k][0][r] * velocity_at_cell_faces[i][j][k][0][r];
               if constexpr( CC::DIM() != Dimension::One )   energy_flux_y += tau[i][j][k][1][r] * velocity_at_cell_faces[i][j][k][1][r];
               if constexpr( CC::DIM() == Dimension::Three ) energy_flux_z += tau[i][j][k][2][r] * velocity_at_cell_faces[i][j][k][2][r];
            }

            dissipative_flux_x[ETI(Equation::Energy)][i][j][k] -= energy_flux_x;
            if constexpr( CC::DIM() != Dimension::One )   dissipative_flux_y[ETI(Equation::Energy)][i][j][k] -= energy_flux_y;
            if constexpr( CC::DIM() == Dimension::Three ) dissipative_flux_z[ETI(Equation::Energy)][i][j][k] -= energy_flux_z;

         }
      }
   }

}

/**
 * @brief Computes the gradient of the velocity field at cell centers.
 * @param u The velocity field in x direction.
 * @param v The velocity field in y direction.
 * @param w The velocity field in z direction.
 * @param cell_size The cell size.
 * @param velocity_gradient The velocity gradient field (indirect return parameter).
 */
void ViscousFluxes::ComputeVelocityGradient(double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()], double const cell_size,
   double (&velocity_gradient)[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())][DTI(CC::DIM())]) const {
   /**
    * @brief Offsets in order to also calculate first derivatives in halo cells. This is necessary for the  reconstruction to cell faces.
    */
   constexpr unsigned int offset_x = ViscousFluxesDerivativeStencil::DownstreamStencilSize();
   constexpr unsigned int offset_y = CC::DIM() != Dimension::One ? ViscousFluxesDerivativeStencil::DownstreamStencilSize() : 0;
   constexpr unsigned int offset_z = CC::DIM() == Dimension::Three ? ViscousFluxesDerivativeStencil::DownstreamStencilSize() : 0;

   for( unsigned int i = 0 + offset_x; i < CC::TCX() - offset_x; ++i ) {
      for( unsigned int j = 0 + offset_y; j < CC::TCY() - offset_y; ++j ) {
         for( unsigned int k = 0 + offset_z; k < CC::TCZ() - offset_z; ++k ) {
            std::array< std::array<double, 3>, 3> const single_gradient = DifferentiationUtilities::ComputeVectorGradient<ViscousFluxesDerivativeStencil, double>( u, v, w, i, j, k, cell_size );
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  velocity_gradient[i][j][k][r][s] = single_gradient[r][s];
               }
            }
         }
      }
   }
}

/**
 * @brief Reconstructs the velocity gradient at cell faces.
 * @param velocity_gradient The velocity gradient at cell centers.
 * @param cell_size The cell size.
 * @param u The velocity field in x direction.
 * @param v The velocity field in y direction.
 * @param w The velocity field in z direction.
 * @param velocity_gradient_at_cell_faces The velocity gradient at cell faces (indirect return parameter).
 */
void ViscousFluxes::ComputeVelocityGradientAtCellFaces( double const (&velocity_gradient)[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())][DTI(CC::DIM())], double const cell_size,
   double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double (&velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())] ) const {

   std::vector<double> interpolation_array;
   interpolation_array.reserve(ViscousFluxesReconstructionStencil::StencilSize());

   for( unsigned int i = CC::FICX() - 1; i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {

                  if( s != 0 ) {
                     interpolation_array.clear();
                     for( unsigned int ii = 0; ii < ViscousFluxesReconstructionStencil::StencilSize(); ++ii ) {
                        interpolation_array.push_back(velocity_gradient[i + (ii - ViscousFluxesReconstructionStencil::DownstreamStencilSize())][j][k][r][s]);
                     }

                     velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][r][s] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, double>(interpolation_array, cell_size);
                  } else {
                     switch( r ) {
                        case 0:
                           velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::X, double>(u, i, j, k, cell_size);
                           break;
                        case 1:
                           velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::X, double>(v, i, j, k, cell_size);
                           break;
                        case 2:
                           velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::X, double>(w, i, j, k, cell_size);
                           break;
                        default:
                           break;
                     }
                  }

               }
            }
         }
      }
   }

   if constexpr( CC::DIM() != Dimension::One ) {
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY() - 1; j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

               for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
                  for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {

                     if( s != 1 ) {
                        interpolation_array.clear();
                        for( unsigned int jj = 0; jj < ViscousFluxesReconstructionStencil::StencilSize(); ++jj ) {
                           interpolation_array.push_back(velocity_gradient[i][j + (jj - ViscousFluxesReconstructionStencil::DownstreamStencilSize())][k][r][s]);
                        }

                        velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][r][s] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, double>(interpolation_array, cell_size);
                     } else {
                        switch( r ) {
                           case 0:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Y, double>(u, i, j, k, cell_size);
                              break;
                           case 1:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Y, double>(v, i, j, k, cell_size);
                              break;
                           case 2:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Y, double>(w, i, j, k, cell_size);
                              break;
                           default:
                              break;
                        }
                     }

                  }
               }
            }
         }
      }
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ() - 1; k <= CC::LICZ(); ++k ) {

               for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
                  for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {

                     if( s != 2 ) {
                        interpolation_array.clear();
                        for( unsigned int kk = 0; kk < ViscousFluxesReconstructionStencil::StencilSize(); ++kk ) {
                           interpolation_array.push_back(velocity_gradient[i][j][k + (kk - ViscousFluxesReconstructionStencil::DownstreamStencilSize())][r][s]);
                        }

                        velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][r][s] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, double>(interpolation_array, cell_size);
                     } else {
                        switch( r ) {
                           case 0:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Z, double>(u, i, j, k, cell_size);
                              break;
                           case 1:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Z, double>(v, i, j, k, cell_size);
                              break;
                           case 2:
                              velocity_gradient_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][r][s] = DifferentiationUtilities::ApplyStencil<FourthOrderCellFace, StencilProperty::Central, Direction::Z, double>(w, i, j, k, cell_size);
                              break;
                           default:
                              break;
                        }
                     }

                  }
               }
            }
         }
      }
   }
}

/**
 * @brief Computes tau, the viscous part of the stress tensor.
 * @param velocity_gradient_at_cell_faces The velocity gradient field at cell faces.
 * @param viscosity A vector containing the shear and bulk viscosity.
 * @param tau The viscous part of the stress tensor.
 */
void ViscousFluxes::ComputeTau( double const (&velocity_gradient_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())],
   std::vector<double> const viscosity, double (&tau)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const {

   double const mu_1 = viscosity[0];
   double const mu_2 = viscosity[1] - 2.0 * viscosity[0] / 3.0;

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               double volumetric_part = 0.0;
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  tau[i][j][k][r][s] = mu_1 * (velocity_gradient_at_cell_faces[i][j][k][r][r][s] + velocity_gradient_at_cell_faces[i][j][k][r][s][r]);
                  volumetric_part += velocity_gradient_at_cell_faces[i][j][k][r][s][s];
               }
               tau[i][j][k][r][r] += volumetric_part * mu_2;
            }
         }
      }
   }
}

/**
 * @brief Reconstruct the velocity at cell faces.
 * @param u The velocity in x direction.
 * @param v The velocity in y direction.
 * @param w The velocity in z direction.
 * @param cell_size The cell size.
 * @param velocity_at_cell_faces The velocity field at cell faces (indirect return parameter).
 */
void ViscousFluxes::ComputeVelocityAtCellFaces( double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
   double const cell_size, double (&velocity_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const {

   for( unsigned int i = CC::FICX()-1; i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][0] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::X, double>(u, i, j, k, cell_size);
            if constexpr( CC::DIM() != Dimension::One )   velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][1] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::X, double>(v, i, j, k, cell_size);
            if constexpr( CC::DIM() == Dimension::Three ) velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][0][2] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::X, double>(w, i, j, k, cell_size);

         }
      }
   }

   if constexpr( CC::DIM() != Dimension::One ) {
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY()-1; j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][0] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Y, double>(u, i, j, k, cell_size);
               if constexpr( CC::DIM() != Dimension::One )   velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][1] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Y, double>(v, i, j, k, cell_size);
               if constexpr( CC::DIM() == Dimension::Three ) velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][1][2] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Y, double>(w, i, j, k, cell_size);
            }
         }
      }
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ()-1; k <= CC::LICZ(); ++k ) {
               velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][0] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Z, double>(u, i, j, k, cell_size);
               if constexpr( CC::DIM() != Dimension::One )   velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][1] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Z, double>(v, i, j, k, cell_size);
               if constexpr( CC::DIM() == Dimension::Three ) velocity_at_cell_faces[BIT::T2FX(i)][BIT::T2FY(j)][BIT::T2FZ(k)][2][2] = DifferentiationUtilities::ApplyStencil<ViscousFluxesReconstructionStencil, StencilProperty::Central, Direction::Z, double>(w, i, j, k, cell_size);
            }
         }
      }
   }
}
