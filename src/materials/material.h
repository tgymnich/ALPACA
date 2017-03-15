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
#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

/**
 * @brief The Material class defines an interface for different materials. Materials may be equations of state or discretizations of
 *        them with e.g. fixed parameters to model a certain fluid.
 */
class Material {

   /**
   * @brief See public function GetPressure.
   */
   virtual double DoGetPressure   ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const = 0;
   /**
   * @brief See public function GetEnthalpy.
   */
   virtual double DoGetEnthalpy   ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const = 0;
   /**
   * @brief See public function GetEnergy.
   */
   virtual double DoGetEnergy     ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const = 0;
   /**
   * @brief See public function GetTemperature.
   */
   virtual double DoGetTemperature( double const, double const, double const, double const, double const ) const {
      return -1.0;
   }
   /**
   * @brief See public function GetGruneisen.
   */
   virtual double DoGetGruneisen() const = 0;
   /**
   * @brief See public function GetGruneisen.
   */
   virtual double DoGetGruneisen( double const ) const {
      return DoGetGruneisen();
   }
   /**
   * @brief See public function GetPsi.
   */
   virtual double DoGetPsi( double const pressure, double const one_density ) const = 0;
   /**
   * @brief See public function GetGamma.
   */
   virtual double DoGetGamma() const {
      return -1.0;
   }
   /**
   * @brief See public function GetB.
   */
   virtual double DoGetB() const {
      return -1.0;
   }
   /**
   * @brief See public function GetSpeedOfSound.
   */
   virtual double DoGetSpeedOfSound( double const density, double const pressure ) const = 0;
   /**
   * @brief See public function GetViscosity.
   */
   virtual std::vector<double> DoGetViscosity() const = 0;
   /**
   * @brief See public function GetThermalConductivity.
   */
   virtual double DoGetThermalConductivity() const {
      return 0.0;
   }
   /**
   * @brief See public function GetSpecificHeat.
   */
   virtual double DoGetSpecificHeat() const {
      return -1.0;
   }

public:
   explicit Material() = default;
   virtual ~Material() = default;
   Material( Material const& ) = delete;
   Material& operator=( Material const& ) = delete;
   Material( Material&& ) = delete;
   Material& operator=( Material&& ) = delete;

   /**
    * @brief Computes the pressure based on the given input of arbitrary density, momenta and energy according to the applied material equation of state.
    * @param density .
    * @param momentum_x .
    * @param momentum_y .
    * @param momentum_z .
    * @param energy .
    * @return Pressure for the state imposed by the inputs of the implemented material.
    */
   double GetPressure( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
      return DoGetPressure(density, momentum_x, momentum_y, momentum_z, energy);
   }

   /**
    * @brief Computes the enthalpy based on the given input of arbitrary density, momenta and energy according to the applied material equation of state.
    * @param density .
    * @param momentum_x .
    * @param momentum_y .
    * @param momentum_z .
    * @param energy .
    * @return Enthalphy for the state imposed by the inputs of the implemented material.
    */
   double GetEnthalpy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
      return DoGetEnthalpy(density, momentum_x, momentum_y, momentum_z, energy);
   }

   /**
    * @brief Computes the temperature based on the given input of arbitrary density, momenta and energy according to the applied material equation of state.
    * @param density .
    * @param momentum_x .
    * @param momentum_y .
    * @param momentum_z .
    * @param energy .
    * @return Temperature for the state imposed by the inputs of the implemented material.
    * @note Returns -1.0 if the equation of state cannot supply a temperature calculation.
    */
   double GetTemperature( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
      return DoGetTemperature(density, momentum_x, momentum_y, momentum_z, energy);
   }

   /**
    * @brief Computes the energy based on the given input of arbitrary density, momenta and pressure according to the applied material equation of state.
    * @param density .
    * @param momentum_x .
    * @param momentum_y .
    * @param momentum_z .
    * @param pressure .
    * @return Energy for the state imposed by the inputs of the implemented material.
    */
   double GetEnergy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const {
      return DoGetEnergy(density, momentum_x, momentum_y, momentum_z, pressure);
   }

   /**
    * @brief Computes the Grueneisen coefficient according to the applied material equation of state in generalized form indpendetend of any conservative state.
    * @return Grueneisen coefficient for the implemented material.
    */
   double GetGruneisen() const {
      return DoGetGruneisen();
   }

   /**
    * @brief Computes the Grueneisen coefficient based on the given input of an arbitrary density according to the applied material equation of state in generalized form.
    * @param density .
    * @return Grueneisen coefficient for the implemented material.
    */
   double GetGruneisen( double const density ) const {
      return DoGetGruneisen(density);
   }

   /**
    * @brief Computes psi based on the given input of arbitrary pressure and density according to the applied material equation of state in generalized form.
    * @param pressure .
    * @param one over density (saves costs) .
    * @return Psi for the state imposed by the inputs of the implemented material.
    */
   double GetPsi( double const pressure, double const one_density ) const {
      return DoGetPsi(pressure, one_density);
   }

   /**
    * @brief Returns the isentropic exponent gamma of the material.
    * @return Gamma value of the implemented material.
    */
   double GetGamma() const {
      return DoGetGamma();
   }

   /**
    * @brief Returns the B of the material.
    * @return B value of the implemented material.
    */
   double GetB() const {
      return DoGetB();
   }

   /**
    * @brief Computes the speed of sound based on the given input of arbitrary density and pressure according to the applied material equation of state
    * @param density .
    * @param pressure .
    * @return Speed of sound for the state imposed by the inputs of the implemented material.
    */
   double GetSpeedOfSound( double const density, double const pressure ) const {
      return DoGetSpeedOfSound(density, pressure);
   }

   /**
    * @brief Provides the viscosity parameters of the material.
    * @return Viscosity Parameters %Currently mu_shear and mu_bulk, might change in future Versions%.
    */
   std::vector<double> GetViscosity() const {
      return DoGetViscosity();
   }

   /**
    * @brief Returns the thermal conductivity of the material.
    * @return Thermal conductivity.
    */
   double GetThermalConductivity() const {
      return DoGetThermalConductivity();
   }

   /**
    * @brief Returns the thermal conductivity of the material.
    * @return Thermal conductivity.
    */
   double GetSpecificHeat() const {
      return DoGetSpecificHeat();
   }
};

#endif //MATERIAL_H
