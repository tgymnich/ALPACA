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
#ifndef UNIT_HANDLER
#define UNIT_HANDLER

#include "input_output/inputfile_parser.h"
#include "log_writer.h"
#include <algorithm>
#include <cmath>

/**
 * @brief The UnitHandler class takes care of the conversions between unitless and non-unitless representations of quantities.
 *        Within the Simulation Kernel only non-dimensional quantities are used. In the user input and output,
 *        however, values are given in unit repesentation.
 */
class UnitHandler{

   const double length_reference_;
   const double velocity_reference_;
   const double density_reference_;
   const double temperature_reference_;
   const double mu_reference_;
   const double thermal_conductivity_reference_ ;
   const double diffusivity_reference_;
   const double surface_tension_coefficient_reference_;
   const double gravity_reference_;

public:
   UnitHandler() = delete;
   explicit UnitHandler(const InputFileParser& parser);
   ~UnitHandler() = default;
   UnitHandler( UnitHandler const& ) = delete;
   UnitHandler& operator=( UnitHandler const& ) = delete;
   UnitHandler( UnitHandler&& ) = delete;
   UnitHandler& operator=( UnitHandler&& ) = delete;

   /**
    * @brief Translates a unit based density value into a unitless one, using the provided reference density.
    * @param rho_dimensional Unit based density.
    * @return Unitless density value.
    */
   inline double NonDimensionalizeDensity(const double density_dimensional) const {return density_dimensional/density_reference_;}

   /**
    * @brief Translates a unit based length/size/position value into a unitless one, using the provided reference length.
    * @param length_dimensional Unit based length.
    * @return Unitless length value.
    */
   inline double NonDimensionalizeLength(const double length_dimensional) const {return length_dimensional/length_reference_;}

   /**
    * @brief Translates a unit based velocity value into a unitless one, using the provided reference velocity.
    * @param velocity_dimensional Unit based velocity.
    * @return Unitless velocity value.
    */
   inline double NonDimensionalizeVelocity(double velocity_dimensional) const {return velocity_dimensional/velocity_reference_;}

   /**
    * @brief Translates a unit based momentum value into a unitless one, using the provided reference density and velocity.
    * @param momentum_dimensional Unit based momentum.
    * @return Unitless momentum value.
    */
   inline double NonDimensionalizeMomentum(double momentum_dimensional) const {return momentum_dimensional/(density_reference_*velocity_reference_);}

   /**
    * @brief Translates a unit based energy value into a unitless one, using the provided reference density and velocity.
    * @param energy_dimensional Unit based energy.
    * @return Unitless energy value.
    */
   inline double NonDimensionalizeEnergy(double energy_dimensional) const {return energy_dimensional/(density_reference_*velocity_reference_*velocity_reference_);}

   /**
    * @brief Translates a unit based temperature value into a unitless one, using the provided reference temperature.
    * @param temperature_dimensional Unit based temperature.
    * @return Unitless temperature value.
    */
   inline double NonDimensionalizeTemperature(double temperature_dimensional) const {return temperature_dimensional/temperature_reference_;}

   /**
    * @brief Translates a unit based viscosity value into a unitless one, using the provided reference mu.
    * @param mu_dimensional Unit based viscosity.
    * @return Unitless viscosity value.
    */
   inline double NonDimensionalizeViscosity(double mu_dimensional) const {return mu_dimensional/mu_reference_;}

   /**
    * @brief Translates a unit based diffusivity value into a unitless one, using the provided reference diffusivity.
    * @param diffusivity_dimensional Unit based diffusivity.
    * @return Unitless diffusivity value.
    */
   inline double NonDimensionalizeDiffusivity(const double diffusivity_dimensional) const {return diffusivity_dimensional/diffusivity_reference_;}

   /**
    * @brief Translates a unit based surface-tension-coefficient value into a unitless one, using the provided reference surface tension coefficient.
    * @param surface_tension_coefficients_dimensional Unit based surface tension coefficients.
    * @return Unitless surface tension coefficients.
    */
   inline std::vector<double> NonDimensionalizeSurfaceTensionCoefficients(const std::vector<double> surface_tension_coefficients_dimensional) const {

      double one_surface_tension_coefficient_reference = 1.0/surface_tension_coefficient_reference_;
      std::vector<double> nondimensional_surface_tension_coefficients = surface_tension_coefficients_dimensional;

      //non-dimensionalization
      std::for_each( nondimensional_surface_tension_coefficients.begin(), nondimensional_surface_tension_coefficients.end(), [ &one_surface_tension_coefficient_reference ]( double& element ){ element*=one_surface_tension_coefficient_reference; } );

      return nondimensional_surface_tension_coefficients;
   }

   /**
    * @brief Translates a unit based time value into a unitless one, using the provided reference time.
    * @param time_dimensional Unit based time.
    * @return Unitless time value.
    */
   inline double NonDimensionalizeTime(const double time_dimensional) const {return time_dimensional*velocity_reference_/length_reference_;}

   /**
    * @brief Translates a unit based thermal conductivity value into a unitless one, using the provided reference thermal conductivity.
    * @param thermal_conductivity_dimensional Unit based thermal conductivity.
    * @return Unitless thermal conductivity value.
    */
   inline double NonDimensionalizeThermalConductivity(const double thermal_conductivity_dimensional) const {
      return thermal_conductivity_dimensional/thermal_conductivity_reference_;
   }

   /**
    * @brief Translates a unit based pressure value into a unitless one, using the provided reference parameters: velocity and density.
    * @param pressure_dimensional Unit pressure density.
    * @return Unitless pressure value.
    */
   inline double NonDimensionalizePressure(const double pressure_dimensional) const {
      return pressure_dimensional/(velocity_reference_*velocity_reference_*density_reference_);
   }

   /**
    * @brief Translates a unit based gravity value into a unitless one, using the provided reference gravity.
    * @param gravity Unit based gravity.
    * @return Unitless gravity value.
    */
   inline std::array<double,3> NonDimensionalizeGravity(const std::array<double,3> gravity) const {
      return {gravity[0]/gravity_reference_,gravity[1]/gravity_reference_,gravity[2]/gravity_reference_};
   }

   /**
    * @brief Nondimensionalzes the material properties of the fluids.
    * @param fluids Vector of fluids with unit based material properties.
    * @return Vector of fluids with unitless material properties
    */
   inline std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>> NonDimensionalizeMaterialProperties(std::vector<std::pair<MaterialName, std::unordered_map<std::string, double>>> fluids) const {
      for(auto& fluid : fluids) {
         auto& material_properties = fluid.second;

         //depending on the EOS, the parameters have different units, thus need different nondimensionilization
         if( material_properties.find("A") != material_properties.end() ) {
            if(fluid.first == MaterialName::WaterlikeFluid) {
               material_properties["A"] = NonDimensionalizePressure(material_properties["A"]); //A used for Tait EOS
            }
            else {
               material_properties["A"] = material_properties["A"]/(velocity_reference_*velocity_reference_); //A aka energy_translation_factor used by StiffenedGasComplete EOS
            }
         }
         if( material_properties.find("B") != material_properties.end() ) {
            material_properties["B"] = NonDimensionalizePressure(material_properties["B"]); //B background pressure used for both StiffenedGas and Tait
         }
         if( material_properties.find("C") != material_properties.end() ) {
            const double gamma = material_properties["gamma"];
            material_properties["C"] = material_properties["C"]*std::pow(density_reference_,gamma)/(velocity_reference_*velocity_reference_); //C aka thermal energy factor used by StiffenedGasComplete EOS
         }
         if( material_properties.find("rho0") != material_properties.end() ) {
            material_properties["rho0"] = NonDimensionalizeDensity(material_properties["rho0"]); //rho0  used by Tait EOS

         }
         if( material_properties.find("specificGasConstant") != material_properties.end() ) {
            material_properties["specificGasConstant"] = material_properties["specificGasConstant"]*temperature_reference_/(velocity_reference_*velocity_reference_); //specific gas constant
         }
         if( material_properties.find("thermalConductivity") != material_properties.end() ) {
            material_properties["thermalConductivity"] = NonDimensionalizeThermalConductivity(material_properties["thermalConductivity"]); //thermal conductivity
         }
         if( material_properties.find("dynamicShear") != material_properties.end() ) {
            material_properties["dynamicShear"] = NonDimensionalizeViscosity(material_properties["dynamicShear"]); //shear viscosity
         }
         if( material_properties.find("dynamicBulk") != material_properties.end() ) {
            material_properties["dynamicBulk"] = NonDimensionalizeViscosity(material_properties["dynamicBulk"]); //bulk viscosity
         }
      }
      return fluids;
   }

   /**
    * @brief Translates a unitless A (energy translation factor) value into a unit based one, unsing the provided reference parameters.
    * @param A_unitless Unitless A value.
    * @return Unit based A value.
    */
   inline double DimensionalizeA(const double A_unitless) const {return A_unitless*velocity_reference_*velocity_reference_;}

   /**
    * @brief Translates a unitless C (thermal energy factor) value into a unit based one, unsing the provided reference parameters.
    * @param C_unitless Unitless C value.
    * @return Unit based C value.
    */
   inline double DimensionalizeC(const double C_unitless, const double gamma) const {return C_unitless*velocity_reference_*velocity_reference_/std::pow(density_reference_,gamma);}


   /**
    * @brief Translates a unitless density value into a unit based one, unsing the provided reference density.
    * @param density_unitless Unitless density.
    * @return Unit based density value.
    */
   inline double DimensionalizeDensity(const double density_unitless) const {return density_unitless*density_reference_;}

   /**
    * @brief Translates a unitless momentum value into a unit based one, unsing the provided reference density and velocity.
    * @param momentum_unitless Unitless momentum.
    * @return Unit based momentum value.
    */
   inline double DimensionalizeMomentum(const double momentum_unitless) const {return momentum_unitless*density_reference_*velocity_reference_;}

   /**
    * @brief Translates a unitless velocity value into a unit based one, unsing the provided reference velocity.
    * @param velocity_unitless Unitless velocity.
    * @return Unit based velocity value.
    */
   inline double DimensionalizeVelocity(const double velocity_unitless) const {return velocity_unitless*velocity_reference_;}

   /**
    * @brief Translates a unitless length/size/position value into a unit based one, unsing the provided reference length.
    * @param length_unitless Unitless length/size/position.
    * @return Unit based length/size/position value.
    */
   inline double DimensionalizeLength(const double length_unitless) const {return length_unitless*length_reference_;}

   /**
    * @brief Translates a unitless time value into a unit based one, unsing the provided reference time.
    * @param time_unitless Unitless time.
    * @return Unit based time value.
    */
   inline double DimensionalizeTime(const double time_unitless) const {return time_unitless*length_reference_/velocity_reference_;}

   /**
    * @brief Translates a unitless energy value into a unit based one, unsing the provided reference energy.
    * @param energy_unitless Unitless energy.
    * @return Unit based energy value.
    */
   inline double DimensionalizeEnergy(const double energy_unitless) const {
      return energy_unitless*density_reference_*velocity_reference_*velocity_reference_;
   }

   /**
    * @brief Translates a unitless temperature value into a unit based one, unsing the provided reference parameters.
    * @param temperature_unitless Unitless temperature.
    * @return Unit based temperature value.
    */
   inline double DimensionalizeTemperature(const double temperature_unitless) const {
      return temperature_unitless*temperature_reference_;
   }

   /**
    * @brief Translates a unitless pressure value into a unit based one, unsing the provided reference parameters.
    * @param pressure_unitless Unitless pressure.
    * @return Unit based pressure value.
    */
   inline double DimensionalizePressure(const double pressure_unitless) const {
      return pressure_unitless*density_reference_*velocity_reference_*velocity_reference_;
   }

   /**
    * @brief Translates a unitless pressure value into a unit based one, unsing the provided reference parameters.
    * @param mu_unitless Unitless viscosity.
    * @return Unit based viscosity value.
    */
   inline double DimensionalizeViscosity(const double mu_unitless) const {
      return mu_unitless*mu_reference_;
   }

   /**
    * @brief Translates a unitless thermal conductivity value into a unit based one, unsing the provided reference parameters.
    * @param thermal_conductivity_unitless Unitless thermal conductivity.
    * @return Unit based thermal conductivity value.
    */
   inline double DimensionalizeThermalConductivity(const double thermal_conductivity_unitless) const {
      return thermal_conductivity_unitless*thermal_conductivity_reference_;
   }

   /**
    * @brief Translates a unitless specific gas constant value into a unit based one, unsing the provided reference parameters.
    * @param specific_gas_constant_unitless Unitless thermal conductivity.
    * @return Unit based thermal specific gas constant value.
    */
   inline double DimensionalizeSpecificGasConstant(const double specific_gas_constant_unitless) const {
      return specific_gas_constant_unitless*velocity_reference_*velocity_reference_/temperature_reference_;
   }

   /**
    * @brief Translates a unitless surface-tension-coefficient value into a unit based one, using the provided reference surface tension coefficient.
    * @param surface_tension_coefficients_unitless Unitless surface tension coefficients.
    * @return Unit based surface tension coefficients.
    */
   inline std::vector<double> DimensionalizeSurfaceTensionCoefficients(const std::vector<double> surface_tension_coefficients_unitless) const {

      std::vector<double> surface_tension_coefficients_dimensional = surface_tension_coefficients_unitless;

      //non-dimensionalization
      std::for_each( surface_tension_coefficients_dimensional.begin(), surface_tension_coefficients_dimensional.end(), [ this ]( double& element){ element*=surface_tension_coefficient_reference_; } );

      return surface_tension_coefficients_dimensional;
   }


   /**
    * @brief Gives the velocity reference
    * @return velocity reference
    */
   inline double GetVelocityReference() const {
      return velocity_reference_;
   }

   /**
    * @brief Gives the density reference
    * @return density reference
    */
   inline double GetDensityReference() const {
      return density_reference_;
   }

   /**
    * @brief Gives the length reference
    * @return length reference
    */
   inline double GetLengthReference() const {
      return length_reference_;
   }

   /**
    * @brief Gives the temperature reference
    * @return temperature reference
    */
   inline double GetTemperatureReference() const {
      return temperature_reference_;
   }

};

#endif // UNIT_HANDLER
