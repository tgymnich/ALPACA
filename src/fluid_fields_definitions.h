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
#ifndef FLUID_FIELDS_DEFINITIONS_H
#define FLUID_FIELDS_DEFINITIONS_H

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>
#include "user_specifications/compile_time_constants.h"
#include "helper_functions.h"
#include "enums/unit_type.h"

/**
 * @brief Unique identifier for the conservative buffer type, i.e. the average, right-hand side, or initial buffer.
 */
enum class ConservativeBufferType : unsigned short {Average = 0, RightHandSide = 1, Initial = 2};

/**
 * @brief Converts an conservative buffer identifier to a (C++11 standard compliant, i. e. positive) array index. "CTTI = Conservative type to index"
 * @param c The conservative buffer identifier
 * @return The index.
 */
constexpr std::underlying_type<ConservativeBufferType>::type CTTI( const ConservativeBufferType c ) {return static_cast<typename std::underlying_type<ConservativeBufferType>::type>( c );}

/**
 * @brief Unique Identifier for the Field type, i.e. a conservative or a prime-state field.
 */
enum class FluidFieldType : unsigned short {Conservatives = 0, PrimeStates = 1};

/**
 * @brief Converts an fluid field identifier to a (C++11 standard compliant, i. e. positive) array index. "FTTI = Field type to Index"
 * @param f The fluid field identifier
 * @return The index.
 */
constexpr std::underlying_type<FluidFieldType>::type FTTI( const FluidFieldType f ) {return static_cast<typename std::underlying_type<FluidFieldType>::type>( f );}


namespace FieldDetails {

   // NH: NEVER EVER change the underlying type.
   /**
   * @brief Unique identifier for all possible conservative equations in arbitrary order.
   * @note  Every member has to be added to the Equation enumeration as well.
   */
   enum class EquationPool : unsigned int {
      // mandatory equations
      Mass, Energy, MomentumX, MomentumY, MomentumZ,
      // optional equations
      // Example
   };

   // NH: NEVER EVER change the underlying type.
   /**
   * @brief Unique Identifier for all possible prime states in arbitrary order.
   * @note  Every member has to be added to the PrimeState enumeration as well.
   */
   enum class PrimeStatePool : unsigned int {
      // mandatory prime states
      Density, Pressure, VelocityX, VelocityY, VelocityZ,
      // optional prime states
      Temperature, // Example
   };

   // NH: NEVER EVER change the underlying type.
   /**
   * @brief Unique Identifier for all possible interface quantities in arbitrary order.
   * @note  Every member has to be added to the InterfaceQuantity enumeration as well.
   */
   enum class InterfaceQuantityPool : unsigned short {
      // mandatory interface quantities
      Velocity,
      // optional interface quantities
      PressurePositive, PressureNegative
   };

   /**
   * @brief Bundles relevant configuration data for active fluid fields (i.e. conservative equations, prime states, interface quantities).
   * @tparam DerivedActiveFieldsDefinition CRTP template parameter
   * @tparam FieldEnum Enumeration type of the underlying fluid field
   */
   template<typename DerivedActiveFieldsDefinition, typename FieldEnum>
   struct ActiveFieldsDefinition {
      struct FieldInfo {
         FieldEnum const Field;
         UnitType const Unit;
         std::string_view const OutputName;
         std::string_view const InputName = "";
      };

      /**
      * @brief Value used to offset inactive from active fluid fields in the later generated enumeration.
      * @note IMPORTANT: Has to be larger than the number of active fields.
      */
      static constexpr unsigned int InactiveFieldOffset = 100;

      /**
      * @brief Gives the number of active fields.
      */
      static constexpr unsigned int Count = DerivedActiveFieldsDefinition::Definition.size();

      /**
      * @brief Gives the index of the given field. Should only be used to construct the resulting enumeration.
      * @param field Field whose index should be found.
      * @return Index of the field.
      */
      static constexpr unsigned int FieldIndex( FieldEnum const field ) {
         for( unsigned int i = 0; i < Count; ++i ) {
            if( DerivedActiveFieldsDefinition::Definition[i].Field == field ) return i;
         }
         return static_cast<unsigned int>( field ) + InactiveFieldOffset;
      }

      /**
      * @brief Gives the unit of the field at the given index.
      * @param index Index of the field.
      * @return UnitType representing the unit of the field.
      */
      static constexpr auto GetUnit( unsigned int const index ) {
         return DerivedActiveFieldsDefinition::Definition[index].Unit;
      }

      /**
      * @brief Gives the name of the field as it should be used in ouput writing.
      * @param index Index of the field.
      * @return Name used for output. Empty if the field should not be used in the output.
      */
      static constexpr auto GetOutputName( unsigned int const index ) {
         return DerivedActiveFieldsDefinition::Definition[index].OutputName;
      }

      /**
      * @brief Gives the name of the field as it appears in input files.
      * @param index Index of the field.
      * @return Name used in input files. Empty if the field is not used in input files.
      * @note Only prime states are currently used in input files.
      */
      static constexpr auto GetInputName( unsigned int const index ) {
         return DerivedActiveFieldsDefinition::Definition[index].InputName;
      }
   };

   /**
   * @brief Configuration of the active conservative equations.
   */
   struct ActiveEquations : public ActiveFieldsDefinition<ActiveEquations, EquationPool> {
      static constexpr auto Definition = MakeArray(
         //           conservative                  unit                    output name*
         // mandatory equations (DO NOT CHANGE):
           FieldInfo{ EquationPool::Mass,           UnitType::Density,      "mass" }
         , FieldInfo{ EquationPool::Energy,         UnitType::Energy,       "energy" }
         , FieldInfo{ EquationPool::MomentumX,      UnitType::Momentum,     "momentumX" }
   #if DIMENSION != 1
         , FieldInfo{ EquationPool::MomentumY,      UnitType::Momentum,     "momentumY" }
   #endif
   #if DIMENSION == 3
         , FieldInfo{ EquationPool::MomentumZ,      UnitType::Momentum,     "momentumZ" }
   #endif
         // optional equations:
      // , FieldInfo{ EquationPool::Example,        UnitType::Example,      "example" }
      );
      // notes: * if output name is empty the field will not be considered for output; conservatives are only considered in debug output
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active equations! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active prime states.
   */
   struct ActivePrimeStates : public ActiveFieldsDefinition<ActivePrimeStates, PrimeStatePool> {
      static constexpr auto Definition = MakeArray(
         //           prime state                      unit                       output name*      input name**
         // mandatory prime states (DO NOT CHANGE):
           FieldInfo{ PrimeStatePool::Density,         UnitType::Density,         "density",        "DENSITY" }
         , FieldInfo{ PrimeStatePool::Pressure,        UnitType::Pressure,        "pressure",       "PRESSURE" }
         , FieldInfo{ PrimeStatePool::VelocityX,       UnitType::Velocity,        "velocityX",      "VELOCITY_X" }
   #if DIMENSION != 1
         , FieldInfo{ PrimeStatePool::VelocityY,       UnitType::Velocity,        "velocityY",      "VELOCITY_Y" }
   #endif
   #if DIMENSION == 3
         , FieldInfo{ PrimeStatePool::VelocityZ,       UnitType::Velocity,        "velocityZ",      "VELOCITY_Z" }
   #endif
         // optional prime states:
      // , FieldInfo{ PrimeStatePool::Temperature,     UnitType::Temperature,     "temperature",    "" }
      // , FieldInfo{ PrimeStatePool::Example,         UnitType::Example          "example",        "EXAMPLE" }
      );
      // notes: *  if output name is empty the field will not be considered for output
      //        ** if input name is empty the field will not be read from input files
      // attention: the OUTPUT name is also used as input name for fixed value boundary conditions (TODO-19 unify them)
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active prime states! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active interface quantities.
   */
   struct ActiveInterfaceQuantities : public ActiveFieldsDefinition<ActiveInterfaceQuantities, InterfaceQuantityPool> {
      static constexpr auto Definition = MakeArray(
         //           interface quantity                          unit                  output name*
         // mandatory interface quantities (DO NOT CHANGE):
           FieldInfo{ InterfaceQuantityPool::Velocity,            UnitType::Velocity,   "interface_velocity" }
         // optional interface quantities:
         , FieldInfo{ InterfaceQuantityPool::PressurePositive,    UnitType::Pressure,   "" }
         , FieldInfo{ InterfaceQuantityPool::PressureNegative,    UnitType::Pressure,   "" }
      // , FieldInfo{ InterfaceQuantityPool::Example,             UnitType::Example,    "example" }
      );
      // notes: * if output name is empty the field will not be considered for output
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active interface quantities! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };
}

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the conservative equation. All active equations are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of EquationPool has to be added here as well.
 */
enum class Equation : unsigned int {
   // mandatory equations:
   Mass           = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Mass ),
   Energy         = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Energy ),
   MomentumX      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumX ),
   MomentumY      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumY ),
   MomentumZ      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumZ ),
   // optional equations
   // Example     = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Example ),
};
/**
 * @brief Converts an equation identifier to a (C++11 standard compliant, i. e. positive) array index. "ETI = Equation to Index"
 * @param e The equation identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Equation>::type ETI( const Equation e ) {return static_cast<typename std::underlying_type<Equation>::type>( e );}

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the prime states. All active prime states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of PrimeStatePool has to be added here as well.
 */
enum class PrimeState : unsigned int {
   // mandatory prime states
   Density        = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Density ),
   Pressure       = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Pressure ),
   VelocityX      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityX ),
   VelocityY      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityY ),
   VelocityZ      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityZ ),
   // optional prime states
   Temperature    = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Temperature ),
   // Example     = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Example ),
};
/**
 * @brief Converts a prime state identifier to a (C++11 standard compliant, i. e. positive) array index. "PTI = Prime state to Index"
 * @param p The prime state identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<PrimeState>::type PTI( const PrimeState p ) {return static_cast<typename std::underlying_type<PrimeState>::type>( p );}

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the interface quantities. All active interface quantities are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceQuantityPool has to be added here as well.
 */
enum class InterfaceQuantity : unsigned int {
   // mandatory interface quantities
   Velocity             = FieldDetails::ActiveInterfaceQuantities::FieldIndex( FieldDetails::InterfaceQuantityPool::Velocity ),
   // optional interface quantities
   PressurePositive     = FieldDetails::ActiveInterfaceQuantities::FieldIndex( FieldDetails::InterfaceQuantityPool::PressurePositive ),
   PressureNegative     = FieldDetails::ActiveInterfaceQuantities::FieldIndex( FieldDetails::InterfaceQuantityPool::PressureNegative ),
   // Example           = FieldDetails::ActiveInterfaceQuantities::FieldIndex( FieldDetails::InterfaceQuantityPool::Example ),
};
/**
 * @brief Converts a interface quantity identifier to a (C++11 standard compliant, i. e. positive) array index. "IQTI = Interface quantity to index"
 * @param iq The interface quantity identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceQuantity>::type IQTI( const InterfaceQuantity iq ) {return static_cast<typename std::underlying_type<InterfaceQuantity>::type>( iq );}


class FluidFieldsDefinitions {

   // get arrays of the consecutively ordered active fields
   static constexpr auto active_equations_            = IndexSequenceToEnumArray<Equation>( std::make_index_sequence<FieldDetails::ActiveEquations::Count>{} );
   static constexpr auto active_prime_states_         = IndexSequenceToEnumArray<PrimeState>( std::make_index_sequence<FieldDetails::ActivePrimeStates::Count>{} );
   static constexpr auto active_interface_quantities_ = IndexSequenceToEnumArray<InterfaceQuantity>( std::make_index_sequence<FieldDetails::ActiveInterfaceQuantities::Count>{} );

   static constexpr const std::array<Equation,DTI( CC::DIM())> active_momenta_ = {Equation::MomentumX
#if DIMENSION!=1
      ,Equation::MomentumY
#endif
#if DIMENSION==3
      ,Equation::MomentumZ
#endif
   };

   static constexpr const std::array<PrimeState,DTI( CC::DIM())> active_velocities_ = {PrimeState::VelocityX
#if DIMENSION!=1
      ,PrimeState::VelocityY
#endif
#if DIMENSION==3
      ,PrimeState::VelocityZ
#endif
   };
   static constexpr std::array<Equation, 2> wavelet_analysis_equations_ {Equation::Mass, Equation::Energy};

   static constexpr std::array<InterfaceQuantity, 1> interface_quantities_to_extend_= {InterfaceQuantity::Velocity};

public:
   FluidFieldsDefinitions() = delete;
   ~FluidFieldsDefinitions() = default;
   FluidFieldsDefinitions( FluidFieldsDefinitions const& ) = delete;
   FluidFieldsDefinitions& operator=( FluidFieldsDefinitions const& ) = delete;
   FluidFieldsDefinitions( FluidFieldsDefinitions&& ) = delete;
   FluidFieldsDefinitions& operator=( FluidFieldsDefinitions&& ) = delete;

   /**
    * @brief Gives whether the given conservative equation is active.
    * @param e Equation identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsEquationActive( Equation const e ) {
      return static_cast<unsigned int>( e ) < FieldDetails::ActiveEquations::InactiveFieldOffset;
   }

   /**
    * @brief Gives whether the given prime state is active.
    * @param p PrimeState identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsPrimeStateActive( PrimeState const p ) {
      return static_cast<unsigned int>( p ) < FieldDetails::ActivePrimeStates::InactiveFieldOffset;
   }

   /**
    * @brief Gives whether the given interface quantity is active.
    * @param iq InterfaceQuantity identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsInterfaceQuantityActive( InterfaceQuantity const iq ) {
      return static_cast<unsigned int>( iq ) < FieldDetails::ActiveInterfaceQuantities::InactiveFieldOffset;
   }

   /**
    * @brief Gives the name used in inputfiles for the given field.
    */
   static constexpr auto FieldInputName( PrimeState const p ) {
      return FieldDetails::ActivePrimeStates::GetInputName( PTI(p) );
   }

   /**
    * @brief Gives the name used in the output for the given field.
    */
   static constexpr auto FieldOutputName( Equation const e ) {
      return FieldDetails::ActiveEquations::GetOutputName( ETI(e) );
   }

   /**
    * @brief Gives the name used in the output for the given field.
    */
   static constexpr auto FieldOutputName( PrimeState const p ) {
      return FieldDetails::ActivePrimeStates::GetOutputName( PTI(p) );
   }

   /**
    * @brief Gives the name used in the output for the given field.
    */
   static constexpr auto FieldOutputName( InterfaceQuantity const iq ) {
      return FieldDetails::ActiveInterfaceQuantities::GetOutputName( IQTI(iq) );
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    */
   static constexpr auto FieldUnit( Equation const e ) {
      return FieldDetails::ActiveEquations::GetUnit( ETI(e) );
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    */
   static constexpr auto FieldUnit( PrimeState const p ) {
      return FieldDetails::ActivePrimeStates::GetUnit( PTI(p) );
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    */
   static constexpr auto FieldUnit( InterfaceQuantity const iq ) {
      return FieldDetails::ActiveInterfaceQuantities::GetUnit( IQTI(iq) );
   }

   /**
    * @brief Gives the number of equations considered in the simulation, i.e. Euler Equations.
    * @return Number of Equations = 5 (mass, energy, X-,Y-,Z-momentum) for 3D
    * @return Number of Equations = 4 (mass, energy, X-,Y-momentum) for 2D
    * @return Number of Equations = 3 (mass, energy, X-momentum) for 1D
    */
   static constexpr unsigned int ANOE() {return FieldDetails::ActiveEquations::Count;}

   /**
    * @brief Gives the number of prime states considered in the simulation.
    * @return Number of prime states = 6 (density, pressure, temperature, X-, Y-, Z-velocity) for 3D
    * @return Number of prime states = 6 (density, pressure, temperature, X-, Y-velocity) for 2D
    * @return Number of prime states = 4 (density, pressure, temperature, X-velocity) for 1D
    */
   static constexpr unsigned int ANOP() {return FieldDetails::ActivePrimeStates::Count;}

   /**
    * @brief Gives the number of interface quantities considered in the simulation.
    * @return Number of active interface quantities
    */
   static constexpr unsigned int ANIQ() {return FieldDetails::ActiveInterfaceQuantities::Count;}

   /**
    * @brief Gives the number of active fields for the given field type, i.e. conservatives or prime states.
    * @param field_type The fluid-field type.
    * @return Number of active fields.
    */
   static constexpr unsigned int ANOF( const FluidFieldType field_type ) {
      switch( field_type ) {
         case FluidFieldType::Conservatives:
            return ANOE();
         default: // FluidFieldType::PrimeStates
            return ANOP();
      }
   }

   /**
    * @brief Gives the set of equations which are worked on in this simulation. I. e. varies with dimension of the simulation.
    * "ASOE = Active Set of Equations".
    * @return Set of equations. E.g. Rho, Energy, X-Momentum for a 1D pure fluid simulation.
    */
   static constexpr auto ASOE() {return active_equations_;}

   /**
    * @brief Gives the active momentum equations, i.e. varies with dimension of the simulation. "AME = Active Momentum Equations".
    * @return Set of momentum equations.
    */
   static constexpr auto AME() {return active_momenta_;}

   /**
    * @brief Gives the set of prime states which are worked on in this simulation. I. e. varies with dimension of the simulation.
    * "ASOP = Active Set of Prime states".
    * @return Set of prime states. E.g. Rho, Pressure, X-Velocity for a 1D pure fluid simulation.
    */
   static constexpr auto ASOP() {return active_prime_states_;}

   /**
    * @brief Gives the active velocity prime states, i.e. varies with dimension of the simulation. "AV = Active Velocities".
    * @return Set of velocity prime states.
    */
   static constexpr auto AV() {return active_velocities_;}

   /**
    * @brief Gives the equations considered for the coarsening/refinement decision. "EWA = Equations for Wavelet-Analysis".
    * @return List of equations.
    */
   static constexpr auto EWA() {return wavelet_analysis_equations_;}

   /**
    * @brief Gives the set of interface quantities which are used in this simulation. "ASIQ = Active set of interface quantities".
    * @return Set of interface quantities.
    */
   static constexpr auto ASIQ() {return active_interface_quantities_;}

   /**
    * @brief Gives the set of interface quantities which have to be extended at the interface. "IQTE = Interface quantities to extend".
    * @return Set of interface quantities which have to be extended at the interface.
    */
   static constexpr auto IQTE() {return interface_quantities_to_extend_;}

   /**
    * @brief Gives the number of interface quantities which have to be extended at the interface. "NOIQTE = Number of interface quantities to extend".
    * @return Number of interface quantities which have to be extended at the interface.
    */
   static constexpr unsigned int NOIQTE() {return interface_quantities_to_extend_.size();}
};

using FF = FluidFieldsDefinitions;


static_assert( std::make_pair( FF::AME()[0], FF::AV()[0] ) == std::make_pair( Equation::MomentumX, PrimeState::VelocityX ), "FF::AME()[0] and FF::AV()[0] have to consistently return MomentumX and VelocityX" );
#if DIMENSION != 1
static_assert( std::make_pair( FF::AME()[1], FF::AV()[1] ) == std::make_pair( Equation::MomentumY, PrimeState::VelocityY ), "FF::AME()[1] and FF::AV()[1] have to consistently return MomentumY and VelocityY" );
#endif
#if DIMENSION == 3
static_assert( std::make_pair( FF::AME()[2], FF::AV()[2] ) == std::make_pair( Equation::MomentumZ, PrimeState::VelocityZ ), "FF::AME()[2] and FF::AV()[2] have to consistently return MomentumZ and VelocityZ" );
#endif

#endif // FLUID_FIELDS_DEFINITIONS_H
