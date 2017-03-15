#!/usr/bin/env python3

import os.path
import subprocess as sp
from argparse import ArgumentParser as ArgParser

def SetupArgumentParser() :
    parser = ArgParser( prog = "Prepare Alpaca user settings for build",
                        description = "Modify the compile-time user specifications in the code ALPACA" )
    parser.add_argument( "basefolder", help = "The path where the src folder is located" )
    parser.add_argument( "--ic", type = int, help = "Number of internal cells", default = 16 )
    parser.add_argument( "--hs", type = int, help = "Number of halo cells", default = 4 )
    parser.add_argument( "--timeintegration", help = "Name of the time integration scheme.", default = "RK3" )
    parser.add_argument( "--riemann", help = "Name of the Riemann solver. So far either Roe or Hllc", default = "Roe" )
    parser.add_argument( "--fluxsplitting", help = "The flux splitting", default = "Roe" )
    parser.add_argument( "--reconstructionstencil", help = "The stencil that should be used for cell-face state reconstruction", default = "WENO5" )
    parser.add_argument( "--derivativestencil", help = "The stencil that should be used for derivatives", default = "HOUC5" )
    parser.add_argument( "--reinit", help = "The re-initialization method that should be used", default = "Weno" )
    parser.add_argument( "--lsadvect", help = "Indicate which level-set advection scheme should be used", default = "ReconstructionStencil" )
    parser.add_argument( "--interfaceriemann", help = "Choose which Riemann solver should be used at the interface", default = "Linearized" )
    parser.add_argument( "--axisymmetric", help = "If flag is given, axisymmetry will be enabled", action = "store_true" )
    parser.add_argument( "--quite", help = "Disable output of this script", dest = "verbose", action = "store_false"  )
    parser.add_argument( "--limit-end-time", help = "If flag is given, the end time of the simulation will match the given ones in the inputfile exactly"\
                         , dest = "limit_end_time", action = "store_true")
    parser.add_argument( "--gruneisen-density", help = "If flag is given, the Gruneisen parmeter becomes a function of temperature", dest = "gruneisen_density"\
                         , action = "store_true" )
    return parser

def ParseArguments() :
    parser = SetupArgumentParser()
    return parser.parse_args()

def main( arguments ):

    base_path = os.path.abspath( os.path.realpath( arguments.basefolder ) )
    user_specification_path = os.path.join( base_path, "src/user_specifications" )
    compiletime_constants_file = os.path.join( user_specification_path, "compile_time_constants.h" )
    numerical_setup_file = os.path.join( user_specification_path, "numerical_setup.h" )
    riemann_settings_file = os.path.join( user_specification_path, "riemann_solver_settings.h" )
    stencil_setup_file = os.path.join( user_specification_path, "stencil_setup.h" )

    if( arguments.verbose ) :
        print( "Preparing user specifications" )
        print( "Src Path:    " + user_specification_path )
        print( "CC File:     " + compiletime_constants_file )
        print( "Setup File:  " + numerical_setup_file )
        print( "Setting ... to:" )
        print( "  IC                          " + str( arguments.ic ) )
        print( "  HS                          " + str( arguments.hs ) )
        print( "  Integrator                  " + str( arguments.timeintegration ) )
        print( "  Riemann                     " + str( arguments.riemann ) )
        print( "  Fluxsplitting               " + str( arguments.fluxsplitting ) )
        print( "  Reconstruction Stencil      " + str( arguments.reconstructionstencil ) )
        print( "  Derivative Stencil          " + str( arguments.derivativestencil ) )
        print( "  Reinitalizer                " + str( arguments.reinit ) )
        print( "  Advector                    " + str( arguments.lsadvect ) )
        print( "  Interface Riemann           " + str( arguments.interfaceriemann ) )
        print( "  Axissymmetry                " + str( arguments.axisymmetric ) )
        print( "  Limit End Time              " + str( arguments.limit_end_time ) )
        print( "  Gruneisen Density Dependent " + str( arguments.gruneisen_density ) )

    command = ["sed", "-i", "s@_and_dimension_\ =\ [0-9]*;@_and_dimension_\ =\ " + str( arguments.ic ) + ";@", compiletime_constants_file]
    sp.run( command )
    command = ["sed", "-i", "s@halo_width_\ =\ [0-9]*;@halo_width_\ =\ " + str( arguments.hs ) + ";@", compiletime_constants_file]
    sp.run( command )
    command = ["sed", "-i", "s@\ riemann_solver\ =\ .*;@\ riemann_solver\ =\ RiemannSolvers::" + str( arguments.riemann ) + ";@", riemann_settings_file]
    sp.run( command )
    command = ["sed", "-i", "s@x_splitting_scheme\ =\ FluxSplitting::.*;@x_splitting_scheme\ =\ FluxSplitting::" + str( arguments.fluxsplitting ) + ";@", riemann_settings_file]
    sp.run( command )
    command = ["sed", "-i", "s@time_integrator\ =\ .*;@time_integrator\ =\ TimeIntegrators::" + str( arguments.timeintegration ) + ";@", numerical_setup_file]
    sp.run( command )
    command = ["sed", "-i", "s@\ reconstruction_stencil\ =\ .*;@\ reconstruction_stencil\ =\ ReconstructionStencils::" + str( arguments.reconstructionstencil ) + ";@", stencil_setup_file]
    sp.run( command )
    command = ["sed", "-i", "s@\ derivative_stencil\ =\ .*;@\ derivative_stencil\ =\ DerivativeStencils::" + str( arguments.derivativestencil ) + ";@", stencil_setup_file]
    sp.run( command )
    command = ["sed", "-i", "s@interface_riemann_solver\ =\ .*;@interface_riemann_solver\ =\ InterfaceRiemannSolvers::" + str( arguments.interfaceriemann ) + ";@", numerical_setup_file]
    sp.run( command )
    command = ["sed", "-i", "s@levelset_advector\ =\ .*;@levelset_advector\ =\ LevelsetAdvectors::" + str( arguments.lsadvect ) + ";@", numerical_setup_file]
    sp.run( command )
    command = ["sed", "-i", "s@levelset_reinitializer\ =\ .*;@levelset_reinitializer\ =\ LevelsetReinitializers::" + str( arguments.reinit ) + ";@", numerical_setup_file]
    sp.run( command )

    if( arguments.axisymmetric ) :
        command = ["sed", "-i", "s@bool\ axisymmetric_\ .*;@bool\ axisymmetric_\ =\ true;@", compiletime_constants_file]
    else :
        command = ["sed", "-i", "s@bool\ axisymmetric_\ .*;@bool\ axisymmetric_\ =\ false;@", compiletime_constants_file]
    sp.run( command )

    if( arguments.limit_end_time ) :
        command = ["sed", "-i", "s@bool\ limit_end_time_\ .*;@bool\ limit_end_time_\ =\ true;@", compiletime_constants_file]
    else :
        command = ["sed", "-i", "s@bool\ limit_end_time_\ .*;@bool\ limit_end_time_\ =\ false;@", compiletime_constants_file]
    sp.run( command )

    if( arguments.gruneisen_density ) :
        command = ["sed", "-i", "s@bool\ gruneisen_density_dependent_\ .*;@bool\ gruneisen_density_dependent_\ =\ true;@", compiletime_constants_file]
    else :
        command = ["sed", "-i", "s@bool\ gruneisen_density_dependent_\ .*;@bool\ gruneisen_density_dependent_\ =\ false;@", compiletime_constants_file]
    sp.run( command )

    if arguments.timeintegration == 'RK2' :
        command = ["sed", "-i", "s@_time_discretization_order_\ =\ [0-9]*;@_time_discretization_order_\ =\ " + str( 2 ) + ";@", compiletime_constants_file]
    elif arguments.timeintegration == 'RK3' :
        command = ["sed", "-i", "s@_time_discretization_order_\ =\ [0-9]*;@_time_discretization_order_\ =\ " + str( 3 ) + ";@", compiletime_constants_file]
    sp.run( command )


def ModifyCompiletimeSettings( root_directory_of_alpaca, number_of_internal_cells, number_of_halo_cells, time_integrator, riemann_solver,
                               flux_splitting_scheme, reconstruction_stencil, derivative_stencil, levelset_reinitializer, levelset_advector,
                               interfaceriemann, quite = False, limit_end_time = False, gruneisen_density = False, axisymmetric = False ) :

    argument_list = [root_directory_of_alpaca, "--ic", str( number_of_internal_cells ), "--hs", str( number_of_halo_cells ), "--timeintegration",
                     time_integrator, "--riemann", riemann_solver, "--fluxsplitting", flux_splitting_scheme, "--reconstructionstencil",
                     reconstruction_stencil, "--derivativestencil", derivative_stencil, "--reinit", levelset_reinitializer, "--lsadvect",
                     levelset_advector, "--interfaceriemann", interfaceriemann]

    if( axisymmetric ) :
        argument_list.append( "--axisymmetric" )

    if( limit_end_time ) :
        argument_list.append( "--limit-end-time" )

    if( quite ) :
        argument_list.append( "--quite" )

    parser = SetupArgumentParser()
    main( parser.parse_args( argument_list ) )

if __name__ == "__main__":
    main( ParseArguments() )
