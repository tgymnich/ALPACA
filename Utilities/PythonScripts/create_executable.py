#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from argparse import ArgumentParser as ArgParser
from prepare_user_settings_for_build import ModifyCompiletimeSettings

def SetupArgumentParser() :
    parser = ArgParser( prog = "Create an Alpaca executable", description = "Generates an executable with the provided compile-time setting" )
    parser.add_argument( "--dim", type = int, help = "The dimension that should be simulated", default = 3 )
    parser.add_argument( "--ic", type = int, help = "Number of internal cells", default = 16 )
    parser.add_argument( "--hs", help = "Number of halo cells", type = int, default = 4 )
    parser.add_argument( "--time-integrator", help = "Name of the time integration scheme.", default = "RK3" )
    parser.add_argument( "--riemann-solver", help = "Name of the Riemann solver.", default = "Roe" )
    parser.add_argument( "--flux-splitting", help = "The flux splitting", default = "Roe" )
    parser.add_argument( "--reconstruction-stencil", help = "The stencil that should be used for cell-face state reconstruction", default = "WENO5" )
    parser.add_argument( "--derivative-stencil", help = "The stencil that should be used for derivatives", default = "HOUC5" )
    parser.add_argument( "--levelset-reinitializer", help = "The re-initialization method that should be used", default = "Min" )
    parser.add_argument( "--levelset-advector", help = "Indicate which level-set advection scheme should be used", default = "ReconstructionStencil" )
    parser.add_argument( "--interface-riemann", help = "Choose which Riemann solver should be used at the interface", default = "Linearized" )
    parser.add_argument( "--physics-test", action = "store_true", help = "Indicates whether the compiled executable is a physics test setup." )
    parser.add_argument( "--use-ninja", action = "store_true", help = "Indicates the cmake generator which is either make (default) or ninja" )
    parser.add_argument( "--compile-core-number", help = "Indicates the number of cores to be used during the compilation", default = 1 )
    parser.add_argument( "--build-root-dir", help = "The directory from where the build process is started", default = "", )
    parser.add_argument( "--executable-name", help = "Name of the executable to be build", default = "ALPACA" )
    parser.add_argument( "--exec-dir", help = "The directory where the executable will be placed. Directory will be created if not existing",
                         default = "Utilities/Executables", dest = "exec_directory" )
    parser.add_argument( "--log-and-err-file", help = "Whether a log and error file of the build process should be generated", action = "store_true" )
    parser.add_argument( "basefolder", help = "The path where the src directory is located")
    return parser

def RefineArguments( arguments ) :
    refined_args = arguments
    refined_args.basepath = os.path.abspath( os.path.realpath( refined_args.basefolder ) )
    del refined_args.basefolder

    if( os.path.isabs( arguments.exec_directory ) ) :
        refined_args.exec_path = refined_args.exec_directory
    else :
        refined_args.exec_path = os.path.join( refined_args.basepath, refined_args.exec_directory  )
    del refined_args.exec_directory

    if( os.path.isabs( refined_args.build_root_dir ) ) :
        refined_args.build_base_path = refined_args.build_root_dir
    else :
        refined_args.build_base_path = os.path.join( refined_args.basepath, refined_args.build_root_dir )
    del refined_args.build_root_dir

    if( refined_args.use_ninja ) :
        refined_args.g_ninja = "-GNinja"
        refined_args.make_command = "ninja"
    else :
        refined_args.g_ninja = ""
        refined_args.make_command = "make"
    del refined_args.use_ninja

    return refined_args

def ParseArguments():
    parser = SetupArgumentParser()
    return RefineArguments( parser.parse_args() )

def PrintStreamsToFiles( executable_name, out_stream, error_stream, write_mode ) :
    print( out_stream, file = open( executable_name + ".compile", write_mode ) )
    print( error_stream, file = open( executable_name + ".err", write_mode ) )


def main( options ):

    ModifyCompiletimeSettings( options.basepath, options.ic, options.hs, options.time_integrator, options.riemann_solver, options.flux_splitting,
                               options.reconstruction_stencil, options.derivative_stencil, options.levelset_reinitializer, options.levelset_advector,
                               options.interface_riemann, True )

    if not os.path.exists( options.exec_path ) :
        os.makedirs( options.exec_path, exist_ok = True )

    os.chdir(options.exec_path)
    os.makedirs( "1D", exist_ok = True )
    os.makedirs( "2D", exist_ok = True )
    os.makedirs( "3D", exist_ok = True )

    log = ["PHYSICSTEST: " + str( options.physics_test )]
    if options.physics_test:
        os.makedirs( "2D/PhysicsTest", exist_ok = True )

    if options.physics_test:
       build_path = os.path.join( options.build_base_path, "Build/Physicstests")
    else:
       build_path = os.path.join( options.build_base_path, "Build/" + str( options.dim ) + "D" )
    os.makedirs( build_path, exist_ok = True )
    os.chdir( build_path )

    cmake_process = sp.run( ["cmake", options.g_ninja, "-DCMAKE_BUILD_TYPE=Release", "-DDIM:STRING=" + str( options.dim ), options.basepath],
                            encoding = "utf-8", stdout = sp.PIPE, stderr = sp.PIPE )
    if( options.log_and_err_file ) :
        PrintStreamsToFiles( options.executable_name, cmake_process.stdout, cmake_process.stderr, "w" )

    if( cmake_process.returncode != 0 ) :
        sys.exit( cmake_process.returncode )
    log.append( cmake_process.stdout )

    make_process = sp.run( [options.make_command, "-j", str( options.compile_core_number )], encoding = "utf-8", stdout = sp.PIPE, stderr = sp.PIPE )
    if( options.log_and_err_file ) :
        PrintStreamsToFiles( options.executable_name, make_process.stdout, make_process.stderr, "a" )
    if( make_process.returncode != 0 ):
        sys.exit( make_process.returncode )
    log.append( make_process.stdout )

    if options.physics_test:
        move_command = ["mv", "ALPACA", str(options.exec_path) + "/" + str(options.dim) + "D/" + "PhysicsTest/" + options.executable_name]
    else:
        move_command = ["mv", "ALPACA", str(options.exec_path) + "/" + str(options.dim) + "D/" + options.executable_name]

    move_process = sp.run( move_command, encoding = "utf-8", stdout = sp.PIPE, stderr = sp.PIPE )
    if( options.log_and_err_file ) :
        PrintStreamsToFiles( options.executable_name, move_process.stdout, move_process.stderr, "a" )
    from_dir = str( os.getcwd() )
    log.append( "Moved executable from " + from_dir  + " to " + move_command[2] )

    return log

def CreateExecutable( root_directory_of_alpaca, time_integrator, riemann_solver, flux_splitting_scheme, reconstruction_stencil, derivative_stencil,
                      levelset_reinitializer, levelset_advector, interface_riemann_solver, build_root_directory, executables_directory,
                      dimension = 3, number_of_internal_cells = 16, number_of_halo_cells = 4, executable_name = "ALPACA", use_ninja = False,
                      number_of_cores_for_compilation = 1, log_and_error_file = False, build_physics_test = False ) :

    args_list = ["--dim", str( dimension ), "--ic", str( number_of_internal_cells ), "--hs", str( number_of_halo_cells ), "--time-integrator",
                 time_integrator, "--riemann-solver", riemann_solver, "--flux-splitting", flux_splitting_scheme, "--reconstruction-stencil",
                 reconstruction_stencil, "--derivative-stencil", derivative_stencil, "--levelset-reinitializer", levelset_reinitializer, "--levelset-advector",
                 levelset_advector, "--interface-riemann", interface_riemann_solver, "--build-root-dir",
                 build_root_directory, "--executable-name", str( executable_name ), "--compile-core-number", str( number_of_cores_for_compilation ),
                 "--exec-dir", executables_directory, str( root_directory_of_alpaca ) ]

    if( use_ninja ) :
        args_list.append( "--use-ninja" )

    if( build_physics_test ) :
        args_list.append( "--physics-test" )

    if( log_and_error_file ) :
        args_list.append( "--log-and-err-file" )

    parser = SetupArgumentParser()
    arguments = RefineArguments( parser.parse_args( args_list ) )
    return main( arguments )

if __name__ == "__main__" :
    log = main( ParseArguments() )
    for el in log :
        print( el )