#!/usr/bin/env python3

import os
import xml.etree.ElementTree as et
from argparse import ArgumentParser as ArgParser

def SetupArgumentParser() :
    parser = ArgParser(prog = "Modify Alpaca xml input file", description = "Modify a specified inputfile for the Code ALPACA" )
    parser.add_argument( "--blocksize", type = float, help = "The block size of a single block", default = None )
    parser.add_argument( "--blockratio_x", type = int, help = "Number of blocks in x-direction", default = None )
    parser.add_argument( "--blockratio_y", type = int, help = "Number of blocks in y-direction", default = None )
    parser.add_argument( "--blockratio_z", type = int, help = "Number of blocks in z-direction", default = None )
    parser.add_argument( "--maximumlevel", type = int, help = "The maximum level of the multi-resolution scheme.", default = None )
    parser.add_argument( "--epsilonreference", type = float, help = "The epsilon for the multi-resolution scheme.", default = None )
    parser.add_argument( "--referencelevel", type = int, help = "The reference level for the multi-resolution scheme.", default = None )
    parser.add_argument( "--restartsimulation", action = "store_true", help = "Indicates whether simulation should be restarted or not." )
    parser.add_argument( "--restartfilename", help = "The name of the restart file that should be used.", default = None )
    parser.add_argument( "--outputname", help = "Name of the modified inputfile. If not set, take the name of the original inputfile.", default = None )
    parser.add_argument( "inputfile", help = "Original inputfile that should be modified" )
    return parser

def ParseArguments() :
    parser = SetupArgumentParser()
    return parser.parse_args()

def main( options ):

    if options.outputname is None:
        options.outputname = options.inputfile

    options.outputname = os.path.abspath(os.path.realpath(options.outputname))

    options.inputfile = os.path.abspath(os.path.realpath(options.inputfile))

    if not os.path.isdir(os.path.dirname(options.outputname)):
        os.makedirs(os.path.dirname(options.outputname))

    tree = et.parse(options.inputfile)
    root = tree.getroot()

    for domain in root.findall("domain"):
        if options.blocksize is not None:
            for child in domain.findall("blockSize"):
                child.text = " "+str(options.blocksize)+" "
        for blockratio in domain.findall("blockRatio"):
            if options.blockratio_x is not None:
                for child in blockratio.findall("x"):
                    child.text = " "+str(options.blockratio_x)+" "
            if options.blockratio_y is not None:
                for child in blockratio.findall("y"):
                    child.text = " "+str(options.blockratio_y)+" "
            if options.blockratio_z is not None:
                for child in blockratio.findall("z"):
                    child.text = " "+str(options.blockratio_z)+" "

    for multiresolution in root.findall("multiResolution"):
        if options.maximumlevel is not None:
            for child in multiresolution.findall("maximumLevel"):
                child.text = " "+str(options.maximumlevel)+" "
            for refinementcriterion in multiresolution.findall("refinementCriterion"):
                if options.epsilonreference is not None:
                    for child in refinementcriterion.findall("epsilonReference"):
                        child.text = " "+str(options.epsilonreference)+" "
                if options.referencelevel is not None:
                    for child in refinementcriterion.findall("levelOfEpsilonReference"):
                        child.text = " "+str(options.referencelevel)+" "

    for restart in root.findall("restart"):
        for child in restart.findall("restoreMode"):
            if options.restartsimulation:
                child.text = " 1 "
            else:
                child.text = " 0 "
        if options.restartfilename is not None:
            for child in restart.findall("restoreFileName"):
                child.text = " "+options.restartfilename+" "

    tree.write(options.outputname)

def ModifyXmlInputfile( base_inputfile, block_size = None, block_ratio_x = None, block_ratio_y = None, block_ratio_z = None, maximum_level = None,
                        reference_epsilon = None, reference_level = None, output_file_name = None, restart_simulation = False, restart_filename = None ) :

    argument_list = [str( base_inputfile )]
    if( block_size is not None ) :
        argument_list.append( "--blocksize" )
        argument_list.append( str( block_size ) )
    if( block_ratio_x is not None ) :
        argument_list.append( "--blockratio_x" )
        argument_list.append( str( block_ratio_x ) )
    if( block_ratio_y is not None ) :
        argument_list.append( "--blockratio_y" )
        argument_list.append( str( block_ratio_y ) )
    if( block_ratio_z is not None ) :
        argument_list.append( "--blockratio_z" )
        argument_list.append( str( block_ratio_z ) )
    if( maximum_level is not None ) :
        argument_list.append( "--maximumlevel" )
        argument_list.append( str( maximum_level ) )
    if( reference_epsilon is not None ) :
        argument_list.append( "--epsilonreference" )
        argument_list.append( str( reference_epsilon ) )
    if( reference_level is not None ) :
        argument_list.append( "--referencelevel" )
        argument_list.append( str( reference_level ) )
    if( output_file_name is not None ) :
        argument_list.append( "--outputname" )
        argument_list.append( str( output_file_name ) )
    if( restart_filename is not None ) :
        argument_list.append( "--restartfilename" )
        argument_list.append( str( restart_filename ) )
    if( restart_simulation ) :
        argument_list.append( "--restartsimulation" )

    parser = SetupArgumentParser()
    main( parser.parse_args( argument_list ) )

if __name__ == "__main__":
    main( ParseArguments() )