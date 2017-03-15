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
#include "xdmf_file_writer.h"

#include <fstream>
#include <iostream>
#include "user_specifications/compile_time_constants.h"
#include "fluid_fields_definitions.h"
#include "helper_functions.h"
#include "xdmf_utilities.h"
#include "textbased_file_writer.h"

namespace {

   /**
    * @brief Give a properly formatted string for a XDMF file header.
    * @return Header of the xdmf file.
    */
   std::string HeaderInformation() {
      return "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf Version=\"2.0\">\n";
   }

   /**
    * @brief Gives a properly formatted string for the domain, time, topology and geometry nodes in a XDMF file.
    * @param filename The name of the hdf5 file (without path).
    * @param number_of_cells The number of cells in the mesh.
    * @param number_of_vertices The respective number of vertices in the mesh (Note: Due to MPI overlapp this number is not directly linked to
    *        the number of cells.
    * @param time The simulation time of the given mesh (snapshot).
    * @return domain information string.
    */
   std::string DomainInformation( std::string const& filename,  unsigned int const number_of_cells, unsigned int const number_of_vertices, double const time ) {

      return " <Domain>\n"
             "  <Grid Name=\"TimeStep\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
             "   <Time TimeType=\"List\">\n"
             "    <DataItem Format=\"XML\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"1\">\n"
             "     " + ToScientificNotationString( time ) + "\n"
             "    </DataItem>\n"
             "   </Time>\n"
             "   <Grid Name=\"SpatialData\" GridType=\"Uniform\">\n"
             "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" + std::to_string( number_of_cells ) + "\">\n"
             "     <DataItem NumberType=\"Int\" Format=\"HDF\" Dimensions=\"" + std::to_string( number_of_cells ) + " 8\">" + filename + ":/domain/cell_vertices</DataItem>\n"
             "    </Topology>\n"
             "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" + std::to_string( number_of_vertices ) + "\">\n"
             "     <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"" + std::to_string( number_of_vertices ) + " 3\">" + filename + ":/domain/vertex_coordinates</DataItem>\n"
             "    </Geometry>"
             "\n";
   }

   /**
    * @brief Gives a properly formatted string for the partition attribute node in the XDMF file.
    * @param filename The name of the hdf5 file (without path).
    * @param number_of_values The number of values in the dataset (number of cells).
    * @return partition information string.
    */
   std::string PartitionInformation( std::string const& filename, unsigned int const number_of_values ) {
      std::string const data_item("     <DataItem NumberType=\"UInt\" Format=\"HDF\" Dimensions=\"" + std::to_string( number_of_values )
                                  + "\">" + filename + ":/simulation/partition</DataItem>");
      return XdmfUtilities::AttributeString( "partition", data_item, 4) + "\n" ;
   }

   /**
    * @brief Gives a properly formatted string for the levelset attribute node in the XDMF file.
    * @param filename The name of the hdf5 file (without path).
    * @param number_of_values The number of values in the dataset (number of cells).
    * @return levelset information string.
    */
   std::string LevelsetInformation( std::string const& filename, unsigned int const number_of_values ) {
      return XdmfUtilities::AttributeString( "levelset", XdmfUtilities::DataItemString( filename, "simulation/levelset", number_of_values, 5 ), 4) + "\n";
   }

   /**
    * @brief Gives a properly formatted string for the interface quantities attribute nodes in the XDMF file.
    * @param filename The name of the hdf5 file (without path).
    * @param number_of_values The number of values in the dataset (number of cells).
    * @return interface quantity information string.
    */
   std::string InterfaceQuantityInformation( std::string const& filename, unsigned int const number_of_values ) {
      std::string information;
      for( InterfaceQuantity const iq : FF::ASIQ() ) {
         std::string const name( FF::FieldOutputName( iq ) );
         if( !name.empty() ) { // output might be disabled for certain interface quantities (by empty name)
            information += XdmfUtilities::AttributeString( name, XdmfUtilities::DataItemString( filename, "simulation/" + name, number_of_values, 5 ), 4) + "\n";
         }
      }
      return information;
   }

   /**
    * @brief Gives a properly formatted string for the prime states (density, pressure, temperature and velocity) attribute node in the XDMF file.
    * @param filename The name of the hdf5 file (without path).
    * @param number_of_values The number of values in the dataset (number of cells).
    * @return prime state information string.
    */
   std::string PrimeStateInformation( std::string const& filename, unsigned int const number_of_values ) {

      std::string prime_states;
      for( PrimeState const prime : FF::ASOP() ) {
         // skip velocities as they are treated separately
         if( prime == PrimeState::VelocityX || prime == PrimeState::VelocityY || prime == PrimeState::VelocityZ ) continue;
         std::string const name( FF::FieldOutputName( prime ) );
         if( !name.empty() ) { // output might be disabled for certain prime states (by empty name)
            prime_states += XdmfUtilities::AttributeString( name, XdmfUtilities::DataItemString( filename, "simulation/" + name, number_of_values, 5 ), 4) + "\n";
         }
      }

      prime_states += "    <Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Cell\">\n";
      // Join functions are used to mimic a 3D vector for lower dimensions (otherwise vector cannot be visualized in ParaView)
      if constexpr( CC::DIM() == Dimension::One ) {
         prime_states += "     <DataItem ItemType=\"Function\" Function=\"JOIN($0,0*$0,0*$0)\" NumberType=\"Float\" Precision=\"8\""
                                " Dimensions=\"" + std::to_string( number_of_values ) + " 3\">\n";
      } else if constexpr( CC::DIM() == Dimension::Two ) {
         prime_states += "     <DataItem ItemType=\"Function\" Function=\"JOIN($0,  $1,0*$0)\" NumberType=\"Float\" Precision=\"8\""
                                " Dimensions=\"" + std::to_string( number_of_values ) + " 3\">\n";
      } else {
         prime_states += "     <DataItem ItemType=\"Function\" Function=\"JOIN($0,  $1,  $2)\" NumberType=\"Float\" Precision=\"8\""
                                " Dimensions=\"" + std::to_string( number_of_values ) + " 3\">\n";
      }
      for( PrimeState const prime : FF::AV() ) {
         std::string const name( FF::FieldOutputName( prime ) );
         prime_states += XdmfUtilities::DataItemString( filename, "simulation/" + name, number_of_values, 6 ) + "\n";
      }
      prime_states += "     </DataItem>\n    </Attribute>\n";

      return prime_states;
   }

   /**
    * @brief Give a properly formatted string for a XDMF file footer.
    * @return Footer of the xdmf file.
    */
   std::string FooterInformation() {
      return "   </Grid>\n  </Grid>\n </Domain>\n</Xdmf>";
   }
}

// Specification of the static variable output_times_ that it exists
std::vector<double> XdmfFileWriter::output_times_;

/**
 * @brief Default constructor.
 * @param filename The name (including the path) of the xdmf file.
 * @param mpi_rank The MPI rank on which this instance is created.
 * @param print_time_series Flag indicating whether a time series XDMF file should be written.
 * @param number_of_fluids The number of fluids used in the simulation.
 */
XdmfFileWriter::XdmfFileWriter(std::string const filename, int const mpi_rank, bool const print_time_series, unsigned int const number_of_fluids) :
   series_filename_(filename),
   mpi_rank_(mpi_rank),
   number_of_fluids_(number_of_fluids),
   print_time_series_(print_time_series)
{
   // Creates the file where the xdmf information of the full time series is written into
   if(mpi_rank_ == 0 && print_time_series_) {
      std::ofstream output_stream(series_filename_, std::ios::app);

      //We print the header to the time series file
      output_stream << "<?xml version=\"1.0\" ?>" << std::endl;
      output_stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
      output_stream << "<Xdmf Version=\"2.0\">" << std::endl;
      output_stream << " <Domain>" << std::endl;
      output_stream << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

      output_stream.flush();
      output_stream.close();
   }
}

/**
 * @brief The destructor writes the closing entries into the xdmf file before it destroys itself.
 */
XdmfFileWriter::~XdmfFileWriter() {
   if(mpi_rank_ == 0 && print_time_series_) {
      std::ofstream output_stream(series_filename_, std::ios::app);

      //We print the final clauses for the time series file
      output_stream << "   <Time TimeType=\"List\">" << std::endl;
      output_stream << "    <DataItem Format=\"XML\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""<< output_times_.size() <<"\">" << std::endl;
      output_stream << "    ";
      for(auto const& time : output_times_) {
         output_stream << ToScientificNotationString( time ) << " ";
      }
      output_stream << std::endl;
      output_stream << "    </DataItem>" << std::endl;
      output_stream << "   </Time>" << std::endl;
      output_stream << "  </Grid>" << std::endl;
      output_stream << " </Domain>" << std::endl;
      output_stream << "</Xdmf>" << std::endl;

      output_stream.flush();
      output_stream.close();
   }
}

/**
 * @brief Appends all necessary data for the current time step to the xdmf file for the series representation.
 * @param hdf5_filename The name of the accompanying HDF5 file (written in another place).
 * @param global_number_of_cells The global number of cells.
 * @param global_number_of_vertices The global number of vertices.
 */
void XdmfFileWriter::AppendTimestep(std::string const hdf5_filename, unsigned int const global_number_of_cells, unsigned int const global_number_of_vertices) const {

   std::string short_file_name = hdf5_filename;
   short_file_name.erase(0,short_file_name.find_last_of("/")+1);

   std::ofstream output_stream(series_filename_, std::ios::app);

   /*** DOMAIN ***/
   output_stream << "   <Grid Name=\"step_" << output_times_.size() << "\" GridType=\"Uniform\">" << std::endl;
   output_stream << "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << global_number_of_cells << "\">" << std::endl;
   output_stream << "     <DataItem NumberType=\"Int\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << " 8\">" << short_file_name << ":/domain/cell_vertices</DataItem>" << std::endl;
   output_stream << "    </Topology>" << std::endl;
   output_stream << "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" << global_number_of_vertices << "\">" << std::endl;
   output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_vertices << " 3\">" << short_file_name << ":/domain/vertex_coordinates</DataItem>" << std::endl;
   output_stream << "    </Geometry>" << std::endl;

   // /*** PRIMES ***/
   output_stream << PrimeStateInformation(short_file_name, global_number_of_cells);

   /*** PARTITION ***/
   output_stream << PartitionInformation(short_file_name, global_number_of_cells);

   /*** LEVELSET ***/
   output_stream << LevelsetInformation(short_file_name, global_number_of_cells);

   /*** INTERFACE QUANTITIES ***/
   output_stream << InterfaceQuantityInformation(short_file_name, global_number_of_cells);

   /*** FOOTER ***/
   output_stream << "   </Grid>" << std::endl;

   output_stream.flush();
   output_stream.close();
}

/**
 * @brief Create a a full XDMF file for the current time step, and the current time step only.
 * @param hdf5_filename Filename of the hdf5 file (without path).
 * @param global_number_of_cells The The number of cells in the output.
 * @param global_number_of_vertices The number of vertices in the output.
 * @param time Output time for which the xdmf file should be created.
 */
void XdmfFileWriter::CreateTimestepFile( std::string const& hdf5_file_with_path, unsigned int const global_number_of_cells,
                                         unsigned int const global_number_of_vertices, double const time ) const {

   std::string const hdf5_filename( RemoveFilePath( hdf5_file_with_path ) );

   std::string const content(
      HeaderInformation() +
      DomainInformation( hdf5_filename, global_number_of_cells, global_number_of_vertices, time ) +
      PrimeStateInformation( hdf5_filename, global_number_of_cells ) +
      PartitionInformation( hdf5_filename, global_number_of_cells ) +
      LevelsetInformation( hdf5_filename, global_number_of_cells ) +
      InterfaceQuantityInformation( hdf5_filename, global_number_of_cells ) +
      FooterInformation()
   );

   // Write the xdmf file to hard disk (with same name as hdf5 but with xdmf extension)
   SimpleFileWriter::WriteTextBasedFile( ChangeFileExtension( hdf5_file_with_path, ".xdmf" ), content );
}

/**
 * @brief Writes a single XDMF file for the current time step and appends the series-XDMF with the current step.
 * @param timestep The time step of the current output.
 * @param hdf5_filename The name of the externally (= on another class) created (heavy data) h5 file.
 * @param global_number_of_cells The The number of cells in the output.
 * @param global_number_of_vertices The number of vertices in the output.
 */
void XdmfFileWriter::WriteXdmfForTimestep(const double timestep, std::string const hdf5_file_with_path, unsigned int const global_number_of_cells, unsigned int const global_number_of_vertices) const {
   if(mpi_rank_ == 0) {
      output_times_.push_back(timestep);
      AppendTimestep( RemoveFilePath( hdf5_file_with_path ),global_number_of_cells,global_number_of_vertices);
      CreateTimestepFile( hdf5_file_with_path,global_number_of_cells,global_number_of_vertices, timestep );
   }
}

/**
 * @brief Creates a single XDMF file for the current debug output.
 * @param timestep The time of the current output.
 * @param hdf5_filename The name of the accompanying HDF5 file (written in another place, without path).
 * @param global_number_of_cells The number of cells in this output (across all ranks).
 * @param global_number_of_vertices The number of vertices in this output (across all ranks).
 */
void XdmfFileWriter::WriteDebugXdmfFile(const double timestep, std::string const hdf5_filename, unsigned int const global_number_of_cells, unsigned int const global_number_of_vertices) const {

   std::string const fluid = "fluid";
   std::vector<std::string> equations = {{ "density", "pressure", "temperature", "velocityX", "mass_avg", "mass_rhs", "energy_avg", "energy_rhs", "momentumX_avg", "momentumX_rhs" }};
   if constexpr(CC::DIM() != Dimension::One) {
      equations.push_back("momentumY_avg");
      equations.push_back("momentumY_rhs");
      equations.push_back("velocityY");
   }
   if constexpr(CC::DIM() == Dimension::Three) {
      equations.push_back("momentumZ_avg");
      equations.push_back("momentumZ_rhs");
      equations.push_back("velocityZ");
   }

   if(mpi_rank_ == 0) {
      std::string short_file_name = hdf5_filename;
      short_file_name.erase(0,short_file_name.find_last_of("/")+1);
      std::string xdmf_filename = hdf5_filename.substr(0,hdf5_filename.size()-3) + ".xdmf"; //Cut ".h5" = three caracters.

      std::ofstream output_stream(xdmf_filename, std::ios::trunc);
      /*** HEADER ***/
      output_stream << "<?xml version=\"1.0\" ?>" << std::endl;
      output_stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
      output_stream << "<Xdmf Version=\"2.0\">" << std::endl;
      output_stream << " <Domain>" << std::endl;
      output_stream << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
      output_stream << "   <Grid Name=\"step_" << timestep << "\" GridType=\"Uniform\">" << std::endl;
      output_stream << "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << global_number_of_cells << "\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Int\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << " 8\">" << short_file_name << ":/domain/cell_vertices</DataItem>" << std::endl;
      output_stream << "    </Topology>" << std::endl;
      output_stream << "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" << global_number_of_vertices << "\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_vertices << " 3\">" << short_file_name << ":/domain/vertex_coordinates</DataItem>" << std::endl;
      output_stream << "    </Geometry>" << std::endl;

      /*** Simulation Values ***/
      std::vector<std::string> names_of_quantities;
      for( Equation const eq : FF::ASOE() ) {
         std::string const equation_name( FF::FieldOutputName( eq ) );
         if( !equation_name.empty() ) {
            names_of_quantities.emplace_back(equation_name + "_avg");
            names_of_quantities.emplace_back(equation_name + "_rhs");
         }
      }
      for( PrimeState const ps : FF::ASOP() ) {
         if( !FF::FieldOutputName( ps ).empty() ) {
            names_of_quantities.emplace_back( FF::FieldOutputName( ps ) );
         }
      }

      for(unsigned int i = 0; i < number_of_fluids_; ++i) {
         std::string const fluid_name = fluid + "_" + std::to_string(i) + "_";
         for(std::string const& quantity_name : names_of_quantities) {
            output_stream << "    <Attribute Name=\"" + fluid_name + quantity_name + "\" Center=\"Cell\">" << std::endl;
            output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/" + fluid_name + quantity_name  + "</DataItem>" << std::endl;
            output_stream << "    </Attribute>" << std::endl;
         }
      }

      output_stream << "    <Attribute Name=\"partition\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"UInt\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells <<"\">" << short_file_name << ":/simulation/partition</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;

      /*** LEVELSET ***/
      output_stream << "    <Attribute Name=\"interface_tags\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/interface_tags</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;
      output_stream << "    <Attribute Name=\"levelset\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/levelset</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;
      output_stream << "    <Attribute Name=\"levelset_rhs\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/levelset_rhs</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;
      output_stream << "    <Attribute Name=\"levelset_reinitialized\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/levelset_reinitialized</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;
      output_stream << "    <Attribute Name=\"volume_fraction\" Center=\"Cell\">" << std::endl;
      output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/volume_fraction</DataItem>" << std::endl;
      output_stream << "    </Attribute>" << std::endl;
      for( InterfaceQuantity const iq : FF::ASIQ() ) {
         std::string const name( FF::FieldOutputName( iq ) );
         if( !name.empty() ) {
            output_stream << "    <Attribute Name=\"" + name +"\" Center=\"Cell\">" << std::endl;
            output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/simulation/" + name +"</DataItem>" << std::endl;
            output_stream << "    </Attribute>" << std::endl;
         }
      }

      /*** FOOTER ***/
      output_stream << "   </Grid>" << std::endl;
      output_stream << "  </Grid>" << std::endl;
      output_stream << " </Domain>" << std::endl;
      output_stream << "</Xdmf>" << std::endl;
      output_stream.flush();
      output_stream.close();
   }
}
