if [ -z ${HDF5DIR} ]; then
  echo "\$HDF5DIR environment variable is not set! Aborting installation of OpenMPI."
  exit 1
else
  echo "Installing HDF5 into $HDF5DIR"
fi

# create temporary directory
mkdir hdf5-1.8.21_temp && cd hdf5-1.8.21_temp

# download and extract installer
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz
tar -xf hdf5-1.8.21.tar.gz && cd hdf5-1.8.21

# configure and set installation directory
mkdir $HDF5DIR

export CXX=mpic++
export CC=mpicc
./configure --prefix=$HDF5DIR --enable-cxx --enable-parallel --enable-unsupported
make -j 6 && make install

# clean up
cd ../../ && rm -rf hdf5-1.8.21_temp

echo "==================="
echo "Before usage you need to add the binaries to your \$PATH and \$LD_LIBRARY_PATH environmental variable!"
echo ">> export PATH=\"$HDF5DIR/bin:\$PATH\""
echo ">> export LD_LIBRARY_PATH=\"$HDF5DIR/include:\$LD_LIBRARY_PATH\""
echo ">> export LD_LIBRARY_PATH=\"$HDF5DIR/lib:\$LD_LIBRARY_PATH\""
echo "==================="