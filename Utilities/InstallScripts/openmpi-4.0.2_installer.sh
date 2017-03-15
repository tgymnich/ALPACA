if [ -z ${OPENMPIDIR} ]; then
  echo "\$OPENMPIDIR environment variable is not set! Aborting installation of OpenMPI."
  exit 1
else
  if [ -z ${CUDADIR} ]; then
    echo "Installing OpenMPI into $OPENMPIDIR without CUDA-awareness!"
  else
    echo "Installing OpenMPI into $OPENMPIDIR with CUDA-awareness!"
  fi
fi

# create temporary directory
mkdir openmpi-4.0.2_temp && cd openmpi-4.0.2_temp

# download and extract installer
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz
tar -xf openmpi-4.0.2.tar.gz && cd openmpi-4.0.2

# configure and set installation directory
mkdir $OPENMPIDIR
if [ -z ${CUDADIR} ]; then
  ./configure --prefix=$OPENMPIDIR
else
  ./configure --prefix=$OPENMPIDIR --with-cuda=$CUDADIR
fi
make -j 6 && make install

# clean up
cd ../../ && rm -rf openmpi-4.0.2_temp

echo "==================="
echo "Before usage you need to add the binaries to your \$PATH and \$LD_LIBRARY_PATH environmental variable!"
echo ">> export PATH=\"$OPENMPIDIR/bin:\$PATH\""
echo ">> export LD_LIBRARY_PATH=\"$OPENMPIDIR/lib:\$LD_LIBRARY_PATH\""
echo "==================="