if [ -z ${MPICHDIR} ]; then
  echo "\$MPICHDIR environment variable is not set! Aborting installation of MPICH."
  exit 1
else
  if [ -z ${CUDADIR} ]; then
    echo "Installing MPICH into $MPICHDIR without CUDA-awareness!"
  else
    echo "Installing MPICH into $MPICHDIR with CUDA-awareness!"
  fi
fi

# create temporary directory
mkdir mpich-3.3.2_temp && cd mpich-3.3.2_temp

# download and extract installer
wget http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
tar -xf mpich-3.3.2.tar.gz && cd mpich-3.3.2

# configure and set installation directory
mkdir $MPICHDIR
if [ -z ${CUDADIR} ]; then
  ./configure --prefix=$MPICHDIR --enable-fast=all,O3 --enable-g=dbg
else
  ./configure --prefix=$MPICHDIR --enable-fast=all,O3 --enable-g=dbg --with-cuda=$CUDADIR
fi
make -j 6 && make install

# clean up
cd ../../ && rm -rf mpich-3.3.2_temp

echo "==================="
echo "Before usage you need to add the binaries to your \$PATH and \$LD_LIBRARY_PATH environmental variable!"
echo ">> export PATH=\"$MPICHDIR/bin:\$PATH\""
echo ">> export LD_LIBRARY_PATH=\"$MPICHDIR/include:\$LD_LIBRARY_PATH\""
echo ">> export LD_LIBRARY_PATH=\"$MPICHDIR/lib:\$LD_LIBRARY_PATH\""
echo "==================="