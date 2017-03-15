if [ -z ${CMAKEDIR} ]; then
  echo "\$CMAKEDIR environment variable is not set! Aborting installation of CMake."
  exit 1
else
  echo "Installing CMake into $CMAKEDIR"
fi

# create temporary directory
mkdir cmake-3.16.0_temp && cd cmake-3.16.0_temp

# download and extract installer
wget https://github.com/Kitware/CMake/releases/download/v3.16.0-rc3/cmake-3.16.0-rc3.tar.gz
tar -xf cmake-3.16.0-rc3.tar.gz && cd cmake-3.16.0-rc3

# configure and set installation directory
mkdir $CMAKEDIR
./configure --prefix=$CMAKEDIR
make && make install

# clean up
cd ../../ && rm -rf cmake-3.16.0_temp

echo "==================="
echo "Before usage you need to add the binaries to your \$PATH environmental variable!"
echo ">> export PATH=\"$CMAKEDIR/bin:\$PATH\""
echo "==================="