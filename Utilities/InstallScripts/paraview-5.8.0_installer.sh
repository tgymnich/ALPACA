if [ -z ${PARAVIEWDIR} ]; then
  echo "\$PARAVIEWDIR environment variable is not set! Aborting installation of Paraview."
  exit 1
else
  echo "Installing Paraview into $PARAVIEWDIR"
fi

# Sanity checks
if [ -d "$PARAVIEWDIR" ] 
then
  echo "$PARAVIEWDIR already exists. Please choose another folder to install Paraview!"
  exit 1
fi

# download and extract installer
wget "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.8&type=binary&os=Linux&downloadFile=ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz"
tar xf "download.php?submit=Download&version=v5.8&type=binary&os=Linux&downloadFile=ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz"
mv ParaView-5.8.0-MPI-Linux-Python3.7-64bit $PARAVIEWDIR

# clean up
rm -rf "download.php?submit=Download&version=v5.8&type=binary&os=Linux&downloadFile=ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz"

echo "==================="
echo "Before usage you need to add the binaries to your \$PATH environmental variable!"
echo ">> export PATH=\"$PARAVIEWDIR/bin:\$PATH\""
echo "==================="
