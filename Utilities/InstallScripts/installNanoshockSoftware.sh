if [ -z ${ALPACAUTILDIR} ]; then
  echo "\$ALPACAUTILDIR environment variable is not set! Aborting installation of software necessary to compile ALPACA."
  exit 1
else
  echo "Installing software necessary to compile ALPACA into $ALPACAUTILDIR"
fi

#Set installer scripts
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
setupFile="${ALPACAUTILDIR}/setupAlpacaTools.sh"
mpichInstaller="${scriptDir}/mpich-3.3.2_installer.sh"
hdf5Installer="${scriptDir}/hdf5-1.8.21_installer.sh"
cmakeInstaller="${scriptDir}/cmake-3.16.0_installer.sh"
vscodeInstaller="${scriptDir}/vscode-1.44.0_installer.sh"
paraviewInstaller="${scriptDir}/paraview-5.8.0_installer.sh"

# Create directory
mkdir -p ${ALPACAUTILDIR} && cd ${ALPACAUTILDIR}

while true; do
    read -p "Do you wish to install mpich-3.3.2? [Yy, Nn]" yn
    case $yn in
        [Yy]* ) installMpich=true; break;;
        [Nn]* ) installMpich=false; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to install hdf5-1.8.21? [Yy, Nn]" yn
    case $yn in
        [Yy]* ) installHdf5=true; break;;
        [Nn]* ) installHdf5=false; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to install cmake-3.16.0? [Yy, Nn]" yn
    case $yn in
        [Yy]* ) installCmake=true; break;;
        [Nn]* ) installCmake=false; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to install VSCode-1.44.0? [Yy, Nn]" yn
    case $yn in
        [Yy]* ) installVscode=true; break;;
        [Nn]* ) installVscode=false; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to install Paraview-5.8.0? [Yy, Nn]" yn
    case $yn in
        [Yy]* ) installParaview=true; break;;
        [Nn]* ) installParaview=false; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

if $installMpich; then
  # Install MPICH
  export MPICHDIR="${ALPACAUTILDIR}/mpich-3.3.2"
  echo "Install MPICH into ${MPICHDIR}"
  bash ${mpichInstaller}
  # export path to use newly installed MPI compilers
  export PATH="$MPICHDIR/bin:$PATH"
  export LD_LIBRARY_PATH="$MPICHDIR/include:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH="$MPICHDIR/lib:$LD_LIBRARY_PATH"
  echo "# MPICH Compiler" >> ${setupFile}
  echo "export PATH=\"${MPICHDIR}/bin:\$PATH\"" >> ${setupFile}
  echo "export LD_LIBRARY_PATH=\"${MPICHDIR}/include:\$LD_LIBRARY_PATH\"" >> ${setupFile}
  echo "export LD_LIBRARY_PATH=\"${MPICHDIR}/lib:\$LD_LIBRARY_PATH\"" >> ${setupFile}
fi

if $installHdf5; then
  # Install HDF5
  export HDF5DIR="${ALPACAUTILDIR}/hdf5-1.8.21"
  echo "Install HDF5 into ${HDF5DIR}"
  bash ${hdf5Installer}
  # export path to use newly installed HDF5 libraries
  export PATH="$HDF5DIR/bin:$PATH"
  export LD_LIBRARY_PATH="$HDF5DIR/include:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH="$HDF5DIR/lib:$LD_LIBRARY_PATH"
  echo "# HDF5" >> ${setupFile}
  echo "export PATH=\"${HDF5DIR}/bin:\$PATH\"" >> ${setupFile}
  echo "export LD_LIBRARY_PATH=\"${HDF5DIR}/include:\$LD_LIBRARY_PATH\"" >> ${setupFile}
  echo "export LD_LIBRARY_PATH=\"${HDF5DIR}/lib:\$LD_LIBRARY_PATH\"" >> ${setupFile}
fi

if $installCmake; then
  # Install CMake
  export CMAKEDIR="${ALPACAUTILDIR}/cmake-3.16.0"
  echo "Install CMake into ${CMAKEDIR}"
  bash ${cmakeInstaller}
  # export path to use newly installed CMake libraries
  export PATH="$CMAKEDIR/bin:$PATH"
  echo "# CMake" >> ${setupFile}
  echo "export PATH=\"${CMAKEDIR}/bin:\$PATH\"" >> ${setupFile}
fi

if $installVscode; then
  # Install VSCode
  export VSCODEDIR="${ALPACAUTILDIR}/vscode-1.44.0"
  echo "Install VSCode into ${VSCODEDIR}"
  bash ${vscodeInstaller}
  echo "# VSCode" >> ${setupFile}
  echo "alias vscode='${VSCODEDIR}/bin/code --extensions-dir=\"${VSCODEDIR}/.vscode/extensions\" --user-data-dir=\"${VSCODEDIR}/.vscode/settings\"'" >> ${setupFile}
fi

if $installParaview; then
  # Install Paraview 5.8
  export PARAVIEWDIR="${ALPACAUTILDIR}/paraview-5.8.0"
  echo "Install Paraview 5.8 into ${PARAVIEWDIR}"
  bash ${paraviewInstaller}
  # export path to use newly installed Paraview 5.8 libraries
  export PATH="$PARAVIEWDIR/bin:$PATH"
  echo "# Paraview 5.8" >> ${setupFile}
  echo "alias paraview='${PARAVIEWDIR}/bin/paraview'" >> ${setupFile}
  echo "alias pvbatch='${PARAVIEWDIR}/bin/pvbatch'" >> ${setupFile}
  echo "alias pvserver='${PARAVIEWDIR}/bin/pvserver'" >> ${setupFile}
  echo "alias pvpython='${PARAVIEWDIR}/bin/pvpython'" >> ${setupFile}
  echo "" >> ${setupFile}
fi
