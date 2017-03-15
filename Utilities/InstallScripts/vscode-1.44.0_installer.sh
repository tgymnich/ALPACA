if [ -z ${VSCODEDIR} ]; then
  echo "\$VSCODEDIR environment variable is not set! Aborting installation of VSCode."
  exit 1
else
  echo "Installing VSCode into $VSCODEDIR"
fi

# Sanity checks
if [ -d "$VSCODEDIR" ] 
then
  echo "$VSCODEDIR already exists. Please choose another folder to install VSCode!"
  exit 1
fi

# download and extract installer
wget https://update.code.visualstudio.com/1.44.0/linux-x64/stable
tar -xf stable && mv VSCode-linux-x64 $VSCODEDIR

# clean up
rm -rf stable

echo "==================="
echo "Before usage you need to add an alias!"
echo ">> alias vscode='${VSCODEDIR}/bin/code --extensions-dir=\"${VSCODEDIR}/.vscode/extensions\" --user-data-dir=\"${VSCODEDIR}/.vscode/settings\"'"
echo "==================="
