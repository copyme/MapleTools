#!/bin/sh

# Copyright (c) 2017, Kacper Pluta <kacper.pluta@esiee.fr>
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of RigidMotionsMapleTools nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# The list of files to be installed
FILES=("RealAlgebraicNumber.mpl" "RigidMotionsParameterSpaceEvent.mpl"
"third-party/MaplePrimesCode.mpl" "RigidMotionsParameterSpaceCommon.mpl"
"RigidMotionsParameterSpaceRegister.mpl" "RigidMotionsParameterSpaceDecompositionRecursive.mpl"
"RigidMotionsParameterSpaceDecomposition.mpl" )

# This function checks if all the files exist in the install directory
function Check_Files()
{
  for i in ${FILES[@]}; do
    if [ ! -e "${TOOLS_DIR}/${i}" ]; then
      echo "The file: \"${TOOLS_DIR}/${i}\" does not exists!"
      return 1
    fi
  done
  return 0
}

# The main installation function. Its aim is to add correct paths to the ~/.mapleinit file.
function Install()
{
  if [ -d "${TOOLS_DIR}" ] && Check_Files; then
    for i in ${FILES[@]}; do
      echo "read \"${TOOLS_DIR}/${i}\":" >> "${HOME}/.mapleinit"
    done
  elif [ ! -d "${TOOLS_DIR}" ]; then
    echo "The directory: \"${TOOLS_DIR}\" does not exists!"
  fi
}

# This function simple clean ups ~/.mapleinit from the entries created by the Install() function.
function Uninstall()
{
  if [ -e "${HOME}/.mapleinit" ]; then
    for i in ${FILES[@]}; do
      sed -i "/$( basename $i)/d" "${HOME}/.mapleinit"
    done
  fi;
}

# We check if argument is valid and call a specific function related to it.
function Parse_Arguments()
{
  for i in ${@}; do
  case ${i} in
    -d=*|--directory=*)
      TOOLS_DIR="${i#*=}"
      shift # past argument
      Install
    ;;
    -u|--uninstall)
      shift # past argument
      Uninstall
    ;;
    *)
    ;;
  esac
  done
}

# Check if some input parameters were passed.
case ${@} in
  (*[![:blank:]]*)
     Parse_Arguments ${@}
     ;;
  (*)
     echo 'To install RigidMotionsMapleTools type: Install.sh -d= | --directory=<"/path/to/MapleTools">'
     echo 'To clean up ~/.mapleinit type: Install.sh -u | --uninstall'
     exit
    ;;
esac


