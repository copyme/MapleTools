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

# Get rid of national specifications like "," instead of "." for numbers.
export LC_ALL=C

BUILD_DIR="./Installer"
INSTALL_SCRIPT="Install.sh"
NODE_RUNNER_SCRIPT="NodeRunner.sh"
EXEC_SCRIPT="Exec.sh"
ARCHIVE_NAME="DATA"

SOURCE_DIR="../.."
SOURCE_FILES=("RealAlgebraicNumber.mpl" "RigidMotionsParameterSpaceEvent.mpl"
"RigidMotionsParameterSpaceCommon.mpl"
"RigidMotionsParameterSpaceRegister.mpl" "RigidMotionsParameterSpaceDecompositionRecursive.mpl"
"RigidMotionsParameterSpaceDecomposition.mpl" )


THIRD_PARTY_SOURCE_DIR="../../third-party"
THIRD_PARTY_SOURCE=("MaplePrimesCode.mpl")

function Create_Dirs()
{
  # Set the directory in which we create splits of the database.
  if [ ! -d "${BUILD_DIR}" ]; then
    mkdir  "${BUILD_DIR}"
    mkdir  "${BUILD_DIR}/DB"
    mkdir "${BUILD_DIR}/third-party"
  else
    rm -rf "${BUILD_DIR}/"*
    mkdir "${BUILD_DIR}/DB"
    mkdir "${BUILD_DIR}/third-party"
  fi
  rm -f "${ARCHIVE_NAME}.tar"
  rm -f "${ARCHIVE_NAME}.tar.gz"
}


# This function splits the main database into chunks which size depends on the number of nodes in
# the cluster
function Split_DB()
{

  if ! hash sqlite3 2>/dev/null; then
    echo "This script relay on SQLite3! You need to install it."
    exit 1
  fi
  if [ ! -e "${DB_FILE}" ]; then
    echo "The database file: \"${DB_FILE}\" does not exist!"
    exit 1
  fi

  # Set the directory in which we create splits of the database.
  if [ ! -d "${BUILD_DIR}/DB" ]; then
    echo "The database directory: \"${BUILD_DIR}/DB\" does not exists!"
    exit 1
  else
    cd  "${BUILD_DIR}/DB"
  fi 

  EVENTS_PER_NODE=`echo -e "SELECT ROUND(COUNT(*)/${NO_NODES}) FROM RealAlgebraicNumber;" | sqlite3 "${DB_FILE}"`
  # Convert to integer
  EVENTS_PER_NODE=`printf '%.*f' 0 ${EVENTS_PER_NODE}`

  for i in $( seq 0 $(( NO_NODES - 1 )) ); do
    # Recreate schemes
    sqlite3 "${DB_FILE}" '.schema SamplePoint' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema ComputedNumbers' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema Events' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema Quadric' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema RealAlgebraicNumber' | sqlite3 "${i}.db"
    echo "PRAGMA synchronous = NORMAL;" | sqlite3 "${i}.db" > /dev/null
    echo "PRAGMA journal_mode = WAL;" | sqlite3 "${i}.db" > /dev/null
    # Copy data
    sqlite3 "${DB_FILE}" "SELECT * FROM Quadric;" | sqlite3 "${i}.db" ".import /dev/stdin Quadric"
    sqlite3 "${DB_FILE}" "SELECT * FROM RealAlgebraicNumber WHERE ID BETWEEN $(( i * EVENTS_PER_NODE
    )) AND $(( ( i + 1 ) * (EVENTS_PER_NODE + 1) ));" | sqlite3 "${i}.db" ".import /dev/stdin RealAlgebraicNumber"
    sqlite3 "${DB_FILE}" "SELECT * FROM Events WHERE RANumID BETWEEN $(( i * EVENTS_PER_NODE
    )) AND $(( ( i + 1 ) * (EVENTS_PER_NODE + 1) ));" | sqlite3 "${i}.db" ".import /dev/stdin Events"
  done
  # Restore path
  cd "../.."
}

# This function put all necessary files into the BUILD_DIR and compress them.
function Prepare_Files()
{
  if [ ! -e "../${INSTALL_SCRIPT}" ]; then
    echo "The \"${INSTALL_SCRIPT}\" file does not exist!"
    exit 1
  fi
  cp "../${INSTALL_SCRIPT}" "${BUILD_DIR}" 
  chmod +x "${BUILD_DIR}/${INSTALL_SCRIPT}"

  if [ ! -e "./${NODE_RUNNER_SCRIPT}" ]; then
    echo "The \"${NODE_RUNNER_SCRIPT}\" file does not exist!"
    exit 1
  fi
  cp "./${NODE_RUNNER_SCRIPT}" "${BUILD_DIR}" 
  chmod +x "${BUILD_DIR}/${NODE_RUNNER_SCRIPT}"
  # Copy Maple script
  cp "${MAPLE_FILE}" "${BUILD_DIR}"

  # Copy source files
  for i in ${SOURCE_FILES[@]}; do
    if [ ! -e "${SOURCE_DIR}/${i}" ]; then
      echo "The file: \"${SOURCE_DIR}/${i}\" does not exist!"
      exit 1
    else
      cp "${SOURCE_DIR}/${i}" "${BUILD_DIR}"
    fi
  done

  # Third party source code
  for i in ${THIRD_PARTY_SOURCE[@]}; do
    if [ ! -e "${THIRD_PARTY_SOURCE_DIR}/${i}" ]; then
      echo "The file: \"${THIRD_PARTY_SOURCE_DIR}/${i}\" does not exist!"
      exit 1
    else
      cp "${THIRD_PARTY_SOURCE_DIR}/${i}" "${BUILD_DIR}/third-party"
    fi 
  done

  # Compress files
  cd "${BUILD_DIR}"
  tar cf "../${ARCHIVE_NAME}.tar" ./*
  cd ..
  if [ -e "./${ARCHIVE_NAME}.tar" ]; then
    gzip "${ARCHIVE_NAME}.tar"
  else
    echo  "The file: \"${ARCHIVE_NAME}.tar\" does not exist!"
    exit 1
  fi

  # Set variables in the EXEC_SCRIPT
  if [ ! -e "./${EXEC_SCRIPT}" ]; then
    echo "The \"${EXEC_SCRIPT}\" file does not exist!"
    exit 1
  else
    TMP_EXEC=$(mktemp -u "Exec.XXXXXX.sh")
    cp "./${EXEC_SCRIPT}" "${TMP_EXEC}"
    sed -i "s~__REPLACE_MPL__~\"$( basename ${MAPLE_FILE} )\"~g" "${TMP_EXEC}"
  fi

  # Create a self executable archive
  if [ -e "${ARCHIVE_NAME}.tar.gz" ]; then
    cat "${TMP_EXEC}" "${ARCHIVE_NAME}.tar.gz" > "${OUTPUT}"
    rm "${ARCHIVE_NAME}.tar.gz"
    rm "${TMP_EXEC}"
  else
    echo "The file: \"${ARCHIVE_NAME}.tar.gz\" does not exist!"
    exit 1
  fi

  # Test if the archive was created.
  if [ -e "${OUTPUT}" ]; then
    chmod +x "${OUTPUT}"
    MD5=($( md5sum "${OUTPUT}" )) 
    echo "${MD5} $( basename ${OUTPUT} )" > "${OUTPUT}.md5sum"
  else
    echo  "The file: \"${OUTPUT}\" does not exist!"
  fi
}

function Print_ControlSum()
{
  echo  "The file: \"${OUTPUT}\" was successfully created. Its control md5sum is: \"${MD5}\"."
}

# We check if argument is valid and call a specific function related to it.
function Parse_Arguments()
{
  for i in ${@}; do
  case ${i} in
    -d=*|--database=*)
      DB_FILE="${i#*=}"
      shift # past argument
    ;;
    -n=*|--nodes=*)
      NO_NODES="${i#*=}"
      shift # past argument
    ;;
    -o=*|--output=*)
      OUTPUT="${i#*=}"
      shift # past argument
    ;;
    -m=*|--maple=*)
      MAPLE_FILE="${i#*=}"
      shift # past argument
    ;;
    *)
    ;;
  esac
  done
  if [ -z "${DB_FILE}" ]; then
    echo "You have to use the parameter -d=<file.db> to provide the path to the database file!"
    exit 1
  fi
  if [ -z "${NO_NODES}" ]; then
    echo "You have to use the parameter -n=<integer> to provide the number of nodes in the Sun Grid Engine!"
    exit 1
  fi
  if [ -z "${OUTPUT}" ]; then
    echo "You have to use the parameter -o=<file.shx> to provide the name of the output archive!"
    exit 1
  fi
  if [ -e "${OUTPUT}" ]; then
    echo "The file: \"${OUTPUT}\" already exists!"
    exit 1
  fi
  if [ -z "${MAPLE_FILE}" ]; then
    echo "You have to use the parameter -m=<file.mpl> to provide the path to the maple code! This code will be run on each node."
    exit 1
  fi
  if [ ! -e "${MAPLE_FILE}" ]; then
    echo "The file: \"${MAPLE_FILE}\" does not exist!"
    exit 1
  fi
}

# Check if some input parameters were passed.
case ${@} in
  (*[![:blank:]]*)
     Parse_Arguments ${@}
     Create_Dirs
     Split_DB
     Prepare_Files
     Print_ControlSum
     ;;
  (*)
    echo "You have to use the parameter -d=<file.db> to provide the path to the database file!"
    echo "You have to use the parameter -n=<integer> to provide the number of nodes in the Sun Grid Engine!"
    echo "You have to use the parameter -o=<file.shx> to provide the name of the output archive!"
    echo "You have to use the parameter -m=<file.mpl> to provide the path to the maple code! This code will be run on each node."
    exit
    ;;
esac




