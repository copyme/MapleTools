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

BUILD_DIR="./DB"

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
  EVENTS_PER_NODE=`echo -e "SELECT ROUND(COUNT(*)/${NO_NODES}) FROM RealAlgebraicNumber;" | sqlite3 "${DB_FILE}"`
  # Convert to integer
  EVENTS_PER_NODE=`printf '%.*f' 0 ${EVENTS_PER_NODE}`

  # Set the directory in which we create splits of the database.
  if [ ! -d "${BUILD_DIR}" ]; then
    mkdir  "${BUILD_DIR}"
    cd "${BUILD_DIR}"
  else
    cd  "${BUILD_DIR}"
    rm -f *
  fi

  for i in $( seq 0 $(( NO_NODES - 1 )) ); do
    # Recreate schemes
    sqlite3 "${DB_FILE}" '.schema SamplePoint' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema ComputedNumbers' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema Events' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema Quadric' | sqlite3 "${i}.db"
    sqlite3 "${DB_FILE}" '.schema RealAlgebraicNumber' | sqlite3 "${i}.db"
    # Copy data
    sqlite3 "${DB_FILE}" "SELECT * FROM Quadric;" | sqlite3 "${i}.db" ".import /dev/stdin Quadric"
    sqlite3 "${DB_FILE}" "SELECT * FROM RealAlgebraicNumber WHERE ID BETWEEN $(( i * EVENTS_PER_NODE
    )) AND $(( ( i + 1 ) * EVENTS_PER_NODE ));" | sqlite3 "${i}.db" ".import /dev/stdin RealAlgebraicNumber"
    sqlite3 "${DB_FILE}" "SELECT * FROM Events WHERE RANumID BETWEEN $(( i * EVENTS_PER_NODE
    )) AND $(( ( i + 1 ) * EVENTS_PER_NODE ));" | sqlite3 "${i}.db" ".import /dev/stdin Events"
  done
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
    *)
    ;;
  esac
  done
}

# Check if some input parameters were passed.
case ${@} in
  (*[![:blank:]]*)
     Parse_Arguments ${@}
     Split_DB
     ;;
  (*)
     echo 'TODO'
     exit
    ;;
esac




