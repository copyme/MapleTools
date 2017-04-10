#!/bin/sh

# The script is inspired by http://www.linuxjournal.com/node/1005818

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

# This variable is set by the Build.sh script.
INSTALL_SCRIPT="Install.sh"
MAPLE_FILE=__REPLACE_MPL__
NODE_RUNNER_SCRIPT="NodeRunner.sh"

function Init()
{
  if ! hash qsub 2>/dev/null; then
    echo "This script relay on Sun Engine Grid tools! You need to install them."
    exit 1
  fi 
}

# We check if argument is valid and call a specific function related to it.
function Parse_Arguments()
{
  for i in ${@}; do
  case ${i} in
    -d=*|--dir=*)
      SHARED_DIR="${i#*=}"
      shift # past argument
    ;;
    -b=*|--begin=*)
      RANGE_BEGIN="${i#*=}"
      shift # past argument
    ;;
    -e=*|--end=*)
      RANGE_END="${i#*=}"
      shift # past argument
    ;;
    *)
    ;;
  esac
  done
  if [ -z "${SHARED_DIR}" ]; then
    echo "You have to provide a path to the shared directory! Please, type: ${0} -d=</path/to/shared/directory> -b=<an optional begin of the range> -e=<an optional end of the range> <optional arguments to qsub>"
    exit 1
  fi
  if [ ! -e "${SHARED_DIR}" ]; then
    echo "Provided shared directory: \"${SHARED_DIR}\" does not exist!"
    exit 1
  fi
  if [ -z "${RANGE_BEGIN}" ]; then
    RANGE_BEGIN=0
  fi
  if [ -z "${RANGE_END}" ]; then
    RANGE_END=`ls ${TMP_DIR}/DB/*.db`
  fi
  QSUB_ARGS=${@}
}

# Check if some input parameters were passed.
case ${@} in
  (*[![:blank:]]*)
     Init
     Parse_Arguments ${@}
     ;;
  (*)
    echo "Please, type: ${0} -d=</path/to/shared/directory> -b=<an optional begin of the range> -e=<an optional end of the range> <optional arguments to qsub>"
    exit
    ;;
esac

export TMP_DIR=`mktemp -d ${SHARED_DIR}/selfextract.XXXXXX`

ARCHIVE=`awk '/^__ARCHIVE_BELOW__/ {print NR + 1; exit 0; }' "${0}"`

tail -n+"${ARCHIVE}" "${0}" | tar xzv -C "${TMP_DIR}"

cd ${TMP_DIR}

# Uninstall scripts -- if installed before
./${INSTALL_SCRIPT} -u

# Install Maple scripts
./${INSTALL_SCRIPT} -d="${TMP_DIR}"

for i in `seq ${RANGE_BEGIN} ${RANGE_END}`; do
  qsub "${QSUB_ARGS}" "${TMP_DIR}/${NODE_RUNNER_SCRIPT}" "${TMP_DIR}/DB/${i}.db" "${TMP_DIR}/${MAPLE_FILE}"
done

# Prevent calling the rest of the script
exit 0

# Do not add anything below the next line!!!!
__ARCHIVE_BELOW__
