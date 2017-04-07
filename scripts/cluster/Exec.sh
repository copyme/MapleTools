#!/bin/sh

# The script inspired by http://www.linuxjournal.com/node/1005818

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

# This variable is set by the Build.sh script.
SHARED_DIR=${1}
INSTALL_SCRIPT="Install.sh"
MAPLE_FILE=__REPLACE_MPL__
NODE_RUNNER_SCRIPT="NodeRunner.sh"


if [ -z "${SHARED_DIR}" ]; then
    echo "You have to provide a path to the shared directory! Please, type: ${0}
    </path/to/shared/directory> <optional arguments to qsub>"
    exit 1
fi

if [ ! -e "${SHARED_DIR}" ]; then
    echo "Provided shared directory: \"${SHARED_DIR}\" does not exist!"
    exit 1
fi

if ! hash qsub 2>/dev/null; then
  echo "This script relay on Sun Engine Grid tools! You need to install them."
  exit 1
fi

export TMP_DIR=`mktemp -d ${SHARED_DIR}/selfextract.XXXXXX`

ARCHIVE=`awk '/^__ARCHIVE_BELOW__/ {print NR + 1; exit 0; }' "${0}"`

tail -n+"${ARCHIVE}" "${0}" | tar xzv -C "${TMP_DIR}"

cd ${TMP_DIR}

# Uninstall scripts -- if installed before
./${INSTALL_SCRIPT} -u

# Install Maple scripts
./${INSTALL_SCRIPT} -d="${TMP_DIR}"

for f in ./DB/*.db; do
  qsub ${2} "${TMP_DIR}/${NODE_RUNNER_SCRIPT}" "${TMP_DIR}/DB/$( basename ${f} )" "${TMP_DIR}/${MAPLE_FILE}"
done

exit 0

# Do not add anything below the next line!!!!
__ARCHIVE_BELOW__
