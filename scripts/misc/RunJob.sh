#!/bin/sh


# Get rid of national specifications like "," instead of "." for numbers.
export LC_ALL=C
MAPLE_FILE="${1}"

for f in *.db; do
  TMP_MPL=$(mktemp -u "/tmp/MapleScript.XXXXXX.mpl")
  cp "${MAPLE_FILE}" "${TMP_MPL}"
  sed -i "s~__REPLACE__~${f}~g" "${TMP_MPL}"
  maple "${TMP_MPL}" >> ${f}.log
done

exit 0
