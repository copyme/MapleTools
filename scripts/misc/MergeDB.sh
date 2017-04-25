#!/bin/bash 
  
if [ -z "${1}" ]; then
  echo "You have to provide a name for the output database! Try: ${0} <filename.db>"
  exit 1
fi

if [ ! -e "${1}" ]; then
  echo "CREATE TABLE NMM (ID INTEGER PRIMARY KEY AUTOINCREMENT, NMM TEXT NOT NULL UNIQUE CHECK (length(NMM) > 0), A TEXT NOT NULL CHECK (length(A) > 0), B TEXT NOT NULL CHECK (length(B) > 0), C TEXT NOT NULL CHECK (length(C) > 0), T1 TEXT NOT NULL CHECK (length(T1) > 0), T2 TEXT NOT NULL CHECK (length(T2) > 0), T3 TEXT NOT NULL CHECK (length(T3) > 0));" | sqlite3 "${1}" > /dev/null
fi

for f in ./*.db; do
  VERSION=`echo -e "PRAGMA user_version;" | sqlite3 "${f}"`
  if [ $VERSION = 3 ]; then
    f=$(echo "'$f'")
    sqlite3 $1 << EOF
    attach $f as t0; 
    INSERT OR IGNORE INTO NMM (NMM, A, B, C, T1, T2, T3) SELECT  NMM, A, B, C, T1, T2, T3 FROM t0.NMM, t0.SamplePoint WHERE t0.NMM.SP_ID = t0.SamplePoint.ID;
EOF
  fi
done
