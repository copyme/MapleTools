# Class: ComputationRegister
#
# Description:
#  This class implements a register of the computations (see
#  RigidMotionsParameterSpaceDecompostion) which is based on SQLite.
#  Thanks to the register the computations can be restored after a crash.
#
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  11/12/2016 
#
# License:
#  Simplified BSD License
#
# Copyright (c) 2016, Kacper Pluta
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Kacper Pluta BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
module ComputationRegister()
  option object;
  
  (* Connection object returned by Database[SQLite]:-Open *)
  local connection;

# Method: ModuleCopy
#   Standard constructor / copy constructor
#
# Parameters:
#   self::ComputationRegister      - a new object to be constructed
#   proto::ComputationRegister     - a prototype object from which self is derived
#   dbPath::string                 - a path to a database which is a copy of the file 
#                                    CompRegister.db
#
# Comment:
#   The database is open in such a way that journal_mode is set to WAL and
#   synchronous to NORMAL. These option should not be changed! 
#   See SQLite documentation for more information. Note that second database
#   is created in memory to speed up inserts. In particular, small chunks of
#   data are first inserted into the memory database and then moved in one
#   transaction into the main database.
#
# Output:
#   An object of type ComputationRegister.
#
# Exceptions:
#  "There is no database" occurs when a given file (dbPath) does not exists.
# 
  export ModuleCopy::static := proc( self::ComputationRegister,
                                     proto::ComputationRegister,
                                     dbPath::string, $ )
  local fileStatus := false;
    if _passed = 2 then
      self:-connection := proto:-connection;
    else
      if not FileTools:-Exists(dbPath) then
        error "There is no database %1.", dbPath;
      fi;
      self:-connection := Database[SQLite]:-Open(dbPath, create=false);
      Database[SQLite]:-Attach(self:-connection, ":memory:", "cacheDB");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA synchronous = NORMAL");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA journal_mode = WAL");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.RealAlgebraicNumber " ||
      "(ID INTEGER PRIMARY KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, " ||
      "IntervalR TEXT NOT NULL);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Quadric (ID PRIMARY KEY " ||
      "NOT NULL UNIQUE, polynom TEXT NOT NULL UNIQUE);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Events (RANumID INTEGER " ||
      "REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL, " ||
      "QuadID INTEGER REFERENCES Quadric (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.SamplePoint (" ||
      "A TEXT NOT NULL, B TEXT NOT NULL, C TEXT NOT NULL);");
      fi;

    return self;
  end proc;


# Method: ModulePrint
#   Standard printout of an object of type ComputationRegister.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
  export ModulePrint::static := proc( self::ComputationRegister )
    print(self:-connection);
  end proc;

# Method: InsertQuadric
#   Used to insert quadrics into a database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   id::integer                    - an identifier of a given quadric
#   quadric::polynom               - a second degree polynomial
#
# Comments:
#   Each quadric is inserted into a cache, memory stored, database.
#   SynchronizeQuadrics has to be called to move inserted polynomials
#   into the register.
#
  export InsertQuadric::static := proc(self::ComputationRegister, id::integer, quadric::polynom)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO cacheDB.Quadric(ID, polynom) "
                                            || "VALUES (?, ?);");
    Database[SQLite]:-Bind(stmt, 1, id);
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", quadric));
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: InsertEvent
#   Used to insert quadrics into a database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   idNum::integer                 - an identifier of a given event
#   num::RealAlgebraicNumber       - a real algebraic number representing
#                                    a given event
#   quadrics::list                 - a list of integer ids' of quadrics
#                                    each of which has to be the same id
#                                    as ones used with InsertQuadric
#
# Comments:
#   Each event is inserted into a cache, memory stored, database.
#   SynchronizeEvents has to be called to move inserted events
#   into the register.
#
  export InsertEvent::static := proc(self::ComputationRegister, idNum::integer,
                                     num::RealAlgebraicNumber, quadrics::list)
    local x::integer;
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                                             "cacheDB.RealAlgebraicNumber(polynom, " || 
                                             "IntervalL, IntervalR) VALUES (?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, sprintf("%a", GetPolynomial(num)));
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", GetInterval(num)[1]));
    Database[SQLite]:-Bind(stmt, 3, sprintf("%a", GetInterval(num)[2]));
    Database[SQLite]:-Step(stmt);
    Database[SQLite]:-Finalize(stmt);
    for x in quadrics do
      stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO cacheDB.Events(RANumID, " ||
                                                         "QuadID) VALUES(?, ?);");
      Database[SQLite]:-Bind(stmt, 1, idNum);
      Database[SQLite]:-Bind(stmt, 2, x);
      while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
      Database[SQLite]:-Finalize(stmt);
    od;
  end proc;


# Method: SynchronizeQuadrics
#   Synchronize quadrics between memory, cache, database and a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
  export SynchronizeQuadrics::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO Quadric SELECT * FROM " ||
                                                  "cacheDB.Quadric WHERE NOT EXISTS(SELECT 1 FROM " 
                                                  || "Quadric AS Q, cacheDB.Quadric AS qc WHERE " ||
                                                  "q.ID = qc.ID AND q.POLYNOM = qc.POLYNOM);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: SynchronizeEvents
#   Synchronize events between memory, cache, database and a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
  export SynchronizeEvents::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO RealAlgebraicNumber " ||
                                "SELECT * FROM cacheDB.RealAlgebraicNumber WHERE NOT " ||
                                "EXISTS(SELECT 1 FROM RealAlgebraicNumber AS R, " ||
                                "cacheDB.RealAlgebraicNumber AS rc WHERE r.ID = rc.ID AND " ||
                                "r.POLYNOM = rc.POLYNOM AND r.INTERVALL = rc.INTERVALL AND " ||
                                "r.INTERVALR = rc.INTERVALR);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
    
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO Events SELECT * FROM " ||
                           "cacheDB.Events WHERE NOT EXISTS( SELECT 1 FROM Events AS E, " ||
                           "cacheDB.Events AS ev WHERE e.RANUMID = ev.RANUMID AND e.QUADID = " ||
                           "ev.QUADID);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: InsertSkippedCluster
#    Inserts ids' of cluster into a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   id::integer                    - an integer which stands for the identifier of a cluster
#
#
  export InsertSkippedCluster::static := proc(self::ComputationRegister, id::integer)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO SkippedCluster " || 
                                          "(clusterID) VALUES (?);");
    Database[SQLite]:-Bind(stmt, 1, id);    
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: InsertSamplePoint
#    Inserts ids' of cluster into a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   samp::list                     - a list of three rational numbers which represent a sample point
#
#
  export InsertSamplePoint::static := proc(self::ComputationRegister, samp::list)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                                 "cacheDB.SamplePoint (A, B, C) VALUES (?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, sprintf("%a", samp[1]));    
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", samp[2]));    
    Database[SQLite]:-Bind(stmt, 3, sprintf("%a", samp[3]));    
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: SynchronizeSamplePoints
#    Synchronize sample points between cache, memory, database and a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
  export SynchronizeSamplePoints::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO SamplePoint SELECT * " ||
                                 "FROM cacheDB.SamplePoint WHERE NOT EXISTS(SELECT 1 FROM " ||
                                 "SamplePoint AS s, cacheDB.SamplePoint AS sc WHERE s.A = sc.A " ||
                                 "AND s.B = sc.B AND s.C = sc.C);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

# Method: FetchSkippedClusters
#    Fetch skipped clusters from the database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
# Output:
#   A list of integers.
#
  export FetchSkippedClusters::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT * FROM SkippedCluster;");
    local result := convert(Database[SQLite]:-FetchAll(stmt), list);
    Database[SQLite]:-Finalize(stmt);
    return result;
  end proc;


# Method: FetchQuadrics
#    Fetch quadrics from the database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
# Output:
#   A list of second degree polynomials.
#
  export FetchQuadrics::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT polynom FROM Quadric " ||
                                            "ORDER BY ID;");
    local result := map(parse, convert(Database[SQLite]:-FetchAll(stmt), list));
    Database[SQLite]:-Finalize(stmt);
    return result;
  end proc;


# Method: FetchEvents
#    Fetch events from the database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
# Output:
#   A list of pairs: [RealAlgebraicNumber, [quadrics' ids]].
#
  export FetchEvents::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT ID, polynom, IntervalL, " ||
    "IntervalR, group_concat(QuadID) AS quads  FROM RealAlgebraicNumber JOIN EVENTS ON " ||
    "ID = RANUMID GROUP BY ID;");
    local allRows, s;
    allRows := Database[SQLite]:-FetchAll(stmt);
    Database[SQLite]:-Finalize(stmt);
    s := proc(i::integer, allRows)
      local stmp, rowAlg, quads;
      rowAlg := allRows[i];
      return [Object(RealAlgebraicNumber, parse(rowAlg[2]), parse(rowAlg[3]), parse(rowAlg[4])),
                     [parse(rowAlg[5])]];
    end proc;
    return [Threads:-Seq(s(i, allRows),i=1..upperbound(allRows)[1])];
  end proc;

end module;
