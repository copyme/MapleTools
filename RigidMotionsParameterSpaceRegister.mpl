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
  local version;

# Method: ModuleCopy
#   Standard constructor / copy constructor
#
# Parameters:
#   self::ComputationRegister      - a new object to be constructed
#   proto::ComputationRegister     - a prototype object from which self is derived
#   dbPath::string                 - a path to a database
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
  local fileStatus := false, version, stmt;
    if _passed = 2 then
      self:-connection := proto:-connection;
      self:-version := proto:-version;
    else
      fileStatus:=FileTools:-Exists(dbPath);
      self:-connection := Database[SQLite]:-Open(dbPath);
      Database[SQLite]:-Attach(self:-connection, ":memory:", "cacheDB");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA synchronous = OFF;");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA journal_mode = WAL;");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA cacheDB.auto_vacuum = FULL;");

      #Create tables
      if not fileStatus then
        Database[SQLite]:-Execute(self:-connection,"CREATE TABLE Quadric (ID PRIMARY KEY NOT " ||
        "NULL UNIQUE, polynom TEXT NOT NULL UNIQUE);");
        Database[SQLite]:-Execute(self:-connection,"CREATE TABLE RealAlgebraicNumber (ID " ||
        "INTEGER PRIMARY KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, IntervalR " ||
        "TEXT NOT NULL);");
        Database[SQLite]:-Execute(self:-connection,"CREATE TABLE Events (RANumID INTEGER " ||
        "REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL, " ||
        "QuadID INTEGER REFERENCES Quadric (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL);");
        Database[SQLite]:-Execute(self:-connection,"CREATE TABLE SamplePoint (A TEXT NOT NULL, " ||
        "B TEXT NOT NULL, C TEXT NOT NULL);");
        Database[SQLite]:-Execute(self:-connection, "CREATE TABLE ComputedNumbers (RANumID " ||
        "INTEGER REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT " ||
        "NULL);");
      fi;

      stmt := Database[SQLite]:-Prepare(self:-connection, "PRAGMA user_version;");
      while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
      self:-version := Database[SQLite]:-Fetch(stmt, 0);
      Database:-Finalize(stmt);

      if self:-version = 0 then
        Database[SQLite]:-Execute(self:-connection,"CREATE TABLE cacheDB.RealAlgebraicNumber (ID " ||
        "INTEGER PRIMARY KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, IntervalR " ||
        "TEXT NOT NULL);");
        Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Quadric (ID PRIMARY KEY " ||
        "NOT NULL UNIQUE, polynom TEXT NOT NULL UNIQUE);");
        Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Events (RANumID INTEGER " ||
        "REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL, " ||
        "QuadID INTEGER REFERENCES Quadric (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL);");
        Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.SamplePoint (" ||
        "A TEXT NOT NULL, B TEXT NOT NULL, C TEXT NOT NULL);");
      else
         Database[SQLite]:-Execute(self:-connection,"CREATE TABLE cacheDB.RealAlgebraicNumber (ID " ||
        "INTEGER PRIMARY KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, IntervalR " ||
        "TEXT NOT NULL);");
      fi;

    return self;
  end proc;


  export Close::static := proc( self::ComputationRegister )
    Database[SQLite]:-Close(self:-connection);
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
    local stmt;
    if self:-version > 0 then
      error "Adding new quadrics is blocked! Re-run computations with a new database.";
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                                 "cacheDB.Quadric(ID, polynom) VALUES (?, ?);");
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
#   event::EventType               - a given event
# Comments:
#   Each event is inserted into a cache, memory stored, database.
#   SynchronizeEvents has to be called to move inserted events
#   into the register.
#
  export InsertEvent::static := proc(self::ComputationRegister, idNum::integer,
                                     event::EventType)
    local x::integer, num, quadrics, stmt;
    if self:-version > 0 then
      error "Adding new events is blocked! Re-run computations with a new database.";
    fi;
    num := GetRealAlgebraicNumber(event); quadrics := GetQuadrics(event);
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                                             "cacheDB.RealAlgebraicNumber(ID, polynom, " || 
                                             "IntervalL, IntervalR) " ||
                                             "VALUES (?, ?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, idNum);
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", GetPolynomial(num)));
    Database[SQLite]:-Bind(stmt, 3, sprintf("%a", GetInterval(num)[1]));
    Database[SQLite]:-Bind(stmt, 4, sprintf("%a", GetInterval(num)[2]));
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
    for x in quadrics do
      stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                             "cacheDB.Events(RANumID, QuadID) VALUES(?, ?);");
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
    local stmt;
    if self:-version > 0 then
      return NULL;
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO Quadric SELECT * FROM " ||
                                                  "cacheDB.Quadric WHERE NOT EXISTS(SELECT 1 FROM " 
                                                  || "Quadric AS Q, cacheDB.Quadric AS qc WHERE " ||
                                                  "q.POLYNOM = qc.POLYNOM);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: SynchronizeEvents
#   Synchronize events between memory, cache, database and a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
# Comments:
#   Cache database is cleared. Therefore, the method should be called only when all events were
#   inserted.
  export SynchronizeEvents::static := proc(self::ComputationRegister)
    local stmt;
    if self:-version > 0 then
      return NULL;
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO RealAlgebraicNumber " ||
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

    #clean up cacheDB
    stmt := Database[SQLite]:-Prepare(self:-connection,"DELETE FROM cacheDB.RealAlgebraicNumber;");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: InsertComputedNumber
#    Inserts ids' of computed events into a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   id::integer                    - an integer which stands for the identifier of a cluster
#
#
  export InsertComputedNumber::static := proc(self::ComputationRegister, id::integer)
    local stmt;
    if self:-version > 0 then
      error "Adding new computed events is blocked! Re-run computations with a new database.";
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
                                 "ComputedNumbers (RANumID) VALUES (?);");
    Database[SQLite]:-Bind(stmt, 1, id);    
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: InsertSamplePoint
#    Inserts a sample point into a given database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#   samp::list                     - a list of three rational numbers which represent a sample point
#
#
  export InsertSamplePoint::static := proc(self::ComputationRegister, samp::list)
    local stmt;
    if self:-version > 0 then
      error "Adding new computed sample points is blocked! Re-run computations with a new database.";
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO " ||
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
# Comments:
#   Cache database is cleared.
  export SynchronizeSamplePoints::static := proc(self::ComputationRegister)
    local stmt;
    if self:-version > 0 then
      return NULL;
    fi;
    stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO SamplePoint SELECT * " ||
                                 "FROM cacheDB.SamplePoint WHERE NOT EXISTS(SELECT 1 FROM " ||
                                 "SamplePoint AS s, cacheDB.SamplePoint AS sc WHERE s.A = sc.A " ||
                                 "AND s.B = sc.B AND s.C = sc.C);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);

    #clean up cacheDB
    stmt := Database[SQLite]:-Prepare(self:-connection,"DELETE FROM cacheDB.SamplePoint;");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;


# Method: FetchComputedNumbers
#    Fetch ids' of events which were processed already from the database.
#
# Parameters:
#   self::ComputationRegister      - an instance of ComputationRegister
#
#
# Output:
#   A list of integers.
#
  export FetchComputedNumbers::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT * FROM ComputedNumbers;");
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
#   first::integer                 - an id of the first event from a given range
#   last::integer                  - an id of the last event from a given range
#
#
# Output:
#   An Array of EventType.
#
  export FetchEvents::static := proc(self::ComputationRegister, first::integer, last::integer)
    local row, events:=Array([]);
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT ID, polynom, IntervalL, " ||
    "IntervalR, group_concat(QuadID) FROM RealAlgebraicNumber JOIN EVENTS " ||
    "ON ID = RANUMID WHERE ID BETWEEN ? AND ? GROUP BY ID ORDER BY ID;"); 
    Database[SQLite]:-Bind(stmt, 1, first);
    Database[SQLite]:-Bind(stmt, 2, last);

    #Slow but Fetching all can kill with memory consumption
    while Database[SQLite]:-Step(stmt) = Database[SQLite]:-RESULT_ROW do
      row := Database[SQLite]:-FetchRow(stmt);
      events(row[1]) := EventType(RealAlgebraicNumber(parse(row[2]), parse(row[3]), 
                    parse(row[4])), [parse(row[5])]);
    od;
    Database[SQLite]:-Finalize(stmt);
    return events;
  end proc;

  export NumberOfEvents::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT MAX(ID) " ||
                                            "FROM RealAlgebraicNumber;"); 
    # Fetch() does not work - "no data left"
    local num::integer := Database[SQLite]:-FetchAll(stmt)[1][1];
    Database[SQLite]:-Finalize(stmt);
    return num;
  end proc;
  

  export PrepareSamplePoints::static := proc(self::ComputationRegister)
    local stmt, toCompute;
    stmt := Database[SQLite]:-Prepare(self:-connection,"SELECT COUNT(ID) FROM RealAlgebraicNumber "
    || "WHERE ID NOT IN (SELECT RANUMID FROM ComputedNumbers) AND ID NOT IN (SELECT MAX(ID) " ||
    "FROM RealAlgebraicNumber);"); 
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    toCompute := Database[SQLite]:-Fetch(stmt, 0);
    Database[SQLite]:-Finalize(stmt);

    if self:-version = 0 and toCompute = 0 then
      Database[SQLite]:-Execute(self:-connection, "ALTER TABLE SamplePoint RENAME TO " || 
      "sqlitestudio_temp_table; CREATE TABLE SamplePoint (A TEXT NOT NULL, B TEXT NOT NULL, " ||
      "C TEXT NOT NULL, ID INTEGER PRIMARY KEY AUTOINCREMENT); INSERT INTO SamplePoint (A, B, C) " ||
      "SELECT A, B, C FROM sqlitestudio_temp_table; DROP TABLE sqlitestudio_temp_table;");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA user_version = 1";);
      self:-version = 1;
    elif toCompute <> 0 then
      Close(self);
      error "Before running computation of NMM it is necessary to compute all sample points!
      Please, run first RigidMotionsParameterSpaceDecompostion:-LaunchResumeComputations().";
    fi;
  end proc;


  export NumberOfSamplePoints::static := proc(self::ComputationRegister)
    local stmt, num::integer;
    stmt := Database[SQLite]:-Prepare(self:-connection, "SELECT COUNT(*) FROM SamplePoint;"); 
    num::integer := Database[SQLite]:-FetchAll(stmt)[1][1];
    Database[SQLite]:-Finalize(stmt);
    return num;
  end proc;
end module;
