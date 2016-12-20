module ComputationRegister()
  option object;

  local connection;
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
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.RealAlgebraicNumber (ID INTEGER PRIMARY" ||    
      " KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, IntervalR TEXT NOT NULL);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Quadric (ID PRIMARY KEY NOT NULL UNIQUE, "
      || "polynom TEXT NOT NULL UNIQUE);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.Events (RANumID INTEGER " ||
      "REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL, " ||
      "QuadID INTEGER REFERENCES Quadric (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE cacheDB.SamplePoint (ID INTEGER PRIMARY KEY" ||
      " AUTOINCREMENT, A DOUBLE NOT NULL, B DOUBLE NOT NULL, C DOUBLE NOT NULL);");
      fi;

    return self;
  end proc;

  export ModulePrint::static := proc( self::ComputationRegister )
    print(self:-connection);
  end proc;


  export InsertQuadric::static := proc(self::ComputationRegister, id::integer, quadric::polynom)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO cacheDB.Quadric(ID, polynom) "
                                            || "VALUES (?, ?);");
    Database[SQLite]:-Bind(stmt, 1, id);
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", quadric));
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

  export InsertEvent::static := proc(self::ComputationRegister, idNum::integer,
                                    num::RealAlgebraicNumber, quadrics::list)
    local x::integer;
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO cacheDB.RealAlgebraicNumber(polynom,"
                                             || " IntervalL, IntervalR) VALUES (?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, sprintf("%a", GetPolynomial(num)));
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", GetInterval(num)[1]));
    Database[SQLite]:-Bind(stmt, 3, sprintf("%a", GetInterval(num)[2]));
    Database[SQLite]:-Step(stmt);
    Database[SQLite]:-Finalize(stmt);
    #insert quadrics for event
    for x in quadrics do
      stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO cacheDB.Events(RANumID, QuadID) " ||
                                     "VALUES(?, ?);");
      Database[SQLite]:-Bind(stmt, 1, idNum);
      Database[SQLite]:-Bind(stmt, 2, x);
      while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
      Database[SQLite]:-Finalize(stmt);
    od;
  end proc;

  export SynchronizeAlgebraicNumbers::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO RealAlgebraicNumber SELECT * FROM " ||
                                "cacheDB.RealAlgebraicNumber WHERE NOT EXISTS(SELECT 1 FROM " ||
                                "RealAlgebraicNumber AS R, cacheDB.RealAlgebraicNumber AS rc WHERE r.ID = " ||
                                "rc.ID AND r.POLYNOM = rc.POLYNOM AND r.INTERVALL = " ||
                                "rc.INTERVALL AND r.INTERVALR = rc.INTERVALR);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
    stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO Quadric SELECT * FROM cacheDB.Quadric" ||
                                                  " WHERE NOT EXISTS(SELECT 1 FROM Quadric AS Q, " ||
                                                  "cacheDB.Quadric AS qc WHERE q.ID = qc.ID AND q.POLYNOM = " ||
                                                  "qc.POLYNOM);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
    stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO Events SELECT * FROM cacheDB.Events " ||
                                                  "WHERE NOT EXISTS( SELECT 1 FROM Events AS E, " ||
                                                  "cacheDB.Events AS ev WHERE e.RANUMID = ev.RANUMID AND " ||
                                                  "e.QUADID = ev.QUADID);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

  export InsertSkippedCluster::static := proc(self::ComputationRegister, id::integer)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO SkippedCluster " || 
                                          "(clusterID) VALUES (?);");
    Database[SQLite]:-Bind(stmt, 1, id);    
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

  export InsertSamplePoint::static := proc(self::ComputationRegister, samp::list)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT OR IGNORE INTO cacheDB.SamplePoint " || 
                                          "(A, B, C) VALUES (?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, samp[1]);    
    Database[SQLite]:-Bind(stmt, 2, samp[2]);    
    Database[SQLite]:-Bind(stmt, 3, samp[3]);    
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

  export SynchronizeSamplePoints::static := proc(self::ComputationRegister)
    local stmt := Database[SQLite]:-Prepare(self:-connection, "INSERT INTO SamplePoint SELECT * FROM " ||
                                                "cacheDB.SamplePoint " ||
                                                "WHERE NOT EXISTS(SELECT 1 FROM SamplePoint AS s, " ||
                                                "cacheDB.SamplePoint AS sc WHERE s.A = sc.A AND " ||
                                                "s.B = sc.B AND s.C = sc.C);");
    while Database[SQLite]:-Step(stmt) <> Database[SQLite]:-RESULT_DONE do; od;
    Database[SQLite]:-Finalize(stmt);
  end proc;

end module;
