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
      fileStatus := FileTools:-Exists(dbPath);
      self:-connection := Database[SQLite]:-Open(":memory:");
      Database[SQLite]:-Attach(self:-connection, dbPath, "outDB");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA synchronous = OFF");
      Database[SQLite]:-Execute(self:-connection, "PRAGMA journal_mode = MEMORY");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE RealAlgebraicNumber (ID INTEGER PRIMARY" ||    
      " KEY UNIQUE, polynom TEXT NOT NULL, IntervalL TEXT NOT NULL, IntervalR TEXT NOT NULL);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE Quadric (ID PRIMARY KEY NOT NULL UNIQUE, "
      || "polynom TEXT NOT NULL UNIQUE);");
      Database[SQLite]:-Execute(self:-connection, "CREATE TABLE Events (RANumID INTEGER " ||
      "REFERENCES RealAlgebraicNumber (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL, " ||
      "QuadID INTEGER REFERENCES Quadric (ID) ON DELETE CASCADE ON UPDATE CASCADE NOT NULL);");
      fi;

    return self;
  end proc;

  export ModulePrint::static := proc( self::ComputationRegister )
    print(self:-connection);
  end proc;


  export InsertQuadric::static := proc(self::ComputationRegister, id::integer, quadric::polynom)
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO Quadric(ID, polynom) "
                                            || "VALUES (?, ?);");
    Database[SQLite]:-Bind(stmt, 1, id);
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", quadric));
    Database[SQLite]:-Step(stmt);
    Database[SQLite]:-Finalize(stmt);
  end proc;

  export InsertEvent::static := proc(self::ComputationRegister, idNum::integer,
                                    num::RealAlgebraicNumber, quadrics::list)
    local x::integer;
    local stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO RealAlgebraicNumber(polynom,"
                                             || " IntervalL, IntervalR) VALUES (?, ?, ?);");
    Database[SQLite]:-Bind(stmt, 1, sprintf("%a", GetPolynomial(num)));
    Database[SQLite]:-Bind(stmt, 2, sprintf("%a", GetInterval(num)[1]));
    Database[SQLite]:-Bind(stmt, 3, sprintf("%a", GetInterval(num)[2]));
    Database[SQLite]:-Step(stmt);
    Database[SQLite]:-Finalize(stmt);
    #insert quadrics for event
    for x in quadrics do
      stmt := Database[SQLite]:-Prepare(self:-connection,"INSERT INTO Events(RANumID, QuadID) " ||
                                     "VALUES(?, ?);");
      Database[SQLite]:-Bind(stmt, 1, idNum);
      Database[SQLite]:-Bind(stmt, 2, x);
      Database[SQLite]:-Step(stmt);
      Database[SQLite]:-Finalize(stmt);
    od;
  end proc;

  export Synchronize::static := proc(self::ComputationRegister)
      Database[SQLite]:-Execute(self:-connection, "INSERT INTO outDB.RealAlgebraicNumber SELECT * FROM " ||
                                "RealAlgebraicNumber;");
      Database[SQLite]:-Execute(self:-connection, "INSERT INTO outDB.Quadric SELECT * FROM Quadric;");
      Database[SQLite]:-Execute(self:-connection, "INSERT INTO outDB.Events SELECT * FROM Events;");
  end proc;

end module;
