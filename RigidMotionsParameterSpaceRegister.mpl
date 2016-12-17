module SamplePointsWriter()
  option object;

  local midpoint::rational;
  local path::string;
  local prefix::string;
  export ModuleCopy::static := proc( self::SamplePointsWriter,
                                     proto::SamplePointsWriter,
                                     midpoint::rational,
                                     path::string,
                                     prefix::string, $ )
    if _passed = 2 then
      self:-midpoint := proto:-midpoint;
      self:-path := proto:-path;
      self:-prefix := proto:-prefix;
    else
      self:-midpoint := midpoint;
      self:-path := path;
      self:-prefix := prefix;
    fi;

    return self;
  end proc;

  export ModulePrint::static := proc( self::SamplePointsWriter )
    nprintf("(midpoint: '%a', path: '%a', prefix: '%a')", self:-midpoint, self:-path, self:-prefix);
  end proc:

  export Write::static := proc( self::SamplePointsWriter, samplePoints, id )
    local fileID, midpoint := self:-midpoint;
    if nops(samplePoints) <> 0 then
      fileID := fopen(sprintf("%s/%s%a.tsv", self:-path, self:-prefix, id), APPEND, TEXT):
      writedata(fileID, samplePoints, string, proc (f, x) fprintf(f, "%a\t%a\t%a", midpoint, op(x)) end proc):
      fclose(fileID):
    fi:
  end proc:

end module;
