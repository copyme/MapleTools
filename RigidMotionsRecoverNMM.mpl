# Procedure: GetOrderedCriticalPlanes
#   Compute critical planes in the remainder range
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#   samplePoints     - values of a, b, c (see CayleyTransform)
#
# Output:
#   Returns signature of order and ordered critical planes 
#   in the remainder range for X, Y and Z directions.
GetOrderedCriticalPlanes := proc(vars::list, samplePoints::list, planes::list) 
  local params;
  local sdPlanes := [[],[],[]];
  local xSig, ySig, zSig, Signature;
  
  if nops(vars) <> 3 then
    error "Only 3D arrangement is supported.";
  fi;
  
  sdPlanes[1] := eval(sdPlanes[1], {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  sdPlanes[2] := eval(sdPlanes[2], {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  sdPlanes[3] := eval(sdPlanes[3], {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  
  xSig, sdPlanes[1] := Isort(sdPlanes[1]); 
  ySig, sdPlanes[2] := Isort(sdPlanes[2]); 
  zSig, sdPlanes[3] := Isort(sdPlanes[3]);
  
  # remove all out of the range regions 
  sdPlanes[1] := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, sdPlanes[1]));
  sdPlanes[2] := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, sdPlanes[2]));
  sdPlanes[3] := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, sdPlanes[3]));
  # add region boarders
  sdPlanes[1] := [-1/2,op(sdPlanes[1]),1/2];
  sdPlanes[2] := [-1/2,op(sdPlanes[2]),1/2];
  sdPlanes[3] := [-1/2,op(sdPlanes[3]),1/2]; 

  Signature := cat(op(map(proc (x) sprintf("%d", x) end proc, [op(xSig), op(ySig), op(zSig)])));
  return Signature, sdPlanes;
end proc:


# Procedure: RecoverTranslationSamplePoints
#   Compute midpoints of each frame in the remainder range
#
# Parameters:
#   xPlanes            - ordered X critical planes in the remainder range
#   yPlanes            - ordered Y critical planes in the remainder range
#   zPlanes            - ordered Z critical planes in the remainder range
#
# Output:
#   Returns centers of frames in the remainder range
RecoverTranslationSamplePoints := proc(planes::list) 
  local i, j, k; 
  local samples := []; 
  for i to upperbound(planes[1]) - 1 do
    for j to upperbound(planes[2]) - 1 do 
      for k to upperbound(planes[3]) - 1 do 
        samples := [op(samples), [(1/2)*add(planes[1][i .. i+1]), 
                    (1/2)*add(planes[2][j.. j+1]), (1/2)*add(planes[3][k .. k+1])]];
      end do;
    end do;
  end do; 
  return samples; 
end proc:


# Procedure: Get3DNMM
#   Compute neighbourhood motion maps
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   samplesTrans     - midpoints of frames in the remainder range
#   sampleRot        - values of a, b, c (see CayleyTransform)
#
# Output:
#   Returns 3D neighbourhood motion maps for given a,b,c and corresponding translations.
Get3DNMM := proc(nType::string, samplesTrans::list, sampleRot::list, NMMContainer::uneval) 
  local x;
  local R := eval(CayleyTransform({a,b,c}), {a = sampleRot[1], b = sampleRot[2], c = sampleRot[3]});
  local n := GetNeighborhood(nType);
  local NMM := eval(NMMContainer); 
  NMMContainer := NMM union {seq(map(proc (y) map(round, convert(R.Vector(3, y)+Vector(3,x), list))
  end proc, n),x in samplesTrans)};
end proc;

CriticalPlanes := proc(vars::list, nType::string, kRange::list)
  local R := CayleyTransform(vars);
  # we remove 7th element -- [0, 0, 0]
  local n := subsop(7=NULL, GetNeighborhood(nType)); 
  local T := combinat:-cartprod([n, kRange]);
  local planes := [[],[],[]]; 
  local params;
  while not T[finished] do 
    params := T[nextvalue](); 
    planes[1] := [op(planes[1]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[1]]; 
    planes[2] := [op(planes[2]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[2]]; 
    planes[3] := [op(planes[3]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[3]];
  end do;

 return planes;

end proc;


# Procedure: CalculateNMM
#   Reads data from hard drive and generates NMM
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#   dir              - a directory which contains sample points
#   id               - an id of a file which contains sample points
#
# Output:
#   Set of NMM
#
# TODO:
#   Replace hard coded path by variable
#   Export printing code into a procedure
#   Simplify the code
#   After adding the communication approach the number of NMM is much smaller, check this
CalculateNMM := proc(vars::list, nType::string, kRange::list, dir::string, id::integer) 
  local csvFile, data, s, fileID;
  local i, sdPlanes := [], trans, NMM:={}, rotSamp; 
  local percent := -1, Signatures := {}, sig, msg;
  local n := Grid:-NumNodes();
  local planes := CriticalPlanes(vars, nTypes, kRange);

# read sample points from a file and convert them to appropriate format
  csvFile := FileTools:-JoinPath([sprintf(cat(dir,"/sam_%d.csv"), id)], base = homedir);
  data := ImportMatrix(csvFile, 'source' = 'tsv', 'datatype' ='string');
  data := Threads:-Map(proc (x) op(sscanf(x, %a)) end proc, data); 

  # print a progress message
  for i to upperbound(data)[1] do
    if percent <>  round( i / upperbound(data)[1] * 100 ) then
      percent := round( i / upperbound(data)[1] * 100 );
      if percent mod 10 = 0 then
        printf("[%s]:> Node: %d, finished: %d%%.\n",
               StringTools:-FormatTime("%Y-%m-%d -- %R"),id,percent); 
      end if;
    end if;
    rotSamp := convert(data[i], list); 
    # discard the sample points for which NMM are already computed
    sig, sdPlanes := GetOrderedCriticalPlanes(nType, kRange, rotSamp, planes); 
    if member( sig, Signatures ) then
      next;
    end if;

    trans := RecoverTranslationSamplePoints(sdPlanes); 
    Get3DNMM(nType, trans, rotSamp, NMM);

    Signatures := Signatures union {sig};
  end do;

  # write the result
  fileID := fopen(sprintf("/home/plutak/debugNMM2/NMM_%d.tsv", id), WRITE, TEXT);
  writedata(fileID, [op(NMM)], string, proc (f, x) fprintf(f, "%a", x) end proc);
  fclose(fileID);
end proc:

# Procedure: ParallelCalculateNMM
#   Uses Grid framework to reads data from hard drive and generates NMM
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#   dir              - a directory which contains sample points
#
# Output:
#   Set of NMM
ParallelCalculateNMM := proc(vars::list, nType::string, kRange::list, dir::string) 
  local me;
  me := Grid:-MyNode(); 
  CalculateNMM(vars, nType, kRange, dir, me);
  Grid:-Barrier();
end proc:

# Procedure: LaunchOnGridGetNMM
#   Setup and run computation on a local grid
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#   dir              - a directory which contains sample points
#
# Output:
#   Write a set of NMM to a file given by fileName.
LaunchOnGridGetNMM := proc (nType::string, kRange::list, dir::string, nodes:=20) 
  Grid:-Setup("local", numnodes=nodes); 
  Grid:-Launch(ParallelCalculateNMM, nType, kRange, dir):
end proc:

