# File: RigidMotionsRecoverNMM.mpl 
#
# Description:
#  This file contains functions used to obtain an arrangement 6 dimensional parameter space of 3D
#  digitized rigid motions.
#  This code has been written for research propose and its aim is to calculate a particular
#  arrangement of quadrics. Therefore, it can or it cannot be useful in study of generic
#  arrangements. The final output are sample points of full dimensional open cells.
#
#  The code was written in relation with the paper: Kacper Pluta, Guillaume Moroz, Yukiko
#  Kenmochi, Pascal Romon, Quadric arrangement in classifying rigid motions of a 3D digital image,
#  2016, https://hal.archives-ouvertes.fr/hal-01334257 referred late on as [Quadrics:2016].
#
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  11/12/2015 
#
# License:
#  Simplified BSD License
#
# Copyright (c) 2015, Kacper Pluta
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
#

RigidMotionsRecoverNMM := module() 
  option package;
  uses   RigidMotionsParameterSpaceCommon, RigidMotionsMaplePrimesCode;
  
  (* String which represents a type of the neighborhood i.e. N1, N2 or N3. *)
  global nTypeGlobal;
  (* String which represents the path to the database. *)
  global dbPathGlobal; 
  (* List of integers which are the indices of half-grid planes. *)
  global kRangeGlobal;
  (* List of the variables in which the problem is expressed. *)
  global varsGlobal;

  (*Controls how many sample points should be fetch from the database in one quary.*)
  export BUFFER_SIZE := 1000;

  export ParallelCalculateNMM, Get3DNMM, RecoverTranslationSamplePoints, GetOrderedCriticalPlanes,
  CriticalPlanes, CalculateNMM, LaunchFindDistinctSamplePoints,  LaunchComputeNMM,
  FetchSamplePointsFromDB, FetchTopologicallyDistinctSamplePointsFromDB,
  ParallelFindTopologicallyDistinctSamplePoints;


# Procedure: CriticalPlanes
#   Compute critical planes in the remainder range
#
# Parameters:
#   R                - the rotation matrix obtained from CayleyTransform
#   neighborhood     - a neighborhood for which one wants to compute NMM
#   kRange           - a range of planes to consider
#
# Output :
#   Returns a list of lists each containing critical planes for one direction
CriticalPlanes := proc(R::Matrix, neighborhood::list, kRange::list)
  # we remove the element [0, 0, 0]
  local n := subsop(ListTools:-Search([0,0,0], neighborhood)=NULL, neighborhood); 
  local T := combinat:-cartprod([n, kRange]);
  local planes := [[],[],[]], params;
  while not T[finished] do 
    params := T[nextvalue](); 
    planes[1] := [op(planes[1]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[1]]; 
    planes[2] := [op(planes[2]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[2]]; 
    planes[3] := [op(planes[3]), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[3]];
  end do;
 return planes;
end proc;


# Procedure: Get3DNMM
#   Compute neighbourhood motion maps
#
# Parameters:
#   neighborhood     - a neighborhood for which one wants to compute NMM
#   sampleTrans     - a midpoint of a frame in the remainder range
#   R                - the rotation matrix obtained from CayleyTransform
#
# Output:
#   Returns 3D neighborhood motion map for given vars and corresponding translations.
Get3DNMM := proc(neighborhood::list, sampleTrans::Vector, R::Matrix) 
  local n, NMM := [];
  for n in neighborhood do
    NMM := [op(NMM), convert(map[inplace](round, R.Vector(3, n) + sampleTrans), list)];
  od;
  return NMM;
end proc;


# Procedure: RecoverTranslationSamplePoints
#   Compute midpoints of each frame in the remainder range
#
# Parameters:
#   planes            - ordered list of critical planes in the remainder range
#
# Output:
#   Returns centers of frames in the remainder range
RecoverTranslationSamplePoints := proc(planes::list) 
  local s:
    s := proc(planes::list, i::integer, j::integer, k::integer)
      return [(1/2) * add(planes[1][i .. i+1]), (1/2)*add(planes[2][j.. j+1]), 
                                            (1/2)*add(planes[3][k .. k+1])];
  end proc;
  return [seq(seq(seq(s(planes, i, j, k), k=1..nops(planes[3]-1)),
               j=1..nops(planes[2]-1)), i =1..nops(planes[1])-1)];
end proc;


# Procedure: GetOrderedCriticalPlanes
#   Compute critical planes in the remainder range
#
# Parameters:
#   vars             - a list of variables
#   samplePoint      - a rotational sample point 
#   planes           - a precomputed list of planes in the remainder range
#
# Output:
#   Returns the ordered critical planes in the remainder range for X, Y and Z directions and order
#   signature.
GetOrderedCriticalPlanes := proc(vars::list, samplePoint::list, planes::list) 
  local params, sdPlanes := [[],[],[]];
  local xSig, ySig, zSig, Signature;
  
  if nops(vars) <> 3 then
    error "Only 3D arrangement is supported.";
  fi;
  
  sdPlanes[1] := eval(planes[1], [vars[1] = samplePoint[1], vars[2] = samplePoint[2], vars[3] =
                                                                                  samplePoint[3]]);
  sdPlanes[2] := eval(planes[2], [vars[1] = samplePoint[1], vars[2] = samplePoint[2], vars[3] =
                                                                                  samplePoint[3]]);
  sdPlanes[3] := eval(planes[3], [vars[1] = samplePoint[1], vars[2] = samplePoint[2], vars[3] =
                                                                                  samplePoint[3]]);
  
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


# Procedure: CalculateNMM
#   Reads data from hard drive and generates NMM
#
# Parameters:
#   vars             - list of variables in which the problem is expressed
#   planes           - a precomputed list of planes in the remainder range
#   buffer           - an Array which contains rotational sample points
#   N                - a given neighborhood
#   R                - a Matrix computed with Cayley transform
#   db               - an instance of ComputationRegister with open connection to the database
# Output:
#   A neighborhood motion map is saved in the database.
#
CalculateNMM := proc(vars::list, planes::list, buffer::Array, N::list, R::Matrix,
                                                         db::ComputationRegister) 
  local sdPlanes, trans, x, RR, y, NMM; 
  for x in buffer do
    RR := eval(R, [vars[1] = x[2], vars[2] = x[3], vars[3] = x[4]]);
    sdPlanes := GetOrderedCriticalPlanes(vars, x[2..()], planes)[2]; 
    trans := RecoverTranslationSamplePoints(sdPlanes); 
    for y in trans do
      NMM := Get3DNMM(N, Vector(3, y), RR);
      InsertNMM(db, x[1], NMM, y);
    od;
  end do;
end proc:



# Procedure: FetchTopologicallyDistinctSamplePointsFromDB
#   Used to fetch topologically distinct rotational sample points from the database. The number of 
#   sample points is controlled by BUFFER_SIZE.
#
# Parameters:
#   db::ComputationRegister       - an instance of ComputationRegister with open connection to the
#                                   database.
#   fromID::integer               - an id of the first sample point from a range to be fetched
#   last::integer                 - an id of the very last sample points to be fetched
#
# Output:
#   An Array of lists which of each represent a sample point. Note that the first element of each
#   list is an id of a given sample point.
FetchTopologicallyDistinctSamplePointsFromDB := proc(db::ComputationRegister, fromID::integer, 
                                                                                last::integer)
  local n := fromID + RigidMotionsRecoverNMM:-BUFFER_SIZE - 1;
  if n > last then
    n := last - fromID;
  fi;
  return FetchTopologicallyDistinctSamplePoints(db, fromID, n); 
end proc;


# Procedure: ParallelCalculateNMM
#   Uses Grid framework to generates unique NMM.
#
ParallelCalculateNMM := proc() 
  local db:= Object(ComputationRegister, dbPathGlobal), n;
  local noTPoints;
  local first::integer, last::integer;
  local R := CayleyTransform(varsGlobal), N := GetNeighborhood(nTypeGlobal); 
  local planes := RigidMotionsRecoverNMM:-CriticalPlanes(R, N, kRangeGlobal);
  local i::integer, buffer, samplePoint, sig::string;
  noTPoints := NumberOfTopologicallyDistinctSamplePoints(db);

  n := trunc(noTPoints / Grid:-NumNodes());
  first :=  Grid:-MyNode() * n + 1; last := (Grid:-MyNode() + 1) * n;
  if Grid:-MyNode() = Grid:-NumNodes() - 1 then
    last := noTPoints;
  fi;

  for i from first by RigidMotionsRecoverNMM:-BUFFER_SIZE to last do
    buffer := RigidMotionsRecoverNMM:-FetchTopologicallyDistinctSamplePointsFromDB(db, i, last);
    RigidMotionsRecoverNMM:-CalculateNMM(varsGlobal, planes, buffer, N, R, db);
    SynchronizeNMM(db);
  od;
  Close(db);
  Grid:-Barrier();
end proc:


# Procedure: FetchSamplePointsFromDB
#   Used to fetch rotational sample points from the database. The number of sample points is
#   controlled by BUFFER_SIZE.
#
# Parameters:
#   db::ComputationRegister       - an instance of ComputationRegister with open connection to the
#                                   database.
#   fromID::integer               - an id of the first sample point from a range to be fetched
#   last::integer                 - an id of the very last sample points to be fetched
#
# Output:
#   An Array of lists which of each represent a sample point. Note that the first element of each
#   list is an id of a given sample point.
FetchSamplePointsFromDB := proc(db::ComputationRegister, fromID::integer, last::integer)
  local n := fromID + RigidMotionsRecoverNMM:-BUFFER_SIZE - 1;
  if n > last then
    n := last - fromID;
  fi;
  return FetchSamplePointsWithoutSignature(db, fromID, n ); 
end proc;


# Procedure: ParallelFindTopologicallyDistinctSamplePoints
#   Finds rotational sample points which lead to unique arrangement of the critical planes in the
#   remainder range. 
#
 ParallelFindTopologicallyDistinctSamplePoints := proc() 
  local first::integer, last::integer;
  local R := CayleyTransform(varsGlobal), N := GetNeighborhood(nTypeGlobal); 
  local planes := RigidMotionsRecoverNMM:-CriticalPlanes(R, N, kRangeGlobal);
  local db, n, i::integer, buffer, samplePoint, sig, noTPoints;
  db:= Object(ComputationRegister, dbPathGlobal);
  noTPoints := NumberOfSamplePoints(db);
  n := trunc(noTPoints / Grid:-NumNodes());
  first :=  Grid:-MyNode() * n + 1; last := (Grid:-MyNode() + 1) * n;
  if Grid:-MyNode() = Grid:-NumNodes() - 1 then
    last := noTPoints;
  fi;

  for i from first by RigidMotionsRecoverNMM:-BUFFER_SIZE to last do
    buffer := RigidMotionsRecoverNMM:-FetchSamplePointsFromDB(db, i, last);
    for samplePoint in buffer do
      sig := RigidMotionsRecoverNMM:-GetOrderedCriticalPlanes(varsGlobal,samplePoint[2..()],planes)[1];
      InsertSignature(db, samplePoint[1], sig);
    od;
  od;
  SynchronizeSamplePointsSignatures(db);
  Close(db);
  Grid:-Barrier();
end proc;


LaunchFindDistinctSamplePoints := proc(vars::list, nType::string, kRange::list, dbPath::string, 
                                                        nodes:=kernelopts(numcpus))
  local db:=Object(ComputationRegister, dbPath);
  PrepareSamplePoints(db);
  Close(db);
  nTypeGlobal := nType; kRangeGlobal := kRange; dbPathGlobal := dbPath; varsGlobal := vars;
  
  Grid:-Setup("local"); 
  Grid:-Launch(RigidMotionsRecoverNMM:-ParallelFindTopologicallyDistinctSamplePoints,
  imports=['varsGlobal', 'nTypeGlobal', 'kRangeGlobal', 'dbPathGlobal'], numnodes=nodes,
  allexternal=false); 
  db:=Object(ComputationRegister, dbPath);
  CloseSignaturesAddition(db);
  Close(db);
end proc:


# Procedure: LaunchOnGridGetNMM
#   Setup and run computation on a local grid
#
# Parameters:
#   vars             - list of variables in which the problem is expressed
#   nType            - size of neighborhood i.e N1, N2 and N3. 
#   kRange           - a range of planes of the half-grid
#   dbPath           - a path to a database file.
#   nodes            - number of nodes used in the parallel computations
#
# Output:
#   List of unique neighborhood motion maps
LaunchComputeNMM := proc(vars::list, nType::string, kRange::list, dbPath::string, 
                                                        nodes:=kernelopts(numcpus)) 

  nTypeGlobal := nType; kRangeGlobal := kRange; dbPathGlobal := dbPath; varsGlobal := vars;
  
  Grid:-Setup("local"); 
  Grid:-Launch(RigidMotionsRecoverNMM:-ParallelCalculateNMM, imports=['varsGlobal', 'nTypeGlobal',
  'kRangeGlobal', 'dbPathGlobal'], numnodes=nodes, allexternal=false); 
end proc:

end module:

