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

  global nTypeGlobal::string;
  global dbPathGlobal::string
  global kRangeGlobal::list;
  global varsGlobal::list;


  export ParallelCalculateNMM;

# Procedure: Get3DNMM
#   Compute neighbourhood motion maps
#
# Parameters:
#   neighborhood     - a neighborhood for which one wants to compute NMM
#   samplesTrans     - midpoints of frames in the remainder range
#   sampleRot        - values of in determines (see CayleyTransform)
#   NMMContainer     - a reference to a container which stores the output NMM
#   R                - the rotation matrix obtained from CayleyTransform
#   vars             - a list of variables in which R is expressed
#
# Output:
#   Returns 3D neighborhood motion maps for given vars and corresponding translations.
Get3DNMM := proc(neighborhood::list, samplesTrans::list, sampleRot::list, NMMContainer::uneval, R,
                 vars::list) 
  local RR := eval(R, [vars[1] = sampleRot[1], vars[2] = sampleRot[2], vars[3] = sampleRot[3]]);
  local x, NMM := eval(NMMContainer); 
  NMMContainer := NMM union {seq(map(proc (y) map(round, convert(RR.Vector(3, y)+Vector(3,x), list))
  end proc, neighborhood),x in samplesTrans)};
end proc;


# Procedure: RecoverTranslationSamplePoints
#   Compute midpoints of each frame in the remainder range
#
# Parameters:
#   planes            - ordered list of planes in the remainder range
#
# Output:
#   Returns centers of frames in the remainder range
RecoverTranslationSamplePoints := proc(planes::list) 
  local s:
    s := proc(planes::list, i::integer, j::integer, k::integer)
      return [(1/2)*add(planes[1][i .. i+1]), (1/2)*add(planes[2][j.. j+1]), 
              (1/2)*add(planes[3][k .. k+1])];
  end proc:
  return [Threads:-Seq( seq( seq( s(planes, i, j, k), k=1..nops(planes[3]-1) ),
  j=1..nops(planes[2]-1) ), i =1..nops(planes[1])-1 )]:
end proc:


# Procedure: GetOrderedCriticalPlanes
#   Compute critical planes in the remainder range
#
# Parameters:
#   vars             - a list of variables
#   samplePoints     - values of vars[1], vars[2], vars[3] (see CayleyTransform)
#   planes           - a precomputed list of planes in the remainder range
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
  
  sdPlanes[1] := eval(sdPlanes[1], {vars[1] = samplePoints[1], vars[2] = samplePoints[2], vars[3] = samplePoints[3]});
  sdPlanes[2] := eval(sdPlanes[2], {vars[1] = samplePoints[1], vars[2] = samplePoints[2], vars[3] = samplePoints[3]});
  sdPlanes[3] := eval(sdPlanes[3], {vars[1] = samplePoints[1], vars[2] = samplePoints[2], vars[3] = samplePoints[3]});
  
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


# Procedure: CriticalPlanes
#   Compute critical planes in the remainder range
#
# Parameters:
#   R                - the rotation matrix obtained from CayleyTransform
#   neighborhood     - a neighborhood for which one wants to compute NMM
#   kRange           - a range of planes to consider
#
# Output:
#   Returns a list of lists each containing critical planes for one direction
CriticalPlanes := proc(R::Matrix, neighborhood::list, kRange::list)
  # we remove 7th element -- [0, 0, 0]
  local n := subsop(7=NULL, neighborhood); 
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
#   dirIn            - a directory which contains sample points
#   dirOut           - a directory to which save NMM
#   prefixIn         - a prefix of files which contains sample points
#   prefixOut        - a prefix of files which will contain NMM
#   id               - an id of a file which contains sample points
#
# Output:
#   Set of NMM
#
# TODO:
#   Export printing code into a procedure
#   Simplify the code
CalculateNMM := proc(nType::string, kRange::list, dirIn::string, dirOut::string, prefixIn::string,
                     prefixOut::string, id::integer) 
  local data, s, fileID, vars := [a,b,c];
  local R := CayleyTransform(vars);
  local neighborhood := GetNeighborhood(nType); 
  local i, sdPlanes := [], trans, NMM:={}, rotSamp; 
  local percent := -1, Signatures := {}, sig, msg;
  local n := Grid:-NumNodes();
  local planes := CriticalPlanes(R, neighborhood, kRange);

  # read sample points from a file and convert them to appropriate format
  data := ImportMatrix(sprintf("%s/%s%d.tsv", dirIn, prefixIN, id), 'source' = 'tsv', 'datatype'='string'); 
  data := Threads:-Map(proc (x) op(sscanf(x, %a)) end proc, data); 

  # print a progress message
  for i to upperbound(data)[1] do
    if percent <>  round( i / upperbound(data)[1] * 100 ) then
      percent := round( i / upperbound(data)[1] * 100 );
      if percent mod 10 = 0 then
        printf("[%s]:> Node: %d, finished: %d%%.\n",
               StringTools:-FormatTime("%Y-%m-%d -- %R"),id, percent); 
      end if;
    end if;
    rotSamp := convert(data[i], list); 
    # discard the sample points for which NMM are already computed
    sig, sdPlanes := GetOrderedCriticalPlanes(vars, rotSamp, planes); 
    if member( sig, Signatures ) then
      next;
    end if;

    trans := RecoverTranslationSamplePoints(sdPlanes); 
    Get3DNMM(neighborhood, trans, rotSamp, NMM, R, vars);

    Signatures := Signatures union {sig};
  end do;

  # write the result
  fileID := fopen(sprintf("%s/%s%d.tsv", dirOur, prefixOut, id), WRITE, TEXT);
  writedata(fileID, [op(NMM)], string, proc (f, x) fprintf(f, "%a", x) end proc);
  fclose(fileID);
end proc:

# Procedure: ParallelCalculateNMM
#   Uses Grid framework to reads data from hard drive and generates NMM
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#   dirIn            - a directory which contains sample points
#   dirOut           - a directory to which save NMM
#   prefixIn         - a prefix of files which contains sample points
#   prefixOut        - a prefix of files which will contain NMM
#
# Output:
#   Set of NMM
ParallelCalculateNMM := proc() 
  local db:=Object(ComputationRegister, dbPathGlobal), n;
  # events-1 because the last event is a copy of events[-2]
  n := trunc((NumberOfSamplePoints(db)-1) / Grid:-NumNodes());
  CalculateNMM(varsGlobal, nTypeGlobal, kRangeGlobal, db, me);
  Grid:-Barrier();
end proc:

# Procedure: LaunchOnGridGetNMM
#   Setup and run computation on a local grid
#
# Parameters:
#   nType            - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange           - a range of planes to consider
#
# Output:
#   Write a set of NMM to a file given by dirOut/prefixOut(id).tsv
LaunchOnGridGetNMM := proc(vars::list, nType::string, kRange::list, dbPath::string, 
                                                        nodes:=kernelopts(numcpus)) 

  local db:=Object(ComputationRegister, dbPathGlobal);
  PrepareSamplePoints(db);
  Close(db);
  Grid:-Setup("local"); 
  nTypeGlobal := nType; kRangeGlobal := kRange; dbPathGlobal := dbPath; varsGlobal := vars;
  Grid:-Launch(RigidMotionsRecoverNMM:-ParallelCalculateNMM, 
               imports=['varsGlobal', 'nTypeGlobal', 'kRangeGlobal', 'dbPathGlobal'], 
                                                                     numnodes=nodes);
end proc:


end module;
