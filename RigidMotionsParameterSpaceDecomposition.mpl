# File: RigidMotionsParameterSpaceDecompostion.mpl  
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

RigidMotionsParameterSpaceDecompostion := module() 
  option package;
  uses   RigidMotionsParameterSpaceDecompostionRecursive, RigidMotionsParameterSpaceCommon;

  (*Threshold which controls when to synchronize databases.*)
  local RECORDS_TO_SYNCH := 1000;

  local  GetQuadric, IsMonotonic, ComputeSetOfQuadrics,
         ComputeEventsATypeGrid, ComputeEventsBTypeGrid, ComputeEventsCTypeGrid,
         ComputeAsymptoticABEventsGrid, ComputeAsymptoticAAEvents, 
         ComputeEventsFromAlgebraicNumbers, IsAsymptotic;
         
  export IsAsymptoticIntersection, LaunchComputeSamplePoints, 
         LaunchResumeComputations, ComputeSamplePoints, ParallelComputeSamplePoints,
         ParallelComputeSamplePointsResume;

  #Variables shared by grid nodes;
  global Q, events, vars, dbPath, skipped;

# Procedure: GetQuadric
#   Compute a quadric.
#
# Parameters:
#   R          - a 3x3 rotation matrix obtained from Cayley transform
#   neighbor   - a vector which points on a neighbor
#   hGridPlane - a half-grid plane i.e. k_dim + 1/2
#   axis       - a given row of R matrix
#   
# Output:
#   Multivariate polynomial of degree 2 related to the changed of configuration of image patch and
#   given in [Quadrics:2016] as equation (7).
GetQuadric := proc(R::~Matrix, neighbor::~Vector, hGridPlane::~Vector, axis::integer)
  local r := R . neighbor - hGridPlane;
  if axis > 0 and axis < 4 then
    return simplify( ( 2 * denom( R[1][1] ) ) * r[axis] );
  else
    error "Wrong dimension: dim must be in [1,3].";
  end if:
end proc:


# Procedure: IsMonotonic
#   Check if given polynomial of degree 2 is non-positive or non-negative
#
# Parameters:
#   x      - a polynomial
#
# Output:
#   true if polynomial of degree 2 is non-positive and false otherwise
IsMonotonic := proc( x::~polynom )
  local homo, hessian, signmap, clean:
  homo := Groebner:-Homogenize( x, v ):
  hessian := 1 / 2 * VectorCalculus:-Hessian( homo, [ op( indets( x ) ), v ] ):
  signmap := map( signum, LinearAlgebra:-Eigenvalues( hessian ) ):
  signmap := convert( signmap, list ):
  clean := remove( `=`, signmap , 0 ):
  return ( andmap( `=`, clean, 1 ) or andmap( `=`, clean, -1 ) ):
end proc:


# Procedure: ComputeSetOfQuadrics
#   Compute a set of quadrics which are reduced by duplicated and these ones
#   which are strictly positive or negative.
#
# Parameters:
#   R          - a rotation matrix obtained from Cayley transform
#   nType      - size of neighborhood i.e. N_1, N_2, N_3. 
#   axis       - a given axis i.e. 1 = x, 2 = y and 3 = z 
#   kRange     - a range of planes to consider
#
# Output:
#   List of quadrics of form f - g reduced by duplicated---up to constant
#   factor---and these ones which are strictly positive or negative.
ComputeSetOfQuadrics := proc( R::~Matrix,
                               nType::string,
                               axis::~integer,
                               kRange::~list )
  local neighborhood := GetNeighborhood( nType ):
  local quadrics := Array([]):
  local f::polynom, g::polynom, quadric::polynom:
  local indexes := []:
  local T := combinat:-cartprod( [ neighborhood, neighborhood, kRange, kRange ] ):

  if LinearAlgebra:-Determinant(R) <> 1 or [upperbound(R)] <> [3,3] or nops(indets(R)) = 0 then
    error "Used matrix is not a correct 3 x 3 rotation matrix!";
  fi;

  if axis < 1 or axis > 3 then
    error "Axis is not from range [1,3].";
  fi:

  while not T[ finished ] do
    indexes := T[ nextvalue ]():
    f := GetQuadric( R, Vector( indexes[1] ), Vector( 3, indexes[3] + 1/2 ), axis ):
    g := GetQuadric( R, Vector( indexes[2] ), Vector( 3, indexes[4] + 1/2 ), axis ):
    quadric := f - g:
    if quadric <> 0 then
      quadric := quadric / lcoeff( quadric ):
      ArrayTools:-Append(quadrics,quadric):
    end if:
  end do:
  quadrics := convert(quadrics,set):
  quadrics := remove( IsMonotonic, quadrics ):
  quadrics := Threads:-Map(proc(x) normal(x) * denom(x) end proc, quadrics):
  return quadrics:
end proc:


# Procedure: IsAsymptotic
#   Checks if given quadric has an asymptotic critical value. For this moment a direction is fixed
#   to a.
#
# Parameters:
#   x          - a quadric in three variables
#   vars       - a list of variables
#
# Output:
#   List of solution for which partial derivatives in vars[2] and vars[3] are collinear.
IsAsymptotic := proc(x::polynom, vars::list)
  local vec, Vb, Vc, VV, sols;
  if nops(vars) < 3 then
    error "Expected at least three variables!";
  fi;
  vec := VectorCalculus:-Gradient(x, vars);
  Vb := Vector(3, [coeff(vec[2],vars[2]),coeff(vec[2],vars[3]),eval(vec[2],[vars[-2]=0,vars[-1]=0])]);
  Vc := Vector(3, [coeff(vec[3],vars[2]),coeff(vec[3],vars[3]),eval(vec[3],[vars[-2]=0,vars[-1]=0])]);
  VV := LinearAlgebra:-CrossProduct(Vb, Vc);
  if norm(VV,1) = 0 then
    return {{vars[1]=0}};
  fi;
  sols := solve({VV[1] = 0, VV[2] = 0,VV[3] = 0,vars[1]>=0}, [vars[1]]);
  if sols = NULL then
    return {};
  else
    return sols;
  fi;
end proc:


# Procedure: IsAsymptoticIntersection
#   Checks if intersection of two quadrics has an asymptotic critical value. For this moment a
#   direction is fixed to a.
#
# Parameters:
#   p          - a quadric in three variables
#   q          - a quadric in three variables
#   vars       - a list of variables
#
# Output:
#   List of solution for which intersection of p and q has an asymptotic intersection
# Comment:
#   - Since Groebner package seems to have memory leak I should rather replace 
#     PolynomialIdeals:-EliminationIdeal by resultant elimination similarly to what I did with
#     univariate polynomials.
# TODO:
#   - Allow user to chose a direction. 
#   - if there is no intersection between quadrics then skip it.
IsAsymptoticIntersection := proc( p::polynom, q::polynom, vars::list )
  local J := PolynomialIdeals:-`<,>`(p,q);
  local Pb, Pc, Cb, Cc, sols, univ;
  if nops(vars) < 3 then
    error "Expected at least three variables!";
  fi;
  Pb := PolynomialIdeals:-EliminationIdeal(J,{op(vars[1..2])}):
  Pb := PolynomialIdeals:-IdealInfo:-Generators(Pb)[1];
  Pc := PolynomialIdeals:-EliminationIdeal(J,{op(vars[[1,3]])});
  Pc := PolynomialIdeals:-IdealInfo:-Generators(Pc)[1];
  Cb := lcoeff(Pb, vars[2]);
  Cc := lcoeff(Pc, vars[3]);
  return gcd(Cb,Cc);
end proc:


# Procedure: ComputeEventsATypeGrid
#   Compute events such that a sweep plane is tangent to a quadric. 
#
# Parameters:
#   Q          - a set of quadrics or conics
#   dim        - a list of indexes of variables used to calculate partial derivatives
#   vars       - a list of variables
#
# Output:
#   List of ranges which contains roots of a system(q, d/db q, d/dc q).
#
# Comment:
#  - only the first direction is supported since EliminationResultant is used
ComputeEventsATypeGrid := proc( Q, dim::list, vars::list )
  local s, result := Array([]);
  if nops(vars) < 3 then
    error "Expected at least three variables!";
  fi;
  s := proc(i::integer, vars::list)
   local sys, univ;
   local q := Q[i];
   sys := { q, diff( q, vars[ dim[1] ] ), diff( q, vars[ dim[2] ] ) };
   univ := EliminationResultant(sys, vars);
   return SerializeEvents(GenerateEvents(univ, [i]));
  end proc:
  map[inplace](proc(x) ArrayTools:-Extend(result, x, inplace=true) end proc, 
                                             [Grid:-Seq(s(i, vars), i=1..nops(Q))]);
  return ReconstructEvents(result);
end proc:


# Procedure: ComputeEventsBType
#   Compute events such that intersection of two quadrics is tangent to sweeping plane.
#
# Parameters:
#   dir        - a direction of a gradient product it should
#                be the same as director of sweep 
#   Q          - a set of conics
#   vars       - a list of variables
#
# Output:
#   Indexes of quadrics which intersect and a component of a vector product of 
#   their gradients in given direction have a common root.
ComputeEventsBTypeGrid := proc( Q, dir::integer, vars::list )
  local s, result := Array([]);
  if nops(vars) < 3 then
    error "Expected at least three a variables!";
  fi;
  s := proc(i, j, vars::list)
    local p, prod, univ, sys:
    prod := LinearAlgebra:-CrossProduct( VectorCalculus:-Gradient( Q[i], vars ),
                                    VectorCalculus:-Gradient( Q[j], vars ) )[dir]:
    sys := { Q[i], Q[j], prod };
    univ := EliminationResultant(sys, vars);
    return SerializeEvents(GenerateEvents(univ, [i, j]));
  end proc:
  map[inplace](proc(x) ArrayTools:-Extend(result, x, inplace=true) end proc, 
                                   [Grid:-Seq(seq(s(i, j, vars), j=i+1..nops(Q)), i=1..nops(Q))]);
  return ReconstructEvents(result);
end proc;


# Procedure: ComputeEventsCTypeGrid
#   Compute events such that three quadrics intersects in a point. 
#
# Parameters:
#   Q          - a set of quadrics
#   vars       - a list of variables
#
# Output:
#   Indexes of quadrics which intersect in a point.
ComputeEventsCTypeGrid := proc( Q, vars::list )
  uses ArrayTools;
  local result := Array([]), s;
  if nops(vars) < 3 then
    error "Expected at least three a variables!";
  fi;
  s := proc (i, j, k, vars::list)
    local univ, sys;
    sys := { Q[i], Q[j], Q[k] }:
    univ := EliminationResultant(sys, vars);
    return SerializeEvents(GenerateEvents(univ, [i, j, k]));
  end proc;
  map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [Grid:-Seq(seq(seq(s(i, j, k, vars), 
                                       k=j+1..nops(Q)),j=i+1..nops(Q)),i=1..nops(Q))]);
  return ReconstructEvents(result);
end proc:


# Procedure: ComputeAsymptoticAAEvents
#   Compute real algebraic numbers which corresponds to 
#   asymptotic cases given by one quadrics.
#
# Parameters:
#   Q          - a set of quadrics
#   vars       - a list of variables
#
# Output:
#   A list of real algebraic numbers and indexes of quadrics -- events.
ComputeAsymptoticAAEvents:=proc(Q, vars::list)
  uses ArrayTools;
  local result := Array([]), s;
  if nops(vars) < 3 then
    error "Expected at least three a variables!";
  fi;
  s:=proc(i::integer, vars::list)
    local events := Array([]), sol, rootsF, tmp:
    rootsF := IsAsymptotic(Q[i], vars):
    for sol in rootsF do
      sol := op(sol);
      if not type(rhs(sol), rational) then
        error "Irrational asymptotic case! Are you sure the input is a set of quadrics?"
      fi:
      ArrayTools:-Append(events, EventType(RealAlgebraicNumber(lhs(sol) * denom(rhs(sol)) -
      numer(rhs(sol)), rhs(sol), rhs(sol)), [i]));
    od:
    return events;
  end proc;
  map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [seq(s(i, vars), i=1..nops(Q))]);
  return result;
end proc;


# Procedure: ComputeAsymptoticABEvents
#   Compute real algebraic numbers which corresponds to 
#   asymptotic cases given by one quadrics.
#
# Parameters:
#   Q          - a set of quadrics
#   vars       - a list of variables
#
# Output:
#   A list of real algebraic numbers and indexes of quadrics -- events.
ComputeAsymptoticABEventsGrid:=proc(Q, vars::list)
  uses ArrayTools;
  local result := Array([]), s;
  if nops(vars) < 3 then
    error "Expected at least three a variables!";
  fi;
  s:=proc(i::integer, j::integer, vars::list)
   local poly;
   poly := RigidMotionsParameterSpaceDecompostion:-IsAsymptoticIntersection(Q[i], Q[j], vars);
   if poly = NULL or nops(poly) = 0 then
     return [];
   fi;
   return SerializeEvents(GenerateEvents(poly, [i, j]));
  end proc;
  map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [Grid:-Seq(seq(s(i, j, vars),
  j=i+1..nops(Q)), i=1..nops(Q))]);
  return ReconstructEvents(result);
end proc;


# Procedure: ComputeEventsFromAlgebraicNumbers
#   Compute and sort events
#
# Parameters:
#   Q     - set of quadrics or conics
#   vars       - a list of variables
# Output:
#   Sorted set of events
ComputeEventsFromAlgebraicNumbers := proc( Q, vars::list )
  local events := Array([]);
  ArrayTools:-Extend(events, ComputeEventsATypeGrid( Q, [2, 3], vars ), inplace=true);
  ArrayTools:-Extend(events, ComputeEventsBTypeGrid( Q, 1, vars ), inplace=true);
  ArrayTools:-Extend(events, ComputeEventsCTypeGrid( Q, vars ), inplace=true);
  ArrayTools:-Extend(events, ComputeAsymptoticAAEvents(Q, vars), inplace=true);
  ArrayTools:-Extend(events, ComputeAsymptoticABEventsGrid(Q, vars), inplace=true);
  return AlgebraicSort(events);
end proc:


# Procedure: ComputeSamplePoints
#   Computes sample points for rotational part of rigid motions
#
#
# Parameters:
#   Q                  - list of quadrics
#   events             - an array of events
#   first              - integer value which indicates a first event to proceed.
#   last               - integer value which indicates a last event to proceed.
#   vars               - list of variables in which conics are expressed
#   db                 - an instance of the class ComputationRegister
#   skipped            - a list of the events' indices to be skipped
#
# Output:
#   It populates a database, given by databasePath, with sample points.
ComputeSamplePoints := proc (Q, events::Array, first::integer, last::integer, 
                             vars::list, db::ComputationRegister, skipped::list:=[]) 
local i, x, midpoint, sys, samplePoints, disjointEvent:=[], ranumI, ranumJ;
  if first < 0 or last < 0 or last < first or upperbound(events) <= last then 
    error "Bounds of the list of the events range are incorrect.": 
  end if:
  for i from first to last do 
    if i in skipped then
      next;
    fi:
    sys := Q[GetQuadrics(events[i])];
    ranumI := GetRealAlgebraicNumber(events[i]);
    ranumJ := GetRealAlgebraicNumber(events[i+1]);
    disjointEvent:=DisjointRanges(ranumI, ranumJ);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
    # never call eval with sets!!
    sys := eval(sys, vars[1] = midpoint);
    LaunchComputeSamplePoints2D(sys, midpoint, 1, false, vars[2..], db);
    SynchronizeSamplePoints(db);
    InsertSkippedCluster(db, i);
  end do;
end proc:


# Procedure: ParallelComputeSamplePoints
#   Computes sample points for rotational part of rigid motions. It should be call via Grid
#   framework.
ParallelComputeSamplePoints := proc()
  local me, numNodes, n;
  local db:=Object(ComputationRegister, dbPath);
  me := Grid:-MyNode();
  numNodes := Grid:-NumNodes();
  # events-1 because the last event is a copy of events[-2]
  n := trunc((upperbound(events)-1)/numNodes);
  # recreate events
  ReconstructEvents(events);
  RigidMotionsParameterSpaceDecompostion:-ComputeSamplePoints(Q, events, me*n+1,(me+1)*n, vars, 
                                                                                  db, []);
  Close(db);
  Grid:-Barrier();
end proc:


# Procedure: LaunchOnComputeSamplePoints
#   Computes sample points for rotational part of rigid motions using the grid framework
#
#
# Parameters:
#   variables     - list of variables in which the problem is expressed
#   databasePath  - a path to a copy of the database CompRegister.db
#   nType         - neighborhood type: N1, N2 or N3.
#   kRange        - range of grid lines passed as a list
#   nodes         - number of nodes used in the parallel computations
# Output:
#   It populates a database given by databasePath.
LaunchComputeSamplePoints := proc(variables::list, databasePath::string, nType::string, 
                                  kRange::list, nodes:=kernelopts(numcpus)) 
  local events, lastEvent, R, boundTmp, i, j, k, mesg;
  local db:=Object(ComputationRegister, databasePath);
  vars:=variables;
  dbPath:=databasePath;
  mesg:=kernelopts(printbytes=false);
  R := CayleyTransform(variables);
  Q := ListTools:-MakeUnique([op(ComputeSetOfQuadrics(R, nType, 1, kRange)), 
       op(ComputeSetOfQuadrics(R, nType, 2, kRange)),
       op(ComputeSetOfQuadrics(R, nType, 3, kRange)), op(variables)]);
  for i from 1 to nops(Q) do
    InsertQuadric(db, i, Q[i]);
  od;
  SynchronizeQuadrics(db);
  events := ComputeEventsFromAlgebraicNumbers(Q, variables);
  events := select[flatten](proc(x) evalb(GetInterval(GetRealAlgebraicNumber(x))[2] >= 0) end proc,
                           events);
  events := ReduceEvents(events);
  AdjustEvents(events, upperbound(Q), variables);
  #Insert events into the register
  for i from 1 to upperbound(events) do
    InsertEvent(db, i, events[i]);
  od;
  SynchronizeEvents(db);
  Close(db);

  if nodes > 1 then
    SerializeEvents(events);
    Grid:-Setup("local"):
    Grid:-Launch(RigidMotionsParameterSpaceDecompostion:-ParallelComputeSamplePoints, 
                 imports=['Q', 'events', 'vars', 'dbPath'], numnodes=nodes);
  else
    ComputeSamplePoints(Q, events, 1, upperbound(events) - 1, variables, db, skipped);             
  fi;
  mesg:=kernelopts(printbytes=mesg);
end proc:


# Procedure: ParallelComputeSamplePointsResume
#   Computes sample points for rotational part of rigid motions. It should be call via Grid
#   framework.
ParallelComputeSamplePointsResume := proc()
  local me, numNodes, n, events, skipped;
  local db:=Object(ComputationRegister, dbPath);
  me := Grid:-MyNode();
  numNodes := Grid:-NumNodes();
  skipped := FetchSkippedClusters(db);
  # events-1 because the last event is a copy of events[-2]
  n := trunc((NumberOfClusters(db)-1)/numNodes);
  # recreate events
  skipped := FetchSkippedClusters(db);
  events := FetchClusters(db, me* n+1,(me+1)*n+1); 
  RigidMotionsParameterSpaceDecompostion:-ComputeSamplePoints(Q, events, me*n+1,(me+1)*n, vars, 
                                                                                  db, skipped);
  Close(db);
  Grid:-Barrier();
end proc:

# Procedure: LaunchResumeComputations
#   Resumes computations of sample points for rotational part of rigid motions using
#   the grid framework.
#
#
# Parameters:
#   variables     - list of variables in which the problem is expressed
#   databasePath  - a path to a database file. If file does not exist it will be crated.
#   nType         - neighborhood type: N1, N2 or N3.
#   kRange        - range of grid lines passed as a list
#   nodes         - number of nodes used in the parallel computations
# Output:
#   It populates a database given by databasePath.
LaunchResumeComputations := proc(variables::list, databasePath::string, nType::string, 
                                 nodes:=kernelopts(numcpus))
  local events, firstEvent, rootTmp, i, mesg;
  local db:=Object(ComputationRegister, databasePath);
  vars:=variables;
  dbPath:=databasePath;
  mesg:=kernelopts(printbytes=false):
  Q := FetchQuadrics(db);
  Close(db);
  Grid:-Setup("local");
  Grid:-Launch(RigidMotionsParameterSpaceDecompostion:-ParallelComputeSamplePointsResume, 
               imports=['Q', 'vars', 'dbPath'], numnodes=nodes);
  mesg:=kernelopts(printbytes=mesg);
end proc;

end module:
