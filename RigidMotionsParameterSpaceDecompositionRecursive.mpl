# File: RigidMotionsParameterSpaceDecompostionRecursive.mpl  
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
#  12/24/2016 
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
#
RigidMotionsParameterSpaceDecompostionRecursive := module() 
  option package;

  (*Threshold which controls when to synchronize databases.*)
  local RECORDS_TO_SYNCH := 1000;

  uses   RigidMotionsParameterSpaceCommon;
  local  ComputeEventsAType2D, ComputeEventsBType2D, ComputeEventsAType1D,
         ComputeAsymptoticAAEvents2D, ComputeEventsAlgebraicNumbers2D, ComputeSamplePoints2D,
         ComputeSamplePoints2DBounded;
         
  export IsAsymptotic2D, LaunchComputeSamplePoints2D;

# Procedure: ComputeEventsAType2D
#   Compute events such that a sweep line is tangent to a conic with respect to the first provided
#   variable. 
#
# Parameters:
#   Q2D2       - a set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars2D       - a list of the variables
#
# Output:
#   list of couples: univariate polynomial and index of the generating conic: [poly, [index]]
#
# Comment/Limitations:
#  - Univariate polynomials are expressed in the first variable.
ComputeEventsAType2D := proc(Q2D2, grid::boolean, vars2D::list)
  uses ArrayTools;
  local result := Array([]), s;
  if nops(vars2D) < 2 then
    error "Only systems in at least two variables are supported. The variables are: %1.", vars2D;
  fi:
  s := proc(i::integer, vars2D::list)
    local sys, univ, sol:
    local q := Q2D2[i];
    sys := { q, diff( q, vars2D[2] ) };
    univ := UnivariatePolynomial(sys, vars2D):
    return SerializeEvents(GenerateEvents(univ, [i], {"A"}));
  end proc;
  if grid then
    map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [Grid:-Seq( s( i, vars2D ), 
                                                                   i=1..nops( Q2D2 ) )]);
  else 
    map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [seq( s( i, vars2D ), 
                                                             i=1..nops( Q2D2 ) )]);
  fi;
  return ReconstructEvents(result);
end proc:

# Procedure: ComputeEventsBType2D
#   Compute events such that two conics intersects in a point.
#
# Parameters:
#   Q2D2       - a set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars2D       - a list of the variables
#
# Output:
#   list of couples: univariate polynomial and indices of the generating conics: [poly,
#   [index1,index2]]
#
# Comment/Limitations:
#  - Univariate polynomials are expressed in the first variable.
ComputeEventsBType2D := proc(Q2D2, grid::boolean, vars2D::list)
  uses ArrayTools;
  local result := Array([]), s:
  s := proc (i, j, vars2D::list)
    local p, sol, univ, sys;
    sys := {Q2D2[i], Q2D2[j]}:
    univ := UnivariatePolynomial(sys, vars2D);
    return SerializeEvents(GenerateEvents(univ, [i, j], {"B"}));
  end proc;
  if grid then
    map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [Grid:-Seq( seq( s(i, j, vars2D), 
                                                       j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]);
  else 
    map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [seq( seq( s(i, j, vars2D), 
                                                 j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]);
  fi;
  return ReconstructEvents(result);
end proc:

# Procedure: IsAsymptotic2D
#   Checks if a curve is an asymptote.
#
# Parameters:
#   p          - a curve given as a polynomial
#   var        - variable to be reduced
#
# Output:
#   Monomial such that there exists a line tangent to the curve at infinity represented in a
#   variable different than 'var' or a constant.
IsAsymptotic2D := proc(p::polynom, var)
  return lcoeff(p, var);
end proc:


# Procedure: ComputeAsymptoticAAEvents2D
#   Compute real algebraic numbers which corresponds to asymptotic cases given by one curve.
#
# Parameters:
#   Q2D2          - a set of conics
#   vars2D       - a list of the variables
#
# Output: Real
#   A list of Events (see EventType).
ComputeAsymptoticAAEvents2D := proc(Q2D2, vars2D::list)
  uses ArrayTools;
  local result := Array([]), s;
  s:=proc(i::integer, vars2D::list)
    local asy;
    asy := RigidMotionsParameterSpaceDecompostionRecursive:-IsAsymptotic2D(Q2D2[i], vars2D[-1]);
    return GenerateEvents(asy, [i], {"AA"});
  end proc:
  map[inplace](proc(x) Extend(result, x, inplace=true) end proc, [seq(s(i, vars2D),
                                                                i=1..nops(Q2D2))]);
  return result;
end:


# Procedure: ComputeEventsAlgebraicNumbers2D
#   Compute and sort events
#
# Parameters:
#   Q2D2     - set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars2D       - a list of the variables
#
# Output:
#   Sorted Array of real algebraic numbers
ComputeEventsAlgebraicNumbers2D := proc(Q2D2, grid::boolean, vars2D::list)
  uses ArrayTools;
  local events2D := Array([]);
  Extend(events2D, ComputeEventsAType2D(Q2D2, grid, vars2D), inplace=true);
  Extend(events2D, ComputeEventsBType2D(Q2D2, grid, vars2D), inplace=true);
  Extend(events2D, ComputeAsymptoticAAEvents2D(Q2D2, vars2D), inplace=true);
  return AlgebraicSort(events2D);
end proc:


# Procedure: ComputeEventsAType1D
#   Compute and sort events as algebraic numbers 
#
# Parameters:
#   Q2D2     - a list of conics
#
# Output:
#   Sorted list of real algebraic numbers without repetitions.
ComputeEventsAType1D := proc(Q2D2::list)
  local q, factored, sqrFree, rootsF, rf, numbers := Array([]);
  for q in Q2D2 do
    if RootFinding:-HasRealRoots(q) then
      factored := factors( q )[2,..,1];
      for sqrFree in factored do
        rootsF := RootFinding:-Isolate(sqrFree, output='interval');
        for rf in rootsF do
          ArrayTools:-Append(numbers, RealAlgebraicNumber(sqrFree, op(rf)[2][1],
                             op(rf)[2][2]));
        od;
      od;
   fi;
  od;
  numbers := AlgebraicSort(numbers);
  return ListTools:-MakeUnique(convert(numbers, list), 1, proc(a,b) evalb(Compare(a,b) = 0) end proc);
end proc:


# Procedure: ComputeSamplePoints2D
#   Computes sample points for rotational part of rigid motions
#
# Parameters:
#   Q2D         - a list of conics
#   events2D    - each element contains a list of algebraically unique events
#   first       - integer value which indicates the first event to proceed.
#   last        - integer value which indicates the last event to proceed.
#   vars2D      - list of variables in which conics are expressed
#   db          - an instance of the class ComputationRegister
#
# Output:
#   It populates a database, given by databasePath, with sample points.
ComputeSamplePoints2D := proc(Q2D, events2D::Array, first::integer, last::integer,
                              vars2D::list, amid::rational, db::ComputationRegister, BOUNDED:=false)
  local i::integer, j::integer, x::list, midpoint::rational, sys::list, records := 0;
  local disjointEvent::list, oneD::list, oneDNeg::list, ranumI, ranumJ;
  if first < 0 or last < 0 or last < first or upperbound(events2D) <= last then 
    error "Bounds of the array range are incorrect.": 
  end if:
  for i from first to last do 
    sys := Q2D[GetQuadrics(events2D[i])];
    ranumI := GetRealAlgebraicNumber(events2D[i]);
    ranumJ := GetRealAlgebraicNumber(events2D[i+1]);
    disjointEvent:=DisjointRanges(ranumI, ranumJ);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
    # never call eval with sets!
    sys := eval(sys, vars2D[1] = midpoint);
    oneD := ComputeEventsAType1D(sys);
    oneD := select(proc(x) return evalb(GetInterval(x)[2] >= 0); end proc, oneD);

    if upperbound(oneD) > 0 then
      for j from 1 to upperbound(oneD) - 1 do
        disjointEvent := DisjointRanges(oneD[j],oneD[j+1]);
        InsertSamplePoint(db, [amid, midpoint, (GetInterval(disjointEvent[1])[2] +
                                            GetInterval(disjointEvent[2])[1])/2]);
        records := records + 1;
        if records mod RECORDS_TO_SYNCH = 0 then
          SynchronizeSamplePoints(db);
        fi;
      od:
      if not BOUNDED then
        InsertSamplePoint(db, [amid, midpoint, GetInterval(oneD[-1])[2] + 1/2]);
      fi;
      records := records + 1;
      if records mod RECORDS_TO_SYNCH = 0 then
        SynchronizeSamplePoints(db);
      fi;
    fi;
 od:
end proc:


# Procedure: ComputeSamplePoints2D
#   Computes sample points for rotational part of rigid motions using the grid framework
#
#
# Parameters:
#   s         - a list of conics
#   midpoint  - the first dimensional midpoint obtained from the 3D decomposition
#   nodes     - number of nodes used in the parallel computations
#   grid      - a control variable for parallel computations. If true and additional conditions on
#               the size of the problem are fulfilled the problem is solved in the grid framework
#   variables - list of variables in which the problem is expressed
#   db        - an instance of ComputationRegister which provides interface to the database
# Output:
#   It populates a database with sample points.
LaunchComputeSamplePoints2D := proc (s::list, midpoint::rational, nodes::integer,
                                           grid::boolean, variables::list, db::ComputationRegister,
                                           boundedOnly := false) 
  local events2D, Q2D := ListTools:-MakeUnique([op(variables), op(s)]);
  if grid and nops(s) > 20 then
     events2D := ComputeEventsAlgebraicNumbers2D(Q2D, true, variables);
  else
     events2D := ComputeEventsAlgebraicNumbers2D(Q2D, false, variables);
  fi;
  events2D := remove[flatten](proc(x) evalb(GetInterval(GetRealAlgebraicNumber(x))[2] < 0) end proc, 
                            events2D);
  if upperbound(events2D) = 0 then
    return NULL;
  fi;
  events2D := ReduceEvents(events2D);
  AdjustEvents(events2D, upperbound(Q2D), variables);
  if not boundedOnly then
    events2D := events2D[..-2];
  fi;
  ComputeSamplePoints2D(Q2D, events2D, 1, upperbound(events2D) - 1, variables, midpoint, db,
                                                                               boundedOnly);
end proc:

end module;
