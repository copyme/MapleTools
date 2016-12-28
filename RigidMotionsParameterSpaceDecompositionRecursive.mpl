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
  local  ComputeEventsAType2D, ComputeEventsBType2D, ComputeEventsAType1D, IsAsymptotic2D,
         ComputeAsymptoticAAEvents2D, ComputeEventsAlgebraicNumbers2D, ComputeSamplePoints2D;
         
  export LaunchComputeSamplePoints2D;

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
  local s:
  if nops(vars2D) < 2 then
    error "Only systems in at least two variables are supported. The variables are: %1.", vars2D;
  fi:
  s := proc(i::integer, vars2D::list)
    local sys, univ, sol:
    local q := Q2D2[i];
    sys := { q, diff( q, vars2D[2] ) };
    univ := EliminationResultant(sys, vars2D):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, vars2D[1..1]):
      if nops(sol) > 0 then
        return [univ, [i]];
      else
        return NULL;
      fi;
    end if;
  end proc;
  if grid then
    return [Grid:-Seq( s( i, vars2D ), i=1..nops( Q2D2 ) )]:
  else 
    return [seq( s( i, vars2D ), i=1..nops( Q2D2 ) )]:
  fi;
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
  local s:
  s := proc (i, j, vars2D::list)
    local p, sol, univ, sys;
    sys := {Q2D2[i], Q2D2[j]}:
    univ := EliminationResultant(sys, vars2D):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, vars2D[1..1]):
      if nops(sol) > 0 then
        return [ univ, [i,j] ]:
      fi:
    fi;
    return NULL:
   end proc:
   if grid then
     return [Grid:-Seq( seq( s(i, j, vars2D), j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]:
   else
     return [seq( seq( s(i, j, vars2D), j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]:
   fi;
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
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#
# Output: Real
#   A list of real algebraic numbers and indexes of conics: [ RealAlgebraicNumber, [index]].
ComputeAsymptoticAAEvents2D := proc(Q2D2, vars2D::list, grid::boolean)
  local out := [], s, factored, sqrFree;
  s:=proc(i::integer, vars2D::list)
    local rf, rootsF;
    local numbers := [], asy:
     asy := IsAsymptotic2D(Q2D2[i], vars2D[-1]);
     if not type(asy, constant) then
       factored := factors(asy)[2,..,1];
       for sqrFree in factored do
         rootsF := RootFinding:-Isolate(sqrFree, vars2D[1], output='interval');
         for rf in rootsF do
           numbers:=[op(numbers), [Object(RealAlgebraicNumber, sqrFree, op(rf)[2][1],
           op(rf)[2][2]), [i]]];
         od;
       od;
     fi:
    return numbers;
  end proc:
  if grid then
    out:=select(proc(x) return evalb(x<>[]) end, [Grid:-Seq(s(i, vars2D),i=1..nops(Q2D2))]);
  else
    out:=select(proc(x) return evalb(x<>[]) end, [seq(s(i, vars2D),i=1..nops(Q2D2))]);
  fi;
  return ListTools:-Flatten(out, 1);
end:


# Procedure: ComputeEventsAlgebraicNumbers2D
#   Compute and sort events as algebraic numbers 
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
  local events, rootsF, rf, poly;
  local numbers := Array([]);
  local factored, sqrFree;

  events:= [op(ComputeEventsAType2D(Q2D2, grid, vars2D)), op(ComputeEventsBType2D(Q2D2, grid, 
            vars2D))];
  for poly in events do
    factored := factors(poly[1])[2,..,1];
    for sqrFree in factored do
      rootsF := RootFinding:-Isolate(sqrFree, output='interval'):
      for rf in rootsF do
        ArrayTools:-Append(numbers, [Object(RealAlgebraicNumber, sqrFree, op(rf)[2][1],
        op(rf)[2][2] ), poly[2]]);
      od;
    od;
  od;
  ArrayTools:-Concatenate(2, numbers, Vector[row]([ComputeAsymptoticAAEvents2D(Q2D2, vars2D)]));
  return SortEvents(numbers);
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
          ArrayTools:-Append(numbers, Object(RealAlgebraicNumber, sqrFree, op(rf)[2][1],
                             op(rf)[2][2]));
        od;
      od;
   fi;
  od;
    numbers := convert(numbers, list); 
    numbers := SortAlgebraicNumbers(numbers);
    return ListTools:-MakeUnique(numbers, 1, proc(a,b) evalb(Compare(a,b) = 0) end proc);
end proc:


# Procedure: ComputeSamplePoints2D
#   Computes sample points for rotational part of rigid motions
#
# Parameters:
#   Q2D         - a list of conics
#   cluster2D   - each element contains a list of equal real algebraic number and related
#                 conics.
#   first       - integer value which indicates a first cluster to proceed.
#   last        - integer value which indicates a last cluster to proceed.
#   vars2D      - list of variables in which conics are expressed
#   db          - an instance of the class ComputationRegister
#
# Output:
#   It populates a database, given by databasePath, with sample points.
ComputeSamplePoints2D := proc(Q2D, cluster2D::list, first::integer, last::integer,
                              vars2D::list, amid::rational, db::ComputationRegister)
  local i::integer, j::integer, x::list, midpoint::rational, sys::list, samplePoints::list;
  local disjointEvent::list, oneD::list, oneDNeg::list;
  if first < 0 or last < 0 or last < first or upperbound(cluster2D) <= last then 
    error "Bounds of the cluster range are incorrect.": 
  end if:
  for i from first to last do 
    sys := []: 
    for x in cluster2D[i] do 
      sys := [op(sys), op(Q2D[x[2]])]:
    end do:
    print("in 2");
    sys := ListTools:-MakeUnique(sys);
    disjointEvent := DisjointRanges(cluster2D[i][1][1],cluster2D[i+1][1][1]);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
   
    # intersection of a line with  conics
    # never call eval with sets!
    print("sys 2", sys);
    sys := eval(sys, vars2D[1] = midpoint):
    oneD := ComputeEventsAType1D(sys);
    if oneD = NULL then
      next;
    fi:

    oneDNeg, oneD := selectremove(proc(x) return evalb(GetInterval(x)[2] < 0); end proc, oneD):
    if nops(oneDNeg) <> 0 then
      oneD := [oneDNeg[-1], op(oneD)];
    fi:

    if upperbound(oneD) > 0 then
      samplePoints := [];
      for j from 1 to upperbound(oneD) - 1 do
        disjointEvent := DisjointRanges(oneD[j],oneD[j+1]);
        samplePoints := [op(samplePoints), [amid, midpoint, (GetInterval(disjointEvent[1])[2] +
                                            GetInterval(disjointEvent[2])[1])/2]];
      od:
      samplePoints := [op(samplePoints), [amid, midpoint, GetInterval(oneD[-1])[2] + 1/2]];
      map(proc(x) InsertSamplePoint(db, x) end proc, samplePoints);
    fi:
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
                                           grid::boolean, variables::list, db::ComputationRegister) 
  local events, firstEvent, rootTmp;
  local Q2D := ListTools:-MakeUnique([op(variables), op(s)]), cluster2D;
  if grid and nops(s) > 20 then
     events := convert(ComputeEventsAlgebraicNumbers2D(Q2D, true, variables), list);
  else
     events := convert(ComputeEventsAlgebraicNumbers2D(Q2D, false, variables), list);
  fi;
  events := remove(proc(x) return evalb(GetInterval(x[1])[2] < 0); end proc, events);
  if upperbound(events) = 0 then
    return NULL;
  fi;
  cluster2D := ClusterEvents(events);

  # assign all conics to the first event
  cluster2D := [[[cluster2D[1][1][1], [seq(1..nops(Q2D))]]], op(cluster2D[2..])]:
  rootTmp:= GetInterval(cluster2D[-1][1][1])[2]+1;
  firstEvent := Object(RealAlgebraicNumber, denom(rootTmp)*variables[1]-numer(rootTmp), rootTmp, 
                       rootTmp);
  cluster2D := [op(cluster2D), [[firstEvent, cluster2D[-1][1][2]]]];
  print("in 1");

  ComputeSamplePoints2D(Q2D, cluster2D, 1, nops(cluster2D) - 1, variables, midpoint, db);
end proc:

end module;
