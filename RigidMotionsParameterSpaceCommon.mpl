# File: RigidMotionsParameterSpaceCommon.mpl  
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
#  Guillaume Moroz - guillaume.moroz@inria.fr 
#  INRIA Nancy, France
#
# Date:
#  11/12/2015 
#
# License:
#  Simplified BSD License
#
# Copyright (c) 2015, Kacper Pluta, Guillaume Moroz
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
# DISCLAIMED. IN NO EVENT SHALL Kacper Pluta and Guillaume Moroz BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#

RigidMotionsParameterSpaceCommon := module() 
  option package;
  uses   RigidMotionsMaplePrimesCode;
  local EliminationResultant, RemoveExponants, OneVariableElimination, EliminationGroebner;
  export CayleyTransform, GetNeighborhood, AlgebraicSort, ReduceEvents, AdjustEvents, GenerateEvents,
  SerializeEvents, ReconstructEvents, UnivariatePolynomial;


# Procedure: CayleyTransform
#   Compute Cayley transform for a 3x3 skew-symmetric matrix.
#
# Parameters:
#   vars   - list of variables
#
# Output:
#   3x3 (or 2x2) rotation matrix
# 
# Links:
#   https://en.wikipedia.org/wiki/Cayley_transform
CayleyTransform := proc( vars::list )
  local A::Matrix, QLSide::Matrix, QRSide::Matrix, dim:
  dim := nops(vars);
  if dim = 1 then
   A := Matrix( [ [ 0, vars[1] ], [ -vars[1], 0 ] ] ):
  elif dim = 2 then
   A := Matrix( [ [ 0, vars[2]/vars[1] ], [ -vars[2]/vars[1], 0 ] ] );
  elif dim = 3 then
   A := Matrix( [ [ 0, vars[1], vars[2] ], [ -vars[1], 0, vars[3] ], 
                                      [ -vars[2], -vars[3], 0 ] ] ):
  else
   error "Unsupported dimension! Check 1, 2 or 3":
  end if:
  dim := upperbound(A)[1];
  QLSide := Matrix( dim, shape = identity ) - A:
  QRSide := LinearAlgebra:-MatrixInverse( Matrix( dim, shape = identity ) + A ):
  return simplify( QLSide . QRSide ):
end proc:

# Procedure: GetNeighborhood
#   Compute a neighborhood
#
# Parameters:
#   nType      - size of neighborhood i.e. N1, N2, N3. 
#
# Output:
#   List of vectors
GetNeighborhood := proc( nType::string )
  local n6 := [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1], [0, 0, 0]]:
  local n18 := [[1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, -1, 0], [-1, 0, -1], [0, -1, -1], 
                [-1, 1, 0], [-1, 0, 1], [1, -1, 0], [1, 0, -1], [0, -1, 1], [0, 1, -1]]:
  local n26 := [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], [-1, 1, -1], 
                [-1, -1, -1], [1, -1, -1]]:

  if nType = "N1" then
    return n6:
  elif nType = "N2" then
    return( [ op( n6 ), op( n18 ) ] ):
  elif nType = "N3" then
    return( [ op( n6 ), op( n18 ), op(n26) ] ):
  else 
    error "Not supported type. Try N1, N2, N3.":
  end if:
end proc:


# Procedure: EliminationGroebner
#   Computes univariate polynomial using Groebner basis and FGb library.
#
# Parameters:
#   S          - a set of quadrics
#   vars       - a list of variables
#
# Comments:
#   To use this procedure you need to install Jean-Charles Faugère's FGb library
#   --http://www-polsys.lip6.fr/~jcf/FGb/FGb/index.html 
#
# Output:
#   Univariate polynomial obtained from S.
#
EliminationGroebner := proc(S::list, vars::list)
  option cache;
  local univ := op(FGb:-fgb_gbasis_elim(S, 0, vars[2..], vars[1]));
  if univ = NULL then
    univ := 0;
  fi;
  return univ;
end proc;

# Procedure: EliminationResultant
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials in at least two variables
#   vars       - list of all the variables where the output univariate polynomial is represented in
#                the first one -- the first element of the list
#
# Output:
#   Univariate polynomial obtained from S in the first variable.
EliminationResultant := proc( S::list, vars::list )
  option cache;
  local p, var, res := S, permm;
  if nops(S) < 2 then
    error "Wrong size of the input set: %1. Expected size is at least 2.", S;
  fi;
  if not andmap(type, S, `polynom`) then
    error "Wrong type of elements. Expected argument is a set of polynomials but received  %1."; S; 
  fi;
  if nops(vars) < 1 then
    error "Wrong number of indeterminates. It should be at least 2.";
  fi;
  if nops(vars) = 1 then
    permm := combinat:-permute(res, 2);
    res := [];
    for p in permm do
      res := [op(res), OneVariableElimination(p[1], p[2], vars[1])];
    od;
    return foldl( gcd, op(res) );
  fi;

   for var in vars[2..] do
       permm := combinat:-permute(res, 2);
       res := [];
     for p in permm do
       res := [op(res), OneVariableElimination(p[1], p[2], var)];
     od;
   od;

  return foldl( gcd, op(res) );
end proc:

# Procedure: RemoveExponants
#    Removes exponants in an expression
#
# Parameters:
#    r - expression, the expression to simplify
#
# Output:
#    An arithmetic expression that has the same square free part as r.
RemoveExponants := proc(r)
        local remove_exponant, sqrr, result;
        remove_exponant := e -> if type(e,`^`) then op(1,e) else e end if;
        result := remove_exponant(r);
        if type(r,`*`) then
            result := map(remove_exponant, result);
        end if;
        return result;
end proc;

OneVariableElimination := proc( p, q, v)
    local r;
    if degree(p,v)>0 or degree(q,v)>0 then
        r := resultant(p, q, v); 
        r := RemoveExponants(r);
        return r;
    elif nops(indets({p,q}))=1 then
        return gcd(p, q);
    else
        return p;
    end if;
end proc;


# Procedure: UnivariatePolynomial
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of quadrics
#   vars       - a list of variables
#
# Comments:
#   If FGb is installed then fgb_gbasis_elim() is used and EliminationResultant, otherwise.
#   Note that, build-in solutions like Groebner:-UnivariatePolynomial are not used because of
#   potential memory explosion or other problems which make them practically useless.
#
# Output:
#   Univariate polynomial obtained from S.
#
UnivariatePolynomial := proc(S::set, vars::list)
  if type(FGb, package) then
    return EliminationGroebner([op(S)], vars);
  fi;
    return EliminationResultant([op(S)], vars);
end proc;


# Procedure: GenerateEvents
#   Generates events from a univariate polynomial (see EventType)
#
# Parameters:
#    x::polynom           - a univariate polynomial
#    quads::              - a list of quadrics
#    eventType::string    - Type of events: A, B, C, AA, AB
#
# Output:
#    An Array of events generated from a given univariate polynomial.
GenerateEvents := proc(x::polynom, quads::list, eventType::set)
  uses ArrayTools;
  local events := Array([]), factored, rootsF, sqrFree, rf;
  if nops(indets(x)) = 1 then
    factored := factors(x)[2,..,1];
    for sqrFree in factored do
      rootsF := RootFinding:-Isolate(sqrFree, output='interval');
      for rf in rootsF do
        Append(events, EventType(RealAlgebraicNumber(sqrFree, op(rf)[2][1], op(rf)[2][2]), quads,
        eventType));
      od;
    od;
  fi;
  return events;
end proc;


# Procedure: SerializeEvents
#   Changes (inplace) an Array of events into an Array of strings of unevaluated calls to
#   the constructor of EventType.
#
# Parameters:
#    events::Array     - an Array of events (see EventType)
#
# Output:
#    An Array of strings to unevaluated calls of the constructor of EventType.
SerializeEvents := proc(events::Array)
  return map[inplace](proc(x) sprintf("%a", x) end proc, events);
end proc;


# Procedure: ReconstructEvents
#   Changes (inplace) an Array of strings of unevaluated calls to the constructor of EventType into
#   the proper objects.
#
# Parameters:
#    events::Array     - an Array of events (see SerializeEvents)
#
# Output:
#    An Array of EventType.
ReconstructEvents := proc(events::Array)
  return map[inplace](proc(x) eval(parse(x)) end proc, events);
end proc;


# Procedure: ReduceEvents
#   Returns an array of events such that they are distinct. Each pair of events which are equal are
#   merged in such a way that a list of quadrics of the second on from a pair is merged with the
#   list of the first event.
#
# Parameters:
#    L::Array     - an Array of events
#
# Output:
#    An Array of EventType.
ReduceEvents := proc(L::Array) 
  local R, k, j, last, x, quadrics := [], eventTypes := {};
  R := Array([]); k := 0; last := 1; 
  for j from 2 to upperbound(L) do 
    if Compare(L[j-1], L[j], _rest) <> 0 then 
      k := k+1; 
      quadrics := ListTools:-MakeUnique([seq(op(GetQuadrics(x)), x=L[last .. j-1])]);
      eventTypes := {seq(op(GetEventTypes(x)), x=L[last .. j-1])};
      R(k) := EventType(GetRealAlgebraicNumber(L[last]), quadrics, eventTypes);
      last := j;
    end if;
  end do;
  if last <> j then
    quadrics := ListTools:-MakeUnique([seq(op(GetQuadrics(x)), x=L[last .. ()])]);
    eventTypes := {seq(op(GetEventTypes(x)), x=L[last .. ()])};
    ArrayTools:-Append(R, EventType(GetRealAlgebraicNumber(L[last]), quadrics, eventTypes));
  fi;
  return R;
end proc;


# Procedure: AdjustEvents
#   Returns an array of events such that the first event will contains all the quadris and the last
#   event will be duplicated with changed interval of the underlaying real algebraic number.
#
# Parameters:
#    L::Array     - an Array of events
#
# Output:
#    An Array of EventType.
AdjustEvents := proc(events::Array, quadNum::integer, variables::list)
  local boundTmp, lastEvent;
  # assign all quadrics to the first even
  events[1] := EventType(GetRealAlgebraicNumber(events[1]), [seq(1..quadNum)], {"A", "B", "C", "AA", "AB"});
  # add the last slice twice but shifted to calculate correctly last quadrics
  boundTmp:= GetInterval(GetRealAlgebraicNumber(events[-1]))[2]+1;
  lastEvent := RealAlgebraicNumber(denom(boundTmp)*variables[1]-numer(boundTmp), boundTmp, boundTmp);
  ArrayTools:-Append(events, EventType(lastEvent, GetQuadrics(events[-1]), GetEventTypes(events[-1])), inplace=true)
end proc;


# Procedure: AlgebraicSort
#   Sorts [if possible in place] RealAlgebraicNumbers or Events. (see types: RealAlgebraicNumber 
#   and EventType).
#
# Parameters:
#   events  - a list or Array of RealAlgebraicNumbers or Events
#
# Output:
#   A increasingly sorted list or Array.
AlgebraicSort :=proc(events)
  # In Maple  <2016 there is a bug which causes: stack limit reached if sorting an empty Array
  if upperbound(events) <> 0 then
     return sort['inplace'](events, 
                           proc( l, r ) 
                             if Compare( l, r ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  );
  fi;
  return events;
end proc;

end module;

