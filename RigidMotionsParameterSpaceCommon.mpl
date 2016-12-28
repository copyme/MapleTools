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
  export CayleyTransform, GetNeighborhood, EliminationResultant, RemoveExponants,
  OneVariableElimination, ClusterEvents, SortAlgebraicNumbers, SortEvents;


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
#   nType      - size of neighborhood i.e. N1, N2, N3, N1_2D. 
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

# Procedure: EliminationResultant
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials in at least two variables
#   vars       - variables to be eliminated
#
# Output:
#   Univariate polynomial obtained from S in the first variable.
EliminationResultant := proc( S::~set, vars::~list )
  option cache:
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
      res := [op(res), OneVariableElimination(p[1], p[2], vars[1])]:
    od;
    return foldl( gcd, op(res) ):
  fi;

   for var in vars[2..] do
       permm := combinat:-permute(res, 2);
       res := [];
     for p in permm do
       res := [op(res), OneVariableElimination(p[1], p[2], var)]:
     od;
   od;

  return foldl( gcd, op(res) ):
end proc:

# Procedure: RemoveExponants
#    Removes exponants in an expression
#
# Parameters:
#    r - expression, the expression to simplify
#
# Output:
#    An arithmetic expression that has the same squarefree part as r
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


# Procedure: ClusterEvents
#   Compute sublists which contains equal real algebraic numbers
#
# Comments:
#   Used to split calculations over a grid of nodes for parallel computations.
#
# Parameters:
#   nType            - list of real algebraic numbers and quadrics related eg numbers[i] =
#   [RealAlgebraicNumber,[Quadrics involved]]
#
# Output:
#   A list of sublists where each of them contains equal real algebraic numbers
ClusterEvents := proc(numbers::list)
  return SplitScan(
            proc (x, y) 
              if not type(x[1], RealAlgebraicNumber) or not type(y[1], RealAlgebraicNumber) then
                error "Elements are not of type RealAlgebraicNumber" 
              end if;
            return evalb(Compare(x[1], y[1]) <> 0) 
            end proc, numbers)
end proc:


SortAlgebraicNumbers := proc (numbers::list)
  return sort( numbers, 
        proc( l, r ) 
        if Compare( l, r ) = -1 then
          return true:
        else 
          return false:
        fi;
        end proc);
end proc;


SortEvents :=proc(events)
  # In maple 2015.2 there is a bug which causes: stack limit reached if sorting an empty Array
  if not StringTools:-Has(kernelopts(version), "Maple 2016") and upperbound(events) <> 0 then
      events := sort(events, 
                           proc( l, r ) 
                             if Compare( l[1], r[1] ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  );
  elif StringTools:-Has(kernelopts(version), "Maple 2016") then
      sort['inplace'](events, 
                           proc( l, r ) 
                             if Compare( l[1], r[1] ) = -1 then
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
