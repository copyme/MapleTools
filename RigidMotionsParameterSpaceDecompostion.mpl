# File: RigidMotionsParameterSpaceDecompostion.mpl  
#
# Description:
#  This file contains functions used to obtain an arrangement 6 dimensional parameter space of 3D
#  digitized rigid motions.
#  This code has been written for research propose and its aim is to calculate a particular
#  arrangement of quadrics. Therefore, it can or it cannot be useful in study of generic
#  arrangements. The final output are sample points of full dimensional open cells.
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
  export LaunchOnGridComputeSamplePoints, LaunchOnGridGetNMM, CalculateNMM,
    CayleyTransform, GetQuadric, GetNeighborhood, EliminationResultant, IsMonotonic, 
         ComputeSetOfQuadrics, IsAsymptotic, IsAsymptoticIntersection, ComputeEventsATypeGrid,
         ComputeEventsBTypeGrid, ComputeEventsCTypeGrid, ComputeEventsAlgebraicNumbers, SplitScan,
         ClusterEvents, ComputeSamplePoints, ParallelComputeSamplePoints, Isort, 
         GetOrderedCriticalPlanes, RecoverTranslationSamplePoints, Get3DNMM,
         ParallelCalculateNMM:


# Procedure: CayleyTransform
#   Compute Cayley transform for a 3x3 skew-symmetric matrix.
#
# Parameters:
#   vars   - set of variables
#
# Output:
#   3x3 (or 2x2) rotation matrix
# 
# Links:
#   https://en.wikipedia.org/wiki/Cayley_transform
CayleyTransform := proc( vars::~set )
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


# Procedure: GetQuadric
#   Compute a quadric
#
# Parameters:
#   R          - a rotation matrix obtained from Cayley transform
#   neighbor   - a vector which points on a neighbor
#   hGridPlane - a half-grid plane i.e. k_dim + 1/2
#   axis       - a given axis i.e. 1 = x, 2 = y and 3 = z 
#   
#
# Output:
#   Multivariate polynomial
GetQuadric := proc( R::~Matrix,
                     neighbor::~Vector,
                     hGridPlane::~Vector,
                     axis::integer )
  local r := R . neighbor - hGridPlane;
  if axis > 0 and axis < 4 then
    return simplify( ( 2 * denom( R[1][1] ) ) * r[axis] );
  else
    error "Wrong dimension: dim must be in [1,3].";
  end if:
end proc:


# Procedure: EliminationResultant
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials in three variables
#
# Output:
#   Univariate polynomial obtained from S in the first variable.
EliminationResultant := proc( S::~set )
  option cache;
  local r1, r2, r3, rr1, rr2, rr3, vars;
  vars := indets(S);
  if nops(S) <> 3 then
    error "Wrong size of the input set. Expected size is 3.";
  fi;
  if not type(S[1], polynom) or not type(S[2], polynom) or not type(S[3], polynom) then
    error "Wrong type of elements. Expected argument is a set of polynomials!"; 
  fi;
  if nops(vars) <> 3 then
    error "Wrong number of indeterminates. It should be 3.";
  fi;
  r1 := OneVariableElimination( S[1], S[2], vars[3] ):
  r2 := OneVariableElimination( S[1], S[3], vars[3] ):
  r3 := OneVariableElimination( S[2], S[3], vars[3] ):
  rr1 := OneVariableElimination( r1, r2, vars[2] ):
  rr2 := OneVariableElimination( r1, r3, vars[2] ):
  rr3 := OneVariableElimination( r2, r3, vars[2] ):
  return foldl( gcd, rr1, rr2, rr3 ):
end proc:

# Procedure: RemoveExponant
#    Removes exponants in an expression
#
# Parameters:
#    r - expression, the expression to simplify
#
# Output:
#    An arithmetic expression that has the same squarefree part as r
RemoveExponants := proc(r)
        local remove_exponant, sqrr;
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
    elif  nops(indets({p,q}))=1 then
        return gcd(p, q);
    else
        return p;
    end if;
end proc;

# Procedure: EliminationResultant2
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials of degree 2 in three variables
#
# Output:
#   Univariate polynomial obtained from S in the first variable, with
#   formula from Chapter 3, p. 89 of "Using Algebraic Geometry" from Cox,
#   Little, O'Shea.
EliminationResultant2 := proc( S::~set )
  #option cache;
  local vars, monomials, L, J, J1, J2, dJ, result;
  vars := indets(S);
  if nops(S) <> 3 then
    error "Wrong size of the input set. Expected size is 3.";
  fi;
  if not type(S[1], polynom) or not type(S[2], polynom) or not type(S[3], polynom) then
    error "Wrong type of elements. Expected argument is a set of polynomials!"; 
  fi;
  if nops(vars) <> 3 then
    error "Wrong number of indeterminates. It should be 3.";
  fi;
  vars := [op(vars)][2..3];
  monomials := map2(map,`-`,[op(combinat:-composition(5,3))], 1)[..,2..3];
  L := map(p->[diff(p,vars[1]), diff(p,vars[2]),p], [op(S)]);
  J := LinearAlgebra:-Determinant(L);
  J1 := diff(J,vars[1]);
  J2 := diff(J,vars[2]);
  dJ := [J1, J2, degree(J,vars)*J-vars[1]*J1 - vars[2]*J2];
  L := map2(map2,(p,d)-> coeftayl(p, vars=[0,0], d),
                 [op(S),op(dJ)], monomials);
  result := -LinearAlgebra:-Determinant(L)/64;
  return result;
end proc:

# Procedure: IsMonotonic
#   Check if given polynomial is strictly positive or negative
#
# Parameters:
#   x      - a polynomial
#
# Output:
#   true if polynomial is positive or negative and false otherwise.
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
#
# Output:
#   List of solution for which partial derivatives in b and c are collinear
IsAsymptotic := proc(x::polynom)
  local vars, vec, Vb, Vc, VV, sols;
  vars := [op(indets(x))];
  vec := VectorCalculus:-Gradient(x, vars);
  Vb := Vector(3, [coeff(vec[2],vars[2]),coeff(vec[2],vars[3]),eval(vec[2],[vars[-2]=0,vars[-1]=0])]);
  Vc := Vector(3, [coeff(vec[3],vars[2]),coeff(vec[3],vars[3]),eval(vec[3],[vars[-2]=0,vars[-1]=0])]);
  VV := LinearAlgebra:-CrossProduct(Vb, Vc);
  if norm(VV,1) = 0 then
    return {a=0};
  fi;
  sols := solve({VV[1] = 0, VV[2] = 0,VV[3] = 0,vars[1]>=0});
  if sols = NULL or sols = {} then
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
IsAsymptoticIntersection := proc( p::polynom, q::polynom  )
  local J := PolynomialIdeals:-`<,>`(p,q), vars := indets([p,q]);
  local Pb, Pc, Cb, Cc, sols, univ:

  Pb := PolynomialIdeals:-EliminationIdeal(J,vars[1..2]):
  Pb := PolynomialIdeals:-IdealInfo:-Generators(Pb)[1]:
  Pc := PolynomialIdeals:-EliminationIdeal(J,vars[[1,3]]):
  Pc := PolynomialIdeals:-IdealInfo:-Generators(Pc)[1]:
  Cb := lcoeff(Pb, vars[2]):
  Cc := lcoeff(Pc, vars[3]):
  return gcd(Cb,Cc);
end proc:


# Procedure: ComputeEventsATypeGrid
#   Compute events such that a sweep plane is tangent to a quadric. 
#
# Parameters:
#   Q          - a set of quadrics or conics
#   dim        - a list of indexes of variables used to calculate partial derivatives
#
# Output:
#   List of ranges which contains roots of a system(q, d/db q, d/dc q).
#
# Comment:
#  - only the first direction is supported since EliminationResultant is used
ComputeEventsATypeGrid := proc( Q, dim::list )
  local s:
   s := proc(i::integer)
    local sys, univ, sol, vars:
    local q := Q[i];
    vars := [ op( indets( q ) ) ];
    if nops(vars) = 3 then
      sys := { q, diff( q, vars[ dim[1] ] ), diff( q, vars[ dim[2] ] ) };
    elif nops(vars) = 2 then 
      sys := [ q, diff( q, vars[ dim[1] ] ) ];
    fi:
    univ := RigidMotionsParameterSpaceDecompostion:-EliminationResultant(sys):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, [ op( indets(univ ) ) ]):
      sol := nops(select(e -> rhs(e) >= 0, sol)):
      if sol > 0 then
        return [univ, [i]]:
      else
        return NULL:
      fi:
    end if:
  end proc:
  return [Grid:-Seq(s(i),i=1..nops(Q))]:
end proc:


# Procedure: ComputeEventsBType
#   Compute events such that intersection of two quadrics is tangent to sweeping plane.
#
# Parameters:
#   dir        - a direction of a gradient product it should
#                be the same as director of sweep 
#   Q          - a set of quadrics or conics
#
# Output:
#   Indexes of quadrics which intersect and a component of a vector product of 
#   their gradients in given direction have a common root.
ComputeEventsBTypeGrid := proc( Q, dir::integer )
  local s:
  s := proc(i, j)
      local p, prod, ivars, jvars, univ, sys, vars, sol:
      if i + j mod 100 = 0 then
        #(*Clear remember tables*)
        seq(forget(p, forgetpermanent = true), p in {anames('procedure')}):
      fi:
      vars := indets ( [ Q[i], Q[j] ] ):
      if nops ( vars ) = 3 then
        ivars := [ op( indets( [ Q[i] ] ) ) ]:
        jvars := [ op( indets( [ Q[j] ] ) ) ]:
        prod := LinearAlgebra:-CrossProduct( VectorCalculus:-Gradient( Q[i], ivars ),
                                    VectorCalculus:-Gradient( Q[j], jvars ) )[dir]:
        sys := { Q[i], Q[j], prod }:
      elif nops(vars) = 2 then
        sys := [ Q[i], Q[j] ]:
      fi:
      univ := RigidMotionsParameterSpaceDecompostion:-EliminationResultant(sys):
      if not type( univ, constant ) then
        sol := RootFinding:-Isolate( univ, [ op( indets(univ ) ) ]):
        sol := nops(select(e -> rhs(e) >= 0, sol)):
       if sol > 0 then
         return [univ, [i,j]]:
       else
         return NULL:
       fi;
      fi:
  end proc:
  return [Grid:-Seq(seq(s(i,j),j=i+1..nops(Q)),i=1..nops(Q))]:
end proc:


# Procedure: ComputeEventsCTypeGrid
#   Compute events such that three quadrics intersects in a point. 
#
# Parameters:
#   Q          - a set of quadrics
#
# Output:
#   Indexes of quadrics which intersect in a point.
ComputeEventsCTypeGrid := proc( Q )
  local s:
  s := proc (i, j, k)
    local p, sol, univ, sys;
    if i + j + k mod 100 = 0 then
      #(*Clear remember tables*)
      seq(forget(p, forgetpermanent = true), p in {anames('procedure')}):
    fi:
    sys := { Q[i], Q[j], Q[k] }:
    univ := RigidMotionsParameterSpaceDecompostion:-EliminationResultant(sys):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, [ op( indets(univ ) ) ]):
      sol := nops(select(e -> rhs(e) >= 0, sol)):
      if sol > 0 then
        return [ univ, [i,j,k] ]:
      fi:
    fi;
    return NULL:
   end proc:
   return [Grid:-Seq(seq(seq(s(i,j,k),k=j+1..nops(Q)),j=i+1..nops(Q)),i=1..nops(Q))]:
end proc:


# Procedure: ComputeEventsAlgebraicNumbers
#   Compute and sort events as algebraic numbers 
#
# Parameters:
#   Q     - set of quadrics or conics
# Output:
#   Sorted set of real algebraic numbers
#
# TODO:
#  - Rewrite it in a way which will allow to apply recursive strategy
ComputeEventsAlgebraicNumbers := proc( Q::~set )
  local events, rootsF, rf, poly, i, j, s, sol:
  local numbers := Array([]):
  local factored, sqrFree:

    events:= {op(ComputeEventsATypeGrid( Q, [2, 3] )), op(ComputeEventsBTypeGrid( Q, 1 )),
                                                         op(ComputeEventsCTypeGrid( Q ))}:
  for poly in events do
    factored := factors( poly[1] )[2,..,1]: 
    for sqrFree in factored do
      rootsF := RootFinding:-Isolate( sqrFree, output='interval' ):
      for rf in rootsF do
        ArrayTools:-Append(numbers, [ Object( RealAlgebraicNumber, sqrFree, op(rf)[2][1],
        op(rf)[2][2] ), poly[2]]):
      od:
    od:
  od:

  # asymptotic case occurred for one quadric
  for i from 1 to nops(Q) do
    rootsF := IsAsymptotic(Q[i]):
    for sol in rootsF do
      rootsF := rhs(sol):
      if not type(rootsF, rational) then
        error "Irrational asymptotic case! Are you sure the input is a set of quadrics?"
      fi:
      ArrayTools:-Append(numbers,[Object(RealAlgebraicNumber, lhs(sol) * denom(rootsF) -
      numer(rootsF), rootsF, rootsF), [i]]): 
    od:
  od:
  print("asymptotic case AA computed");

  for i from 1 to nops(Q) do
    for j from i + 1 to nops(Q) do
      poly := IsAsymptoticIntersection(Q[i], Q[j]):
      if poly = NULL then
        next;
      fi:
      factored := sqrfree( factor( poly ) )[2,..,1]: 
      for sqrFree in factored do
        rootsF := RootFinding:-Isolate( sqrFree, op( indets( sqrFree ) ), output='interval'):
        for rf in rootsF do
          ArrayTools:-Append(numbers, [ Object( RealAlgebraicNumber, sqrFree, op(rf)[2][1],
          op(rf)[2][2] ), [i,j]]): 
        od:
      od:
    od:
  od:
  print("asymptotic case AB computed");

  numbers := sort( numbers, 
                           proc( l, r ) 
                             if Compare( l[1], r[1] ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  ):
  print("events sorted");
  return numbers:
end proc:


# Procedure: SplitScan
#   Can split a list into sublists of the same elements. Code written by Carl Love.
#   See http://www.mapleprimes.com/questions/205030-Split-The-List-Into-Sublists-With-Identical
#
# Comments:
#   The input list has to be sorted.
#
# Parameters:
#   f            - a procedure used to check if two elements are different eg `<>`
#   L            - a list to split
#
# Output:
#   Returns a list of sublists where each of them contains elements which are same.
# Call example:
#   SplitList(`<>`,[1, 1, 2, 3, 3, 4, 9, 9, 0, 11]);
SplitScan := proc(f, L::list) 
  local R, k, j, last; R := Vector(); k := 0; last := 1; 
  for j from 2 to nops(L) do 
    if f(L[j-1], L[j], _rest) then 
      k := k+1; R(k) := L[last .. j-1]; 
      last := j 
    end if 
    end do;
  [seq(k, k = R), L[last .. ()]]
end proc:


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


# Procedure: ComputeSamplePoints
#   Computes sample points for rotational part of rigid motions
#
#
# Parameters:
#   cluster            - each element contains a list of equal real algebraic number and quadrics
#                        related.
#   first              - integer value which indicates a first cluster to proceed.
#   last               - integer value which indicates a last cluster to proceed.
#   id                 - id which indicates a node
#
# Output:
#   Writes a list of sample points into a file "sam_id.csv"
#
# TODO:
#  Split and allow recursive calls
ComputeSamplePoints := proc (Q::~set, cluster::list, first::integer, last::integer, id::integer)
  local i, x, midpoint, sys, samplePoints, fileID, vars:
  if first < 0 or last < 0 or last < first or upperbound(cluster) < last then 
    error "Bounds of cluster rangers are incorrect": 
  end if:
  for i from first to last do 
    if i mod 100 = 0 then
      #(*Clear remember tables*)
      seq(forget(p, forgetpermanent = true), p in {anames('procedure')}):
    fi:
    sys := {}: 
    for x in cluster[i] do 
      sys := sys union Q[x[2]]:
    end do:
    vars := indets(sys):
    midpoint := (GetInterval(cluster[i][1][1])[2] + GetInterval(cluster[i+1][1][1])[1])/2:
    sys := Threads:-Map(proc (x) return eval(x, vars[1] = midpoint) <> 0 end proc, sys):
    # Some how negative sample points are still generated.
    sys := [op(sys), 0 < vars[2], 0 < vars[3]]:
    # this has to be replaced by recursive approach
    samplePoints := RootFinding[Parametric]:-CellDecomposition(sys, [], [op(vars[2..3])],
                                      'output'='witnesspoints'):-SamplePoints:
    
    samplePoints := remove(proc(x) if rhs(x[1]) < 0 or rhs(x[2]) < 0 then return true: else return
                                                                     false: fi: end, samplePoints):

    fileID := fopen(sprintf("sam_%d.csv", id), APPEND, TEXT):
    writedata(fileID, Threads:-Map(proc (x) return [midpoint, 
               rhs(x[1]), rhs(x[2])] end proc, samplePoints),
              string, proc (f, x) fprintf(f, %a, x) end proc):
    fclose(fileID): 
  end do:
  return NULL:
end proc:


# Procedure: ParallelComputeSamplePoints
#   Computes sample points for rotational part of rigid motions. It should be call via Grid
#   framework.
#
ParallelComputeSamplePoints := proc () 
  local me, numNodes, n;
  me := Grid:-MyNode();
  numNodes := Grid:-NumNodes();
  # cluster-1 because the last cluster is a doubled cluster[-2]
  n := trunc((upperbound(cluster)-1)/numNodes);
  ComputeSamplePoints(Q, cluster, me*n+1,(me+1)*n, me);
  Grid:-Barrier();
end proc:

# Procedure: ComputeSamplePoints
#   Computes sample points for rotational part of rigid motions using the grid framework
#
#
# Parameters:
#   nType      - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange     - a range of planes to consider
# Output:
#   Writes a list of sample points into a file "sam_id.csv" where id corresponds to an id of used
#   thread during computations.
LaunchOnGridComputeSamplePoints := proc (nType::string, kRange::list) 
  local numbers, events, R: 
  global Q, cluster,rootTmp:
  kernelopts(printbytes=false):
  R := CayleyTransform({a,b,c}):
  Q := ComputeSetOfQuadrics(R, nType, 1, kRange) union 
       ComputeSetOfQuadrics(R, nType, 2, kRange) union
       ComputeSetOfQuadrics(R, nType, 3, kRange):
  numbers := ComputeEventsAlgebraicNumbers(Q):
  numbers := convert(numbers,list):
  print("numbers converted");
  numbers := remove( proc(x) return evalb(GetInterval(x[1])[2] < 0); end proc, numbers):
  print("numbers removed");
  cluster := ClusterEvents(numbers):
  print("numbers clustered");
  # Collect all unique quadrics
  events := ListTools:-MakeUnique(Threads:-Map(op, [seq(op(cluster[i][..,2]), i = 1..nops(cluster))])):
  print("all quadrics found in cluster");
  # assign all quadrics to the second event
  cluster := [[[cluster[1][1][1], events]], op(cluster[2..])]:
  # add the last slice twice but shifted to calculate correctly last quadrics
  events := cluster[-1][1][1]:
  rootTmp:= GetInterval(events)[2]+1/4;
  events := Object(RealAlgebraicNumber, denom(rootTmp) * indets(GetPolynomial(events))[1] -
  numer(rootTmp), rootTmp, rootTmp):
  cluster := [op(cluster), [events,cluster[-1][1][2]]]:
  print("all quadrics added into the first event");
  print("cluster size ", upperbound(cluster));
  # The first cluster is heavy so we compute it separately;
  ComputeSamplePoints(Q,cluster[1..2],1,1, kernelopts(numcpus));
  cluster := cluster[2..-1];
  Grid:-Setup("local"):
  # We define printer as a procedure which returns NULL to avoid a memory leak problem while
  # writing to a file from a node. It seems that while fprintf is called it also calls printf
  # which is a default printer function. Therefore, data are returned to node of ID 0. 
  print("check point before run of grid");
  Grid:-Launch(ParallelComputeSamplePoints, imports = ['Q, cluster'],numnodes=kernelopts(numcpus),
  printer=proc(x) return NULL: end proc):
end proc:


# Procedure: Isort
#   Sort elements and return their indices in original order after sorting
#
# Author:
#   Alec Mihailovs - http://www.mapleprimes.com/posts/43507-Sorting-With-Indices#comment80473  
# Parameters:
#   L            - a container to be sorted. 
#
# Output:
#   indices of original order after sorting, sorted collection.
Isort:=proc(L)
  local a;
  a := sort([$1..nops(L)],(i,j)->L[i]<=L[j]); 
  return a, [seq(L[i],i=a)]; 
end:


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
GetOrderedCriticalPlanes := proc(nType::string, kRange::list, samplePoints::list) 
  local params; 
  local R := CayleyTransform({a,b,c}); 
  local n := GetNeighborhood(nType)[1 .. 6]; 
  local T := combinat:-cartprod([n, kRange]); 
  local xPlanes := []; 
  local yPlanes := []; 
  local zPlanes := []; 
  local xSig, ySig, zSig, Signature;
  
  while not T[finished] do 
    params := T[nextvalue](); 
    xPlanes := [op(xPlanes), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[1]]; 
    yPlanes := [op(yPlanes), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[2]]; 
    zPlanes := [op(zPlanes), (R . Vector(3, params[1]) - Vector(3, params[2]-1/2))[3]];
  end do; 
  xPlanes := eval(xPlanes, {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  yPlanes := eval(yPlanes, {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  zPlanes := eval(zPlanes, {a = samplePoints[1], b = samplePoints[2], c = samplePoints[3]});
  
  xSig, xPlanes := Isort(xPlanes); 
  ySig, yPlanes := Isort(yPlanes); 
  zSig, zPlanes := Isort(zPlanes);
  
  # remove all out of the range regions 
  xPlanes := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, xPlanes));
  yPlanes := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, yPlanes));
  zPlanes := remove(proc(x) return evalb(x <= -1/2) end proc,
  remove(proc(x) return evalb(x >= 1/2) end proc, zPlanes));
  # add region boarders
  xPlanes := [-1/2,op(xPlanes),1/2];
  yPlanes := [-1/2,op(yPlanes),1/2];
  zPlanes := [-1/2,op(zPlanes),1/2]; 

  Signature := cat(op(map(proc (x) sprintf("%d", x) end proc, [op(xSig), op(ySig), op(zSig)])));
  return Signature, xPlanes, yPlanes, zPlanes;
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
RecoverTranslationSamplePoints := proc(xPlanes::list, yPlanes::list, zPlanes::list) 
  local i, j, k; 
  local samples := []; 
  for i to upperbound(xPlanes) - 1 do
    for j to upperbound(yPlanes) - 1 do 
      for k to upperbound(zPlanes) - 1 do 
        samples := [op(samples), [(1/2)*add(xPlanes[i .. i+1]), 
                    (1/2)*add(yPlanes[j.. j+1]), (1/2)*add(zPlanes[k .. k+1])]];
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
CalculateNMM := proc(nType::string, kRange::list, dir::string, id::integer) 
  local csvFile, data, i, xPlanes, yPlanes, zPlanes, trans, NMM:={}, rotSamp; 
  local percent := -1, Signatures := {}, sig,fileID, msg;
  local n := Grid:-NumNodes();

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
    sig, xPlanes, yPlanes, zPlanes := GetOrderedCriticalPlanes(nType, kRange, rotSamp); 
    if member( sig, Signatures ) then
      next;
    end if;

    trans := RecoverTranslationSamplePoints(xPlanes,yPlanes, zPlanes); 
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
ParallelCalculateNMM := proc(nType::string, kRange::list, dir::string) 
  local me;
  me := Grid:-MyNode(); 
  RigidMotionsParameterSpaceDecompostion:-CalculateNMM(nType, kRange, dir, me);
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
LaunchOnGridGetNMM := proc (nType::string, kRange::list, dir::string) 
  Grid:-Setup("local"); 
  Grid:-Launch(ParallelCalculateNMM, nType, kRange, dir, numnodes =
  kernelopts(numcpus)):
end proc:

end module:
