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
    error "Wrong size of the input set. Expected size is at least 2.";
  fi;
  if not andmap(type, S, `polynom`) then
    error "Wrong type of elements. Expected argument is a set of polynomials!"; 
  fi;
  if nops(vars) <= 1 then
    error "Wrong number of indeterminates. It should be at least 2.";
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

# Procedure: RemoveExponant
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
    elif  nops(indets({p,q}))=1 then
        return gcd(p, q);
    else
        return p;
    end if;
end proc;


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

SortAlgebraicNumbers := proc (numbers::list)
  return sort( numbers, 
        proc( l, r ) 
        if Compare( l, r ) = -1 then
          return true:
        else 
          return false:
        fi;
        end proc);
end proc:

