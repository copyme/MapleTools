
# Procedure: ComputeEventsAType2D
#   Compute events such that a sweep line is tangent to a conic. 
#
# Parameters:
#   Q          - a set of conics
#   dim        - a list of indexes of variables used to calculate partial derivatives
#
# Output:
#   List of ranges which contains roots of a system(q, d/db q, d/dc q).
#
# Comment:
#  - only the first direction is supported since EliminationResultant is used
ComputeEventsAType2D := proc( Q )
  local s:
   s := proc(i::integer)
    local sys, univ, sol, vars:
    local q := Q[i];
    vars := [ op( indets( q ) ) ];
    if nops(vars) = 2 then
      sys := { q, diff( q, vars[2] ) };
    else
      error "Only system in two variables is supported!";
    fi:
    univ := EliminationResultant(sys, [ op( indets(sys ) ) ]):
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
  return [seq(s(i),i=1..nops(Q))]:
end proc:

# Procedure: ComputeEventsBType2D
#   Compute events such that two conics intersects
#
# Parameters:
#   dir        - a direction of a gradient product it should
#                be the same as director of sweep 
#   Q          - a set of conics
#
# Output:
#   Indexes of quadrics which intersect and a component of a vector product of 
#   their gradients in given direction have a common root.
ComputeEventsBType2D := proc( Q )
  local s:
  s := proc (i, j)
    local p, sol, univ, sys;
    sys := {Q[i], Q[j]}:
    univ := EliminationResultant(sys, [ op( indets(sys) ) ]):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, [ op( indets(univ ) ) ]):
      sol := nops(select(e -> rhs(e) >= 0, sol)):
      if sol > 0 then
        return [ univ, [i,j] ]:
      fi:
    fi;
    return NULL:
   end proc:
   return [seq(seq(s(i,j),j=i+1..nops(Q)),i=1..nops(Q))]:
end proc:

# Procedure: ComputeEventsAlgebraicNumbers2D
#   Compute and sort events as algebraic numbers 
#
# Parameters:
#   Q     - set of conics
# Output:
#   Sorted set of real algebraic numbers
ComputeEventsAlgebraicNumbers2D := proc( Q::~set )
  local events, rootsF, rf, poly:
  local numbers := Array([]):
  local factored, sqrFree:

  events:= {op(ComputeEventsAType2D( Q )), op(ComputeEventsBType2D( Q ))}:
  for poly in events do
    factored := factors( poly[1] )[2,..,1]: 
    for sqrFree in factored do
      rootsF := RootFinding:-Isolate(sqrFree, output='interval'):
      for rf in rootsF do
        ArrayTools:-Append(numbers, [ Object( RealAlgebraicNumber, sqrFree, op(rf)[2][1],
        op(rf)[2][2] ), poly[2]]):
      od:
    od:
  od:
 
  numbers := sort(numbers, 
                           proc( l, r ) 
                             if Compare( l[1], r[1] ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  ):
  return numbers:
end proc:


