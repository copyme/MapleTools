
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
      error "Only system in two variables is supported! (%1)", q;
    fi:
    #if not RootFinding:-HasRealRoots(sys) then
      #return NULL;
    #fi;
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
    #if not RootFinding:-HasRealRoots(sys) then
      #return NULL;
    #fi;
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

ComputeEventsAType1D := proc( Q )
  local q, factored, sqrFree, rootsF, rf, numbers := Array([]);
  for q in Q do
    if RootFinding:-HasRealRoots(q) then
      factored := factors( q )[2,..,1]: 
      for sqrFree in factored do
        rootsF := RootFinding:-Isolate(sqrFree, output='interval'):
        for rf in rootsF do
          ArrayTools:-Append(numbers, Object(RealAlgebraicNumber, sqrFree, op(rf)[2][1], op(rf)[2][2])):
        od:
      od:
   fi;
  od:
  return numbers;
end proc:

# Procedure: ComputeSamplePoints2D
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
#   Writes a list of sample points into a file "sam_id.csv". Note that all sample points are
#   positive since other variation are same up to some similarities (reflections and rotations).
#
ComputeSamplePoints2D := proc (Q::~set, cluster::list)
  local i, j, x, midpoint, sys, samplePoints := [], fileID, vars, disjointEvent:=[]:
  local oneD, tmp;
  for i from 1 to upperbound(cluster) - 1 do 
    sys := {}: 
    for x in cluster[i] do 
      sys := sys union Q[x[2]]:
    end do:

    vars := indets(sys):
    disjointEvent := DisjointRanges(cluster[i][1][1],cluster[i+1][1][1]);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
   
   # intersection of a line with  conics
    sys := eval(sys, vars[1] = midpoint):
    oneD := ComputeEventsAType1D(sys);

    if oneD = NULL then
      next;
    fi:
    oneD := convert(oneD, list);
    oneD := sort(oneD, 
                           proc( l, r ) 
                             if Compare( l, r ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  ):
    tmp, oneD := selectremove(proc(x) return evalb(GetInterval(x)[2] < 0); end proc, oneD):
    if nops(tmp) <> 0 then
      oneD := [tmp[-1], op(oneD)];
    fi:
    if upperbound(oneD) > 0 then
      for j from 1 to upperbound(oneD) - 1 do
        disjointEvent := DisjointRanges(oneD[j],oneD[j+1]);
        samplePoints := [op(samplePoints), [midpoint, (GetInterval(disjointEvent[1])[2] +
                                                       GetInterval(disjointEvent[2])[1])/2]];
      od:
      samplePoints := [op(samplePoints), [midpoint, GetInterval(oneD[-1])[2] + 1/2 ]];
    fi:


    #samplePoints := remove(proc(x) if rhs(x[1]) < 0 or rhs(x[2]) < 0 then return true: else return
                                                                     #false: fi: end, samplePoints):

    #fileID := fopen(sprintf("sam_%d.csv", id), APPEND, TEXT):
    #writedata(fileID, Threads:-Map(proc (x) return [midpoint, 
               #rhs(x[1]), rhs(x[2])] end proc, samplePoints),
              #string, proc (f, x) fprintf(f, %a, x) end proc):
    #fclose(fileID):
  od:

  return samplePoints:
end proc:


