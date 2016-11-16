
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
ComputeEventsAType2D := proc( Q2D2, grid::boolean )
  local s:
   s := proc(i::integer)
    local sys, univ, sol, vars:
    local q := Q2D2[i];
    vars := [ op( indets( q ) ) ];
    if nops(vars) = 2 then
      sys := { q, diff( q, vars[2] ) };
    else
      error "Only system in two variables is supported! (%1)", q;
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
  if grid then
    return [Grid:-Seq(s(i),i=1..nops(Q2D2))]:
  else 
    return [seq(s(i),i=1..nops(Q2D2))]:
  fi;
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
ComputeEventsBType2D := proc( Q2D2, grid::boolean )
  local s:
  s := proc (i, j)
    local p, sol, univ, sys;
    sys := {Q2D2[i], Q2D2[j]}:
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
   if grid then
     return [Grid:-Seq(seq(s(i,j),j=i+1..nops(Q2D2)),i=1..nops(Q2D2))]:
   else
     return [seq(seq(s(i,j),j=i+1..nops(Q2D2)),i=1..nops(Q2D2))]:
   fi;
end proc:

# Procedure: IsAsymptotic2D
#   Checks if intersection of two quadrics has an asymptotic critical value. For this moment a
#   direction is fixed to a.
#
# Parameters:
#   p          - a quadric in three variables
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
IsAsymptotic2D := proc( p::polynom, var  )
  return lcoeff(p, var);
end proc:


# Procedure: ComputeAsymptoticAAEvents2D
#   Compute real algebraic numbers which corresponds to 
#   asymptotic cases given by one quadrics.
#
# Parameters:
#   Q          - a set of quadrics
#
# Output:
#   A list of real algebraic numbers and indexes of quadrics
#   which corresponds to them.
ComputeAsymptoticAAEvents2DGrid:=proc(Q2D2)
  local list := [], s;
  s:=proc(i::integer)
    local rf, rootsF;
    local numbers := [], sol:
     sol := IsAsymptotic2D(Q2D2[i], indets(Q2D2[i])[-1]):
     if sol <> NULL and nops(sol) <> 0 then
       rootsF := RootFinding:-Isolate(sol, op(indets(sol)), output='interval');
       for rf in rootsF do
         numbers:=[op(numbers), [Object(RealAlgebraicNumber, sol, op(rf)[2][1],
         op(rf)[2][2]), [i]]]: 
       od:
      numbers:=[op(numbers), [Object(RealAlgebraicNumber, sol * denom(sol) -
      numer(sol), sol, sol), [i]]]:
     fi:
    return numbers;
  end proc:
  return [Grid:-Seq(s(i),i=1..nops(Q))];
end:


# Procedure: ComputeEventsAlgebraicNumbers2D
#   Compute and sort events as algebraic numbers 
#
# Parameters:
#   Q     - set of conics
# Output:
#   Sorted set of real algebraic numbers
ComputeEventsAlgebraicNumbers2D := proc( Q2D2, grid::boolean )
  local events, rootsF, rf, poly:
  local numbers := Array([]):
  local factored, sqrFree:

  events:= {op(ComputeEventsAType2D( Q2D2, grid )), op(ComputeEventsBType2D( Q2D2, grid )),
  op(ComputeAsymptoticAAEvents2DGrid(Q2D2))}:
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

ComputeEventsAType1D := proc( Q2D2 )
  local q, factored, sqrFree, rootsF, rf, numbers := Array([]);
  for q in Q2D2 do
  print(q);
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


# Procedure: ParallelComputeSamplePoints2D
#   Computes sample points for rotational part of rigid motions. It should be call via Grid
#   framework.
#
ParallelComputeSamplePoints2D := proc () 
  local me, numNodes, n;
  me := Grid:-MyNode();
  numNodes := Grid:-NumNodes();
  # cluster-1 because the last cluster is a doubled cluster[-2]
  n := trunc((upperbound(cluster2D)-1)/numNodes);
  ComputeSamplePoints2D(Q2D, cluster2D, me*n+1,(me+1)*n, me, aValue);
  Grid:-Barrier();
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
ComputeSamplePoints2D := proc(Q2D, cluster2D::list, first::integer,
                             last::integer, id::integer, aValue)
  local i, j, x, midpoint, sys, samplePoints := [], fileID, vars, disjointEvent:=[]:
  local oneD, tmp;
  if first < 0 or last < 0 or last < first or upperbound(cluster2D) <= last then 
    error "Bounds of the cluster range are incorrect.": 
  end if:
  for i from first to last do 

    sys := []: 
    for x in cluster2D[i] do 
      sys := [op(sys), Q2D[x[2]]]:
    end do:
    sys := ListTools:-Flatten(sys);
    sys := ListTools:-MakeUnique(sys);
    vars := indets(sys):
    disjointEvent := DisjointRanges(cluster2D[i][1][1],cluster2D[i+1][1][1]);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
   
   # intersection of a line with  conics
    # never call eval with sets!
    sys := eval(sys, vars[1] = midpoint):
    oneD := ComputeEventsAType1D(sys);

    if oneD = NULL then
      next;
    fi:
    oneD := convert(oneD, list);
    oneD := SortAlgebraicNumbers(oneD);
    tmp, oneD := selectremove(proc(x) return evalb(GetInterval(x)[2] < 0); end proc, oneD):
    if nops(tmp) <> 0 then
      oneD := [tmp[-1], op(oneD)];
    fi:
    if upperbound(oneD) > 0 then
      for j from 1 to upperbound(oneD) - 1 do
        disjointEvent := DisjointRanges(oneD[j],oneD[j+1]);
        samplePoints := [op(samplePoints), [aValue, midpoint, (GetInterval(disjointEvent[1])[2] +
                                                       GetInterval(disjointEvent[2])[1])/2]];
      od:
      samplePoints := [op(samplePoints), [aValue, midpoint, GetInterval(oneD[-1])[2] + 1/2]];
    fi:
 od:
 if nops(samplePoints) <> 0 then
   fileID := fopen(sprintf("sam_%d.csv", id), APPEND, TEXT):
   writedata(fileID, samplePoints, string, proc (f, x) fprintf(f, %a, x) end proc):
   fclose(fileID):
 fi:
 return NULL;
end proc:


# Procedure: ComputeSamplePoints2D
#   Computes sample points for rotational part of rigid motions using the grid framework
#
#
# Parameters:
#   nType      - size of neighborhood i.e. N_1, N_2, N_3. 
#   kRange     - a range of planes to consider
# Output:
#   Writes a list of sample points into a file "sam_id.csv" where id corresponds to an id of used
#   thread during computations.
LaunchOnGridComputeSamplePoints2D := proc (s::list, midpoint, nodes::integer, grid::boolean, id::integer) 

  local numbers, events, R, rootTmp, n := nodes: 
  global Q2D := s, cluster2D, aValue := midpoint:
  if nodes > 1 then
     numbers := convert(ComputeEventsAlgebraicNumbers2D(Q2D, true), list);
  else
     numbers := convert(ComputeEventsAlgebraicNumbers2D(Q2D, false), list);
  fi;
  numbers := ThreadsRemove(proc(x) return evalb(GetInterval(x[1])[2] < 0); end proc, numbers):
  if upperbound(numbers) = 0 then
    return NULL;
  fi;
  cluster2D := ClusterEvents(numbers):
  if upperbound(cluster2D) < nodes then
    n := upperbound(cluster2D);
  fi;
  cluster2D := [[[cluster2D[1][1][1], convert(Q2D, list)]], op(cluster2D[2..])]:
  # add the last slice twice but shifted to calculate correctly last quadrics
  events := cluster2D[-1][1][1]:
  rootTmp:= GetInterval(events)[2]+1/4;
  events := Object(RealAlgebraicNumber, denom(rootTmp) * indets(GetPolynomial(events))[1] -
  numer(rootTmp), rootTmp, rootTmp):
  cluster2D := [op(cluster2D), [[events, cluster2D[-1][1][2]]]]:
  # The first cluster2D is heavy so we compute it separately;
  # We define printer as a procedure which returns NULL to avoid a memory leak problem while
  # writing to a file from a node. It seems that while fprintf is called it also calls printf
  # which is a default printer function. Therefore, data are returned to node of ID 0. 
  if grid then
    Grid:-Launch(ParallelComputeSamplePoints2D,
                 imports = ['Q2D, cluster2D, aValue'], numnodes=n, printer=proc(x) return NULL: end proc
                 ):
  else
    ComputeSamplePoints2D(Q2D, cluster2D, 1, nops(cluster2D) - 1, id, aValue);
  fi;
end proc:

