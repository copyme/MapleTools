
# Procedure: ComputeEventsAType2D
#   Compute events such that a sweep line is tangent to a conic with respect to the first provided
#   variable. 
#
# Parameters:
#   Q2D2       - a set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars       - a list of the variables
#
# Output:
#   list of couples: univariate polynomial and index of the generating conic: [poly, [index]]
#
# Comment/Limitations:
#  - Univariate polynomials are expressed in the first variable.
ComputeEventsAType2D := proc(Q2D2, grid::boolean, vars::list)
  local s:
  if nops(vars) < 2 then
    error "ComputeEventsAType2D: Only systems in at least two variables are supported.";
  fi:
  s := proc(i::integer, vars::list)
    local sys, univ, sol:
    local q := Q2D2[i];
    sys := { q, diff( q, vars[2] ) };
    univ := EliminationResultant(sys, vars):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, vars[1..1]):
      if nops(sol) > 0 then
        return [univ, [i]]:
      else
        return NULL:
      fi:
    end if:
  end proc:
  if grid then
    return [Grid:-Seq( s( i, vars ), i=1..nops( Q2D2 ) )]:
  else 
    return [seq( s( i, vars ), i=1..nops( Q2D2 ) )]:
  fi;
end proc:

# Procedure: ComputeEventsBType2D
#   Compute events such that two conics intersects in a point.
#
# Parameters:
#   Q2D2       - a set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars       - a list of the variables
#
# Output:
#   list of couples: univariate polynomial and indices of the generating conics: [poly,
#   [index1,index2]]
#
# Comment/Limitations:
#  - Univariate polynomials are expressed in the first variable.
ComputeEventsBType2D := proc(Q2D2, grid::boolean, vars::list)
  local s:
  s := proc (i, j, vars::list)
    local p, sol, univ, sys;
    sys := {Q2D2[i], Q2D2[j]}:
    univ := EliminationResultant(sys, vars):
    if not type( univ, constant ) then
      sol := RootFinding:-Isolate( univ, vars[1..1]):
      if nops(sol) > 0 then
        return [ univ, [i,j] ]:
      fi:
    fi;
    return NULL:
   end proc:
   if grid then
     return [Grid:-Seq( seq( s( i, j, vars ), j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]:
   else
     return [seq( seq( s(i, j, vars ), j=i+1..nops( Q2D2 ) ), i=1..nops( Q2D2 ) )]:
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
#   vars       - a list of the variables
#
# Output: Real
#   A list of real algebraic numbers and indexes of conics: [ RealAlgebraicNumber, [index]].
ComputeAsymptoticAAEvents2DGrid := proc(Q2D2, vars::list)
  local out := [], s;
  s:=proc(i::integer, vars::list)
    local rf, rootsF;
    local numbers := [], asy:
     asy := IsAsymptotic2D(Q2D2[i], vars[-1]);
     if not type(asy, constant) then
       rootsF := RootFinding:-Isolate(asy, vars[1], output='interval');
       for rf in rootsF do
         numbers:=[op(numbers), [Object(RealAlgebraicNumber, asy, op(rf)[2][1],
         op(rf)[2][2]), [i]]];
       od:
     fi:
    return numbers;
  end proc:
  out:=select(proc(x) return evalb(x<>[]) end, [seq(s(i, vars),i=1..nops(Q2D2))]);
  return ListTools:-Flatten(out, 1);
end:


# Procedure: ComputeEventsAlgebraicNumbers2D
#   Compute and sort events as algebraic numbers 
#
# Parameters:
#   Q2D2     - set of conics
#   grid       - if true then computations are performed in the parallel computation framework
#                called Grid (see Maple documentation)
#   vars       - a list of the variables
#
# Output:
#   Sorted Array of real algebraic numbers
ComputeEventsAlgebraicNumbers2D := proc(Q2D2, grid::boolean, vars::list)
  local events, rootsF, rf, poly;
  local numbers := Array([]);
  local factored, sqrFree;

  events:= [op(ComputeEventsAType2D(Q2D2, grid, vars)), op(ComputeEventsBType2D(Q2D2, grid, vars))];
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
  ArrayTools:-Concatenate(2, numbers, Vector[row]([ComputeAsymptoticAAEvents2DGrid(Q2D2, vars)]));
  
  # In maple 2015.2 there is a bug which causes: stack limit reached
  if upperbound(numbers) <> 0 then
    numbers := sort(numbers, 
                           proc( l, r ) 
                             if Compare( l[1], r[1] ) = -1 then
                               return true:
                             else 
                               return false:
                             fi:
                           end proc
                  ):
  fi:
  return numbers;
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
  ComputeSamplePoints2D(Q2D, cluster2D, me*n+1,(me+1)*n, me, vars, writer);
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
ComputeSamplePoints2D := proc(Q2D, cluster2D::list, first::integer, last::integer, id::integer,
vars::list, writer)
  local i::integer, j::integer, x::list, midpoint::rational, sys::list, samplePoints::list;
  local disjointEvent::list, oneD::list, oneDNeg::list;
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
    disjointEvent := DisjointRanges(cluster2D[i][1][1],cluster2D[i+1][1][1]);
    midpoint := (GetInterval(disjointEvent[1])[2] + GetInterval(disjointEvent[2])[1])/2:
   
    # intersection of a line with  conics
    # never call eval with sets!
    sys := eval(sys, vars[1] = midpoint):
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
        samplePoints := [op(samplePoints), [[midpoint, (GetInterval(disjointEvent[1])[2] +
                                                       GetInterval(disjointEvent[2])[1])/2]]];
      od:
      samplePoints := [op(samplePoints), [[midpoint, GetInterval(oneD[-1])[2] + 1/2]]];
      Write(writer, samplePoints, id);
    fi:
 od:
 return NULL;
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
#               the size of the problem are fullfiled the problem is solved in the grid framework.
#   id        - id of a file used when grid computations are set to false
#   variables - list of variables in which the problem is expressed
#   path      - directory in which the output is going to be saved
#   prefix    - file name prefix
# Output:
#   The output file(s) are saved into a files: path/prefix(id).tsv
LaunchOnGridComputeSamplePoints2D := proc (s::list, midpoint::rational, nodes::integer,
grid::boolean, id::integer, variables::list, path::string, prefix::string) 
  local numbers, firstEvent, R, rootTmp, n := nodes;
  global Q2D := ListTools:-MakeUnique([op(variables),op(s)]), cluster2D, writer, vars := variables;
  if grid and nops(s) > 20 then
     numbers := convert(ComputeEventsAlgebraicNumbers2D(Q2D, true, vars), list);
  else
     numbers := convert(ComputeEventsAlgebraicNumbers2D(Q2D, false, vars), list);
  fi;
  numbers := remove(proc(x) return evalb(GetInterval(x[1])[2] < 0); end proc, numbers):
  if upperbound(numbers) = 0 then
    return NULL;
  fi;
  cluster2D := ClusterEvents(numbers):
  if upperbound(cluster2D) < nodes then
    n := upperbound(cluster2D);
  fi;

  # assign all conics to the first event
  cluster2D := [[[cluster2D[1][1][1], [seq(1..nops(Q2D))]]], op(cluster2D[2..])]:
  rootTmp:= GetInterval(cluster2D[-1][1][1])[2]+1;
  firstEvent := Object(RealAlgebraicNumber, denom(rootTmp) * indets(GetPolynomial(firstEvent))[1] -
  numer(rootTmp), rootTmp, rootTmp):
  cluster2D := [op(cluster2D), [[firstEvent, cluster2D[-1][1][2]]]];

  writer := Object(SamplePointsWriter, midpoint, path, prefix);

  # The first cluster2D is heavy so we compute it separately;
  # We define printer as a procedure which returns NULL to avoid a memory leak problem while
  # writing to a file from a node. It seems that while fprintf is called it also calls printf
  # which is a default printer function. Therefore, data are returned to node of ID 0. 
  if grid and nodes > 1 then
    Grid:-Launch(ParallelComputeSamplePoints2D, imports = ['Q2D, cluster2D, vars, writer'], numnodes=n,
                 printer=proc(x) return NULL: end proc);
  else
    ComputeSamplePoints2D(Q2D, cluster2D, 1, nops(cluster2D) - 1, id, vars, writer);
  fi;
end proc:

