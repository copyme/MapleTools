# File: MapleCodePrime.mpl  
#
# Description:
#  This file contains code written by external developers -- users of mapleprimes.com
# Author:
#  Carl Love, Alec, Mihailovs
# 
#
# Date:
#  27/12/2016
#
# License:
#  Contact the authors.
#
# Copyright (c) Carl Love, Alec Mihailovs
# All rights reserved.


RigidMotionsMaplePrimesCode := module() 
  option package;
  export SplitScan, Isort;
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

end module;
