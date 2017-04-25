# File: MapleCodePrime.mpl  
#
# Description:
#  This file contains code written by external developers -- users of mapleprimes.com
# Author:
#  Alec, Mihailovs
# 
#
# Date:
#  27/12/2016
#
# License:
#  Contact the authors.
#
# Copyright (c) Alec Mihailovs
# All rights reserved.


RigidMotionsMaplePrimesCode := module() 
  option package;
  
  export Isort;

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
