# File: TestRigidMotionsParameterSpaceDecompositionRecursive.mpl  
#
# Description:
#  This file contains tests for the module RigidMotionsParameterSpaceDecompositionRecursive.
#
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  05/05/2016
#
# License:
#  Simplified BSD License
#
# Copyright (c) 2016, Kacper Pluta
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


TestRigidMotionsParameterSpaceDecompositionRecursive := module()
  option package;
  local init, Test, CheckInit;
  export InitModule, TestSpecialCases, TestAll:


  InitModule := proc()
    local real::RealAlgebraicNumber, writer::SamplePointsWriter;
    Test := CodeTools:-Test;
    init := true:
  end proc:

  CheckInit := proc()
    if init <> true then
      error "You have to first call InitModule!";
    fi:
  end proc;

  TestSpecialCases :=proc()
    local quadrics := Array([]);
    local i::integer;
    CheckInit();
    quadrics(1) := [x*y-1];
    quadrics(2) := [x+y-1/4];
    quadrics(3) := [x+y-1/9];
    quadrics(4) := [x^2+y^2-1];
    quadrics(5) := [y-(x-1/9)^2+1/100];
    quadrics(6) := [x*y-3*x-3*y+8, x-3-y/10^50];
    quadrics(7) := [x*y-3*x-3*y+8, y-3-x/10^50];
    quadrics(8) := [x*y-1, x, x-1]; # error expected
    quadrics(9) := [x-1, y-1];
    quadrics(10) := [x^2-2*x-2, y^2-2*y-2];
    quadrics(11) := [x*y-3*x-3*y+8, x-3, y-5/2];
    quadrics(12) := [x[1]^2+x[2]^2-1];

    # Fix the seed for reproducible random generation
    RandomTools:-MersenneTwister:-SetState(state=123456789);
    quadrics(13) := [randpoly $ 60]([x,y], dense, degree=2, coeffs=rand(0..1));
    quadrics(13) := eval(quadrics(13), [x = -x, y = -y]);

    for i from 1 to upperbound(quadrics) do
    try
       LaunchOnGridComputeSamplePoints2D(quadrics[i], 0, 20, true, i, [op(indets(quadrics(i)))],
       "/tmp", "sam_" );
    catch:
        print("Error", lastexception);
    end try;
    od;
  end proc;

  TestAll := proc()
    CheckInit();
    TestSpecialCases();
  end proc;

end module;
