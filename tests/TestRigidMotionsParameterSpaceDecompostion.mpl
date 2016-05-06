# File: TestRigidMotionsParameterSpaceDecompostion.mpl  
#
# Description:
#  This file contains tests for the module RigidMotionsParameterSpaceDecompostion.
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


TestRigidMotionsParameterSpaceDecompostion := module()
  option package;
  local init, Test, CheckInit;
  export InitModule, TestCayleyTransform, TestGetNeighborhood, TestGetQuadric,
  TestEliminationResultant, TestIsMonotonic, TestComputeSetOfQuadrics, TestAll:


  InitModule := proc(path::string)
    local mplFile;
    mplFile :=
    FileTools:-JoinPath([sprintf("%s/RigidMotionsParameterSpaceDecompostion.mpl",path)],base=homedir):
    read mplFile: 
    Test := CodeTools:-Test;
    init := true:
  end proc:

  CheckInit := proc()
    if init <> true then
      error "You have to first call InitModule!";
    fi:
  end proc;

   
  TestCayleyTransform :=proc()
    local CayleyTransform;
    CheckInit();
    CayleyTransform := RigidMotionsParameterSpaceDecompostion:-CayleyTransform;
    Test( indets( CayleyTransform( {a} ) ), {a}, label="CayleyTransform: Test Pass 1");
    Test( whattype( CayleyTransform( {a} ) ), Matrix, label="CayleyTransform: Test Pass 2");
    Test( {upperbound( CayleyTransform( {a} ) )}, {2, 2}, label="CayleyTransform: Test Pass 3");
    Test( indets( CayleyTransform( {a, b} ) ), {a, b}, label="CayleyTransform: Test Pass 4");
    Test( whattype( CayleyTransform( {a, b} ) ), Matrix, label="CayleyTransform: Test Pass 5");
    Test( {upperbound( CayleyTransform( {a, b} ) )}, {2, 2}, label="CayleyTransform: Test Pass 6");
    Test( indets( CayleyTransform( {a, b, c} ) ), {a, b, c}, label="CayleyTransform: Test Pass 7");
    Test( whattype( CayleyTransform( {a, b, c} ) ), Matrix, label="CayleyTransform: Test Pass 8");
    Test( {upperbound( CayleyTransform( {a, b, c} ) )}, {3, 3}, label="CayleyTransform: Test Pass 9");
    Test( CayleyTransform( {a, b, c, d} ), "Unsupported dimension", testerror,
    label="CayleyTransform: Test Pass 10"); 
  end proc;


  TestGetNeighborhood := proc()
    local GetNeighborhood;
    CheckInit();
    GetNeighborhood := RigidMotionsParameterSpaceDecompostion:-GetNeighborhood;
    Test( whattype( GetNeighborhood("N1") ), list, label="GetNeighborhood: Test Pass 1");
    Test( nops( ListTools:-MakeUnique( GetNeighborhood("N1") ) ), 7, label="GetNeighborhood: Test Pass 2");
    Test( whattype( GetNeighborhood("N2") ), list, label="GetNeighborhood: Test Pass 3");
    Test( nops( ListTools:-MakeUnique( GetNeighborhood("N2") ) ), 19, label="GetNeighborhood: Test Pass 4");
    Test( whattype( GetNeighborhood("N3") ), list, label="GetNeighborhood: Test Pass 5");
    Test( nops( ListTools:-MakeUnique( GetNeighborhood("N3") ) ), 27, label="GetNeighborhood: Test Pass 6");
    Test( nops( ListTools:-MakeUnique( GetNeighborhood("TEST") ) ), "Not supported type", testerror,
    label="GetNeighborhood: Test Pass 7");
  end proc;
  

  TestGetQuadric := proc()
    local GetQuadric, R, N1, N2, N3;
    CheckInit();
    GetQuadric := RigidMotionsParameterSpaceDecompostion:-GetQuadric;
    R := RigidMotionsParameterSpaceDecompostion:-CayleyTransform( {a, b, c} );
    N1 := RigidMotionsParameterSpaceDecompostion:-GetNeighborhood( "N1" );
    N2 := RigidMotionsParameterSpaceDecompostion:-GetNeighborhood( "N2" );
    N3 := RigidMotionsParameterSpaceDecompostion:-GetNeighborhood( "N3" );

    Test( type( GetQuadric( R, Vector( N1[1] ), 
          Vector([1/2, 1/2, 1/2]), 1 ), polynom ),
          true, label="GetQuadric: Test Pass 1" );
    Test( type( GetQuadric( R, Vector( N2[8] ),
          Vector([1/2, 1/2, 1/2]), 1 ), polynom ),
          true, label="GetQuadric: Test Pass 2" );
    Test( type( GetQuadric( R, Vector( N3[27] ),
          Vector([1/2, 1/2, 1/2]), 1 ), polynom ),
          true, label="GetQuadric: Test Pass 3" );
    Test( type( GetQuadric( R, Vector( N1[7] ),
          Vector([1/2, 1/2, 1/2]), 1 ), polynom ),
          true, label="GetQuadric: Test Pass 4" );
    Test( GetQuadric( R, Vector( N1[7] ),
          Vector([1/2, 1/2, 1/2]), 4 ),
          "Wrong dimension", testerror, 
          label="GetQuadric: Test Pass 5" );
    Test( GetQuadric( R, Vector( N1[7] ),
          Vector([1/2, 1/2, 1/2]), -1 ),
          "Wrong dimension", testerror, 
          label="GetQuadric: Test Pass 6" );
  end proc;


  TestEliminationResultant := proc()
    local p1, p2, p3, EliminationResultant;
    CheckInit();
    EliminationResultant := RigidMotionsParameterSpaceDecompostion:-EliminationResultant;
    p1 := -15*a^4*b-59*a*b*c^3+30*a^2*c^2-27*a*b^3+16*c^2-28*b;  
    p2 := -48*a^3*b^2+53*a^2*b*c^2-91*a^2*c-88*b*c^2+92*c^2+43*c;  
    p3 := 9*a^3*c^2-60*a*b*c^3-83*b*c^4+83*a^2*c^2+16*b^2*c+71*b^2;

    Test( EliminationResultant( {p1, p2, p3} ), 0, label="EliminationResultant: Test Pass 1" );
    Test( EliminationResultant( {p1, p2} ), "Wrong size of the input", testerror,
          label="EliminationResultant: Test Pass 2" );
    Test( EliminationResultant( {"test", p2, p3} ), "Wrong type of", testerror,
          label="EliminationResultant: Test Pass 3" );
    Test( EliminationResultant( {p1, "test", p3} ), "Wrong type of", testerror,
          label="EliminationResultant: Test Pass 4" );
    Test( EliminationResultant( {p1, p2, "test"} ), "Wrong type of", testerror,
          label="EliminationResultant: Test Pass 5" );
    Test( EliminationResultant( {x, x^2, 2*x^4} ), "Wrong number of", testerror,
          label="EliminationResultant: Test Pass 6" );
    Test( EliminationResultant( {x, x^2, 2*y^4} ), "Wrong number of", testerror,
          label="EliminationResultant: Test Pass 7" );
    Test( EliminationResultant( {x*z*d, x^2, 2*y^4} ), "Wrong number of", testerror,
          label="EliminationResultant: Test Pass 8" );
    Test( type(EliminationResultant( {2*b, -6*c, a^2+b^2-3*c^2-3} ), polynom ), true,
          label="EliminationResultant: Test Pass 9" );
  end proc;

 TestIsMonotonic := proc()
 local p1, p2, IsMonotonic;
    CheckInit();
    IsMonotonic := RigidMotionsParameterSpaceDecompostion:-IsMonotonic;
    p1 := a^2+b^2-3*c^2-3;
    p2 := x^2 + 1;
    Test( IsMonotonic( p1 ), false, label="IsMonotonic: Test Pass 1" );
    Test( IsMonotonic( p2 ), true, label="IsMonotonic: Test Pass 2" );
    Test( IsMonotonic( -p2 ), true, label="IsMonotonic: Test Pass 3" );
    Test( IsMonotonic( p2 - 1 ), true, label="IsMonotonic: Test Pass 4" );
    Test( IsMonotonic( p2 - y ), false, label="IsMonotonic: Test Pass 5" );
 end proc;


 TestComputeSetOfQuadrics := proc()
 end proc;

  TestAll := proc()
    CheckInit();
    TestCayleyTransform();
    TestGetNeighborhood();
    TestGetQuadric();
    TestEliminationResultant();
    TestIsMonotonic();
  end proc;

end module;
