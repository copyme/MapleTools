# Class: TestRealAlgebraicNumbersComparison.mpl
#
# Description:
#  Tests for comparisons of objects implemented by the class RealAlgebraicNumber.mpl
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  14/04/2016
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


# Procedure: TestRealAlgebraicNumbersComparisonLower
#   Some test of implementation of real algebraic numbers.
#
# Parameters:
#   loops::integer             - a number of random polynomials to consider,
#                                 default value is 10
#   deg::integer               - a maximal possible degree of random polynomial
#
TestRealAlgebraicNumbersComparison := proc( loops::integer := 10, deg::integer := 80 )
  local numA := Object( RealAlgebraicNumber, 3 * x - 1, 91625968981/274877906944, 45812984491/137438953472 );
  local numB := Object( RealAlgebraicNumber, x^2 - 3, -238051250353/137438953472, -14878203147/8589934592 );
  local numC := Object( RealAlgebraicNumber, x - 1/2, -1, 1 );
  local numD := Object( RealAlgebraicNumber, 3*x - 1/2, 0, 1 );
  local numE := Object( RealAlgebraicNumber, x, 0, 1 );
  local numF := Object( RealAlgebraicNumber, x, -1, 0 );
  local f::polynom, g::polynom, rootsF := [], rootsG := [], numbers := [];
  local i, rf, rg;
  print( "numA : ", numA );
  print( "numB : ", numB );
  kernelopts( assertlevel = 1 );
    (* Test for A rational and B non-rational, or equal. *)
    ASSERT( evalb( numB  < numA ) = true, " ( numA < numB ) gives false should true" );
    ASSERT( evalb( numA  < numB ) = false, " ( numA < numB ) gives true should false" );
    ASSERT( evalb( numA  < numA ) = false, " ( numA < numA ) gives true should false" );
    ASSERT( evalb( numB  > numA ) = false, " ( numA > numB ) gives false should true" );
    ASSERT( evalb( numA  > numB ) = true, " ( numA > numB ) gives false should true" );
    ASSERT( evalb( numA  > numA ) = false, " ( numA > numA ) gives true should false" );
    ASSERT( evalb( numB  <= numA ) = true, " ( numA <= numB ) gives false should true" );
    ASSERT( evalb( numA  <= numB ) = false, " ( numA <= numB ) gives true should false" );
    ASSERT( evalb( numA  <= numA ) = true, " ( numA <= numA ) gives false should true" );
    ASSERT( evalb( numB  >= numA ) = false, " ( numA >= numB ) gives true should false" );
    ASSERT( evalb( numA  >= numB ) = true, " ( numA >= numB ) gives false should true" );
    ASSERT( evalb( numA  >= numA ) = true, " ( numA >= numA ) gives false should true" );
    ASSERT( evalb( numA  = numA ) = true, " ( numA = numA ) gives false should true" );
    ASSERT( evalb( numB  = numB ) = true, " ( numB = numB ) gives false should true" );
    ASSERT( evalb( numA  = numB ) = false, " ( numrA = numB ) gives true should false" );
    ASSERT( evalb( numB  = numA ) = false, " ( numrB = numA ) gives true should false" );
    (*Test for numbers with ill ranges*)
    ASSERT( evalb( numC  = numD ) = false, " ( numrC = numD ) gives true should false" );
    ASSERT( evalb( numE  = numF ) = true, " ( numrE = numF ) gives false should true" );
    ASSERT( evalb( numC  > numE ) = true, " ( numrC > numE ) gives false should true" );
    ASSERT( evalb( numC  > numD ) = true, " ( numrC > numD ) gives false should true" );
    ASSERT( evalb( numC  < numD ) = false, " ( numrC < numD ) gives true should false" );
    ASSERT( evalb( numC  = numF ) = false, " ( numrC = numF ) gives true should false" );
  kernelopts( assertlevel = 0 );
    for i from 1 to loops do
       f := randpoly( x, degree=MapleTA[Builtin][rint]( 1, deg ) );
       f := sqrfree( factor( f ) )[2,..,1][1];
       g := randpoly( x, degree=MapleTA[Builtin][rint]( 1, deg )  );
       g := sqrfree( factor( g ) )[2,..,1][1];
       rootsF := RootFinding:-Isolate ( f, x, output='interval' );
       rootsG := RootFinding:-Isolate ( g, x, output='interval' );
       for rf in rootsF do
        numbers := [ op( numbers ), Object( RealAlgebraicNumber, f, op( rf )[2][1], op( rf )[2][2] ) ];
       end do;
       for rg in rootsG do
        numbers := [ op( numbers ), Object( RealAlgebraicNumber, g, op( rg )[2][1], op( rg )[2][2] ) ];
       end do;
    end do;
    for numA in numbers do
      for numB in numbers do
        i := Compare( numA, numB );
        if i = -1 then 
          print( numA, " < ", numB  );
        elif i = 0 then
          print( numA, " = ", numB  );
        elif i = 1 then
          print( numA, " > ", numB  );
        else
          error "Two numbers was not compared correctly.";
        end if;
       end do;
    end do;
    print( "Use sort() to sort the list of algebraic numbers." );
    numbers := sort( numbers, proc(l,r) 
                     if Compare( l, r ) = -1 then
                       return true;
                     else 
                       return false;
                     end if;
                  end proc );
    print( numbers );
end proc:

