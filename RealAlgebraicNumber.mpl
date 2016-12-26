# Class: RealAlgebraicNumber
#
# Description:
#  Implementation of real algebraic numbers together with their comparison.
#  This implementation is inspired by the implementation in CGAL 4.7.
#
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  11/12/2015 
#
# License:
#  Simplified BSD License
#
# Copyright (c) 2015, Kacper Pluta
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
module RealAlgebraicNumber()
  option object;
  (* Univariete polynomaial *)
  local poly::polynom;
  (* Lower bound of the range in which exists only one real root of poly. *)
  local a::rational;
  (* Upper bound of the range in which exists only one real root of poly. *)
  local b::rational;
  (* Real algebraic number is rational when a = b and sign of poly at a/b is 0. *)
  local isRational_;
  (* Note that isolating interval has to be open iff a real algebraic number is not rational and
  closed, a = b, otherwise.*)

# Method: ModuleCopy
#   Standard constructor / copy constructor
#
# Parameters:
#   self::RealAlgebraicNumber     - a new object to be constructed
#   proto::RealAlgebraicNumber    - a prototype object from which self is derived
#   poly::polynom                 - a univariate polynomial
#   a::rational                   - a lower bound of the range in which exists only one real root of poly 
#   b::rational                   - an upper bound of the range in which exists only one real root of poly 
#
# Output:
#   An object of type RealAlgebraicNumber.
#
# Exceptions:
#  "Invalid range. A range is valid when: a <= b."
#  "Degree of %1 is invalid."
#
  export ModuleCopy::static := proc( self::RealAlgebraicNumber,
                                     proto::RealAlgebraicNumber,
                                     poly::polynom,
                                     a::rational,
                                     b::rational, $ )
    local signAtA, signAtB;
    if _passed = 2 then
      self:-poly := proto:-poly;
      self:-a := proto:-a;
      self:-b := proto:-b;
      self:-isRational_ := proto:-isRational_;
    else
      if upperbound( indets( poly ) ) > 1 then
        error "%1 is not univariate!", poly;
      end if;
      if a > b then
        error "Invalid range. A range is valid when: a <= b.";
      end if;
      if gcd( poly, diff( poly, op( indets( poly ) ) ) ) <> 1 then
        error "Polynomial: %1 is not square-free.", poly;
      end if;
      self:-poly := poly;
      if degree(poly) >= 1 then
        signAtA := signum( eval( poly, indets( poly )[1] = a ) );
        signAtB := signum( eval( poly, indets( poly )[1] = b ) );
        self:-a := a;
        self:-b := b;
        self:-isRational_ := evalb( self:-a = self:-b and signAtA = 0 );
        if signAtA = 0 and signAtB <> 0 then
          WARNING("Incorrect interval. Sign of univariate polynomial on one side of the interval"
          " evaluated to zero but not on the another. Interval fixed.");
          self:-b := self:-a;
          self:-isRational_ := true:
        elif signAtB = 0 and signAtA <> 0 then 
          WARNING("Incorrect interval. Sign of univariate polynomial on one side of the interval"
          " evaluated to zero but not on the another. Interval fixed.");
          self:-a := self:-b;
          self:-isRational_ := true:
        elif signAtA = signAtB and self:-a <> self:-b then
          error "Interval incorrect! No root in the interval: (%1, %2), for %3 .", self:-a, self:-b,
          self:-poly;
        fi:
      elif degree( poly ) = 0 then
        self:- denom( poly ) * 'a'  - numer( poly );
        self:-a := poly;
        self:-b := poly;
        self:-isRational_ := true;
      else
        error "Degree of %1 is invalid.", poly;
      end if;
    end if;
    return self;
  end proc:

# Method: ModulePrint
#   Standard printout of an object of type RealAlgebraicNumber.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#
  export ModulePrint::static := proc( self::RealAlgebraicNumber )
    if(self:-a = self:-b) then
        nprintf( "( %a, [%a, %a] )", self:-poly, self:-a, self:-b );
    else
        nprintf( "( %a, ]%a, %a[ )", self:-poly, self:-a, self:-b );
    end if;
  end proc:

# Method: GetPolynomial
#   A getter method to access the univariate polynomial of RealAlgebraicNumber.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   Univariate polynomial stored in self.
#
  export GetPolynomial::static := proc( self::RealAlgebraicNumber )
    return self:-poly;
  end proc:

# Method: GetInterval
#   A getter method to access the range isolating a root of univariate polynomial.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   The range isolating a root of univariate polynomial -- self:-poly.
#   Upper and lower bounds of a range are rationals. When lower = upper
#   then a real algebraic number is rational.
#
  export GetInterval::static := proc( self::RealAlgebraicNumber )
    return [ self:-a, self:-b ];
  end proc:

# Method: IsRational
#   A method to check if a real algebraic number is rational.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   True when a real algebraic number is rational, false
#   otherwise. A real algebraic number is meant as rational when a = b and
#   when sign of poly at a/b is zero.
#
  export IsRational::static := proc( self::RealAlgebraicNumber )
    return self:-isRational_;
  end proc:

# Method: CompareRational
#   A method used to compare a real algebraic number with a
#   rational number.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#   m::rational                    - a rational number
#
# Output:
#   -1 when a real algebraic number is smaller than a rational
#    number, 0 when they are equal and 1 when a real algebraic
#    number is bigger than a rational.
#
  local CompareRational::static := proc( self::RealAlgebraicNumber, m::rational )
    local refined:
    refined := StrongRefineAt(self,m);
    if evalb( refined:-a < m ) then
      return -1;
    elif evalb( refined:-a > m ) then
      return 1;
    elif evalb( signum( eval( refined:-poly, op( indets( refined:-poly ) ) = m ) ) = 0 ) then
      return 0;
    end if;
  end proc:

# Method: RefineAt
#   A method used to refine a real algebraic number using a rational
#   number for adaptation of a range isolating a root of poly.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#   m::rational                    - a rational number
#
# Output:
#   A RealAlgebraicNumber obtaind from self refined at m.
#
# Comment:
#   Not that the type can change to rational.
#
  local RefineAt::static := proc( self::RealAlgebraicNumber, m::rational )
    local signAtM, f::polynom, g::polynom;
    local var := op( indets( self:-poly ) );
    if self:-isRational_ or m <= self:-a or self:-b <= m then
      return self;
    end if;
    signAtM := signum( eval( self:-poly, var = m ) );
    if evalb( signAtM = 0 ) then
     g := denom( m ) * var  - numer( m );
     return Object( RealAlgebraicNumber, g, m, m ); 
    elif evalb( signum( eval( self:-poly, var = self:-a ) ) = signAtM ) then
      return Object( RealAlgebraicNumber, self:-poly, m, self:-b );
    elif evalb( signum( eval( self:-poly, var = self:-b ) ) = signAtM ) then
      return Object( RealAlgebraicNumber, self:-poly, self:-a, m );
    else
      return self;
    end if;
  end proc:


# Method: BisectRange
#   A method used to compare a real algebraic number with a
#   rational number.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   Refine an isolating range at ( self:-a + self:-b ) / 2
#
  local BisectRange::static := proc( self::RealAlgebraicNumber )
    return RefineAt( self, ( self:-a + self:-b ) / 2 );
  end proc:


# Method: StrongRefineAt
#   A method used to refine a real algebraic number using a rational
#   number for adaptation of a range isolating a root of poly.
#
# Parameters:
#   self::RealAlgebraicNumber      - a real algebraic number
#   m::rational                    - a rational number
#
  local StrongRefineAt::static := proc( self::RealAlgebraicNumber, m::rational )
    local refined:
    if self:-isRational_ or signum( eval( self:-poly, indets( self:-poly )[1] = m ) ) = 0 then
      return self;
    fi:
    refined := self;
    while refined:-a <= m and m <= refined:-b do
      refined := BisectRange(refined);
    od:
    return refined:
  end proc:

  export DisjointRanges::static := proc(a::RealAlgebraicNumber, b::RealAlgebraicNumber)
    local ll := a, rr := b, i;

    (* No intersection.*)
    if evalb( ll:-b < rr:-a ) or evalb( ll:-a > rr:-b ) then
      return [a,b];
    fi:

    for i from 1 while 1 = 1 do
      ll := BisectRange( ll ):
      rr := BisectRange( rr ):

      (* No intersection.*)
      if evalb( ll:-b < rr:-a ) or evalb( ll:-a > rr:-b ) then
        return [ll,rr];
      fi:
    od:
  end proc:

# Method: Compare
#   A method used to compare two real algebraic numbers.
#
# Parameters:
#   l::RealAlgebraicNumber      - a real algebraic number
#   r::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   -1 when l is smaller than r, 0 when they are equal and 1 when l is bigger than r.
#
  export Compare::static := proc( l::RealAlgebraicNumber, r::RealAlgebraicNumber, $ )          
    local i::integer, a::rational, b::rational, F1::polynom, F2::polynom, G::polynom;
    local ll::RealAlgebraicNumber, rr::RealAlgebraicNumber;

    if indets(l:-poly) <> indets(r:-poly) then
      error "Univariate polynomials have different variables: %1 and %2.", indets(l:-poly),
      indets(r:-poly);
    fi;

    if evalb( l:-poly = r:-poly and l:-a = r:-a and l:-b = r:-b ) then
        return 0;
    end if; 
 
    (* When rationals *)
    if r:-isRational_ then
      return CompareRational( l, r:-a );
    elif l:-isRational_ then
      return -CompareRational( r, l:-a );
    end if;

    (* Check if there is no intersection of the ranges *)
    if evalb( l:-b < r:-a ) then
      return -1;
    elif evalb( l:-a > r:-b ) then
      return 1;
    end if:

    (* Get the intersecting interval *)
    if evalb( l:-a > r:-a ) then
      a := l:-a;
    else
      a := r:-a;
    end if;
    if evalb( l:-b < r:-b ) then
      b := l:-b;
    else
      b := r:-b;
    end if;

    (* refine at the intersecting interval *)
    ll := RefineAt( l, a ):
    ll := RefineAt( ll, b ):
    rr := RefineAt( r, a ):
    rr := RefineAt( rr, b ):

    (* Refiment can change type to rational. *)
    if rr:-isRational_ then
      return CompareRational( ll, rr:-a );
    elif ll:-isRational_ then
      return -CompareRational( rr, ll:-a );
    end if;

    (* Check if there is no intersection after refiment. *)
    if evalb( ll:-b < rr:-a ) then
      return -1;
    elif evalb( ll:-a > rr:-b ) then
      return 1;
    end if;

    (* The number of roots of the GCD of two polynomials is equal to the number of common roots.
       use this to simplify the problem in the intersecting range.*)
    G := gcd( ll:-poly, rr:-poly );
    F1 := simplify( ll:-poly / G );
    F2 := simplify( rr:-poly / G );
    
    if evalb( signum( eval( G, op( indets( G ) ) = ll:-a ) ) <> signum( eval( G,
      op( indets( G ) ) = ll:-b ) ) ) then
      ll := Object( ll, G, ll:-a, ll:-b ):
    else
      ll := Object( ll, F1, ll:-a, ll:-b ):
    end if:

    if evalb( signum( eval( G, op( indets( G ) ) = rr:-a ) ) <> signum( eval( G,
      op( indets( G ) ) = rr:-b ) ) ) then
      rr := Object( rr, G, rr:-a, rr:-b ):
    else
      rr := Object( rr, F2, rr:-a, rr:-b ):
    end if:

    (* Use of GCD can change type to rational. *)
    if rr:-isRational_ then
      return CompareRational( ll, rr:-a );
    elif ll:-isRational_ then
      return -CompareRational( rr, ll:-a );
    end if;

    (* Check for equality. *)
    if evalb( signum( eval( G, op( indets( G ) ) = a ) ) <> signum( eval( G,
      op( indets( G ) ) = b ) ) ) then
      return 0;
    end if;
 
    (* Refiment until disjoitness. *)
    for i from 1 while 1 = 1 do
      ll := BisectRange( ll ):
      rr := BisectRange( rr ):

      (* Rationals after refiment. *)
      if rr:-isRational_ then
        return CompareRational( ll, rr:-a );
      elif ll:-isRational_ then
        return -CompareRational( rr, ll:-a );
      end if;

      (* No intersection after refiment. *)
      if evalb( ll:-b < rr:-a ) then
        return -1;
      elif evalb( ll:-a > rr:-b ) then
        return 1;
      end if:   
    end do:
  end proc:

# Method: < operator
#   A method used to compare two real algebraic numbers.
#
# Parameters:
#   l::RealAlgebraicNumber      - a real algebraic number
#   r::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   true when l is smaller than r and false otherwise.
#
  export `<`::static := proc( l, r, $ )
   if ( _npassed <> 2 or not l::RealAlgebraicNumber or not r::RealAlgebraicNumber ) then
    return false;          
  end if; 
  if Compare( l, r ) = -1 then
    return true;
  else
    return false;
  end if;
  end proc:

# Method: <= operator
#   A method used to compare two real algebraic numbers.
#
# Parameters:
#   l::RealAlgebraicNumber      - a real algebraic number
#   r::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   true when l is smaller or equal to r and false otherwise.
#
  export `<=`::static := proc( l, r, $ )
   if ( _npassed <> 2 or not l::RealAlgebraicNumber or not r::RealAlgebraicNumber ) then
    return false;          
  end if; 
  if Compare( l, r ) <= 0 then
    return true;
  else
    return false;
  end if;
  end proc:

# Method: = operator
#   A method used to compare two real algebraic numbers.
#
# Parameters:
#   l::RealAlgebraicNumber      - a real algebraic number
#   r::RealAlgebraicNumber      - a real algebraic number
#
# Output:
#   true when l is equal to r and false otherwise.
#
  export `=`::static := proc( l, r, $ )          
    if Compare( l, r ) = 0 then
      return true;
    else
      return false;
    end if;
  end proc:

end module:

