# Class: Event
#
# Description:
#  Implementation of events -- a real algebraic number and list of associated quadrics.
#
# Author:
#  Kacper Pluta - kacper.pluta@esiee.fr
#  Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
#
# Date:
#  29/12/2016
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
Event := module()
  option object;
   
  (*Real Algebraic Number*)
  local ranum::RealAlgebraicNumber;
  
  (*List of indices of associeted quadrics*)
  local quadrics::list;


# Method: ModuleCopy
#   Standard constructor / copy constructor
#
# Parameters:
#   self::Event                   - a new object to be constructed
#   proto::Event                  - a prototype object from which self is derived
#   ranum::RealAlgebraicNumber    - a real algebraic number
#   quadrics::rational            - a list of indices which corresponds to associated quadrics
#
# Output:
#   An object of type RealAlgebraicNumber.
#
# Exceptions:
#  "Invalid range. A range is valid when: a <= b."
#  "Degree of %1 is invalid."
#
  export ModuleCopy::static := proc( self::Event,
                                     proto::Event,
                                     ranum::RealAlgebraicNumber,
                                     quadrics::list, $ )
    if _passed = 2 then
      self:-ranum := proto:-ranum;
      self:-quadrics := proto:-quadrics;
    else
      self:-ranum := ranum;
      self:-quadrics := quadrics;
    fi;
  end proc;


# Method: ModulePrint
#   Standard printout of an object of type Event.
#
# Parameters:
#   self::Event                  - an event
#
  export ModulePrint::static := proc( self::Event )
   nprintf("(%a, %a)",self:-ranum, self:-quadrics); 
  end proc;


# Method: ModuleApply
#   Define standard constructor.
#
  export ModuleApply::static := proc()
   Object(Event, args)
  end proc;


# Method: ModuleDeconstruct
#   Provides information how to recreate an object after being serialized.
#
# Parameters:
#   self::Event                  - an event
#
  export ModuleDeconstruct::static := proc( self::Event )
    ('Event')(('RealAlgebraicNumber')(GetPolynomial(self:-ranum), GetInterval(self:-ranum)[1],
    GetInterval(self:-ranum)[2]), self:-quadrics)
  end proc;

# Method: GetRealAlgebraicNumber
#   Getter of the real algebraic number.
#
# Parameters:
#   self::Event                  - an event
# Output:
#   Associated real algebraic number.
  export GetRealAlgebraicNumber::static := proc(self::Event)
    return self:-ranum;
  end proc;


# Method: GetQuadrics
#   Getter of the list of indices.
#
# Parameters:
#   self::Event                  - an event
# Output:
#   Associated quadrics' indices.
  export GetQuadrics::static := proc(self::Event)
    return self:-quadrics;
  end proc;


# Method: Compare
#   A method used to compare two events.
#
# Parameters:
#   l::Event      - an event
#   r::Event      - an event
#
# Output:
#   -1 when l is smaller than r, 0 when they are equal and 1 when l is bigger than r.
#
  export Compare::static := proc( l::Event, r::Event, $ )
    return Compare(l:-ranum, r:-ranum);
  end proc;

end module;

