# Procedure: EliminationResultant
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials in three variables
#   vars       - variables to be eliminated
#
# Output:
#   Univariate polynomial obtained from S in the first variable.
EliminationResultant := proc( S::~set, vars::~list )
  option cache:
  local p, var, res := S, permm;
  if nops(S) < 2 then
    error "Wrong size of the input set. Expected size is at least 2.";
  fi;
  if not andmap(type, S, `polynom`) then
    error "Wrong type of elements. Expected argument is a set of polynomials!"; 
  fi;
  if nops(vars) <= 1 then
    error "Wrong number of indeterminates. It should be 3.";
  fi;

   for var in vars[2..] do
       permm := combinat:-permute(res, 2);
       res := [];
     for p in permm do
       res := [op(res), OneVariableElimination(p[1], p[2], var)]:
     od;
   od;

  return foldl( gcd, op(res) ):
end proc:

# Procedure: RemoveExponant
#    Removes exponants in an expression
#
# Parameters:
#    r - expression, the expression to simplify
#
# Output:
#    An arithmetic expression that has the same squarefree part as r
RemoveExponants := proc(r)
        local remove_exponant, sqrr, result;
        remove_exponant := e -> if type(e,`^`) then op(1,e) else e end if;
        result := remove_exponant(r);
        if type(r,`*`) then
            result := map(remove_exponant, result);
        end if;
        return result;
end proc;

OneVariableElimination := proc( p, q, v)
    local r;
    if degree(p,v)>0 or degree(q,v)>0 then
        r := resultant(p, q, v); 
        r := RemoveExponants(r);
        return r;
    elif  nops(indets({p,q}))=1 then
        return gcd(p, q);
    else
        return p;
    end if;
end proc;

# Procedure: EliminationResultant2
#   Computes univariate polynomial.
#
# Parameters:
#   S          - a set of polynomials of degree 2 in three variables
#
# Output:
#   Univariate polynomial obtained from S in the first variable, with
#   formula from Chapter 3, p. 89 of "Using Algebraic Geometry" from Cox,
#   Little, O'Shea.
EliminationResultant2 := proc( S::~set )
  #option cache;
  local vars, monomials, L, J, J1, J2, dJ, result;
  vars := indets(S);
  if nops(S) <> 3 then
    error "Wrong size of the input set. Expected size is 3.";
  fi;
  if not type(S[1], polynom) or not type(S[2], polynom) or not type(S[3], polynom) then
    error "Wrong type of elements. Expected argument is a set of polynomials!"; 
  fi;
  if nops(vars) <> 3 then
    error "Wrong number of indeterminates. It should be 3.";
  fi;
  vars := [op(vars)][2..3];
  monomials := map2(map,`-`,[op(combinat:-composition(5,3))], 1)[..,2..3];
  L := map(p->[diff(p,vars[1]), diff(p,vars[2]),p], [op(S)]);
  J := LinearAlgebra:-Determinant(L);
  J1 := diff(J,vars[1]);
  J2 := diff(J,vars[2]);
  dJ := [J1, J2, degree(J,vars)*J-vars[1]*J1 - vars[2]*J2];
  L := map2(map2,(p,d)-> coeftayl(p, vars=[0,0], d),
                 [op(S),op(dJ)], monomials);
  result := -LinearAlgebra:-Determinant(L)/64;
  return result;
end proc:
