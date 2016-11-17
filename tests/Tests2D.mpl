read "../RealAlgebraicNumber.mpl";
read "../RigidMotionsParameterSpaceCommon.mpl";
read "../RigidMotionsParameterSpaceDecompositionRecursive.mpl";

# Table of examples to test
quadrics := table();

quadrics[1] := [x*y-1];
quadrics[2] := [x+y-1/4];
quadrics[3] := [x+y-1/9];
quadrics[4] := [x^2+y^2-1];
quadrics[5] := [x*y-3*x-3*y+8, x-3-y/10^50];
quadrics[6] := [x*y-3*x-3*y+8, y-3-x/10^50];
quadrics[7] := [x*y-1, x, x-1]; # error expected
quadrics[8] := [x-1, y-1];
quadrics[9] := [x^2-2*x-2, y^2-2*y-2];
quadrics[10] := [x*y-3*x-3*y+8, x-3, y-5/2];
quadrics[11] := [x[1]^2+x[2]^2-1];

# Fix the seed for reproductible random generation
RandomTools:-MersenneTwister:-SetState(state=123456789);
quadrics[12] := [randpoly $ 20]([x,y], dense, degree=2, coeffs=rand(0..1));

#stopat(LaunchOnGridComputeSamplePoints2D);
#stopat(ComputeEventsAlgebraicNumbers2D);

for i from 1 to 12 do
    print(i);
    try
        LaunchOnGridComputeSamplePoints2D(quadrics[i], 0, 1, false, i);
    catch "Only system in two variables":
        print(lasterror);
    catch:
        print("Error", lastexception);
    end try;
end do;
