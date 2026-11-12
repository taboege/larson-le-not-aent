-- Run with M2 --script algebraic_dual.m2 (requires the Matroids package).
-- Output files are written to the current directory.
-- Source: Larson--Le, arXiv:2609.02834v1, Example 1.3 (printed page 3).
-- No list of expected bases/nonbases is used in this computation.

needsPackage "Matroids";
print("Macaulay2 version: " | toString version#"VERSION");
print("Matroids version: " | toString (options Matroids)#Version);
started = currentTime();

-- The order of this list is the paper's labeling 1,...,10.
R = (ZZ/2)[t_1..t_6, x_1..x_10];
parameters = take(gens R, 6);
coordinates = drop(gens R, 6);
functions = {
    t_1,
    t_2,
    t_3,
    t_1*t_2*t_3*(t_4+t_6),
    t_1*t_3*(t_4+t_5+t_6),
    t_1^2*t_2*t_3*(t_4+t_5),
    t_1*t_2*(t_5+t_6),
    t_4,
    t_5,
    t_6
};
parametrizationGraph = ideal apply(10, i -> coordinates#i - functions#i);
relationsInR = eliminate(parameters, parametrizationGraph);

-- Move the elimination ideal to the ring of the ten output coordinates.
-- The parameters are sent to zero only AFTER their elimination.
S = (ZZ/2)[y_1..y_10];
toCoordinates = map(S, R, apply(6, i -> 0_S) | gens S);
coordinateIdeal = toCoordinates relationsInR;
assert(dim(S/coordinateIdeal) == 6);
-- An additional transparent check: six outputs are the original parameters,
-- so the full relation ideal should be this four-generator graph ideal.
assert(coordinateIdeal == ideal(
    y_4-y_1*y_2*y_3*(y_8+y_10),
    y_5-y_1*y_3*(y_8+y_9+y_10),
    y_6-y_1^2*y_2*y_3*(y_8+y_9),
    y_7-y_1*y_2*(y_9+y_10)
));
kernelFile = openOut "kernel.m2";
kernelFile << "S = (ZZ/2)[y_1..y_10];" << endl;
kernelFile << "coordinateIdeal = " << toString coordinateIdeal << ";" << endl;
close kernelFile;
print("Eliminated the six parameters. dim(S/coordinateIdeal) = " | toString dim(S/coordinateIdeal));

-- For each six-subset B, eliminate ALL output variables outside B.
-- The resulting ideal is zero iff the six corresponding functions are
-- algebraically independent. This works in characteristic two, including
-- inseparable dependencies; no Jacobian or random specialization is used.
ground = toList(0..9);
csv = L -> concatenate between(",", apply(L, toString));
oneBased = B -> apply(sort toList B, i -> i+1);
sixSubsets = subsets(ground, 6);
eliminationFile = openOut "eliminations.tsv";
eliminationFile << "subset\tindependent\telimination_ideal" << endl;
independentSix = {};
scan(#sixSubsets, k -> (
    B := sixSubsets#k;
    outside := toList(set ground - set B);
    eliminated := eliminate(apply(outside, i -> S_i), coordinateIdeal);
    independent := eliminated == ideal(0_S);
    if independent then independentSix = append(independentSix, B);
    eliminationFile << csv(oneBased B) << "\t" << toString independent
        << "\t" << toString eliminated << endl;
    if (k+1)%30 == 0 then print("Elimination tests: " | toString(k+1) | "/210");
));
close eliminationFile;
assert(#sixSubsets == 210);

-- Matroids stores subsets using indices 0,...,9. Convert to 1,...,10 ONLY
-- when exporting, preserving the paper's labels without an isomorphism test.
N = matroid(10, independentSix);
assert(isWellDefined N);
assert(rank N == 6);
D = dual N;
assert(isWellDefined D);
assert(rank D == 4);
computedBases = sort apply(bases D, oneBased);
assert(#computedBases == #independentSix);
basisFile = openOut "dual_bases.csv";
basisFile << "a,b,c,d" << endl;
scan(computedBases, B -> basisFile << csv B << endl);
close basisFile;
computedNonbases = sort apply(nonbases D, oneBased);
nonbasisFile = openOut "dual_nonbases.csv";
nonbasisFile << "a,b,c,d" << endl;
scan(computedNonbases, B -> nonbasisFile << csv B << endl);
close nonbasisFile;
print("rank N = " | toString(rank N) | "; bases N = " | toString(#(bases N)));
print("rank dual N = " | toString(rank D) | "; nonbases dual N = " | toString(#computedNonbases));
print("Wrote all " | toString(#computedBases) | " dual bases to dual_bases.csv.");
print("Computed dual nonbases (paper labels):");
scan(computedNonbases, B -> print toString B);
print("Elapsed seconds: " | toString(currentTime()-started));
