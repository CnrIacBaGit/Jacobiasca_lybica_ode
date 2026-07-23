function O = observed_2002()
% OBSERVED_2002  Field data for the 2002 Feudo Arancio season
%
% Risk thresholds (nymphs/leaf): low 0.5, medium 1.0, high 2.0.
% Returns observed nymph-risk crossing dates 

O = struct();
O.riskNames = ["low","med","high"];
O.riskVals  = [0.5, 1.0, 2.0];

O.datesN = datetime(2002, [5 5 6 6 6 7 7 8 8 9 9 10 10], ...
                          [5 19 2 16 30 14 28 11 25 8 22 6 20]).';
O.nymphsAllLeaf = [0.00 0.00 0.08 0.12 0.20 0.25 0.27 0.28 0.61 2.00 4.10 2.00 0.43].';

O.datesA = datetime(2002, [5 6 6 7 7 7 8 8 9 9 10 10], ...
                          [22 5 19 3 17 31 14 28 11 25 9 23]).';
O.adultsTrap = [0 0 0 10 35 80 150 140 200 265 230 190].';

% interpolated observed crossings (same rule as the original solver)
O.crossings = struct();
for ir = 1:numel(O.riskVals)
    thr = O.riskVals(ir);
    tc = NaT;
    for k = 1:numel(O.nymphsAllLeaf)-1
        if O.nymphsAllLeaf(k) < thr && O.nymphsAllLeaf(k+1) >= thr
            frac = (thr - O.nymphsAllLeaf(k)) / ...
                   (O.nymphsAllLeaf(k+1) - O.nymphsAllLeaf(k));
            dd = round(frac * days(O.datesN(k+1) - O.datesN(k)));
            tc = O.datesN(k) + days(dd);
            break
        end
    end
    O.crossings.(O.riskNames(ir)) = tc;
end
end
