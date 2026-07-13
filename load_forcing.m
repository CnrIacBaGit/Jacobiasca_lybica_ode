function F = load_forcing(matPath)
% LOAD_FORCING  Load the 2000-2025 daily environmental forcing and clean
% sentinel values. The distributed archive spans 1 Jan 2000 to 1 Dec 2025
% (9467 daily records) and has been cleaned. The conversion/interpolation
% below is retained as a defensive safeguard for replacement NASA POWER files
% that may still contain the standard -999 missing-data code.
%
% OUTPUT F: struct with
%   .dates (datetime vector, day-resolution, starting 2000-01-01)
%   .T .H .R .W  (double column vectors, cleaned)
%   .year_of(k)  helper is not stored; use year(F.dates).

S = load(matPath);
req = {'Tvec','Hvec','Rvec','Wvec'};
for i = 1:numel(req)
    if ~isfield(S, req{i})
        error('load_forcing:missingVar','%s not found in %s.', req{i}, matPath);
    end
end
T = S.Tvec(:); H = S.Hvec(:); R = S.Rvec(:); W = S.Wvec(:);
n = numel(T);
if any([numel(H) numel(R) numel(W)] ~= n)
    error('load_forcing:lengthMismatch','Forcing vectors have unequal length.');
end

dates = (datetime(2000,1,1) + caldays(0:n-1)).';

% sentinel -> NaN (any value <= -900 in any channel), then interpolate
T = clean_channel(T, dates, 'T');
H = clean_channel(H, dates, 'H');
R = clean_channel(R, dates, 'R');
W = clean_channel(W, dates, 'W');

F = struct('dates',dates,'T',T,'H',H,'R',R,'W',W);
fprintf('[load_forcing] %d daily records, %s to %s.\n', n, ...
    datestr(dates(1),'dd-mmm-yyyy'), datestr(dates(end),'dd-mmm-yyyy'));
end

function x = clean_channel(x, dates, name)
bad = x <= -900;
nb = nnz(bad);
if nb > 0
    fprintf('[load_forcing] %s: %d sentinel value(s) -> NaN + interpolate (e.g. %s).\n', ...
        name, nb, datestr(dates(find(bad,1)),'dd-mmm-yyyy'));
    x(bad) = NaN;
    good = ~isnan(x);
    if nnz(good) >= 2
        x(~good) = interp1(find(good), x(good), find(~good), 'linear', 'extrap');
    end
end
end
