function Param = lsNeuroMetric(X, P, Param0)
% least-square fit of neuro-metric funtion by using the cumulative Gaussian
% function

%-----------INPUT ----------
% X: the x-axis data
% P: the empirical cumulative probability, whichi is interpreted as the
% correct ratio
% Param0: the initial parameters.

% ---------OUTPUT -------
% Param: the 1-by-2 vector, 1st element: mu; 2nd element: sigma

% Author: Wenhao Zhang, July-23-2013


% parse the input
if nargin < 3   
    [P, Idx, ~] = unique(P);
    X = X(Idx);
    Param0(1) = interp1(P, X, 0.5);
    Idx = ceil(length(X)/2);
    Idx = (Idx-3: Idx+2);
    Idx(Idx==0) = [];
    tt= polyfit(X(Idx), P(Idx), 1);
    Param0(2) = 1/(2*sqrt(2)*tt(1));
    clear tt
end

if isnan(Param0(2)) || Param0(2)<0
    Param0(2) = 10;
end
options = optimset('TolX', 1e-7, 'display', 'off');
% [Param, err] = fminsearch(@(x) F(x), Param0, options);

A = [0, -1];
b = 0;
[Param, err] = fmincon(@(x) F(x), Param0, A, b, [], [],[],[],[], options);

    function err = F(Param)
        y = normcdf(X, Param(1), Param(2));
        err = sum((y-P).^2);
    end
end