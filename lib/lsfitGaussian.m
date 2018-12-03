function Param = lsfitGaussian(X, P, Param0)
% least-square fit of neuro-metric funtion by using the cumulative Gaussian
% function

%-----------INPUT ----------
% X: the x-axis data
% P: the empirical density function
% Param0: the initial parameters.

% ---------OUTPUT -------
% Param: the 1-by-2 vector, 1st element: mu; 2nd element: sigma

% Author: Wenhao Zhang, July-23-2013


% parse the input
if nargin < 3
Param0(1) = 0;
Param0(2) = 1;
end

options = optimset('TolX', 1e-7, 'display', 'off');
[Param, err] = fminsearch(@(x) F(x), Param0, options);

    function err = F(Param)
        y = normpdf(X, Param(1), Param(2));
        err = sum((y-P).^2);
    end
end