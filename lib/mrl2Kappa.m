function kappa = mrl2Kappa(mrl)
% kappa = mrl2Kappa(mrl)
%  Convert the mean resultant length of a circular distribribution to the
%  concentration parameter kappa of the von-Mises distribution, by using an
%  approximated ML method

%  Input:
%    mrl: mean resultant length

%  Output:
%    kappa: estimated concentration parameter

% Wen-Hao Zhang, Apr-7,2016

kappa = zeros(size(mrl));

Idx = (mrl<0.53);
kappa(Idx) = 2*mrl(Idx) + mrl(Idx).^3 + 5*mrl(Idx).^5/6;

Idx = (mrl>=0.53) & (mrl<0.85);
kappa(Idx) = -.4 + 1.39*mrl(Idx) + 0.43./(1-mrl(Idx));

Idx = (mrl>=0.85);
kappa(Idx) = 1./(mrl(Idx).^3 - 4*mrl(Idx).^2 + 3*mrl(Idx));

% if mrl < 0.53
%   kappa = 2*mrl + mrl.^3 + 5*mrl.^5/6;
% elseif mrl>=0.53 && mrl<0.85
%   kappa = -.4 + 1.39*mrl + 0.43/(1-mrl);
% else
%   kappa = 1/(mrl.^3 - 4*mrl.^2 + 3*mrl);
% end

% if N<15 && N>1
%   if kappa < 2
%     kappa = max(kappa-2*(N*kappa)^-1,0);    
%   else
%     kappa = (N-1)^3*kappa/(N^3+N);
%   end
% end
