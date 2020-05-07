function varargout = runge(fx,meth,n)
%RUNGE Show Interpolation of Runge function
%
% h = RUNGE(FX,METH) show iterpolation by a flash, where METHOD is one of
% 'GL','CHEB','UNI'.They are interpolation of Runge function based on
% Gauss-Lobatto, Chebyshev, and Equispaced points. FX should be a
% function handle.
%
% h = RUNGE(FX,METH, N) give N points interpolation of Runge function. 
%
% h =RUNGE(FX,N) give N points interpolation of Runge function based on Equispaced
% points.
%
%Example:
% 
% h = runge(@(x) sin(pi*x),10);
% h = runge(@(x) sin(pi*x),'cheb');
% h = runge(@(x) 1./(1+25.*x.^2),'uni')
% h = runge(@(x) sign(x), 'cheb')

% Checked: 29-Aug-2017
% $Date: 29-Aug-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com

x = -1:0.01:1;

if nargin < 1
    error('Input not enough.\n')
end

if ( (nargin == 2) && ischar(meth) )
    runge_flash(fx,meth,x)
elseif ((nargin == 2) && isnumeric(meth) )
    n = meth; meth = 'uni';
    runge_interp(fx,meth,n,x)
else
    runge_interp(fx,meth,n,x)
end


varargout = {gcf};

end    



function runge_interp(fx,meth,n,x)

meth = lower(meth);

switch meth
    case 'uni'
    
            xk = linspace(-1,1,n); xk = xk(:);
            vk = baryweights(xk); fvals = fx(xk);
            funi = baryinterp(x, fvals, xk ,vk);
            clf
            plot(x,fx(x),'b-')
            hold on
            plot(x,funi,'r-');
            plot(xk,fvals,'r.');
            title(['Equispaced interpolation of Runge function,',...
                num2str(n), ' points'])

    case 'cheb'
       
            [xk,~,vk] =  chebpts(n); fvals = fx(xk);
            fcheb = baryinterp(x, fvals, xk ,vk);
            
            clf
            plot(x,fx(x),'b-')
            hold on
            plot(x,fcheb,'r-');
            plot(xk,fvals,'r.');
            title(['Chebyshev interpolation of Runge function,',...
                num2str(n), ' points'])
            
            
    case 'gl'
            [xk,~,vk] =  glpts(n); fvals = fx(xk);
            fgl = baryinterp(x, fvals, xk ,vk);
            clf
            plot(x,fx(x),'b-');
            hold on
            plot(x,fgl,'r-');
            plot(xk,fvals,'r.');
            title(['Gauss-Lobatto interpolation of Runge function,',...
                num2str(n), ' points'])
            
            
    case 'gauss'
        [xk,~,vk] =  gausspts(n,1);fvals = fx(xk);
        fgauss = baryinterp(x, fvals, xk ,vk);
        clf
        plot(x,fx(x),'b-');
        hold on
        plot(x,fgauss,'r-');
        plot(xk,fvals,'r.');
        title(['Gauss interpolation of Runge function,',...
            num2str(n), ' points'])
            
 
end

end


function runge_flash(fx,meth,x)

meth = lower(meth);

switch meth
    case 'uni'
        
        for i = 2:40
            xk = linspace(-1,1,i); xk = xk(:);
            vk = baryweights(xk); fvals = fx(xk);
            funi = baryinterp(x, fvals, xk ,vk);
            clf
            plot(x,funi,'r-');
            hold on
            plot(x,fx(x),'b-')
            
            plot(xk,fvals,'r.');
            title('Interpolation on Equispaced points')
            legend(['n=',num2str(i)],'Location','eastoutside')
            
            pause(0.5)
            drawnow
        end
    case 'cheb'
        for i = 2:40
            [xk,~,vk] =  chebpts(i); fvals = fx(xk);
            fcheb = baryinterp(x, fvals, xk ,vk);
            
            clf
            plot(x,fcheb,'r-');
            hold on
            plot(xk,fvals,'r.');
            plot(x,fx(x),'b-')
            title('Interpolation on second Chebyshev points')
            legend(['n=',num2str(i)],'Location','eastoutside')
            
            pause(0.5)
            drawnow
        end
    case 'gl'
        for i = 2:40
            [xk,~,vk] =  glpts(i);fvals = fx(xk);
            fgl = baryinterp(x, fvals, xk ,vk);
            clf
            plot(x,fgl,'r-');
            hold on
            plot(xk,fvals,'r.');
            plot(x,fx(x),'b-')
            title('Interpolation on Gauss-Lobatto points')
            legend(['n=',num2str(i)],'Location','eastoutside')
            
            pause(0.5)
            drawnow
        end      
    case 'gauss'
        for i = 2:40
            [xk,~,vk] =  gausspts(i); fvals = fx(xk);
            fgauss = baryinterp(x, fvals, xk ,vk);
            clf
            plot(x,fgauss,'r-');
            hold on
            plot(xk,fvals,'r.');
            plot(x,fx(x),'b-')
            title('Interpolation on Gauss-Lobatto points')
            legend(['n=',num2str(i)],'Location','eastoutside')
            
            pause(0.5)
            drawnow
        end
end
end

function fx = baryinterp(x, fvals, xk ,vk)
%BARYINTERP Interpolation by barycentric fromular.
%   FX = BARYINTERP(X, FVALK, XK ,VK) returns the interpolation of a
%   function at X, where the baryweights are VK, interpolate points are XK,
%   function values are FVALS.  Note that X, XK and VK should be column
%   vectors, and FVALS, XK, and VK should have the same length.
%
%   If size(FVALS, 2) > 1 then BARYINTERP(X, FVALS) returns values in the
%   form [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
%
%   Example:
%     xk = gausspts(5);
%     vk = baryweights(xk);
%     fvals = eye(5,5);
%     x = -1:0.02:1;
%     fx = baryinterp(x, fvals, xk ,vk);
%     for i = 1:5
%       plot(x,fx(:,i));
%       hold on
%     end
%     plot(xk,zeros(size(xk)),'bo');
%     plot(xk,ones(size(xk)),'rx' );
%




% $Date: 04-Aug-2017$.
% Copyright (c) G.Wang, Email: wangguanjie0@126.com

% This program is modified from chebfun/bary.


% Check that input (X,XK,VK) is a column vector: 

if ( (min(size(xk)) > 1) || (min(size(vk)) > 1) || (min(size(x)) > 1) )
    error('Interpolation:baryinterp', ...
        'XK and VK must be column vector with same length');
end

% get column vector if X, XK and Vk are row vector 
x = x(:); xk = xk(:); vk = vk(:);


if ( length(xk) ~= length(vk) )
    error('Interpolation:baryinterp', ...
        'XK and VK must be column vector with same length');
end

[n, m] = size(fvals);

if ( length(xk) ~= n )
    error('Interpolation:baryinterp', ...
        'FVALS, XK and VK must be column vector with same length');
end

% The function is a constant.
if ( n == 1 )
    fx = repmat(fvals, length(x), 1);
    return
end

% The function is NaN.
if ( any(isnan(fvals)) )
    fx = NaN(length(x), m);
    return
end

% The main loop:
if ( numel(x) < 4*length(xk) )  % Loop over evaluation points
    % Note: The value "4" here was detemined experimentally.

    % Initialise return value:
    fx = zeros(size(x, 1), m);

    % Loop:
    for j = 1:numel(x)
        xx = vk ./ (x(j) - xk);
        fx(j,:) = (xx.'*fvals) / sum(xx);
    end
else                            % Loop over barycentric nodes
    % Initialise:
    num = zeros(size(x, 1), m);
    denom = num;

    % Loop:
    for j = 1:length(xk)
        tmp = (vk(j) ./ (x - xk(j)));
        num = num + bsxfun(@times, tmp, fvals(j,:));
        denom = bsxfun(@plus, denom, tmp);
    end
    fx = num ./ denom;
end

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(x(k) == xk, 1);    % Find the corresponding node
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end


end

function w = baryweights(x)
%BARYWEIGHTS   Barycentric weights.
%   W = BARYWEIGHTS(X) returns scaled barycentric weights for the points in the
%   column vector X. The weights are scaled such that norm(W, inf) == 1.
% 
% This program is adopted from chebfun/baryWeights.
% See http://www.chebfun.org/ for Chebfun information.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% This program is due to Chebfun Developers. 

% Check inout dimension:
[n, m] = size(x);
if ( m > 1 )
    if ( n > 1 )
        error('CHEBFUN:baryWeights:matrix', 'Input must be a vector.')
    else
        % Allow a row vector:
        n = m;
        x = x.';
    end
end

% Capacity:
if ( isreal(x) )
    C = 4/(max(x) - min(x));   % Capacity of interval.
else
    C = 1; % Scaling by capacity doesn't apply for complex nodes.
end

% Compute the weights:
if ( (n < 2001) )              % For small n using matrices is faster.
   V = C*bsxfun(@minus, x, x.');
   V(1:n+1:end) = 1;
   VV = exp(sum(log(abs(V))));
   w = 1./(prod(sign(V)).*VV).';
   
else                           % For large n use a loop.
   w = ones(n,1);
   for j = 1:n
       v = C*(x - x(j)); v(j) = 1;
       vv = exp(sum(log(abs(v))));
       w(j) = 1./(prod(sign(v))*vv);
   end
end

% Scaling:
w = w./norm(w, inf);

end

function [x,w,v] = chebpts(n)
%CHEBPTS compute the second kind chebyshev points.
%   [X,W,V] = CHEBPTS(N) retuns N Chebyshev points of the 1st kind in [-1, 1].
%   X is the second kind Chebshev points, W is the weights of Clenshaw-Curtis
%   quadrature, and V is the for barycentric polynomial interpolation in the
%   Chebyshev points X.
%

% This programe is adapted from Chebfun/chebtech2.chebpts, http://www.chebfun.org/


if n == 1
    x = 0;
    w = 2;
    v = 1;
else
    % Chebyshev points:
    m = n - 1;
    x = sin(pi*(-m:2:m)/(2*m)).';  % (Use of sine enforces symmetry.)
    
    % quadratrue weights
    c = 2./[1, 1-(2:2:(m)).^2];    % Exact integrals of T_k (even)
    c = [c, c(floor(n/2):-1:2)];   % Mirror for DCT via FFT
    w = ifft(c);                   % Interior weights
    w([1,n]) = w(1)/2;             % Boundary weights
    
    % barycentric interpolate weights
    v = [0.5;ones(n-1,1)];        % Note v(1) is positive.
    v(2:2:end) = -1;
    v(end) = .5*v(end);
end

end