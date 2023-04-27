function [bandwidth, density, xmesh, cdf]=kde_DIY(data, prsistntParams, MIN, MAX)
% At least 3 times faster than the original code!

% Reliable and extremely fast kernel density estimator for one-dimensional data;
%        Gaussian kernel is assumed and the bandwidth is chosen automatically;
%        Unlike many other implementations, this one is immune to problems
%        caused by multimodal densities with widely separated modes (see example). The
%        estimation does not deteriorate for multimodal densities, because we never assume
%        a parametric model for the data.
% INPUTS:
%     data    - a vector of data from which the density estimate is constructed;
%          n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12;
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where Range=max(data)-min(data);
% OUTPUTS:
%   bandwidth - the optimal bandwidth (Gaussian kernel assumed);
%     density - column vector of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%             - If no output is requested, then the code automatically plots a graph of
%               the density estimate.
%        cdf  - column vector of length 'n' with the values of the cdf
%  Reference:
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957.

%
%  Example:
%           data=[randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%              kde(data,2^14,min(data)-5,max(data)+5);

% Copyright (c) 2015, Zdravko Botev
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The University of New South Wales nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


data=data(:); %make data a column vector
N=length(unique(data));
if nargin<2 % if n is not supplied switch to the default
    n = 2^14;
    I = (1:n-1).'.^2;
    len = 7;
    L = 2 .* (pi.^2 .* I).^(1:len);
    R = -I.*pi^2;
    tBase = (1 + 0.5.^((1:len) + 0.5)).*cumprod(1:2:2*len - 1)./(3*N*sqrt(0.5*pi));
    tPower = 2./(3 + 2.*(1:len));
    tBaseMin = 2*N*sqrt(pi);
elseif isstruct(prsistntParams)
    n = prsistntParams.n;
    %     I = prsistntParams.I;
    len = prsistntParams.len;
    L = prsistntParams.L;
    R = prsistntParams.R;
    tBase = prsistntParams.tBase;
    tPower = prsistntParams.tPower;
    tBaseMin = prsistntParams.tBaseMin;
else
    n = 2^ceil(log2(prsistntParams));
    I = (1:n-1).'.^2;
    len = 7;
    L = 2 .* (pi.^2 .* I).^(1:len);
    R = -I.*pi^2;
    tBase = (1 + 0.5.^((1:len) + 0.5)).*cumprod(1:2:2*len - 1)./(3*N*sqrt(0.5*pi));
    tPower = 2./(3 + 2.*(1:len));
    tBaseMin = 2*N*sqrt(pi);
end

if nargin<4 %define the default  interval [MIN,MAX]
    minimum=min(data); maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/2; MAX=maximum+Range/2;
end
% set up the grid over which the density estimate is computed;
Range=MAX-MIN; dx=Range/(n-1); xmesh=MIN+(0:dx:Range);
%bin the data uniformly using the grid defined above;
initial_data=histc(data,xmesh)/N;  initial_data=initial_data/sum(initial_data);
a=dct1d(initial_data); % discrete cosine transform of initial data
% now compute the optimal bandwidth^2 using the referenced method
a2=(a(2:end)/2).^2;
L = (L.*a2).';
% use  fzero to solve the equation t=zeta*gamma^[5](t)
t_star=root(@(t)fixed_point(t, len, L, R, tBase, tPower, tBaseMin),N);
% smooth the discrete cosine transform of initial data using t_star
a_t=a.*[1; exp(R*t_star/2)];
% now apply the inverse discrete cosine transform
if (nargout>1)||(nargout==0)
    density=idct1d(a_t)/Range;
end
% take the rescaling of the data into account
bandwidth=sqrt(t_star)*Range;
density(density<0)=eps; % remove negatives due to round-off error
if nargout==0
    figure(1), plot(xmesh,density)
end
% for cdf estimation
if nargout>3
    f=L(1, :)*exp(R*t_star);
    t_cdf=(sqrt(pi)*f*N)^(-2/3);
    % now get values of cdf on grid points using IDCT and cumsum function
    a_cdf=a.*[1; exp(R*t_cdf/2)];
    cdf=cumsum(idct1d(a_cdf))*(dx/Range);
    % take the rescaling into account if the bandwidth value is required
%     bandwidth_cdf=sqrt(t_cdf)*R;
end

end
%################################################################
function  out=fixed_point(t, len, L, R, tBase, tPower, tBaseMin)
% this implements the function t-zeta*gamma^[l](t)

f = L(len, :)*exp(R*t);
for s = len-1:-1:2
    f = L(s, :)*exp(R*((tBase(s)/f)^tPower(s)));
end
out = t - (tBaseMin*f)^(-tPower(1));

end



%##############################################################
function out = idct1d(data)

% computes the inverse discrete cosine transform
nrows=size(data, 1);
% Compute weights
weights = nrows*exp(1i*(0:nrows-1)*pi/(2*nrows)).';
% Compute x tilde using equation (5.93) in Jain
data = real(ifft(weights.*data));
% Re-order elements of each column according to equations (5.93) and
% (5.94) in Jain
out = zeros(nrows,1);
out(1:2:nrows) = data(1:nrows/2);
out(2:2:nrows) = data(nrows:-1:nrows/2+1);
%   Reference:
%      A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
end
%##############################################################

function data=dct1d(data)
% computes the discrete cosine transform of the column vector data
nrows = size(data, 1);
% Compute weights to multiply DFT coefficients
weight = [1;2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
% Re-order the elements of the columns of x
data = [ data(1:2:end,:); data(end:-2:2,:) ];
% Multiply FFT by weights:
data= real(weight.* fft(data));
end

function t=root(f,N)
% try to find smallest root whenever there is more than one
N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
tol=10^-12+0.01*(N-50)/1000;
flag=0;
while flag==0
    try
        t=fzero(f,[0,tol]);
        flag=1;
    catch
        tol=min(tol*2,.1); % double search interval
    end
    if tol==.1 % if all else fails
        t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
    end
end
end





