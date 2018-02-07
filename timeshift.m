function [tshift,CCt] = timeshift( CCforward, data, dt,fmin,fmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    fmin = 1/dt/length(data);
    fmax = 1/dt-1/dt/length(data);
end

npoint = length(data);

win = zeros(npoint,1);
band = fix(fmin * dt * npoint) +1 : fix(fmax * dt * npoint) +1; band = band';
nband = fix(fmax * dt * npoint)-fix(fmin * dt * npoint) +1 ;
win(band) = tukeywin(nband, 0.1);

CCforward =CCforward ./ max(abs(CCforward));%;normc(CCforward)
data = data./ max(abs(data));%;normc(data)
forwf = fft(ifftshift(CCforward));
dataf = fft(ifftshift(data));
CCf = forwf .* conj(dataf);
CCt = real( ifft(CCf .* win));
CCt = fftshift(CCt);
[~,tshift] = max(CCt);
tshift = ( tshift - (npoint+1)/2 ) * dt;
end

