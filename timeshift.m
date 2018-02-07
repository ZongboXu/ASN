function [tshift,CCt] = timeshift( CCforward, data, dt,twin,fwin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 5
    fmin = fwin(1); fmax = fwin(2);
else
    fmin = 1/dt/length(data);
    fmax = 1/dt-1/dt/length(data);
end

if nargin == 4
    tmin = twin(:,1);tmax = twin(:,2);
else
  tmin = -(length(data)-1)/2 * dt * ones(size(data,2));
  tmax = (length(data)-1)/2 * dt * ones(size(data,2));
end


npoint = length(data);

%% time window
ttaper = zeros(size(data));
npt0 = (npoint+1)/2; 
ntrace = size(data,2);

for itrace = 1 : ntrace
   nptmin = tmin(itrace)/dt + npt0;
   nptmax = tmax(itrace)/dt + npt0;
   ttaper(nptmin:nptmax,itrace) = tukeywin(nptmax-nptmin+1,0.1);
end

data = data .* ttaper;
%% frequency window
ftaper = zeros(npoint,1);
band = fix(fmin * dt * npoint) +1 : fix(fmax * dt * npoint) +1; band = band';
nband = fix(fmax * dt * npoint)-fix(fmin * dt * npoint) +1 ;
ftaper(band) = tukeywin(nband, 0.1);

CCforward = CCforward./max(abs(CCforward));
data = data./max(abs(data));
forwf = fft(ifftshift(CCforward));
dataf = fft(ifftshift(data));
CCf = forwf .* conj(dataf);
CCt = real( ifft(CCf .* ftaper));
CCt = fftshift(CCt);
[~,tshift] = max(CCt);
tshift = ( tshift - (npoint+1)/2 ) * dt;
end
