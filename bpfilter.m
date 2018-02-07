function filtsig = bpfilter(signal,flow,fhigh,sampler)

nlen = length(signal);
df = sampler/nlen;

iflow = fix(flow/df) + 1;
ifhigh = fix(fhigh/df)+1;

filter = zeros(nlen,1);
filter(iflow:ifhigh) = tukeywin(ifhigh-iflow+1,0.1);
filter(nlen-ifhigh+1:nlen-iflow+1) = tukeywin(ifhigh-iflow+1,0.1);

filtsig = ifft( fft(signal) .* filter );
filtsig = real(filtsig);

end
