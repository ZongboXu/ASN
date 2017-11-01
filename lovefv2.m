function fv = lovefv2
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu1 = 8E7;
mu2 = 1.28E9;
beta1 = 200; beta2 = 800; H = 40;

dc = 0.00001;
middv = dc/10;
fres = 4;
fv = zeros(length(fres),3);
fv(:,1) = fres';

fun = @(vel,omega) tan(omega*H*sqrt(1/(beta1^2)-1/(vel^2)))-mu2/mu1*sqrt(1/(vel^2)-1/(beta2^2))/sqrt(1/(beta1^2)-1/(vel^2));
tanomega = @(vel,omega) tan(omega*H*sqrt(1/(beta1^2)-1/(vel^2)));

for ifre = 1 : length(fres) 
    omega = 2 * pi *fres(ifre);
    imode = 0;
    for c = 1+beta1 : dc : beta2 - dc
        a = c; b = c + dc;
        funa = fun(a,omega);
        funb = fun(b,omega);
        if funa*funb > 0 
            continue;
        elseif tanomega(b,omega) - tanomega(a,omega)>0
            while abs(a-b) > middv
              midab = (a+b)/2;
              funmid = fun(midab,omega);
              if funa * funmid > 0 
                  a = midab;
              else
                  b = midab;
              end
            end
            fv(ifre,imode+2) = (a+b)/2;
            imode = imode + 1;
            if imode > 2 
                break
            end
            continue;
        end
    end
            
end

end