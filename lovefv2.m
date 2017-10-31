function fv = lovefv2
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu1 = 8E7;
mu2 = 1.28E9;
beta1 = 200; beta2 = 800; H = 40;

dc = 0.9;
middv = 0.1;
fres = 1:1:100;
fv = zeros(length(fres),2);
fv(:,1) = fres';

fun = @(vel,omega) tan(omega*H*sqrt(1/(beta1^2)-1/(vel^2)))-mu2/mu1*sqrt(1/(vel^2)-1/(beta2^2))/sqrt(1/(beta1^2)-1/(vel^2));

for ifre = 1 : length(fres) 
    fre = fres(ifre);
    omega = 2 * pi *fre;
    for c = 1+beta1 : dc : beta2 - dc
        a = c; b = c + dc;
        funa = fun(a,omega);
        funb = fun(b,omega);
        if funa*funb > 0
            continue;
        else
            while abs(a-b) > middv
              midab = (a+b)/2;
              funmid = fun(midab,omega);
              if funa * funmid > 0
                  a = midab;
              else
                  b = midab;
              end
            end
            fv(ifre,2) = (a+b)/2;
            break;
        end
    end
            
end

end

