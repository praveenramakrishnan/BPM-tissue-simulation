%function [g] = sawtooth(x,T)
%
%Calculate a sawtooth function as a function of x, and period 2*T

function [g] = sawtooth_extend_data(x,T)
    y1 = mod(x,T);
    y2 = mod(floor(x/T),2);
    y3 = mod(floor(x/T),2).*floor(x/T);
    
    g1 = y1 - y1.*y2;
    g2 = y1 - g1;
    
    
    g = g1 - g2 + (abs(y3)>0)*T;
  
    inds = find(mod(x,T)==0);
    res = x/T;
    res = res(inds);
    newval = zeros(size(res));
    
    newval(find(mod(res,2)==1)) = T;
    g(inds) = newval;
    
    %g(find(mod(x,T)==0))=T;
