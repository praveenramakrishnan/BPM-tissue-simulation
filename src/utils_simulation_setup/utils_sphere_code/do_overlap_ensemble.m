function [yn] = do_overlap_ensemble(R,x,R_ens,x_ens)
    
%R_ens should be 1 x N
%x_ens should be 3 x N
    
    R_ens=reshape(R_ens,1,numel(R_ens));
    N=numel(R_ens);
    [m,n]=size(x_ens);
    if m~=3
        x_ens=x_ens.';
    end
    x=reshape(x,3,1);
    d_ens=x_ens-x*ones(1,N);
    sep = sqrt(sum(abs(d_ens).^2));
    prod(double(sep>(R+R_ens)));
    yn = prod(double(sep>(R+R_ens)))~=1;
    
