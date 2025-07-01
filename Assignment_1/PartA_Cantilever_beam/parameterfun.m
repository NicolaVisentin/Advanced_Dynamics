function err=parameterfun(X,f_opt,G_exp)
%
%-----------------------------------------------------------------
%    err=parameterfun(X,f_opt,G_exp)
%
% This function takes as inputs the parameter vector X for the optimization,
% the frequency range over where we are performing the optimization and the
% set of experimental data points that are needed to compute the error, and
% it is built for defining an handle function in the script:
%
%   err=@(X) parameterfun(X,f_opt,G_exp_k);
%
% that is a function ONLY of X and that can be given as a input for the
% lsqnonlin function that performs the optimisation finding X.
%-----------------------------------------------------------------

% compute the number of transfer functions (basing on length(X)) and the
% length of frequency vector

R=(length(X)-2)/3;
U=length(f_opt);

% compute the error

Omega=2*pi*f_opt;
err=0;
for uu=1:U
    for rr=1:R
        G_num_iter=(X(rr+2)./(-Omega(uu).^2+1j*2*X(2)*X(1)*Omega(uu)+X(1)^2)) + (X(rr+R+2)./(Omega(uu).^2)) + (X(rr+2*R+2));
        err_iter=real(G_exp(uu,rr)-G_num_iter)^2+imag(G_exp(uu,rr)-G_num_iter)^2;
        err=err+err_iter;
    end
end
end