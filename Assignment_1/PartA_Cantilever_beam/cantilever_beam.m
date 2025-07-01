%% Data definition

clear
close all
clc

% geometrical data [m]

L=1.2;
b=40e-3;
h=8e-3;

% physical data [kg], [m]

rho=2700;
E=68000e6;
xi=0.01;    % adimensional damping ratio

m=rho*b*h; % mass per unit length
J=b*h^3/12; % inertia

% domains definition

f=linspace(0,200,1e5); % frequencies [Hz]
w=2*pi*f;
gamma=(m*w.^2/E/J).^(1/4);

x=linspace(0,L,1001);   % space [m]

%% POINT a) EIGENMODES AND EIGENSHAPES

%% Computation of the eigenmodes (1) (via Matlab)

% define matrix H representing the system of boundary conditions

H=@(gamma) [1 0 1 0;
            0 1 0 1;
            sin(gamma*L) -cos(gamma*L) sinh(gamma*L) cosh(gamma*L);
            -cos(gamma*L) -sin(gamma*L) cosh(gamma*L) sinh(gamma*L)];

% compute the determinant for different values of gamma (w)

detH=zeros(1,length(gamma));
for ii=1:length(gamma)
    gamma_ii=gamma(ii);
    H_ii=H(gamma_ii);
    detH(ii)=det(H_ii);
end

% impose detH=0 and find the N natural frequencies and (allocate them into
% w_0)

[~,indx]=findpeaks(-log(abs(detH)),'MinPeakProminence',2);
w_0=w(indx);
gamma_0=gamma(indx);
f_0=f(indx);
N=length(indx);

fprintf(['Natural frequencies [Hz] are:\n', num2str(f_0,5),'\n\n'])

% plot the determinant as a function of the frequecy (just to visualize
% it)

figure
semilogy(f,abs(detH))
hold on
semilogy(f_0,abs(detH(indx)),'r*')
grid on
xlabel('f [Hz]')
ylabel('|detH|')
title('Determinant of H as a function of the frequency')
legend('','natural frequencies','location','southeast')
hold off

%% Computation of the eigenmodes (2) (determinant computed with Wolfram)

% % insert expression of the determinant (computed via wolfram)
% 
% detH_wolf=@(gamma) -2*(cos(gamma*L).*cosh(gamma*L)+1);
% detH_wolf=detH_wolf(gamma);
% 
% % find where it is zero
% 
% [~,indx]=findpeaks(-log(abs(detH_wolf)));
% w_0=w(indx);
% gamma_0=gamma(indx);
% f_0_wolf=f(indx);
% N=length(indx);
% 
% fprintf(['Natural frequencies [Hz] are:\n', num2str(f_0_wolf,5),'\n\n'])

%% Computation of modeshapes' coefficients A, B, C, D

% actually we cannot find A, B, C, D, but let's impose A=1 and define
% B_hat=B/A, C_hat=C/A, D_hat=D/A (allocated in a vector X_hat).
% solve the reduced system.

coeff=zeros(4,N); % initialise matrix containing the ceoff 1, B_hat, C_hat, D_hat (each column is associated to a mode)
coeff(1,:)=1;
for k=1:N
    gamma_k=gamma_0(k); % compute "natural gamma" for k-th mode
    H_k=H(gamma_k); % compute H matrix for k-th mode
    H_hat_k=H_k(2:end,2:end);   % extract reduced H matrix for k-th mode
    N_k=H_k(2:end,1);   % extract reduced N vector for k-th mode

    X_hat_k=-H_hat_k\N_k;   % compute H_hat for the k-th mode

    coeff(2:end,k)=X_hat_k; % allocate H_hat in the k-th column
end

%% Computation and plot of modeshapes

% allocate modeshapes in a column vector

modeshapes=zeros(length(x),N);
modeshapes_norm=zeros(length(x),N);  
phi=@(b,c,d,g,x) cos(g*x)+b*sin(g*x)+c*cosh(g*x)+d*sinh(g*x);   % general form of a modeshape
for k=1:N
    gamma_k=gamma_0(k); % compute parameters for k-th mode
    B_hat_k=coeff(2,k);
    C_hat_k=coeff(3,k);
    D_hat_k=coeff(4,k);

    phi_k=phi(B_hat_k,C_hat_k,D_hat_k,gamma_k,x); % compute k-th modeshape
    phi_k_norm=phi_k./max(abs(phi_k)); % normalise the modeshapes, for the plot

    modeshapes(:,k)=phi_k;  % allocate the modeshape in the matrix
    modeshapes_norm(:,k)=phi_k_norm;
end

% plot modeshapes

figure

if floor(N/2)==N/2

    for k=1:N
        subplot(N/2,N/2,k)
        plot(x,modeshapes_norm(:,k),'r')
        hold on
        plot(x,zeros(1,length(x)),'k','LineWidth',2)
        axis([0 1.1*x(end) -1.5 1.5])
        grid on
        xlabel('x [m]')
        title(['Modeshape for \omega_',num2str(k),'=',num2str(w_0(k),4),' Hz'])
        hold off
    end

else
    
    for k=1:N
        subplot(ceil(N/2),floor(N/2),k)
        plot(x,modeshapes_norm(:,k),'r')
        hold on
        plot(x,zeros(1,length(x)),'k','LineWidth',2)
        axis([0 1.1*x(end) -1.5 1.5])
        grid on
        xlabel('x [m]')
        title(['Modeshape for \omega_',num2str(k),'=',num2str(w_0(k),4),' Hz'])
        hold off
    end

end

% modeshapes animation

n_modeshapes=size(modeshapes_norm,2);
k=str2double(inputdlg(['Enter the mode shape to display (1 to ',num2str(n_modeshapes),'):']));

figure
hold on
grid on
title(['Modeshape for \omega_',num2str(k),'=',num2str(w_0(k),4),' Hz'])

plot(x,modeshapes_norm(:,k),':k','LineWidth',2)
h1=plot(x,zeros(size(x)),'LineWidth',2);
xlabel('x [m]')
ylim([-1.5 1.5])
xlim([0 1.1*L])

for t=linspace(0,2/f_0(k),200)

    if ishandle(h1)
        w1=real(modeshapes_norm(:,k)*exp(1i*w_0(k)*t));
        h1.YData=w1;
        pause(.03)
    else
        return
    end

end 


%% POINT b) TRANSFER FUNCTIONS

%% Assigning forces and positions

% provide positions of inputs (forces) and outputs (amplitude of 
% oscillation) in order to compute those transfer functions

x_m=[1.2];        % vector of inputs positions [m]
x_i=[0.3 0.95 L]; % vector of outputs positions [m]

n_in=length(x_m);
n_out=length(x_i);

% plotting inputs and outputs positions

figure
hold on
plot(x,zeros(1,length(x)),'k','LineWidth',1.5)
xlim([0 1.1*L])
grid on
xlabel('x [m]')
title('Positions of inputs and outputs')
plot(x_m,zeros(1,n_in),'rx','LineWidth',2)
plot(x_i,zeros(1,n_out),'go','LineWidth',1)
for k=1:N
    plot(x,modeshapes(:,k),'k--')
end
legend('','inputs (forces)','outputs (measurements)','modeshapes','Location','best')
hold off

%% Modal parameters computation

% modal mass matrix

M_q=zeros(N);
for k=1:N
    phi_k=modeshapes(:,k);
    fun_k=m*phi_k.*phi_k;
    m_k=trapz(x,fun_k);

    M_q(k,k)=m_k;
end

% modal stiffness matrix (actually useless)

K_q=diag(diag(M_q)).*diag(w_0.^2);

% natural frequencies and damping coefficients are already defined

%% Matrix of transfer functions (rows referred to inputs, columns to outputs)

% compute it

G_mat=zeros(n_in,n_out,length(f));  % initialise 3D matrix (row --> input, column --> output)

for mm=1:n_in
    for ii=1:n_out
        for ff=1:length(f)
            sum=0;
            sum_k=0;
            for k=1:N
                gamma_k=gamma_0(k);
                w_k=w_0(k);
                B_hat_k=coeff(2,k);
                C_hat_k=coeff(3,k);
                D_hat_k=coeff(4,k);
                phi_k=@(x) cos(gamma_k*x)+B_hat_k*sin(gamma_k*x)+C_hat_k*cosh(gamma_k*x)+D_hat_k*sinh(gamma_k*x);
                phi_km=phi_k(x_m(mm));
                phi_ki=phi_k(x_i(ii));
                m_k=M_q(k,k);
                Omega=2*pi*f(ff);
    
                sum_k=(phi_ki*phi_km/m_k) ./ (-Omega.^2+1j*2*xi*w_k*Omega+w_k^2);
                sum=sum+sum_k;   % this is the specific transfer function between input mm and output ii evaluated for freq ff
            end
            G_mat(mm,ii,ff)=sum;
        end
    end
end

% plot amplitudes

figure
mi=1;

for mm=1:n_in
    for ii=1:n_out

        ampl=zeros(length(f),1);
        for ff=1:length(f)
            ampl(ff)=abs(G_mat(mm,ii,ff));
        end
        
        subplot(n_in,n_out,mi)
        semilogy(f,ampl)
        grid on
        xlabel('f [Hz]')
        if ii==1
            ylabel(['x_i_n=',num2str(x_m(mm)),' m'])
        end
        if mm==1
            title(['x_o_u_t=',num2str(x_i(ii)),' m'])
        end
        mi=mi+1;

    end
end

% plot phases

figure
mi=1;

for mm=1:n_in
    for ii=1:n_out

        phase=zeros(length(f),1);
        for ff=1:length(f)
            phase(ff)=angle(G_mat(mm,ii,ff));
        end
        
        subplot(n_in,n_out,mi)
        plot(f,phase)
        grid on
        xlabel('f [Hz]')
        yticks([-pi -pi/2 0 pi/2 pi])
        yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
        if ii==1
            ylabel(['x_i_n=',num2str(x_m(mm)),' m'])
        end
        if mm==1
            title(['x_o_u_t=',num2str(x_i(ii)),' m'])
        end
        mi=mi+1;

    end
end

% plot a specific tranfer function required by the user

mm=str2double(inputdlg(['Enter the number of the input force (1 to ',num2str(length(x_m)),'):']));
ii=str2double(inputdlg(['Enter the number of the output position(1 to ',num2str(length(x_i)),'):']));

figure

subplot(3,1,1)  % input and output position
hold on
plot(x,zeros(1,length(x)),'k','LineWidth',1.5)
plot(x_m(mm),0,'rx','LineWidth',2)
plot(x_i(ii),0,'go','LineWidth',1)
for k=1:N
    plot(x,modeshapes(:,k),'k--')
end
legend('','inputs (force)','output (measurement)','modeshapes','Location','best')
xlim([0 1.1*L])
grid on
xlabel('x [m]')
title(['Trans. func. between force in x_m= ',num2str(x_m(mm)),' m and the output in x_i= ',num2str(x_i(ii)),' m'])
hold off

subplot(3,1,2)  % amplitude
ampl=zeros(length(f),1);
    for ff=1:length(f)
        ampl(ff)=abs(G_mat(mm,ii,ff));
    end
semilogy(f,ampl)
xlabel('f [Hz]')
ylabel(['|G_',num2str(mm),'_,_',num2str(ii),'| [m/N]'])
grid on

subplot(3,1,3)  % phase
phase=zeros(length(f),1);
    for ff=1:length(f)
        phase(ff)=angle(G_mat(mm,ii,ff));
    end
plot(f,rad2deg(phase))
xlabel('f [Hz]')
ylabel(['∠G_',num2str(mm),'_,_',num2str(ii),' [deg]'])
grid on

%% POINT c) RECONSTRUCTION OF MODAL PARAMETERS
%
% We want to reconstruct modal parameters of our system by looking at 
% experimental data (a matrix of experimental transfer functions). In this
% case we'll take a "ficticious" experimental matrix (the one previously 
% computed, G_mat). 
% This matrix is a A x R matrix, where A is the length of the frequency
% vector, while R is the number of tranfer functions (each r-th column is
% associated to a tranfer function).

% define matrix G_exp (just pass from a 3D matrix to a 2D matrix)

G_exp=zeros(length(f),n_in*n_out);
rr=1;
for mm=1:n_in
    for ii=1:n_out
        G_exp(:,rr)=G_mat(mm,ii,:);
        rr=rr+1;
    end
end

%% Error minimisation procedure 

% Now we need to reconstruct modal parameters by looking at G_exp. So first
% we need to find the N peaks of G_exp.

[~,indx]=findpeaks(log(abs(G_exp(:,1))));   % find natural freq by looking at the first transf funct
N=length(indx);
w_0=w(indx);
f_0=f(indx);
gamma_0=gamma(indx);

opt_length=1000;    % "width" of our approximation (in terms of n° of indices)
f_opt_mat=zeros(opt_length+1,N);   % first output of the for cycle: frequency vectors (as columns) where we perform the optimization
X_mat=zeros(2+3*n_in*n_out,N);   % second output of the for cycle: matrix that contains (each column) the extimated X (containing modal parameters)
for kk=1:N  % we perform the optimization for each mode (resonance peak)

    % compute the frequency range around the k-th mode where we perform the
    % optimization
    
    indx_in=indx(kk)-ceil(opt_length/2);
    indx_end=indx(kk)+floor(opt_length/2);
    f_opt=f(indx_in:indx_end);

    f_opt_mat(:,kk)=f_opt';
    
    % Isolate the experimental transfer funct matrix around the k-th peak.
    % G_exp_k contains U rows and R columns (each column contains a transf
    % funct, defined over the "restricted" U frequencies)
    
    U=length(f_opt);
    R=min(size(G_exp));
    
    G_exp_k=G_exp(indx_in:indx_end,:);
    
    % define the initial guessess vector for the unknown parameters
    
    w_k_guess=w(indx(kk));   % k-th natural pulsation (take it close to the one we know from data)
    
    phase_der=( angle(G_exp(indx(kk),1))-angle(G_exp(indx(kk)-1,1)) ) / ( 2*pi*f(indx(kk))-2*pi*f(indx(kk)-1) ); % adimensional damping estimated with phase derivative method
    xi_k_guess=-1/(w_k_guess*phase_der);
    
    A_k_guess=zeros(1,R);      % vector of coefficients A_k (each element refers to the r-th trans funct):
                               % each element is choosen as minus the imaginary part of the experimental transfer
                               % transf funct evaluated at the k-th natural freq, multiplied by xi_k and w_k
    for rr=1:R
        [~,indx_resonance]=findpeaks(abs(G_exp_k(:,rr)),'NPeaks',1);
        A_k_guess(rr)=-imag(G_exp_k(indx_resonance,rr))*2*xi_k_guess*w_k_guess^2;
    end
    
    RL_k_guess=zeros(1,R)+eps*1j*ones(1,R); % left-residual coefficients
    RH_k_guess=zeros(1,R)+eps*1j*ones(1,R); % right-residual coefficients
    
    X_guess=[w_k_guess,xi_k_guess,A_k_guess,RL_k_guess,RH_k_guess]';   % initial guess parameters vector
    
    % compute the error function that only depends on parameters X to be found
    % by lsqnonlin (see parameterfun for more information)
    
    err=@(X) parameterfun(X,f_opt,G_exp_k);
    
    % minimise with matlab function
    
    options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
    X=lsqnonlin(err,X_guess,[],[],options);
    
    X_mat(:,kk)=X;
end
    
%% Checking the extimated parameters
    
% plot amplitudes

figure
mi=1;

for mm=1:n_in   
    for ii=1:n_out
        
        % plot experimental data

        subplot(n_in,n_out,mi)
        semilogy(f,abs(G_exp(:,mi)))
        hold on

        % plot numerical, reconstructed transfer functions
        
        for kk=1:N
            G_num=@(Omega) (X_mat(mi+2,kk)./(-Omega.^2+1j*2*X_mat(2,kk)*X_mat(1,kk)*Omega+X_mat(1,kk)^2)) + (X_mat(mi+R+2,kk)./(Omega.^2)) + (X_mat(mi+2*R+2,kk));
            ampl=abs(G_num(2*pi*f_opt_mat(:,kk)));
            semilogy(f_opt_mat(:,kk),ampl,'LineWidth',1,'Color','r','Marker','.')           
        end

        % plot grid, axis, etc

        grid on
        xlabel('f [Hz]')
        if ii==1
            ylabel(['x_i_n=',num2str(x_m(mm)),' m'])
        end
        if mm==1
            title(['x_o_u_t=',num2str(x_i(ii)),' m'])
        end

        hold off
        mi=mi+1;

    end
end

% plot phases

figure
mi=1;

for mm=1:n_in 
    for ii=1:n_out

        % plot experimental data

        subplot(n_in,n_out,mi)
        plot(f,angle(G_exp(:,mi)))
        hold on

        % plot numerical, reconstructed transfer functions
        
        for kk=1:N
            G_num=@(Omega) (X_mat(mi+2,kk)./(-Omega.^2+1j*2*X_mat(2,kk)*X_mat(1,kk)*Omega+X_mat(1,kk)^2)) + (X_mat(mi+R+2,kk)./(Omega.^2)) + (X_mat(mi+2*R+2,kk));
            phase=angle(G_num(2*pi*f_opt_mat(:,kk)));
            plot(f_opt_mat(:,kk),phase,'LineWidth',1,'Color','r','Marker','.')           
        end

        % plot grid, axis, etc

        grid on
        xlabel('f [Hz]')
        yticks([-pi -pi/2 0 pi/2 pi])
        yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
        if ii==1
            ylabel(['x_i_n=',num2str(x_m(mm)),' m'])
        end
        if mm==1
            title(['x_o_u_t=',num2str(x_i(ii)),' m'])
        end

        hold off
        mi=mi+1;

    end
end

% plot a specific tranfer function required by the user

mm=str2double(inputdlg(['Enter the number of the input force (1 to ',num2str(length(x_m)),'):']));
ii=str2double(inputdlg(['Enter the number of the output position(1 to ',num2str(length(x_i)),'):']));

mi=ii+(mm-1)*n_out;

figure

subplot(3,1,1)  % input output position
hold on
plot(x,zeros(1,length(x)),'k','LineWidth',1.5)
plot(x_m(mm),0,'rx','LineWidth',2)
plot(x_i(ii),0,'go','LineWidth',1)
for k=1:N
    plot(x,modeshapes(:,k),'k--')
end
legend('','inputs (force)','output (measurement)','modeshapes','Location','best')
xlim([0 1.1*L])
grid on
xlabel('x [m]')
title(['Trans. func. between force in x_m= ',num2str(x_m(mm)),' m and the output in x_i= ',num2str(x_i(ii)),' m'])
hold off

subplot(3,1,2)  % amplitude reconstruction
semilogy(f,abs(G_exp(:,mi)))
hold on
for kk=1:N
    G_num=@(Omega) (X_mat(mi+2,kk)./(-Omega.^2+1j*2*X_mat(2,kk)*X_mat(1,kk)*Omega+X_mat(1,kk)^2)) + (X_mat(mi+R+2,kk)./(Omega.^2)) + (X_mat(mi+2*R+2,kk));
    ampl=abs(G_num(2*pi*f_opt_mat(:,kk)));
    semilogy(f_opt_mat(:,kk),ampl,'LineWidth',1,'Color','r','Marker','.')           
end
grid on
xlabel('f [Hz]')
ylabel(['|G_',num2str(mm),'_,_',num2str(ii),'| [m/N]'])
legend('experimental data','numerical data')
hold off

subplot(3,1,3)  % phase reconstruction
plot(f,rad2deg(angle(G_exp(:,mi))))
hold on
for kk=1:N
    G_num=@(Omega) (X_mat(mi+2,kk)./(-Omega.^2+1j*2*X_mat(2,kk)*X_mat(1,kk)*Omega+X_mat(1,kk)^2)) + (X_mat(mi+R+2,kk)./(Omega.^2)) + (X_mat(mi+2*R+2,kk));
    phase=angle(G_num(2*pi*f_opt_mat(:,kk)));
    plot(f_opt_mat(:,kk),rad2deg(phase),'LineWidth',1,'Color','r','Marker','.')           
end
grid on
xlabel('f [Hz]')
ylabel(['∠G_',num2str(mm),'_,_',num2str(ii),' [deg]']) 
legend('experimental data','numerical data')
hold off

%% Modeshapes reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MUST IMPOSE ONE OUTPUT POINT AT THE TIP IN ORDER FOR THIS SHIT TO WORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if floor(N/2)==N/2
    
    % compute points approximating the modeshape

    reconstructed_modeshapes=zeros(n_out,N);  % pre-allocate matrix containing the 
                                              % points that approximate the modeshapes 
                                              % (N modeshapes as columns, each approximated with n_out points as rows)
    for k=1:N
        points_kk=ones(1,n_out);
        points_kk_norm=zeros(1,n_out);
        for ii=1:n_out
            points_kk(ii)=-imag(G_exp(indx(k),ii));
            points_kk_norm=points_kk/max(abs(points_kk));
            reconstructed_modeshapes(:,k)=points_kk_norm;
        end
    end

    % plot them

    figure
    for k=1:N
        subplot(N/2,N/2,k)
        plot(x,modeshapes_norm(:,k),'r')
        hold on
        plot(x,-modeshapes_norm(:,k),'r--')
        plot(x,zeros(1,length(x)),'k','LineWidth',2)
        plot(x_i,reconstructed_modeshapes(:,k),'b*')
        axis([0 1.1*x(end) -1.5 1.5])
        grid on
        xlabel('x [m]')
        title(['Modeshape for \omega_',num2str(k),'=',num2str(w_0(k),4),' Hz'])
        legend('model','','','reconstructed')
        hold off
    end

else
    
    % compute points approximating the modeshape

    reconstructed_modeshapes=zeros(n_out,N);  % pre-allocate matrix containing the 
                                              % points that approximate the modeshapes 
                                              % (N modeshapes as columns, each approximated with n_out points as rows)
    for k=1:N
        points_kk=imag(G_exp(indx(k),1:n_out));
        points_kk_norm=points_kk/max(abs(points_kk));
        reconstructed_modeshapes(:,k)=points_kk_norm;
    end

    % plot them

    figure
    for k=1:N
        subplot(ceil(N/2),floor(N/2),k)
        plot(x,modeshapes_norm(:,k),'r')
        hold on
        plot(x,-modeshapes_norm(:,k),'r--')
        plot(x,zeros(1,length(x)),'k','LineWidth',2)
        plot(x_i,reconstructed_modeshapes(:,k),'b*')
        axis([0 1.1*x(end) -1.5 1.5])
        grid on
        xlabel('x [m]')
        title(['Modeshape for \omega_',num2str(k),'=',num2str(w_0(k),4),' Hz'])
        legend('model','','','reconstructed')
        hold off
    end

end