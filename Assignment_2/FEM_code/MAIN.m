clear
close all
clc

%% Data of the problem

% structural damping
alpha=0.1;
beta=2.0e-4;

%% 1) BUILD FE MODEL
%% Load input file and assemble structure

[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;
 
%% Draw structure

dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

%% Assemble mass, damping and stiffness matrices

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);
C=alpha.*M+beta.*K;

%% Partitioning matrices M, C and K

M_FF=M(1:ndof,1:ndof);
M_FC=M(1:ndof,ndof+1:end);
M_CF=M(ndof+1:end,1:ndof);
M_CC=M(ndof+1:end,ndof+1:end);

C_FF=C(1:ndof,1:ndof);
C_FC=C(1:ndof,ndof+1:end);
C_CF=C(ndof+1:end,1:ndof);
C_CC=C(ndof+1:end,ndof+1:end);

K_FF=K(1:ndof,1:ndof);
K_FC=K(1:ndof,ndof+1:end);
K_CF=K(ndof+1:end,1:ndof);
K_CC=K(ndof+1:end,ndof+1:end);

%% 2) NATURAL FREQUENCIES AND MODESHAPES
%% Find natural frequencies and modeshapes

[modes,omega2]=eig(M_FF\K_FF);
omega=sqrt(diag(omega2));

% sort frequencies in ascending order 

[omega,i_omega]=sort(omega);
freq0=omega/2/pi;

% sort mode shapes in ascending order 

modes=modes(:,i_omega);

%% Find natural frequencies and modeshapes (also considering damping)

% % state matrix formulation
% 
% A=[-M_FF\C_FF -M_FF\K_FF; eye(size(M_FF)) zeros(size(M_FF))];
% [modes,omega]=eig(A);
% omega=diag(omega);
% modes=modes(floor(end/2)+1:end,:);
% 
% [~,indices]=sort(imag(omega));
% omega=omega(indices);
% omega=omega(floor(end/2)+1:end);
% omega=abs(imag(omega));
% modes=modes(:,indices);
% modes=abs(modes(:,floor(end/2)+1:end));
% 
% % sort natural freq and modeshapes in ascending order
% 
% [~,indices]=sort(omega);
% omega=omega(indices);
% freq0=omega/2/pi;
% modes=modes(:,indices);

%% Displaying natural frequencies and modeshapes

scale_factor=5; % scale factor to display the modeshapes
nmodes=4;   % number of modes we want to consider

% display natural frequencies

fprintf(['\nFirst ',num2str(nmodes)',' natural frequencies are [Hz]: ',num2str(freq0(1:nmodes)'),'\n\n'])

% plot first 3 natural frequencies

for k=1:nmodes
    figure
    diseg2(modes(:,k),scale_factor,incid,l,gamma,posiz,idb,xy)
    title(['Modeshape ',num2str(k),' (',num2str(freq0(k)),' Hz)'])
end

%% 3)
%% Define force vector

F_F0=zeros(ndof,1);
F_F0(idb(9,2))=1;

%% Computing FRFs

% define frequency (pulsation) vector

df=0.01;  % freq resolution
freq_max=8; % max frequency of the input
freq=0:df:freq_max;
om=freq*2*pi;
freq0_of_interest=freq0(freq0<=freq_max);    % natural frequencies inside the frequency range of interest

% computing nodal displacements (actually here they correspond to the FRFs,
% since we set unitary force) and the reactions

x_F0=zeros(ndof,length(om));
for ii=1:length(om) 
    A=-om(ii)^2*M_FF+1j*om(ii)*C_FF+K_FF; 
    x_F0(:,ii)=A\F_F0; 
end

% extract the desired FRFs

input2A=x_F0(idb(9,2),:);
input2B=x_F0(idb(23,1),:);

% plot FRFs

figure
subplot(2,1,1)
semilogy(freq,abs(input2A))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('amplitude')
title('FRF between vertical force in A and vertical displacement in A')
subplot(2,1,2)
plot(freq,unwrap(angle(input2A)))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('phase')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

figure
subplot(2,1,1)
semilogy(freq,abs(input2B))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('amplitude')
title('FRF between vertical force in A and horizontal displacement in B')
subplot(2,1,2)
plot(freq,unwrap(angle(input2B)))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('phase')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

%% 4)
%% Computing FRFs (modal approach: first 2 modes)

considered_modes=2;

% modal matrices

Phi=modes(:,1:considered_modes);
Mq=Phi'*M_FF*Phi; 
Kq=Phi'*K_FF*Phi; 
Cq=Phi'*C_FF*Phi; 

% modal forcing vector

Q=Phi'*F_F0; 

% FRFs in modal superposition approach 

q=zeros(considered_modes,length(om));
for ii=1:length(om) 
    q(:,ii)=(-om(ii)^2*Mq+1j*om(ii)*Cq+Kq)\Q; 
end 
x_mod=Phi*q;

% extract desired FRF

FRF_mod=x_mod(idb(23,1),:);

% plot comparison

figure
subplot(2,1,1)
semilogy(freq,abs(input2B))
hold on
semilogy(freq,abs(FRF_mod))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('amplitude')
legend('FEM','modal superposition')
title('FRF between vertical force in A and horizontal displacement in B')
hold off
subplot(2,1,2)
plot(freq,unwrap(angle(input2B)))
hold on
plot(freq,unwrap(angle(FRF_mod)))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('phase')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legend('FEM','modal superposition')
hold off

%% 5)
%% Define force vector

F_F0=zeros(ndof,1);
F_F0(idb(1,2))=1;

%% Compute nodal displacements

% computing nodal displacements and the reactions

x_F0=zeros(ndof,length(om));
R=zeros(3*nnod-ndof,length(om));
for ii=1:length(om) 
    A=-om(ii)^2*M_FF+1j*om(ii)*C_FF+K_FF; 
    x_F0(:,ii)=A\F_F0; 
    R(:,ii)=(-om(ii)^2*M_CF+1j*om(ii)*C_CF+K_CF)*x_F0(:,ii);
end

%% Compute axial force on element 30

% compute matrix lambda, to go from global to local ref frame

gamma_el=gamma(30);
lambda=[cos(gamma_el) sin(gamma_el) 0; -sin(gamma_el) cos(gamma_el) 0; 0 0 1];

% find nodal displacements of nodes 23 and 25 and convert into local

x_23_G=x_F0(idb(23,:),:);
x_25_G=zeros(3,length(om));
x_25_G(3,:)=x_F0(idb(25,3),:);

x_23_L=lambda*x_23_G;
x_25_L=lambda*x_25_G;

% compute constant axial force

L_el=l(30);
N=EA(30)*(x_25_L(1,:)-x_23_L(1,:))/L_el;

%% Plot FRF

% since we have N for an input force, N is directly the FRF

figure
subplot(2,1,1)
semilogy(freq,abs(N))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('amplitude')
title('FRF between vertical force in C and axial reaction in O_2')
subplot(2,1,2)
plot(freq,unwrap(angle(N)))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('phase')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

%% Comparison with the vertical reaction force in O2

reaction=R(idb(25,2)-ndof,:);

figure
subplot(2,1,1)
semilogy(freq,abs(N))
hold on
semilogy(freq,abs(reaction))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('amplitude')
title('FRF between vertical force in C and axial reaction/vertical reaction in O_2')
legend('axial force','vertical reaction')
hold off
subplot(2,1,2)
plot(freq,unwrap(angle(N)))
hold on
plot(freq,unwrap(angle(reaction)))
xline(freq0_of_interest,'k--')
grid on
axis tight
xlabel('f')
ylabel('phase')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legend('axial force','vertical reaction')
hold off

%% 6) STATIC RESPONSE DUE TO WEIGHT
%% Build the equivalent forcing vector

% APPROACH 1: AS A GENERIC DISTRIBUTED LOAD

    % define the distributed load (in the global reference)

    p0_Ga=[0 -312*9.81]';   % m*g [N/m] for elements from 1 to 8
    p0_Gb=[0 -90*9.81]';   % m*g [N/m] for elements from 9 to 22
    p0_Gc=[0 -200*9.81]';   % m*g [N/m] for elements from 23 to 30

    % compute the equivalent force

    F_0=zeros(3*nnod,1);   % initialise the final equivalent force vector (in the global reference, for all elements, over all nodes (both constrained and free)
    for ii=1:nbeam % cicle over elements involved (all of them)
        gamma_ii=gamma(ii);
        L_k=l(ii);
        if ii<=8
            p0_L=[cos(gamma_ii) sin(gamma_ii);-sin(gamma_ii) cos(gamma_ii)]*p0_Ga;   % distr load in the local reference 
            Fk_L=[L_k/2 0 0 L_k/2 0 0]'*p0_L(1) + [0 L_k/2 L_k^2/12 0 L_k/2 -L_k^2/12]'*p0_L(2);    % equivalent force in the local reference
            Fk_G=[cos(gamma_ii) -sin(gamma_ii) 0 0 0 0;
                  sin(gamma_ii) cos(gamma_ii) 0 0 0 0;
                  0 0 1 0 0 0;
                  0 0 0 cos(gamma_ii) -sin(gamma_ii) 0;
                  0 0 0 sin(gamma_ii) cos(gamma_ii) 0;
                  0 0 0 0 0 1]*Fk_L;    % equivalent force (for k-th element) in the global reference 
            E_k=zeros(6,3*nnod);    % extraction matrix, to place the equivalent force for the k-th element correctly inside the "final" equiv force vector
            for jj=1:6
                E_k(jj,incid(ii,jj))=1;
            end
            F_0=F_0+E_k'*Fk_G;    % add the k-th element contribution to the final equivalent force vector
        end
        if ii>=9 && ii<=22
            p0_L=[cos(gamma_ii) sin(gamma_ii);-sin(gamma_ii) cos(gamma_ii)]*p0_Gb;   % distr load in the local reference 
            Fk_L=[L_k/2 0 0 L_k/2 0 0]'*p0_L(1) + [0 L_k/2 L_k^2/12 0 L_k/2 -L_k^2/12]'*p0_L(2);    % equivalent force in the local reference
            Fk_G=[cos(gamma_ii) -sin(gamma_ii) 0 0 0 0;
                  sin(gamma_ii) cos(gamma_ii) 0 0 0 0;
                  0 0 1 0 0 0;
                  0 0 0 cos(gamma_ii) -sin(gamma_ii) 0;
                  0 0 0 sin(gamma_ii) cos(gamma_ii) 0;
                  0 0 0 0 0 1]*Fk_L;    % equivalent force (for k-th element) in the global reference 
            E_k=zeros(6,3*nnod);    % extraction matrix, to place the equivalent force for the k-th element correctly inside the "final" equiv force vector
            for jj=1:6
                E_k(jj,incid(ii,jj))=1;
            end
            F_0=F_0+E_k'*Fk_G;    % add the k-th element contribution to the final equivalent force vector
        end
        if ii>=23
            p0_L=[cos(gamma_ii) sin(gamma_ii);-sin(gamma_ii) cos(gamma_ii)]*p0_Gc;   % distr load in the local reference 
            Fk_L=[L_k/2 0 0 L_k/2 0 0]'*p0_L(1) + [0 L_k/2 L_k^2/12 0 L_k/2 -L_k^2/12]'*p0_L(2);    % equivalent force in the local reference
            Fk_G=[cos(gamma_ii) -sin(gamma_ii) 0 0 0 0;
                  sin(gamma_ii) cos(gamma_ii) 0 0 0 0;
                  0 0 1 0 0 0;
                  0 0 0 cos(gamma_ii) -sin(gamma_ii) 0;
                  0 0 0 sin(gamma_ii) cos(gamma_ii) 0;
                  0 0 0 0 0 1]*Fk_L;    % equivalent force (for k-th element) in the global reference 
            E_k=zeros(6,3*nnod);    % extraction matrix, to place the equivalent force for the k-th element correctly inside the "final" equiv force vector
            for jj=1:6
                E_k(jj,incid(ii,jj))=1;
            end
            F_0=F_0+E_k'*Fk_G;    % add the k-th element contribution to the final equivalent force vector
        end
    end
    F_F0=F_0(1:ndof);   % extract F_F0 out of F_0=[F_F0; F_C]

    % display the deformed structure

    x_F0=K_FF\F_F0; % vector of nodal displacements
    figure
    diseg2(x_F0,30,incid,l,gamma,posiz,idb,xy)
    title('Deformed structure under its weight - "long" method - x30 scale factor')

% APPROACH 2: USING MASS MATRIX

    % define the g acceleration vector, acting on vertical nodal
    % displacements

    g_vect=zeros(3*nnod,1);
    g_vect(idb(:,2))=-9.81;

    % define the equivalent forcing vector

    F_gravity=M*g_vect;
    F_gravity=F_gravity(1:ndof);

    % compute and plot the deformed structure

    x_gravity=K_FF\F_gravity;
    figure
    diseg2(x_gravity,30,incid,l,gamma,posiz,idb,xy)
    title('Deformed structure under its weight - "short" method - x30 scale factor')

%% Compute the maximum nodal static deflection and display its value

nodal_deflections=sqrt(x_gravity(idb(1:nnod-2,1)).^2+x_gravity(idb(1:nnod-2,2)).^2);
[max_disp,indx]=max(nodal_deflections);
fprintf('Max nodal deflection occurs at node %d and its value is %.1f mm.\n\n',indx,max_disp*1000)

%% Compute maximum displacement in general (also using shapefunctions)
    
% compute vectors with ALL x and y positions of deformed and undeformed structure

    % initialise vectors of deformed and undeformed coordinates
    x_undef=0;
    y_undef=0;
    x_deformed=0;
    y_deformed=0;
    % iterate all elements to compute stuff
    for k=1:nbeam
    % building the nodal displacements vector of each element in the global
    % reference frame
        xkG=zeros(6,1);
        for iri=1:6
            if incid(k,iri) <= ndof
                xkG(iri,1)=x_gravity(incid(k,iri));
            else
                xkG(iri,1)=0.;
            end
        end
        % scale factor (1 in this case)
        xkG=1*xkG;
        % global to Local reference frame rotation
        lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
                  -sin(gamma(k)) cos(gamma(k)) 0.
	                    0.      0.     1. ] ;
        Lambda = [ lambda     zeros(3)
                  zeros(3)   lambda      ] ;
        xkL=Lambda*xkG;
        % computing the axial (u) and transversal (w) displacements by means of
        % shape functions
        csi=l(k)*[0:0.05:1];
        fu=zeros(6,length(csi));
        fu(1,:)=1-csi/l(k);
        fu(4,:)=csi/l(k);
        u=(fu'*xkL)';
        fw=zeros(6,length(csi));
        fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
        fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
        fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
        fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
        w=(fw'*xkL)';
        % local to global transformation of the element's deformation
        xyG=lambda(1:2,1:2)'*[u+csi;w];
        undef=lambda(1:2,1:2)'*[csi;zeros(1,length(csi))];
        % filling the coordinates vectors
        x_undef=[x_undef, undef(1,:)+posiz(k,1)];
        y_undef=[y_undef, undef(2,:)+posiz(k,2)];
        x_deformed=[x_deformed, xyG(1,:)+posiz(k,1)];
        y_deformed=[y_deformed, xyG(2,:)+posiz(k,2)];
    end
    % % plot check
    % figure
    % plot(x_undef,y_undef,'r--')
    % hold on
    % plot(x_deformed,y_deformed,'b')
    % grid on
    % hold off

% compute the displacement for every point

delta_x=x_deformed-x_undef;
delta_y=y_deformed-y_undef;
deflections=sqrt(delta_x.^2+delta_y.^2);
[max_disp,indx]=max(deflections);
fprintf('Max deflection is %.1f mm.\n\n',max_disp*1000)

%% 7) MOVING LOAD
%% Speed and position profiles of the load

% define time vector (from 0 to 20 s)

tspan=linspace(0,20,1001); % [s]

% define position vector, as a function of time (i.e. this vector contains
% the position of the load, starting from point D as the origin, at each 
% time istant)

position=zeros(1,length(tspan));    % position of the load in time
speed=zeros(1,length(tspan)); % speed of the load in time
for ii=1:length(tspan)
    t=tspan(ii);
    if t<=8
        position(ii)=1/8*t^2;
        speed(ii)=1/4*t;
    end
    if t>8 && t<=12
        position(ii)=2*t-8;
        speed(ii)=2;
    end
    if t>12
        position(ii)=5*t-1/8*t^2-26;
        speed(ii)=5-1/4*t;
    end
end

% plot position and speed of the load to visual check

figure

yyaxis left
plot(tspan,position,'b')
ylabel('position [m]')

yyaxis right
plot(tspan,speed,'r')
ylabel('speed [m/s]')

grid on
title('Position and speed profiles of the moving load')
xlabel('t [s]')
axis tight
hold off

%% Build modal matrices

% we are accounting for the first 5 modes

considered_modes=5;
Phi=modes(:,1:considered_modes);    % Phi is a matrix containing the first 5 modes as columns (evaluated on the dofs)

Mq=Phi'*M_FF*Phi;   % these all are 5x5 matrices
Kq=Phi'*K_FF*Phi; 
Cq=Phi'*C_FF*Phi; 

%% Modal forcing vector

% define the load (here we chose 300 000 N, that corresponds more or less to a mass of 30 tons, the mass of a typical loaded naval container)

P=-300000;  % [N]

% filling matrix Phi_time, that now contains as ROWS the modeshapes evaluated at
% "all" points between D and A, as a function of time (i.e. let's take
% first row/modeshape: first element is the vertical displacement of the
% first modeshape at time t=0, second element is the vertical displacement
% of the second modeshape at time t=0+dt, etc)

Phi_time=zeros(considered_modes,length(tspan));

for ii_mode=1:considered_modes % we want to find the first 5 modeshapes in element 7, knowing the nodal diplacements and using the shapefunctions
    mode=modes(:,ii_mode);
    k=7;    % element 7 first
    % building the nodal displacements vector of element 7 in the global reference frame
    for iri=1:6
        if incid(k,iri) <= ndof
            xkG(iri,1)=mode(incid(k,iri));
        else
            xkG(iri,1)=0.;
        end
    end
    % scale factor
    xkG=1*xkG;
    % global to local reference frame rotation (not needed in this case)
    xkL=xkG;
    % computing the transversal (w) displacements by means of shape functions
    w=@(csi) (2*(csi/l(k)).^3-3*(csi/l(k)).^2+1)*xkL(2) + (l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k)))*xkL(3) + (-2*(csi/l(k)).^3+3*(csi/l(k)).^2)*xkL(5) + (l(k)*((csi/l(k)).^3-(csi/l(k)).^2))*xkL(6);
    % filling properly matrix Phi
    indicespos=find(position<12);
    indxpos=indicespos(end);
    Phi_time(ii_mode,1:indxpos)=w(position(1:indxpos));
end

for ii_mode=1:considered_modes % we want to find the first 5 modeshapes in element 8, knowing the nodal diplacements and using the shapefunctions
    mode=modes(:,ii_mode);
    k=8;    % element 8
    % building the nodal displacements vector of element 8 in the global reference frame
    for iri=1:6
        if incid(k,iri) <= ndof
            xkG(iri,1)=mode(incid(k,iri));
        else
            xkG(iri,1)=0.;
        end
    end
    % scale factor
    xkG=1*xkG;
    % global to local reference frame rotation (not needed in this case)
    xkL=xkG;
    % computing the transversal (w) displacements by means of shape functions
    w=@(csi) (2*(csi/l(k)).^3-3*(csi/l(k)).^2+1)*xkL(2) + (l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k)))*xkL(3) + (-2*(csi/l(k)).^3+3*(csi/l(k)).^2)*xkL(5) + (l(k)*((csi/l(k)).^3-(csi/l(k)).^2))*xkL(6);
    % filling properly matrix Phi
    indicespos=find(position>=12);
    indxpos=indicespos(1);
    Phi_time(ii_mode,indxpos:end)=w(-12+position(indxpos:end));
end

% visual check of matrix Phi_time

figure
plot(tspan,Phi_time)
xlabel('t')
ylabel('vertical displacement')
legend('first modesh.','second modesh.','third modesh.','fourth modesh.','fifth modesh.')
title('Temporal vertical component of the modeshapes')

% build the lagrangian component Q(t)

Q=P*Phi_time;

%% Solving equation of motion

% initial condition (static deflection in A in modal coordinates)

q0=Kq\Q(:,1);

% solving the ODE in modal coordinates

q_vect=zeros(considered_modes,length(tspan));   % matrix q of modal coordinates in time (each row is the solution qk(t))

for k=1:considered_modes

    % extract modal parameters
    mk=Mq(k,k); 
    kk=Kq(k,k);
    ck=Cq(k,k);
    Qk=Q(k,:);

    % build state matrix for the state-space problem
    A=[0 1; -kk/mk -ck/mk];

    % build the Qk vector function of time (for each time of tspan we need
    % Qk to change its value)
    struct=spline(tspan,Qk);
    Qk_time=@(t) ppval(struct,t);

    % find the proper initial condition (static deflection, in modal coordinates, due to the load P standing in D)
    x0=[q0(k) 0]';
    %x0=[0 0]';

    % define the problem (state space formulation) and solve it with ode45
    x_dot=@(t,x) A*x+[0 Qk_time(t)/mk]';
    [time,x]=ode45(x_dot,tspan,x0);
    
    % fill solution matrix q_vect
    q_vect(k,:)=x(:,1);
end

% convert from modal coordinates to physical coordinate (vertical displacement in A)

w_A=zeros(1,length(time));
for ii=1:length(tspan)
    for k=1:considered_modes
        Phik_A=Phi_time(k,end);
        qk=q_vect(k,ii);
        w_A(ii)=w_A(ii)+Phik_A*qk;
    end
end

% plot the vertical displacement of A in time

figure
plot(time,w_A)
grid on
xlabel('t [s]')
ylabel('vertical displacement [m]')
title('Vertical displacement of point A - moving load')

% (TO CHECK THE INITIAL CONDITION) Static deflection in D due to load P in D

F_F0=zeros(ndof,1);
F_F0(idb(7,2))=P;

x_F0=K_FF\F_F0;
defl=x_F0(idb(9,2));
fprintf(['Static deflection in A due to load P in D (FEM): ',num2str(defl),' m\n'])
fprintf(['Static deflection in A due to load P in D (modal approach, first 5 modes, at t=0): ',num2str(w_A(1)),' m\n']);