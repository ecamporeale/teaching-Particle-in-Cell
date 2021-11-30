clear all
close all

L=2*pi; %  box length
DT=.5; % time step
NT=200; % number of cycles
NTOUT=25; % 
NG=32; % grid points
N=1000; % number of particles
WP=1;
QM=-1;


V0=0.0; % beam velocity
VT=.0; % thermal velocity
XP1=0; % perturbation in position
V1=0; % perturbation in velocity
mode=1; % perturbed mode
Q=WP^2/(QM*N/L);
rho_back=-Q*N/L;
dx=L/NG; % grid size

%diagnostics
Eg_fft=[];temperature=[];W_k=[];W_E=[];P=[];
position=[];

% initial loading
xp=linspace(0,L-L/N,N)';
vp=VT*(2*rand(N,1)-1.0); % uniform distribution
%vp=VT*randn(N,1); % Maxwellian distribution


% Perturbation
%vp=vp+V1*sin(2*pi*xp/L*mode);
%xp=xp+XP1*sin(2*pi*xp/L*mode);
% apply bc on the particle positions (periodic)
out=(xp<0); xp(out)=xp(out)+L;
out=(xp>=L);xp(out)=xp(out)-L;

% matrix for Poisson equation
p=1:N;p=[p p];
Poisson=spdiags([ones(NG-1,1) -2*ones(NG-1,1) ones(NG-1,1)],[-1 0 1],NG-1,NG-1);

% Main computational cycle
for it=1:NT

  % update xp
  xp=xp+vp*DT;

  % apply bc on the particle positions (periodic)
  out=(xp<0); xp(out)=xp(out)+L;
  out=(xp>=L);xp(out)=xp(out)-L;
  
  % projection p->g
  g1=floor(xp/dx-.5)+1;g=[g1;g1+1];
  fraz1=1-abs(xp/dx-g1+.5);fraz=[fraz1;1-fraz1];
  
  % apply bc on the projection
  out=(g<1);g(out)=g(out)+NG;
  out=(g>NG);g(out)=g(out)-NG;
  mat=sparse(p,g,fraz,N,NG);
  rho=full((Q/dx)*sum(mat))'+rho_back;
  % computing fields
  Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
  Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
  
  % projection q->p and update of vp
  vp=vp+mat*QM*Eg*DT;
  
end


