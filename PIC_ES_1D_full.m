clear all
close all

L=2*pi; %  box length
DT=.5; % time step
NT=500; % number of cycles
NTOUT=25; % 
NG=32; % grid points
N=100000; % number of particles
WP=1;
QM=-1;


V0=0.2; % beam velocity
VT=.05; % thermal velocity
XP1=1e-3; % perturbation in position
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
%vp=VT*(2*rand(N,1)-1.0); % uniform distribution
 vp=VT*randn(N,1); % Maxwellian distribution

pm=[1:N]';pm=1-2*mod(pm,2); vp=vp+pm.*V0; % two stream -instability

% Perturbation
%vp=vp+V1*sin(2*pi*xp/L*mode);
xp=xp+XP1*sin(2*pi*xp/L*mode);
% apply bc on the particle positions (periodic)
out=(xp<0); xp(out)=xp(out)+L;
out=(xp>=L);xp(out)=xp(out)-L;

% matrix for Poisson equation
%p=1:N;p=[p p];
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

%  mat=sparse(p,g,fraz,N,NG);
%  rho=full((Q/dx)*sum(mat))'+rho_back;

  rho=zeros(NG,1);

  for k=1:length(g) % accumulate p->g
    rho(g(k)) = rho(g(k)) + fraz(k)*(Q/dx);
  end
  rho = rho+rho_back;
  
  % computing fields
  Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
  Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
  
  % interpolate g->p
  Ep = zeros(N,1);
  Etemp = [Eg(end); Eg; Eg(1)];
  
  for k=1:length(Ep)
   Ep(k) = fraz1(k)*Etemp(g1(k)+1) + (1-fraz1(k))*Etemp(g1(k)+2);
  end

 
  % update of vp
%  vp=vp+mat*QM*Eg*DT;
   vp = vp + QM*Ep*DT;

  
  set(gcf,'position',[155 1 470 660])
  subplot(3,1,1)
  plot(xp(1:2:end),vp(1:2:end),'.'),hold on
  plot(xp(2:2:end),vp(2:2:end),'.r'),hold off
  axis([0 L -0.5 0.5]); 
  set(gca,'fontsize',16),xlabel('xp'),ylabel('vp')
  title(['T = ' num2str(it)])
 
  subplot(3,1,2)
  
  nXBins = 100; 
  nYBins = 100; 

  vYEdge=linspace(0,L,nYBins);
  vXEdge = linspace(-0.5,0.5,nXBins);
  vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
  vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
  mHist2d = hist2d([xp vp],vYEdge,vXEdge); mHist2d = mHist2d/max(max(mHist2d));
  pcolor(vYLabel,vXLabel, log10(mHist2d')); caxis([-2 0]), shading interp; 
  axis([0 L -0.5 0.5])
  set(gca,'fontsize',16),xlabel('xp'),ylabel('vp')
  
  subplot(3,1,3);
  [m bin]=hist(vp,50);
  bar(bin,m/max(m),1)
  axis([-0.5 0.5 0 1])
  set(gca,'fontsize',16),xlabel('vp'),ylabel('f (v)')
  drawnow
    
  % diagnostics
  E_fft=fft(Eg)/NG;
  Eg_fft = [Eg_fft;abs(E_fft(2))];
  temperature=[temperature; std(vp)];
  W_k = [W_k;sum(vp.^2)];
  W_E = [W_E; sum(Eg.^2)];
  P= [P; sum(vp)];
%  M(it) = getframe(gcf); % only for movie!
end


