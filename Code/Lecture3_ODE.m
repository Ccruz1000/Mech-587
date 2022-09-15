% Mech 587 - CFD
% A sample code for three integration methods

nstep=5;
dt=0.5;
u1=zeros(nstep, 1);
u2=zeros(nstep, 1);
u3=zeros(nstep, 1);
uex=zeros(nstep, 1);
t=zeros(nstep, 1);
t(1) = 0; u1(1)=1; u2(1)=1; u3(1)=1; uex(1)=1;

for i = 2:nstep
  u1(i)=u1(i-1) - dt*u1(i-1); % Forward Euler
  u2(i)=u2(i-1) / (1.0 + dt); % Backward Euler
  u3(i) =u3(i-1) * (1.0-0.5*dt)/(1.0+0.5*dt); %Trapezoidal Rule
  t(i)=t(i-1)+dt;
  uex(i)=exp(-t(i));
endfor
plot(t,u1); hold on; plot(t,u2,'r'); plot(t,u3,'k'); 
plot(t,uex,'b', 'linewidt',3);
set(gca, 'fontsize', 24, 'lindewidt', 2);