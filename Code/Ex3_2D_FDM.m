% MECH 479 - CFD
% EX3: Two-dimensional unsteady diffusion by the FTCS scheme
%------------------------------------------------------------
n=32;
m=32;
nstep=120;
D=0.025;
length=2.0;
h=length/(n-1);
dt=5.0*0.125*h*h/D;
u=zeros(n,m);
uo=zeros(n,m);
time=0.0;
% initial conditions
U=-0.0; V=-1.0; u(12:21,n)=1.0;
for l=1:nstep,l,time;
hold off;mesh(u); axis([0 n 0 m 0 1.5]);pause;
uo=u;
% iterative solver (SOR)
for i=2:n-1, for j=2:m-1
u(i,j)=uo(i,j)-(0.5*dt*U/h)*(uo(i+1,j)-uo(i-1,j))-...
(0.5*dt*V/h)*(uo(i,j+1)-uo(i,j-1))+...
(D*dt/h^2)*(uo(i+1,j)+ uo(i,j+1)+uo(i-1,j)+uo(i,j -1)-4*uo(i,j));
end,end
% boundary conditions
for i=1:n, 
    u(i,1)=u(i,2);
end
for j=1:m,
    u(1,j)=u(2,j);
    u(m,j)=u(m-1,j);
end
time=time+dt;
end