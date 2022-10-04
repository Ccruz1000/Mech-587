% MECH 479 - CFD 
% One-dimensional advection-diffusion by the FTCS scheme
n=81;             % number of mesh points
nstep=100;        % number of time step
length=2.0;       % domain length
h=length/(n-1);   % grid size
dt=0.005;         % time step size
c = 1.0;          % convection speed
D=0.05;           % diffusivity
u=zeros(n,1);     % unknown function
y=zeros(n,1);     % temporary array
ex=zeros(n,1);    % exact solution array
x =zeros(n,1);    % grid location
time=0.0;         % initial time
% Initial conditions
for i=1:n, 
    u(i)=0.5*sin(2*pi*h*(i-1)); 
    x(i)= (i-1)*h;
end; 
% Loop over time
for m=1:nstep, m, time
for i=1:n, 
    ex(i)=exp(-4*pi*pi*D*time)*...
0.5*sin(2*pi*(h*(i-1)-time)); 
end; % exact solution 
hold off; 
plot(x,u,'linewidt',2); axis([0 length -0.5, 0.5]); 
% plot solution
hold on; plot(x,ex,'r','linewidt',2);
legend('numerical', 'exact');
xlabel('location, x');
ylabel('solution, u');
%pause; % plot exact solution
y=u; % store the solution
for i=2:n-1,
u(i)=y(i)-0.5*(dt/h)*c*(y(i+1)-y(i-1))+...
      D*(dt/h^2)*(y(i+1)-2*y(i)+y(i-1)); % advect by centered differences
end;
u(n)=y(n)-0.5*(dt/h)*c*(y(2)-y(n-1))+...
     D*(dt/h^2)*(y(2)-2*y(n)+y(n-1)); % do endpoints for
u(1)=u(n); % periodic boundaries
time=time+dt;
end