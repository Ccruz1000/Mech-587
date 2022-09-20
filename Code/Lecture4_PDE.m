
% PDE demo - advection, diffusion, dissipation
% MECH479 -- CFD
clear; 
clc;

% Discretization
m = 300;
h = 1 / m;
x = h * (1:m)';

% Finite Difference Stencils
A1 = spdiags(ones(m,1) * [1,-4,3] / (2*h), -2:0, m, m);
A1(1,m) = -4/2/h;
A1(1,m-1) = 1/2/h;
A1(2,m) = 1/2/h;

A2 = spdiags(ones(m,1) * [1,-2,1] / h^2, -1:1, m, m);
A2(1,m) = 1/h^2;
A2(m,1) = 1/h^2;

A3 = spdiags(ones(m,1) * [-0.5,1,-1,.5] / h^3, [-2,-1,1,2], m, m);
A3(1,m-1) = -.5/h^3;
A3(1,m) = 1/h^3;
A3(2,m) = -.5/h^3;
A3(m-1:m,1:2) = -A3(1:2,m-1:m)';

I = speye(m, m);
% Solving simple 1D advection/Hyperbolic equation
disp(' Advection, u_t = -u_x ')
u = exp(-(x-.5).^2 / .1^2);
A = -A1;
T = 0.5;
k = .5*h;
for n=1:ceil(T/k)
    
    u = (I - k/2*A) \ ((I + k/2*A)*u);
    if mod(n,10) == 0, plot(x,u, '-o'), axis([0,1,-.1,1.1]), grid on, drawnow, ...
    end
    if mod(n,200) ==0, pause 
        print 'Running step=%d', n
    end 
end
pause
% Solving simple 1D diffusion/parabolic equation
disp(' Diffusion, u_t = u_xx')
u = exp(-(x-.5).^2 / .1^2);
A = A2;
T = 0.01;
k = .5*h^2;
for n=1:ceil(T/k)
    u = (I - k/2*A) \ ((I + k/2*A)*u);
    if mod(n,20) == 0, plot(x,u,'-d'), axis([0,1,-.1,1.1]), grid on, drawnow, end
    if mod(n,900) ==0, pause 
        print 'Running step=%d', n
    end 
end
pause
% Solving advection -diffusion equation
disp(' Advection - Diffusion, u_t = -u_x + 5e-3 * u_xx ')
u = exp(-(x-.5).^2 / .1^2);
A = -A1 + 5e-3*A2;
T = 1.0;
k = .5*h;
for n=1:ceil(T/k)
    u = (I - k/2*A) \ ((I + k/2*A)*u);
    if mod(n,10) == 0, plot(x,u,'-d'), axis([0,1,-.1,1.1]), grid on, drawnow, end
    if mod(n,600) ==0, pause 
        print 'Running step=%d', n
    end 
end
pause
% Solving dispersion equation
disp(' Dispersion, u_t = u_xxx')
u = exp(-(x-.5).^2 / .1^2);
A = A3;
T = 1e-4;
k = 5*h^3;
for n=1:ceil(T/k)
    u = (I - k/2*A) \ ((I + k/2*A)*u);
    if mod(n,20) == 0, plot(x,u,'-d'), axis([0,1,-.5,1.1]), grid on, drawnow, end
     if mod(n,400) ==0, pause 
        print 'Running step=%d', n
    end 
end
pause

 disp(' Advection - Dispersion, u_t = -u_x + 1e-4 u_xxx ')
 u = exp(-(x-.5).^2 / .1^2);
 A = -A1 + 1e-4*A3;
 T = 1.0;
 k = .5*h;
 for n=1:ceil(T/k)
     u = (I - k/2*A) \ ((I + k/2*A)*u);
     if mod(n,10) == 0, plot(x,u,'-d'), axis([0,1,-.5,1.1]), grid on, drawnow, end
     if mod(n,300) ==0, pause 
         print 'Running step=%d', n
     end 
 end
 pause
        
 disp(' Advection - Dispersion, u_t = -u_x - 1e-4 u_xxx ')
 u = exp(-(x-.5).^2 / .1^2);
 A = -A1 - 1e-4*A3;
 T = 1.0;
 k = .5*h;
 for n=1:ceil(T/k)
     u = (I - k/2*A) \ ((I + k/2*A)*u);
     if mod(n,10) == 0, plot(x,u,'-d'), axis([0,1,-.5,1.1]), grid on, drawnow, end
     if mod(n,300) ==0, pause 
         print 'Running step=%d', n
     end 
 end
 pause

 disp(' Advection - Diffusion - Dispersion, u_t = -u_x + 5e-3 u_xx + 5e-4 u_xxx ')
 u = exp(-(x-.5).^2 / .1^2);
 A = -A1 + 5e-3*A2 + 5e-4*A3;
 T = 1.0;
 k = .5*h;
 for n=1:ceil(T/k)
     u = (I - k/2*A) \ ((I + k/2*A)*u);
     if mod(n,10) == 0, plot(x,u,'-d'), axis([0,1,-.5,1.1]), grid on, drawnow, end
     if mod(n,300) ==0, pause 
         print 'Running step=%d', n
     end 
 end
 pause
