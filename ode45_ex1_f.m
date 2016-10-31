function dx = ode45_ex1_f(t,x)
% make sure dx is a column vector
dx = zeros(2,1);

% define inertia matrix and inverse of inertia matrix
I = [3 -0.2;-0.1 4];
invI=inv(I);

%define the input
if t < 1
    u = 0;
else
    u = 1;
end

% define right-hand side of ODEs
dxr(1,1) = -2*x(1)+0.5*x(2)+u;
dxr(2,1) = 0.3*x(1)-5*x(2);

% form function output, dx
dx = invI*dxr
pause

end