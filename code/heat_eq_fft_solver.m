% Program that solves a specific 1-dimensional heat equation IVP
% using properties of the Fourier Transform and the FFT algorithm.
%
% Course: MATH 603
% Author: Aiden Taylor
% UCID: 30092686
% Date: December 5th, 2024

%% Initialize our spatial and frequency domains
L = 1.0;
n = 256;
dx = L/(n-1);
x = (0:dx:L);
omega = (2*pi/n)*(-(n)/2:(n)/2 - 1);
omega = fftshift(omega);

%% Initialize our initial-value function
u0 = 0*x;
u0((n/2 - n/8):(n/2 + n/8)) = 1;

%% Plot the initial-value function
figure
plot((0:255)*dx,u0,"LineWidth",1.5)
axis([0 1 -0.2 1.2])
xlabel("x")
ylabel("u(0,x)")

%% Compute the FFT of the initial-value function
u0_hat = fft_(u0,zeros(1,n),n,exp(2*pi*1i/n));

%% Define our temporal domain
dt = 1;
t = (0:dt:100);

%% Initialize matrices for our Fourier coefficients and unknown function
u_hat = zeros(length(t),n);
u = zeros(length(t),n);

%% Solve each ODE and compute the IFFT of the resulting solutions
for k = 1:length(t)
    u_hat(k,:) = u0_hat.*exp(-omega.^2*t(k));
    u(k,:) = (1/n)*fft_(u_hat(k,:),zeros(1,n),n,exp(-2*pi*1i/n)); % don't need IFFT function, just use FFT function this way
end

%% Plot numerical solutions in 2-dimensions
figure
hold on
for k = 1:10:length(t)
    plot(x,real(u(k,:)),"LineWidth",1.5)
    xlabel("x")
    ylabel("u(t,x)")
end

%% Plot numerical solution in 3-dimensions
figure
hold on, grid on
for k = 1:10:length(t)
    plot3(k*ones(length(x)),x,real(u(k,:)))
    xlabel("t")
    ylabel("x")
    zlabel("u(t,x)")
end
view(3)

%% FFT function
function F = fft_(f,F,n,xi)
    if n == 1
        F(1) = f(1);
    else
        f_even = zeros(1,n/2);
        f_odd = zeros(1,n/2);
        for k = 1:((n/2))
            f_even(k) = f(2*(k-1)+1);
            f_odd(k) = f(2*(k-1)+2);
        end
        F_even = fft_(f_even,zeros(1,n/2),n/2,xi^2);
        F_odd = fft_(f_odd,zeros(1,n/2),n/2,xi^2);
        for k = 1:n
            F(k) = F_even(mod((k-1),n/2)+1) + (xi^(k-1))*F_odd(mod((k-1),n/2)+1);
        end
    end
end
