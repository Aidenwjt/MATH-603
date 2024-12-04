%%% Program to solve numerically the 1D heat equation using the FFT

L = 1.0; % length of spacial domain
n = 256; % number of discrete points
dx = L/n; % spatial step-size
x = (-L/2:dx:L/2 - dx); % spatial domain

zeta = (2*pi/L)*(-(n/2):(n/2)-1);
zeta = fftshift(zeta,1)';

u0 = 0*x;
u0
u0((L/2 - L/5)/dx:dx:(L/2 + L/5)/dx) = 1;
u0_hat = fft(u0');

dt = 0.1;
t = (0.0:dt:10.0);

[t,u_hat] = ode45(@(t,u_hat)rhsHeat(t,u_hat,zeta),t,u0_hat);

u = zeros(length(t),n);
for k = 1:length(t)
    u(k,:) = ifft(u_hat(k,:));
end
u

surf(x,t,u)

%%% Functions
function F = fft_(f,F,n,w)
    if n == 1
        F(1) = f(1);
    else
        f_even = zeros(1,n/2);
        f_odd = zeros(1,n/2);
        for k = 1:((n/2))
            f_even(k) = f(2*(k-1)+1);
            f_odd(k) = f(2*(k-1)+2);
        end
        F_even = fft_(f_even,zeros(1,n/2),n/2,w^2);
        F_odd = fft_(f_odd,zeros(1,n/2),n/2,w^2);
        for k = 1:n
            F(k) = F_even(mod((k-1),n/2)+1) + (w^(k-1))*F_odd(mod((k-1),n/2)+1);
        end
    end
end

function du_hat_dt = rhsHeat(t,u_hat,zeta)
    du_hat_dt = -(zeta.^2).*u_hat;
end