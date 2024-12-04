L = 1.0;
n = 256;
dx = L/(n-1);
x = (0:dx:L);

omega = (2*pi/n)*(-(n)/2:(n)/2 - 1);
omega = fftshift(omega);

u0 = 0*x;
u0((n/2 - n/8):(n/2 + n/8)) = 1;

u0 = fft_(u0,zeros(1,n),n,exp(2*pi*1i/n));
dt = 1;
t = (0:dt:100);
u_hat = zeros(length(t),n);
u = zeros(length(t),n);

for k = 1:length(t)
    u_hat(k,:) = u0_hat.*exp(-omega.^2*t(k));
    u(k,:) = (1/n)*fft_(u_hat(k,:),zeros(1,n),n,exp(-2*pi*1i/n));
end

figure
hold on
for k = 1:10:length(t)
    plot(x,real(u(k,:)),"LineWidth",1.5)
    xlabel("x")
    ylabel("u(t,x)")
end
figure
hold on, grid on
for k = 1:10:length(t)
    plot3(k*ones(length(x)),x,real(u(k,:)))
    xlabel("t")
    ylabel("x")
    zlabel("u(t,x)")
end
view(3)

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