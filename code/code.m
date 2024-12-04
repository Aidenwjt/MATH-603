clear all
a = 1;
L = 100;
N = 256;
dx = L/N;
x = -L/2:dx:L/2-dx;
x

kappa = (2*pi/L)*[-N/2:N/2-1];
kappa = fftshift(kappa);

u0 = x*0;
u0((L/2-L/10)/dx:(L/2 + L/10)) = 1;

t = 0:0.1:20;
[t,uhat]=ode45(@(t,uhat)rhsHeat(t,uhat,kappa,a),t,fft(u0));

for k = 1:length(t)
    u(k,:) = ifft(uhat(k,:));
end
%% Plot solution in time
figure, waterfall((u(1,10:end,:)))
figure, imagesc(flipud(u));

%% Functions
function duhatdt = rhsHeat(t,uhat,kappa,a)
    duhatdt = -a^2*(kappa.^2)'.*uhat;
end