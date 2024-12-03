%%% First let's compare our fft_ function to Matlabs fft function.
for l = 1:10
    n = 2^l;
    f = vec(n);
    matlab_fft = fft(f);
    my_fft = fft_(f,zeros(1,n),n,exp((2*pi*1i)/n));
    fprintf("n = %d: ||matlab_fft - my_fft|| = %f\n",n,norm(matlab_fft - my_fft',2));
    matlab_ifft = ifft(matlab_fft);
    my_ifft = fft_(my_fft,zeros(1,n),n,exp(-(2*pi*1i)/n));
    fprintf("n = %d: ||matlab_ifft - my_ifft|| = %f\n",n,norm(matlab_ifft - (1/n)*my_ifft',2));
end
%%% Next let's solve the 1D heat equation just using the Forward Time
%%% Central Space (FTCS) method
n = 256;
x0 = 0;   x1 = 1.0;   % set up interval for x
t0 = 0;   T = 1.0;   % set up interval for time

M = 100;
k = 1/(2*(n+1)^2);
h = (x1 - x0)/(n+1);
s = k/(h^2);
t = (t0:k:T);
x = (x0:h:T);
assert(s <= 1/2,"Not Stable.");

% TODO:
% - Implement FTCS iterative method
% - Plot results
% - Plot exact solution
% - Implement FFT method
% - Plot against FTCS method

%%% Functions
function f = vec(n)
    f = zeros(n,1);
    for k = 1:n
        f(k,1) = 1/(1 + ((k-1)/n)^2);
    end
end

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