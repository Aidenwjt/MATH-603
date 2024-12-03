for l = 1:10
    n = 2^l;
    f = vec(n);
    matlab_fft = fft(f);
    my_fft = fft_(f,zeros(1,n),n,exp((2*pi*1i)/n));
    dist = norm(matlab_fft - my_fft',2);
    fprintf("n = %d: ||matlab_fft - my_fft|| = %f\n",n,dist);
    matlab_ifft = ifft(matlab_fft);
    my_ifft = fft_(my_fft,zeros(1,n),n,exp(-(2*pi*1i)/n));
    dist = norm(matlab_ifft - (1/n)*my_ifft',2);
    fprintf("n = %d: ||matlab_ifft - my_ifft|| = %f\n",n,dist);
end

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