function sample = mysample(a, b, f, N)
%MYSAMPLE rejection sampler
%   random sample (size N) a distribution function f over domain x = [a,b]
%   Input: x, f, N
%   Output: sample

x = linspace(a, b, 1e6);
fx = f(x);
fmax = max(fx);

sample = zeros(N,1);

cc = 1;
while cc <= N
    xprop = a + (b-a) * rand();
    yprop = fmax * rand();
    if yprop < f(xprop)
        sample(cc) = xprop;
        cc = cc + 1;
    end
end

end

