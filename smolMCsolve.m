function sol = smolMCsolve( a, b, f0, kernel, source, t_out, seed )
%SMOLUCHOWSKIMC solve Smoluchowski equations with MC method
%  Input: domain [a,b], initial distribution function f0, pop size N, kernel, output times, rng seed
%  Output: array of random population corresponding to the times in t_out

% set random number generator seed
rng(seed);

% draw from initial distribution

x = linspace(a, b, 1e6);
N = round(trapz(x,f0(x)));
Mtot = trapz(x,x.*f0(x));
assert(N>1);
sample = mysample(a, b, f0, N);
newsample = zeros(size(sample));

sol = cell(numel(t_out),1);

Proll = rand(N,N);
for i = 1:N
    Proll(i,i) = 0;
    Proll(i,1:i-1) = 1;
end

% start simulation
cc = 1;
t = 0;
t_end = t_out(end);
while t < t_end
    
    % calculate probability matrix of 2 masses interacting
    Pmatrix = kernel(sample,sample')/Mtot;
    assert(numel(Pmatrix) == N^2); % check to see kernel implemented properly
    prob = sum(Pmatrix,2);
    
    % set timestep
    dt = 0.9/max(prob); % timestep is chosen such that any given particle has 90% chance of interacting with something in the timestep
    save_this_turn = 0;
    if t + dt > t_out(cc)
        dt = t_out(cc) - t;
        save_this_turn = 1;
    end
    Pmatrix = Pmatrix * dt;
    
    % simulate random interactions
    % interactMatrix = rand(N,N) < Pmatrix; % the expensive step -- don't reroll all the time, use random matrix + sample permutation
    Proll = Proll(1:N,1:N); %%%
    interactMatrix = Proll < Pmatrix;
    for i = 1:N
        % interactMatrix(i,1:i-1) = 0;
        interactMatrix(i,i) = 1;
    end
    for i = 1:N
        js = interactMatrix(i,:); %%%
        js(i) = 0;
        interactMatrix(js,:) = 0;
        jj = 1:N;
        jj = jj(js);
        for j = jj
            interactMatrix(i+1:end,j) = 0;
        end
    end
    
    % construct new sample
    ss = 1;
    for i = 1:N
        if interactMatrix(i,i) == 1
            newsample(ss) = sum(sample(interactMatrix(i,:))); %%%
            ss = ss + 1;
        end
    end
    newsample = newsample(1:ss-1);
    sample = newsample(randperm(numel(newsample))); % randomly permute sample
    N = numel(sample);
    t = t + dt;
    
    % store solution at requested save times
    if save_this_turn
        sol{cc} = sample;
        sprintf('t=%0.4f, N=%d',t,N)
        cc = cc + 1;
    end
    
end

end





