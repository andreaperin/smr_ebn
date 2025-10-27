function vals = randRange(range, N)
    % Uniform random sampling in [min, max]
    vals = range(1) + (range(2)-range(1)) * rand(1, N);
end