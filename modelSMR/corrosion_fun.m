function C = corrosion_fun(age, Cr, maint_interval, repair_frac)
    % CORROSION_FUN
    % Supports:
    % - age:   vector, Cr: scalar
    % - age:   scalar, Cr: vector
    % - age & Cr: vectors of the same size

    % Broadcast scalar to match the size of the other argument
    if isscalar(age) && ~isscalar(Cr)
        age = ones(size(Cr)) * age;
    elseif ~isscalar(age) && isscalar(Cr)
        Cr = ones(size(age)) * Cr;
    elseif ~isequal(size(age), size(Cr))
        error('age and Cr must be the same size, or one of them must be scalar.');
    end

    C = zeros(size(age));

    for k = 1:numel(age)
        a   = age(k);
        Crk = Cr(k);

        n   = floor(a / maint_interval);   % number of full intervals
        rem = mod(a, maint_interval);      % partial interval

        Ck = 0;

        % full intervals
        for i = 1:n
            new_corr = Crk * maint_interval;          % corrosion in this interval
            Ck = Ck + new_corr * (1 - repair_frac);   % only part of it remains
        end

        % last partial interval: no maintenance yet
        Ck = Ck + Crk * rem;

        C(k) = Ck;
    end
end
