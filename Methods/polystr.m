function s = polystr(p)
% POLYSTR Converts a vector of polynomial coefficients to a character string.
% S is the string representation of P.

if all(p == 0) % All coefficients are 0.
    s = '0';
else
    d = length(p) - 1; % Degree of polynomial.
    s = []; % Initialize s.
    for a = p
        if a ~= 0 % Coefficient is nonzero.
            if ~isempty(s) % String is not empty.
                if a > 0
                    s = [s ' + ']; % Add next term.
                else
                    s = [s ' - ']; % Subtract next term.
                    a = -a; % Value to subtract.
                end
            end
            if a ~= 1 || d == 0 % Add coefficient if it is ~=1 or polynomial is constant.
                s = [s num2str(a)];
                if d > 0 % For nonconstant polynomials, add *.
                    s = [s '*'];
                end
            end
            if d >= 2 % For terms of degree > 1, add power of x.
                s = [s 'x^' int2str(d)];
            elseif d == 1 % No power on x term.
                s = [s 'x'];
            end
        end
        d = d - 1; % Increment loop: Add term of next lowest degree.
    end
end
end

