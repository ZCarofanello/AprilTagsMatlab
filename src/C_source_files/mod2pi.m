function value = mod2pi(value, ref)
if nargin == 1
    TwoPi = 2 * pi;
    TwoPiInv = 1 / TwoPi;
    absVal = abs(value);
    q = absVal*TwoPiInv + 0.5;
    qi = floor(q);
    r = absVal - qi * TwoPi;
    if(value > 0)
        value = r;
    else
        value = -r;
    end
else
    value = mod2pi(value - ref) + ref;
end
end