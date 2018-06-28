function value = mod2pi(value, ref)
if nargin == 1
    value = value - 2*pi*floor( (value+pi)/(2*pi) );
    return;
end
shifted = value - ref;
value = (shifted - 2*pi*floor( (shifted+pi)/(2*pi) ))+ref;
end