function value = mod2pi(value, ref)
if nargin == 1
    value = value - 2*pi*floor( (value+pi)/(2*pi) );
else
    value = mod2pi(value - ref) + ref;
end
end