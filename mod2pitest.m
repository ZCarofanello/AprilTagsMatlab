y = [pi():1/pi():8*pi()]
test = mod2pi(y);
plot(y);
plot(test);




function value = mod2pi(value)
value = value - 2*pi*floor( (value+pi)/(2*pi) );
end