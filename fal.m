function output = fal(e,sigma,delta)
i = size(e,1);
output = zeros(i,1);
for t = 1:i
    if abs(e(t))<=delta
        output(t) = e(t)/delta^(1-sigma);
    else
        output(t) = (abs(e(t))^sigma)*sign(e(t));
    end
end
end