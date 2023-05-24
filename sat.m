function a = sat(b)
i = size(b,1);
a = zeros(i,1);
for t = 1:i
    if b(t,1)>=1
    a(t,1)=1;
    elseif b(t,1)<=-1
    a(t,1)=-1;
    else
    a(t,1)=b(t,1);
    end
end
end
