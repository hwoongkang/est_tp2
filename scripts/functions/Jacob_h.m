function ret=Jacob_h(X)

x=X(1);y=X(2);

ret=[-y/(x^2+y^2),x/(x^2+y^2);
    (50-y)/(x^2+(y-50)^2),x/(x^2+(y-50)^2)];
end