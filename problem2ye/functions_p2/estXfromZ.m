function ret=estXfromZ(z)

ret=[50/(tan(z(1))-tan(z(2))); 50*tan(z(1))/(tan(z(1))-tan(z(2)))];
end