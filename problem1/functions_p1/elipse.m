function ret=elipse(evector,eval,center)

n=20;
ret=zeros(2,n+1);
for i=1:n+1
ret(:,i)=center+evector*[eval(1,1)*cos(2*pi*i/n);eval(2,2)*sin(2*pi*i/n)];

end

end   
