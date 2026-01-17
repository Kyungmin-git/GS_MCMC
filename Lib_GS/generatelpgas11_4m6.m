function y=generatelpgas11_4m6(x)
global stepsize 
%para=rb+1;
randstep=randn(1,9);
%if r == 0
%randstep(9:11) = 0;
%else
%randstep(1:8) = 0;
%end
y=x+stepsize.*randstep;

    if (y(1)<973.15) y(1)=2*973.15-y(1); end 
    if (y(1)>1473.15)  y(1)=2*1473.15-y(1); end    
    if (y(2)<0) y(2)=2*0-y(2); end
    if (y(2)>10000) y(2)=2*10000-y(2); end
    if (y(3)<5) y(3)=2*5-y(3); end
    if (y(3)>40) y(3)=2*40-y(3); end
    if (y(4)<1) y(4)=2*1-y(4); end
    if (y(4)>80) y(4)=2*80-y(4); end
    if (y(5)<(6)) y(5)=2*(6)-y(5); end
    if (y(5)>12) y(5)=2*12-y(5); end
    if (y(6)<-log10(3)) y(6)=2*(-log10(3))-y(6); end
    if (y(6)>2) y(6)=2*2-y(6); end
    if (y(7)<-1) y(7)=2*(-1)-y(7); end
    if (y(7)>(1+log10(3))) y(7)=2*(1+log10(3))-y(7); end
    if (y(8)<2800) y(8)=2*2800-y(8); end
    if (y(8)>3100) y(8)=2*3100-y(8); end
    if (y(9)<20) y(9)=2*20-y(9); end
    if (y(9)>400) y(9)=2*400-y(9); end





end