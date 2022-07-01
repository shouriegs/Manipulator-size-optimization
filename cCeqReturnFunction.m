function [c,ceq]=cCeqReturnFunction(le)

a=le/3; %length OA
b=le; %length OB
h=le*2/3; %length OH
%le=1*150; %length of links
n=100; %length of end-effector

% x=-70:20:70;
% y=(-1.2*le - 70):10:(-1.2*le+70);
% phi=-30:5:30;

wcx=0;
wcy=-1.2*le;
wcr=50;

x=linspace(-70,70,10);
y=linspace((wcy - 70),(wcy + 70),10);
phi=linspace(-30,30,10);

lenx=length(x);
leny=length(y);
lenphi=length(phi);

c=zeros(1,(lenx*leny*lenphi)*2+1);

score=0;
counter=1;
for i=1:lenx
    for j=1:leny
        for k=1:lenphi
            yesInside=inOrOut(x(i),y(j),wcx,wcy,wcr);
            if yesInside==1
                %fprintf("here");
                ICN=findICN(x(i),y(j),phi(k),le);
%                 if abs(ICN(1))<0.1 | abs(ICN(2))<0.1
%                     score=score+1
%                 else
%                     
%                 end
                c(counter)=0.0001-ICN(1);
                counter=counter+1;
                c(counter)=0.0001-ICN(2);
                counter=counter+1;
            end
        end
    end
end
ceq=[];