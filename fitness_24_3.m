function [Z,x]=fitness_24_3(x,lb,ub,data1)
pause(10)
hours=24;
for ww=1:hours
%% initial variables

Dt=data1(ww,1);

%% effect of constraints
% effect of power balance

production=x(ww,1)+x(ww,2)+x(ww,4)+x(ww,5)+x(ww,6)+x(ww,7)+x(ww,8)+x(ww,9);
if production-Dt~=0
if production<Dt
   
        x(ww,3)= Dt-production;
        if x(ww,3)>ub(ww,3)
        penalty1(ww)=(10e16)*abs(Dt-production);
        x(ww,3)=ub(ww,3);
        else
        penalty1(ww)=0;
        end
else
    % dar in ghesmat agar tolid bishtar bashad bayad rahi andishid
    penalty1(ww)=(10e16)*abs(Dt-production);
end
end

% effect of ramp up and down rate
if ww~=1
    if x(ww-1,1:3)-x(ww,1:3)<=[3 5 8] & x(ww-1,1:3)-x(ww,1:3)>=[-3 -5 -8]
        penalty2(ww)=0;
    else
        penalty2(ww)=10*(sum(x(ww-1,1:3)-x(ww,1:3))).^4;
    end
end

D(ww,:)=x(ww,10:12)-([1.079 1.378 1.847]/1000.*x(ww,7:9).^2+[1.32 1.63 1.64]/1000.*x(ww,7:9)-[1.32 1.63 1.64]/1000.*x(ww,7:9).*[0 .45 .9]);

end

%% final cost calculation

A=sum(0.06*x(:,1).^2+0.5*x(:,1)+0.03*x(:,2).^2+0.25*x(:,2)+0.04*x(:,3).^2+0.3*x(:,3)); % conv. gen. cost calculation

gama=1;
B=gama.*sum(abs(x(:,4))); %transmission cost calculation

lambda=data1(:,2);
lambda=repmat(lambda,1,3);
C=sum(-x(:,10:12)+lambda/1000.*x(:,7:9)); % demand response cost calculation
C=sum(C);

penalty3=0;
if sum(D)>=[10 10 10] %sum(sum(D))>=[0 0 0]
else
   penalty3=10e14*abs(sum(sum(D).^2));% benefit must be greater than zero
end

penalty4=0;
if sum(sum(x(:,10:12)))>= 500
    penalty4=10e11*(sum(sum(x(:,10:12)))-500).^6; %pen. for daily budget
end

penalty5=0;
if sum(x(:,7:9)) <= [30 35 40]
else
    penalty5=10e10*(sum(sum(x(:,7:9)))-(30+35+40)).^2; % pen. for daily limit of interruptable power
end

Z = 0.5*(A+B)+(0.5*abs(C))+sum(penalty1)+sum(penalty2)+penalty3+penalty4+penalty5;
end