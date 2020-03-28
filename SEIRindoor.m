function [ dydt ] = SEIRindoor(t,y,alpha,gamma)
% S_LC=y(1),E_LC=y(2),I_LC=y(3),R_LC=y(4)
dydt = zeros(4,1);
dydt(1)=0;
dydt(2)=-(gamma*y(2));
dydt(3)=(gamma*y(2))-(alpha*y(3));
dydt(4)=alpha*y(3);
end

