function [ dydt ] = SEIRoutdoor(t,y,alpha,beta,gamma,N)
% S_LD=y(1),E_LD=y(2),I_LD=y(3),R_LD=y(4),S_FC=y(5),E_FC=y(6),I_FC=y(7),R_FC=y(8),S_FD=y(9),E_FD=y(10),I_FD=y(11),R_FD=y(12)
dydt = zeros(12,1);
%infects=y(2)+y(3)+y(6)+y(10)+y(11);
infects=y(2)+y(6)+y(10);
dydt(1)=-((beta*y(1)*infects)/N);
dydt(2)=((beta*y(1)*infects)/N)-(gamma*y(2));
dydt(3)=(gamma*y(2))-(alpha*y(3));
dydt(4)=alpha*y(3);
dydt(5)=-((beta*y(5)*infects)/N);
dydt(6)=((beta*y(5)*infects)/N)-(gamma*y(6));
dydt(7)=(gamma*y(6))-(alpha*y(7));
dydt(8)=alpha*y(7);
dydt(9)=-((beta*y(9)*infects)/N);
dydt(10)=((beta*y(9)*infects)/N)-(gamma*y(10));
dydt(11)=(gamma*y(10))-(alpha*y(11));
dydt(12)=alpha*y(11);
end

