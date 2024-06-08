function [packetsize,RE,T]=updation(noOfNodes,m)
sor= randi(noOfNodes,1,1); 
des=randi(noOfNodes,1,1);
 if (sor==des);
    sor = sor+1;
 else
[ netXloc,netYloc ] = perf(noOfNodes );
% [X,Y,ert]= networks( noOfNodes, sor,des);
mn=noOfNodes;
A=randi(mn);
for i=1:mn
    for j=1:mn
        if (i==j)
             A(i,j)=0;
        else
%             A(i,j)=A(i,j);
        end
    end
end
tic;
packetsize=5;
rpg=8;
rand('state', 0);
a = 0;
b = 1;
w1 = (b-a).*rand(mn,1) + a;
w2 = (b-a).*rand(mn,1) + a;
message='hai this is to be tested';
a = dec2bin(message,8)';
a=randi(mn,mn);
RE=toc;
for lp1=1:mn
   dt(lp1)=a(lp1).*((a(lp1).^1)); 
end
directtrust=dt;
for lp2=1:mn
   idt(lp2)=sum(directtrust(lp2));
end
indirecttrust=idt;
for lp3=1:mn
    neighbour(lp3)=sum(A(:,lp3));
   tt(lp3)=w1(lp3).*dt(lp3)+w2(lp3).*(idt(lp3)./neighbour(lp3)); 
end
totaltrust=round(tt);
nodes=[netXloc;netYloc]';

for i=1:mn
    distance(:,i)=sqrt((nodes(i,1).^2)-(nodes(i,2)).^2);
end
dis=real(distance);
 end
  pdelivered=sum(a(:));
  T=toc;
 