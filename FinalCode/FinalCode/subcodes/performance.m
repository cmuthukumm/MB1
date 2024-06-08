function [throughput]= performance( packetsize,m,T)

packet=packetsize;
pdr=(packet /packet ) ./((1)^-2/ 5 /packet)*m ;
av = 0;
bv = 1;
rt1 = (bv-av).*rand(m,1) + av;
format short
throughput=packet/T;

end

