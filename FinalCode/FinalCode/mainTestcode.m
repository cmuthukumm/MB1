clc;
clear all;
close all;
warning off all;
restoredefaultpath;
addpath(genpath(pwd));
addpath('subcodes');
global N X Y Xb Yb E zX zY max1

%% PARAMETER INITIALIZATION
P=0.1;
Eo=0.5;
Echeck=Eo;
rmax=5000;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
EDA=5*0.000000001;
do=sqrt(Efs/Emp);
rand('seed',1)
N=100;                                   % Total No. of Nodes
RSU=[];
min1=0;
max11=400;
max1=200;
X = min1+(max11-min1)*rand(1,N);
Y = min1+(max1-min1)*rand(1,N);
Xb = max11+2+rand(1,1);
Yb =(max1+2)+rand(1,1);
maxv=50;
minv=10;
vel=rand(1,N);
E=rand(1,N);
minTh=E/2;    

figure(1)
for i2 = 1:N 
          plot(X(i2),Y(i2),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8'); 
           title('Network Area');
           xlabel('X in m')
          ylabel('Y in m')
          text(X(i2), Y(i2), num2str(i2),'FontSize',10); 
          hold on;
          xlim([0 400])
          ylim([0 200])
end
hold on
         plot(Xb,Yb,'s','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',12'); 
          title('Network Area');
          xlabel('X in m')
          ylabel('Y in m')
          xlim([0 400])
          ylim([0 200])
          text(Xb, Yb, 'Sink','FontSize',10); 
          hold on;
zX=0:100:max11;
zY=0:100:max1;
id=zeros(1,N);
ik1=1;
for ik=1:numel(zX)-1
    for ij=1:numel(zY)-1
    rectangle('Position',[zX(ik) zY(ij) 100 100])
    points =[zX(ik) zY(ij); zX(ik)+100 zY(ij);  zX(ik)+100 zY(ij)+100 ;zX(ik) zY(ij)+100 ;zX(ik) zY(ij)];
    if(~isempty(RSU))
          plot(RSU(1),100,'^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',11'); 
          xlabel('X in m')
          ylabel('Y in m')
          xlim([0 400])
          ylim([0 200])
          text(Xb, Yb, 'Sink','FontSize',10); 
          hold on;
          title('Network Area');
    end
    
    figure(2)
    [in,on]=inpolygon(X,Y,points(:,1),points(:,2));
    ind=find(in==1);
     if(~isempty(RSU))
         if(RSU(1)>zX(ik) && RSU(1)<zX(ik+1))
         ind=[];
         end
     else
    Ec=E(ind);      
    SN1=ind;
    SN(ik1).id=ind;
    L=0.*ind; 
    [val1,ind2]=sort(Ec,'descend');
    MonitorID(ik1)=ind(ind2(1));
    ipd0=find(ind==ind(ind2(1)));
    L(ipd0)=1; 
    plot(X(MonitorID(ik1)),Y(MonitorID(ik1)),'o','LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10'); 
           title('Network Area');xlabel('X in m')
          ylabel('Y in m')
          xlim([0 400])
          ylim([0 200])
          text(Xb, Yb, 'Sink','FontSize',10); 
          hold on;      
    chmtosinkdist=abs(X(MonitorID(ik1))-Xb); 
    pc=distpch(chmtosinkdist);
    s1=ceil((pc/100)*numel(ind));
    SN(ik1).s1=s1;   
    
%% CLUSTER FORMATION & CH SELECTION

    pos1=[X(ind) Y(ind)];
    vel1=vel(ind);
    Bs=E(ind);
    posd=pdist2([X(ind); Y(ind)]',[zX(ik)+50 zY(ij)+50]);
    w1=0.6;w2=0.2;w3=0.2; W=w1.*posd+w2.*vel1'+w3.*Bs';
    rng('default')
    [idx3,C,sumdist3] = kmeans([X(ind)+W'; Y(ind)+W']',s1,'Distance','cityblock','Display','final');
    idx=unique(idx3);
    for ip=1:numel(idx)
       ind3= find(idx3==idx(ip));
       [val,ind2]= min(abs(W(ind3)-C(ip)));
       L((ind3(ind2(1))))=ip+1;     
       ind3(ind2(1))=[];
       L((ind3))=ip+1+0.1;           
    end
    SN(ik1).L=L;
    SN(ik1).Z=(ik1).*ones(1,numel(L));
     end
    L=[];
    ik1=ik1+1;
    hold on
    end
end
L1=[SN.L];
id1=[SN.id];
ind4=find(L1==1);
hold on
plot(X(id1(ind4)),Y(id1(ind4)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8'); 
    title('Network Area');
    xlabel('X in m')
            ylabel('Y in m')
            xlim([0 400])
            ylim([0 200])
            text(Xb, Yb, 'Sink','FontSize',10); 
for iv=2: max([SN.s1])
 ind5=find(round(L1)==iv); 
 ind6=find(L1==iv);  
hold on
plot(X(id1(ind5)),Y(id1(ind5)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8'); 
 title('Network Area');
 xlabel('X in m')
            ylabel('Y in m')    
hold on
plot(X(id1(ind6)),Y(id1(ind6)),'^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',11'); 
 title('Network Area');            xlabel('X in m')
            ylabel('Y in m')            
end
xlim([0 400])
ylim([0 200])
 ProblemFunction = @Ackley;
    DisplayFlag = true;
    RandSeed = round(sum(100*clock));
          
outF=CHselect();
round=200;
delv=[5 10 20 30 40];
R=[80 100 130 160 190];
Az=2;
for iu=1:5
ik=1;
while(ik<=70)
for ik1=1:numel(zX)-1
  for ij=1:numel(zY)-1
     rectangle('Position',[zX(ik1) zY(ij) 100 100])
     points =[zX(ik1) zY(ij); zX(ik1)+100 zY(ij);  zX(ik1)+100 zY(ij)+100 ;zX(ik1) zY(ij)+100 ;zX(ik1) zY(ij)];
  end
end
    z1=[SN.Z];
    L1=[SN.L];
    id1=[SN.id];
    indw=find(mod(z1,2)==0);
    vel1=vel;
    vel1(id1(indw))=-vel(id1(indw));
    if(Az==1)
     X=X+delv(iu).*vel1;
    else
     X=X+delv(2).*vel1;   
    end 
    X(X<-20)=400;
    X(X>420)=0;
    for iv=2: max([SN.s1])
    ind5=find(floor(L1)==iv);   
    ind6=find(L1==iv);  
    hold on
    plot(X(id1(ind5)),Y(id1(ind5)),'o','LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8'); 
      title('Network Area');       xlabel('X in m')
            ylabel('Y in m')  
    hold on
    plot(X(id1(ind6)),Y(id1(ind6)),'^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',11'); 
     title('Network Area');        xlabel('X in m')
            ylabel('Y in m')            
    end
    A=randperm(N);
    path=[];
    L1=[SN.L];
    id1=[SN.id];
    z1=[SN.Z];
    A(1)
    indw=find(id1==A(1));
    Zo=z1(indw);
    Lz=SN(Zo).L;
    indb=find(Lz==floor(L1(indw)));
    path1=SN(Zo).L(indb);
    if(ik>20)
    x=[X' Y' vel' E'];
    tic
    y=evalfis(x,outF);
    tp1=toc
    ind2w=find(uint8(y)==1);
    indCH1=[ind2w' N+1]; 
    else
    tp1=1;
    indCH=(find(L1==1 | L1==2 |  L1==3 | L1 ==4 | L1==5 | L1==6 | L1==7));
    indCH=id1(indCH);
    indCH1=[indCH N+1];
    end
    X1=[X Xb];
    Y1=[Y Yb];
    hold on
    plot(X1(indCH1),Y1(indCH1),'^','LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','r','MarkerSize',11'); 
            xlabel('X in m')
            ylabel('Y in m')  
            hold on
            for icc=1:numel(indCH1)
            text(X1(indCH1(icc)),Y1(indCH1(icc)),'CH','FontSize',10); 
            end
            hold on
            plot(X1(end),Y1(end),'^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',11'); 
            xlabel('X in m')
            ylabel('Y in m')        
   distCHtoCH=pdist2([X1(indCH1);Y1(indCH1)]',[X1(indCH1);Y1(indCH1)]'); 
   trust=distCHtoCH;
   if(Az==2)
   trust(distCHtoCH>(R(iu)))=inf;
   else
   trust(distCHtoCH>(R(3)))=inf;    
   end
   w1=trust;
   [r_path, r_cost] = Predictive(path1,size(trust,1), w1);
   [r_pathE, r_costE] = hopbyhop(path1,size(trust,1), w1);
   indCH1(r_path)
   t_cost(ik)=r_cost;
   t_path{ik}=r_path;
   path=indCH1(r_path);
   for p =1:(length(path)-1) 
   line([X1(path(p)) X1(path(p+1))], [Y1(path(p)) Y1(path(p+1))], 'Color','b','LineWidth',2.5, 'LineStyle','-') 
   arrow([X1(path(p)) Y1(path(p)) ], [X1(path(p+1)) Y1(path(p+1)) ])
   end 
xlim([0 400])
ylim([0 200])

DR(ik)= 0;
if(~isempty(path))
if(path(end)==61)
DR(ik)= 1;
end    
end
E2E(ik)=numel(r_path);
E2Ex(ik)=numel(r_pathE);
CHF(ik)=tp1;
CHFe(ik)=(2).*tp1;     
ik=ik+1; 
pause(0.02);
clf;
end
E2E(isinf(E2E))=0;
E2Ex(isinf(E2Ex))=0;
E2Edelay(iu)=mean(CHFe);
end
[MinC] = MBO(ProblemFunction, DisplayFlag, RandSeed);
[assignment,cost] = munkres(r_cost);

%% PERFORMANCE ANALYSIS

[PDR1,PDR2,PDR3,PDR4,PDR5,RE1,RE2,RE3,throughput1,throughput2,throughput3,RE,throughput,RE4,throughput4]=Analysis(N,P,Eo,ETX,ERX,Efs,Emp,EDA,rmax,do);
results(PDR1,PDR2,PDR3,PDR4,RE1,RE2,RE3,throughput1,throughput2,throughput3,RE,throughput,RE4,throughput4,PDR5);

