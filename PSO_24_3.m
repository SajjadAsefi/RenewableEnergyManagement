%% by Sajjad A
clc
clear
close all
format short g
%% parameters setting

Answer=zeros(24,12);

nvar=12; % number of variable

% load and lambda
data1=[31.83 1.57 
31.40 1.40
31.17 2.20
31.00 3.76
31.17 4.50
32.10 4.70
32.97 5.04
34.10 5.35
37.53 6.70
38.33 6.16
40.03 6.38
41.17 6.82
39.67 7.30
41.70 7.80
42.10 8.50
41.67 7.10
40.70 6.80
40.07 6.30
38.63 5.80
36.40 4.20
34.10 3.80
32.80 3.01
32.50 2.53
32.00 1.42];

lb=[1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0
   1 1 1 -4 0 0 0 0 0 0 0 0]; % lower bound

ub=[4 6 9 4 7.56    0   4 4 4 16 16 16 
    4 6 9 4 7.50    0   4 4 4 16 16 16
    4 6 9 4 8.25    0   4 4 4 16 16 16
    4 6 9 4 8.48    0   4 4 4 16 16 16
    4 6 9 4 8.48    0   4 4 4 16 16 16
    4 6 9 4 9.42    0   4 4 4 16 16 16
    4 6 9 4 9.82    0   4 4 4 16 16 16
    4 6 9 4 10.35 7.99  4 4 4 16 16 16
    4 6 9 4 10.88 10.56 4 4 4 16 16 16
    4 6 9 4 11.01 13.61 4 4 4 16 16 16
    4 6 9 4 10.94 14.97 4 4 4 16 16 16
    4 6 9 4 10.68 15    4 4 4 16 16 16
    4 6 9 4 10.42 14.78 4 4 4 16 16 16
    4 6 9 4 10.15 14.59 4 4 4 16 16 16
    4 6 9 4 9.67  13.56 4 4 4 16 16 16
    4 6 9 4 8.98  11.83 4 4 4 16 16 16
    4 6 9 4 8.37  10.17 4 4 4 16 16 16
    4 6 9 4 7.61  7.66  4 4 4 16 16 16
    4 6 9 4 6.70    0   4 4 4 16 16 16
    4 6 9 4 5.72    0   4 4 4 16 16 16
    4 6 9 4 7.21    0   4 4 4 16 16 16
    4 6 9 4 7.75    0   4 4 4 16 16 16
    4 6 9 4 7.88    0   4 4 4 16 16 16
    4 6 9 4 7.69    0   4 4 4 16 16 16];%upper bound

ub(:,7:9)=8*ones(24,3);
ub(:,10:12)=20*ones(24,3);

NP=20000;              % number particle
T=1000;                  % max of iteration

W=1;
C1=2;
C2=2;

alpha=0.05;

%%% initialization
tic
empty.pos=[];
empty.cost=[];
empty.velocity=[];

%load('result.mat');
particle=repmat(empty,NP,1);
% load('result_main2','gparticle')
% results;

for i=1:NP
particle(i).pos=lb+rand(24,nvar).*(ub-lb);
% particle(i).pos=gparticle.pos;
[particle(i).cost,particle(i).pos]=fitness_24_3(particle(i).pos,lb,ub,data1);
particle(i).velocity=0;
end

bparticle=particle;

[value,index]=min([particle.cost]);

gparticle=particle(index);

% main loop

best=zeros(T,1);
AVR=zeros(T,1);

for t=1:T

     for i=1:NP
         
          particle(i).velocity=W*particle(i).velocity...
                              +C1*rand(24,nvar).*(bparticle(i).pos-particle(i).pos)...
                              +C2*rand(24,nvar).*(gparticle.pos-particle(i).pos);
          
         particle(i).pos=particle(i).pos+particle(i).velocity;
         
         particle(i).pos=min(particle(i).pos,ub);
         particle(i).pos=max(particle(i).pos,lb);
          
         
         [particle(i).cost,particle(i).pos]=fitness_24_3(particle(i).pos,lb,ub,data1);
         
         
         if particle(i).cost<bparticle(i).cost
             bparticle(i)=particle(i);
             
             if bparticle(i).cost<gparticle.cost
                 gparticle=bparticle(i);
             end
         end
         
         

     end
     
     
     
 W=W*(1-alpha);
 
 best(t)=gparticle.cost;
 %AVR(t)=mean([particle.cost]);
 
  disp([ ' t = ' num2str(t)   ' BEST = '  num2str(best(t))]);
 

 
end
%% results
disp('====================================================')
%disp([' BEST solution   =  '  num2str(gparticle.pos)])
disp([' BEST fitness    = '   num2str(gparticle.cost)])
disp(['  Time           =  '  num2str(toc)])


% figure(1)
% plot(best(1:t),'r','LineWidth',2)
% hold on
% plot(AVR(1:t),'b','LineWidth',2)
% 
% xlabel('t')
% ylabel(' fitness ')
% 
% legend('BEST','MEAN')
% 
% title (' PSO ')

Answer=gparticle.pos;
Answer_cost=gparticle.cost;

%save('result','gparticle')
% datax=[5 10 15 30 45 60 120 180 240 300 480 720 960 1200 1440];
% datay=[53.9458525300000,66.0138248800000,68.1278801800000,79.5737327200000,86.3940092200000,105.028801800000,115.358294900000,123.215437800000,127.707373300000,129.505760400000,131.785714300000,131.372695900000,131.474654400000,132.520161300000,132.474078300000];
% p1=gparticle.pos(1);
% p2=gparticle.pos(2);
% semilogx(datax,datay,'og','MarkerSize',10,'MarkerFaceColor','y','LineWidth',2)
% hold on
% semilogx(datax,datay,'r')
% tt=5:0.01:1440;
% qt_new=(p1.*tt.*p2.^2)./(1+p2.*tt);
% semilogx(tt,qt_new)
