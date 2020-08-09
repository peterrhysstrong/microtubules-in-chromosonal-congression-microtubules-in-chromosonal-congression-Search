function param
%parmeters for kinetochore code
NC=46; % number of chromosones
Time=1200;%total time  in seconds
dt=0.01;%length of each timestep in seconds
t=Time/dt;%number of timesteps
Nmt=750;%number of mts
r=8;%radius of cell micro meters
L =[-5,0,0];%postion of left pole on x axis 
R =[5,0,0];%positin of right pole on x axis 
pd = abs(L(1)-R(1));%distance between poles
pr = 1; %pole radius in micrometers
NCGT=999; % number of times the second kinetochore will try to be placed 
critang=asin(pr/pd);%critical angle for interecesting poles for microtubles 
critang2=asin(2*pr/pd); % critial angle for intersecting poles for chromosones

%initial values for mts
mt.phi.L = (2*rand( Nmt, 1)).*pi;%phi angle for left pole
mt.phi.R = (2*rand( Nmt, 1)).*pi;%phi angle for right pole
mt.th.L = acos(-1+ 2*rand(Nmt, 1))- (pi/2) ;%theta angle for left pole
mt.th.R = acos(-1+ 2*rand(Nmt, 1))- (pi/2) ;%theta angle for right pole
mt.xyz.L=zeros(Nmt,3);
mt.dist.L=pr*ones(Nmt,1);
mt.dist.R=pr*ones(Nmt,1);
mt.state.L=ones(Nmt,1);%sets state of mts
mt.state.R=ones(Nmt,1);%""

%properties of the mt
Vpol = 12.8/60; %velocity of polymerization in microns per sec
Vdepol = 14.1/60;%"" depolymerization microns per sec
fcat = (0.058.*dt);%.*ones(Nmt,1);%probability of cat in dt
fres = 0.045*dt;%"" res in dt
resposR=find(mt.state.R==1);%"")
resposL=find(mt.state.L==1);%"")
catposR=[];
catposL=[];
c.catposL=zeros(Nmt,1);
c.catposR=zeros(Nmt,1);
dat=[];


%other values
D =0.001;%%diffusion coefficient
kr=0.15;%kt radius microns
csl=0.8; %chromosone spring length microns
sc=90;%*10^-9;%spring constant 
dc=1.8*10^3;%drag coeficcent 
BC= 1.38064852*10^-11;% boltzmans constant 
a=1687.5;%168*10^-9;%polar ejection force coefficient 
ktbr=0.2+0.15; % kt binding radius+kt radius microns
pefco=168.75; %polar ejection force coefficent 
msc=900;%motor spring coefficent
vmaxde=0.5;% maximum speed for laterally attached detrosinated mt
vmaxtr=0.45;%maximum speed for laterally attached trysoniated mt
fstall=40; % the force needed to cause a stall
dmt=24*10^-3;%diameter of mt

%forces for spring breaking equation 
a1=0.1;%break time 1
a2=10;%break time 2
b1=5;%force for break time 1
b2=70;%force for break time 2
fperb=(b1-b2)/log(a1/a2);%force for perpendicular breaking
r0=a1*exp(-b1/fperb);%breaking of lateral sliding spring

%forces for end on spring breaking
a1=0.1;%break time 1 s^-1
a2=10;%break time 2 s^-1
b1=5;%force for break time 1 pN
b2=100;%force for break time 2 pN
fendonb=(b1-b2)/log(a1/a2);%force for end on breaking
r0eo=a1*exp(-b1/fendonb);%breaking of lateral sliding spring

%forces for end on cat and rescue 
c1=100;%cat rate 1 s^-1
c2=0.015;%cat rate 2 s^-1
d1=0;%force for cat rate 1
d2=5;% force for cat rate 2
eocatratecon=(d1-d2)/log(c1/c2);%force for perpendicular breaking
r0eocatrate=c1*exp(-d1/eocatratecon);%breaking of lateral sliding spring

c1=0.015;%cat rate 1 s^-1
c2=100;%cat rate 2 s^-1
d1=0;%force for cat rate 1
d2=10;% force for cat rate 2
eoresratecon=(d1-d2)/log(c1/c2);%force for perpendicular breaking
r0eoresrate=c1*exp(-d1/eoresratecon);%breaking of lateral sliding spring
vc=190;%vicosity coefficent 
epsi=10^-5;
eodepspeed=0.3;%end on depolyersing speed
eopolspeed=0.5;%end on poltersing speed
save ('param' ,'-v7.3')
end                                                                                                             