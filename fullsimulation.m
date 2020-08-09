function [timeforpast,chrom1]=fullsimulation(xpos,ypos)
param;
load param; % loads in parameters  
%disp(['running kinetochore with ' num2str(t) ' steps, on ' num2str(NC) ' chromosomes'])
cpos=zeros(NC,6); %POSITION OF KTS


NC=1;
%INITIAL CHROMOSONE POSITION GENERATION
%cpos(:,1:3)=randsphere(NC,3,r);%generates the first kinetochroes position randomly in a spher 
cpos(:,1)= -xpos-0.5.*rand(NC,1);
cpos(:,2)= ypos+1/2.*rand(NC,1); %rand(NC,1);
cpos(:,3)= 0;
check1 = find( ((cpos(:,1)-L(1)).^2 + (cpos(:,2)-L(2)).^2+ (cpos(:,3)-L(3)).^2<pr^2) | ((cpos(:,1)-R(1)).^2 + (cpos(:,2)-R(2)).^2+ (cpos(:,3)-R(3)).^2<pr^2)|(cpos(:,1).^2 + cpos(:,2).^2+ cpos(:,3).^2>r^2));%checks the chromosone is indide the sphere and outside of the spindle pole
if ~isempty(check1)
    while isempty(check1)==0
    %cpos(check1,1:3)=randsphere(length(check1),3,r); %if the position is not okay generate it again 
    cpos(check1,1:3)=2*r.*rand(length(check1),3)-r;
    cpos(:,1)= -xpos-0.5.*rand(NC,1);
    cpos(:,2)= ypos+1/2.*rand(NC,1); %rand(NC,1);
    check1 = find( ((cpos(:,1)-L(1)).^2 + (cpos(:,2)-L(2)).^2+ (cpos(:,3)-L(3)).^2<pr^2) | ((cpos(:,1)-R(1)).^2 + (cpos(:,2)-R(2)).^2+ (cpos(:,3)-R(3)).^2<pr^2)|(cpos(:,1).^2 + cpos(:,2).^2+ cpos(:,3).^2>r^2)); %check again
    end
end

polargen4 = 2*rand(NC,1)*pi; %generates the angles for the random vetor
polargen5 = acos(-1+ 2*rand(NC,1));
[cpos(:, 4),cpos(:, 5), cpos(:, 6)]  = sph2cart2(polargen4, polargen5,csl); %converts spherical coridnates to cart
cpos(:,4:6) = cpos(:,4:6) + cpos(:,1:3); %positions the kinetochore next to the other one
check2=find(((cpos(:,4)-L(1)).^2 + (cpos(:,5)-L(2)).^2+ (cpos(:,6)-L(3)).^2<pr^2) | ((cpos(:,4)-R(1)).^2 + (cpos(:,5)-R(2)).^2+ (cpos(:,6)-R(3)).^2<pr^2) | (cpos(:,4).^2 + cpos(:,5).^2+ cpos(:,6).^2>r^2)); %checks the kinethore is inside the sphere and outside of the spindle pole

if ~isempty(check2)
    while isempty(check2)==0 %for kinetochores that dont satifiy the above check it generates them again
    polargen4 = 2*rand(length(check2),1)*pi; 
    polargen5 = acos(2*rand(length(check2),1)-1);
    
    [cpos(check2, 4),cpos(check2, 5), cpos(check2, 6)]  = sph2cart2(polargen4, polargen5,csl); %converts left coridnates to cart
    
    cpos(check2,4:6) = cpos(check2,4:6) + cpos(check2,1:3);
    
    %check2=find(((cpos(check2,4)-L(1)).^2 + (cpos(check2,5)-L(2)).^2+ (cpos(check2,6)-L(3)).^2<pr^2) | ((cpos(check2,4)-R(1)).^2 + (cpos(check2,5)-R(2)).^2+ (cpos(check2,6)-R(3)).^2<pr^2) | (cpos(check2,4).^2 + cpos(check2,5).^2+ cpos(check2,6).^2>r^2))
    check2=find(((cpos(:,4)-L(1)).^2 + (cpos(:,5)-L(2)).^2+ (cpos(:,6)-L(3)).^2<pr^2) | ((cpos(:,4)-R(1)).^2 + (cpos(:,5)-R(2)).^2+ (cpos(:,6)-R(3)).^2<pr^2) | (cpos(:,4).^2 + cpos(:,5).^2+ cpos(:,6).^2>r^2));
    end
end

%sets aside space for lots of variables 
e1=zeros(NC,3);
e2=zeros(NC,3);
%pef.L=zeros(NC,6); 
kt.state=zeros(NC,2);
kt.state2=zeros(NC,2);
chrom1=zeros(t,6);
ktmtvecttime=zeros(t,6);
kt.mt=zeros(NC,2);
pef=zeros(NC,6);
ktmtvect=zeros(NC,6);
eocatrate=zeros(NC,2);
eoresrate=zeros(NC,2);
ktmtbl=zeros(NC,6);%END ON ATTACHMENT!
boundedktmtstate=zeros(NC,2);
mtktforcetime1=zeros(t,NC);
mtktforcetime2=zeros(t,NC);

e1R=zeros(NC,3);
e2R=zeros(NC,3);
%pef.L=zeros(NC,6); 
kt.stateR=zeros(NC,2);
kt.state2R=zeros(NC,2);
kt.mtR=zeros(NC,2);
pefR=zeros(NC,6);
ktmtvectR=zeros(NC,6);
eocatrateR=zeros(NC,2);
eoresrateR=zeros(NC,2);
ktmtblR=zeros(NC,6);%END ON ATTACHMENT!
boundedktmtstateR=zeros(NC,2);
cpos
%allows you to try just one pair of kts in a initial location you want
NC=1;
%cpos=[-7.9,0,0,-8,0,0;-6,1,0,-6,-1,0]
%cpos=[7.9,0,0,-8,0,0;6,1,0,-6,-1,0]
%RUNS SIMULATION IN TIME
timeforpast=t*dt;
for i = 1:t
    
  done=find(((cpos(:,1)+cpos(:,4))./2)>-3);
  if length(done)==NC
        timeforpast=i*dt;
      break
  else    
    ktmtbl;
    kt.state;
    pef;
    
%movements of mts from left spindle pole    
rst=rand(Nmt,1);
    %pol mts
         
         %resposL  
        if (~isempty(resposL)) 
            mt.dist.L(resposL) = mt.dist.L(resposL) + Vpol*dt; %makes pol mts longer
             c.catposL = find(rst(resposL) <= fcat); % sees what mts undergo cat
         if (~isempty(c.catposL))
             mt.state.L(resposL(c.catposL))=-1;
         end
        end
        

     %depol mts   
     
        if (~isempty(catposL))
             mt.dist.L(catposL) = mt.dist.L(catposL) - Vdepol*dt; %same for depol mts
             c.resposL = find(rst(catposL) <= fres);
            if (~isempty(c.resposL))  
             mt.state.L(catposL(c.resposL)) = 1;
            end
        end
   
        out = find(mt.dist.L>=r-abs(L(1))); %sees what mts are outside of the cell
             if (~isempty(out))
                 z=length(out);
                 y=linspace(1,z,z);
                 [mt.xyz.L(out, 1),mt.xyz.L(out, 2), mt.xyz.L(out, 3)]  = sph2cart(mt.phi.L(out), mt.th.L(out), mt.dist.L(out)); %converts left coridnates to cart
                outsidedist(y) = sqrt(sum(abs(mt.xyz.L(out,:)+L).^2,2));
                outside = find(outsidedist(y)>=r);
                if (~isempty(outside))
                    mt.state.L(out(outside)) = -1;
                end
             end

         %if mt is inside the other spindel pole
          otherpole = find(((mt.dist.L>=pd-pr) & ((0<mt.phi.L)& mt.phi.L<critang) |(2*pi-critang<mt.phi.L & mt.phi.L<2*pi) & ((pi/2-critang< mt.th.L & mt.th.L <pi/2+critang))));%can try other restrictions to see if more time effective add and sttatments
          if (~isempty(otherpole))
                [mt.xyz.L(otherpole, 1),mt.xyz.L(otherpole, 2), mt.xyz.L(otherpole, 3)]  = sph2cart(mt.phi.L(otherpole), mt.th.L(otherpole), mt.dist.L(otherpole)); %converts left coridnates to cart
                mt.otherpoledist.L(otherpole) = sqrt(sum(abs(mt.xyz.L(otherpole,:)- R + L).^2,2));
                other = find( mt.otherpoledist.L(otherpole)<=pr);
                if (~isempty(other))
                    otherpole(other);
                    mt.state.L(otherpole(other)) = -1;
                end
          end
         
          % if mt is inside its spindel pole
          mt.new.L = find(mt.dist.L <= pr);
          if(~isempty(mt.new.L))
            mt.phi.L(mt.new.L) = (2*rand*pi);%phi angle for new mt
            mt.th.L(mt.new.L) = acos(-1+ 2*rand);%theta angle for new mt
            mt.dist.L(mt.new.L) = pr;
            mt.state.L(mt.new.L) = 1;
          end
% updates the new states
  resposL = find(mt.state.L == 1);%1
  catposL = find(mt.state.L == -1);
  
[mt.xyz.L(:, 1),mt.xyz.L(:, 2), mt.xyz.L(:, 3)]= sph2cart(mt.phi.L, mt.th.L, mt.dist.L); %gives mt location in cart
mtrxyz.L=mt.xyz.L+L; %mt xyz relative to sphere  
unitmt(:,1:3) = mtrxyz.L./((mtrxyz.L(:,1).^2 + mtrxyz.L(:,2).^2 + mtrxyz.L(:,3).^2).^0.5);%finds the unit vector of the mt from the origin

%movements of mts from right spindle pole    
rst=rand(Nmt,1);
    %pol mts
         
         %resposL  
        if (~isempty(resposR)) 
            mt.dist.R(resposR) = mt.dist.R(resposR) + Vpol*dt; %makes pol mts longer
             c.catposR = find(rst(resposR) <= fcat); % sees what mts undergo cat
         if (~isempty(c.catposR))
             mt.state.R(resposR(c.catposR))=-1;
         end
        end
        

     %depol mts   
     
        if (~isempty(catposR))
             mt.dist.R(catposR) = mt.dist.R(catposR) - Vdepol*dt; %same for depol mts
             c.resposR = find(rst(catposR) <= fres);
            if (~isempty(c.resposR))  
             mt.state.R(catposR(c.resposR)) = 1;
            end
        end
   
        out = find(mt.dist.R>=r-abs(R(1))); %sees what mts are outside of the cell
             if (~isempty(out))
                 z=length(out);
                 y=linspace(1,z,z);
                 [mt.xyz.R(out, 1),mt.xyz.R(out, 2), mt.xyz.R(out, 3)]  = sph2cart(mt.phi.R(out), mt.th.R(out), mt.dist.R(out)); %converts left coridnates to cart
                outsidedist(y) = sqrt(sum(abs(mt.xyz.R(out,:)+R).^2,2));
                outside = find(outsidedist(y)>=r);
                if (~isempty(outside))
                    mt.state.R(out(outside)) = -1;
                end
             end

         %if mt is inside the other spindel pole
          otherpole = find(((mt.dist.R>=pd-pr) & ((pi-critang<mt.phi.R)& mt.phi.R<pi+critang)  & ((pi/2-critang< mt.th.R & mt.th.R <pi/2+critang))));%can try other restrictions to see if more time effective add and sttatments
          if (~isempty(otherpole))
                [mt.xyz.R(otherpole, 1),mt.xyz.R(otherpole, 2), mt.xyz.R(otherpole, 3)]  = sph2cart(mt.phi.R(otherpole), mt.th.R(otherpole), mt.dist.R(otherpole)); %converts left coridnates to cart
                mt.otherpoledist.R(otherpole) = sqrt(sum(abs(mt.xyz.R(otherpole,:)- L + R).^2,2));
                other = find( mt.otherpoledist.R(otherpole)<=pr);
                if (~isempty(other))
                    otherpole(other);
                    mt.state.R(otherpole(other)) = -1;
                end
          end
         
          % if mt is inside its spindel pole
          mt.new.R = find(mt.dist.R <= pr);
          if(~isempty(mt.new.R))
            mt.phi.R(mt.new.R) = (2*rand*pi);%phi angle for new mt
            mt.th.R(mt.new.R) = acos(-1+ 2*rand);%theta angle for new mt
            mt.dist.R(mt.new.R) = pr;
            mt.state.R(mt.new.R) = 1;
          end
% updates the new states
  resposR = find(mt.state.R == 1);%1
  catposR = find(mt.state.R == -1);
  
[mt.xyz.R(:, 1),mt.xyz.R(:, 2), mt.xyz.R(:, 3)]= sph2cart(mt.phi.R, mt.th.R, mt.dist.R); %gives mt location in cart
mtrxyz.R=mt.xyz.R+R; %mt xyz relative to sphere  
unitmtR(:,1:3) = mtrxyz.R./((mtrxyz.R(:,1).^2 + mtrxyz.R(:,2).^2 + mtrxyz.R(:,3).^2).^0.5);%finds the unit vector of the mt from the origin




e=(cpos(:,1:3)-cpos(:,4:6))./((cpos(:,1)-cpos(:,4)).^2+(cpos(:,2)-cpos(:,5)).^2+(cpos(:,3)-cpos(:,6)).^2).^0.5; %finds the vector between the two kts

%Calates PEF
pefv(:,1:3)=(cpos(:,1:3)-L);%./(sqrt(sum((cpos(:,1:3)-L).^2))); FINDS VECTOR FROM KTS TO THE SPINDLE POLE
pefv(:,4:6)=(cpos(:,4:6)-L);%./(sqrt(sum((cpos(:,4:6)-L).^2)));
peflength(:,1)=(pefv(:,1).^2 + pefv(:,2).^2 + pefv(:,3).^2).^0.5; % FIND THE LENGTH OF THIS VECTOR
peflength(:,2)=(pefv(:,4).^2 + pefv(:,5).^2 + pefv(:,6).^2).^0.5;
pefnv(:,1:3)=pefv(:,1:3)./peflength(:,1); %CREATES A UNIT PEF VECTOR
pefnv(:,4:6)=pefv(:,4:6)./peflength(:,2);
pefcostheta(:,1)=abs((e(1).*pefnv(:,1)+e(2).*pefnv(:,2)+e(3).*pefnv(:,3)));%./(peflength(:,1))); FINDS THE ANGLE THAT THIS MAKES WITH THE ORITATION OF THE KTS
pefcostheta(:,2)=abs((e(1).*pefnv(:,4)+e(2).*pefnv(:,5)+e(3).*pefnv(:,6)));%./(peflength(:,2)));
pef(:,1:3)= (pefnv(:,1:3).*abs(pefcostheta(:,1)).*pefco)./(peflength(:,1)).^2;%3; WORKS OUT THE PEF FORCE
pef(:,4:6)= (pefnv(:,4:6).*abs(pefcostheta(:,2)).*pefco)./(peflength(:,2)).^2;%3;

%Calates PEF
pefvR(:,1:3)=(cpos(:,1:3)-R);%./(sqrt(sum((cpos(:,1:3)-L).^2))); FINDS VECTOR FROM KTS TO THE SPINDLE POLE
pefvR(:,4:6)=(cpos(:,4:6)-R);%./(sqrt(sum((cpos(:,4:6)-L).^2)));
peflengthR(:,1)=(pefvR(:,1).^2 + pefvR(:,2).^2 + pefvR(:,3).^2).^0.5; % FIND THE LENGTH OF THIS VECTOR
peflengthR(:,2)=(pefvR(:,4).^2 + pefvR(:,5).^2 + pefvR(:,6).^2).^0.5;
pefnvR(:,1:3)=pefvR(:,1:3)./peflengthR(:,1); %CREATES A UNIT PEF VECTOR
pefnvR(:,4:6)=pefvR(:,4:6)./peflengthR(:,2);
pefcosthetaR(:,1)=abs((e(1).*pefnvR(:,1)+e(2).*pefnvR(:,2)+e(3).*pefnvR(:,3)));%./(peflength(:,1))); FINDS THE ANGLE THAT THIS MAKES WITH THE ORITATION OF THE KTS
pefcosthetaR(:,2)=abs((e(1).*pefnvR(:,4)+e(2).*pefnvR(:,5)+e(3).*pefnvR(:,6)));%./(peflength(:,2)));
pefR(:,1:3)= (pefnvR(:,1:3).*abs(pefcosthetaR(:,1)).*pefco)./(peflengthR(:,1)).^2;%3; WORKS OUT THE PEF FORCE
pefR(:,4:6)= (pefnvR(:,4:6).*abs(pefcosthetaR(:,2)).*pefco)./(peflengthR(:,2)).^2;%3;

disL = cpos;
disL(:,1:3) = disL(:,1:3)-L;
disL(:,4:6) = disL(:,4:6)-L;
[poldisL(:,1),poldisL(:,2), poldisL(:,3)]  = cart2sph(disL(:,1), disL(:,2),disL(:,3));
[poldisL(:,4),poldisL(:,5), poldisL(:,6)]  = cart2sph(disL(:,4), disL(:,5),disL(:,6));

disR = cpos;
disR(:,1:3) = disR(:,1:3)-L;
disR(:,4:6) = disR(:,4:6)-L;
[poldisR(:,1),poldisR(:,2), poldisR(:,3)]  = cart2sph(disR(:,1), disR(:,2),disR(:,3));
[poldisR(:,4),poldisR(:,5), poldisR(:,6)]  = cart2sph(disR(:,4), disR(:,5),disR(:,6));

for w = 1:NC % go through all the kts
    mtktconnection=kt.mt(w,1).*(kt.state(w,1)+kt.state2(w,1));%find the connections
    if mtktconnection~=0
    howmany=find(kt.mt.*(kt.state+kt.state2) == mtktconnection);%find all the connection to the mt
    if length(howmany)==1 %find how many kts are connected to the mt
    rw=((ktmtbl(w,1)-L(1)).^2+(ktmtbl(w,2)).^2+(ktmtbl(w,3)).^2).^0.5; %find the distance to the binding region
    forcew=dt.*msc.*(cpos(w,1:3)-ktmtbl(w,1:3)-ktbr*e1(w,:));%finds the force acting on the binding region
    vdrag=8.*pi.*vc.*mt.dist.L(mtktconnection)./(2.*log(mt.dist.L(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.L(mtktconnection).^3./(3.*log(mt.dist.L(mtktconnection)./dmt)));%calulates rotational drag
    angvel=cross(ktmtvect(w,1:3),rw.*forcew)./(tdrag+vdrag*mt.dist.L(mtktconnection).^2);%calculates the angular velcoity
    ktmtbl(w,1:3)=cross(angvel,(ktmtbl(w,1:3)-L)).*dt+ktmtbl(w,1:3);%find the new location of the binding region
    ktmtbl(w,1:3);
    mt.xyz.L(mtktconnection,:)=cross(angvel,mt.xyz.L(mtktconnection,:)).*dt+mt.xyz.L(mtktconnection,:);
    else
    if howmany(1)==w
        forcesum=zeros(1,3);
        for x=1:length(howmany)
            if howmany(x)<=NC
                rw=((ktmtbl(howmany(x),1)-L(1)).^2+(ktmtbl(howmany(x),2)).^2+(ktmtbl(howmany(x),3)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x),1:3)-ktmtbl(howmany(x),1:3)-ktbr*e1(howmany(x),:));
                forcesum=rw.*forcew;
            else
                rw=((ktmtbl(howmany(x)-NC,4)-L(1)).^2+(ktmtbl(howmany(x)-NC,5)).^2+(ktmtbl(howmany(x)-NC,6)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x)-NC,4:6)-ktmtbl(howmany(x)-NC,4:6)-ktbr*e2(howmany(x)-NC,:));   
                forcesum=rw.*forcew;
            end
        end
        vdrag=8.*pi.*vc.*mt.dist.L(mtktconnection)./(2.*log(mt.dist.L(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.L(mtktconnection).^3./(3.*log(mt.dist.L(mtktconnection)./dmt)));%calulates rotational drag
        angvel=cross(ktmtvect(w,1:3),forcesum)./(tdrag+vdrag*mt.dist.L(mtktconnection).^2);%calculates the angular velcoity
        mt.xyz.L(mtktconnection,:)=cross(angvel,mt.xyz.L(mtktconnection,:)).*dt+mt.xyz.L(mtktconnection,:);%calulates the movement of the end of the mt
        for x=1:length(howmany)
            if howmany(x)<=NC
               ktmtbl(howmany(x),1:3)=cross(angvel,(ktmtbl(howmany(x),1:3)-L)).*dt+ktmtbl(howmany(x),1:3);%find the new location of the binding region
            else
               ktmtbl(howmany(x)-NC,4:6)=cross(angvel,(ktmtbl(howmany(x)-NC,4:6)-L)).*dt+ktmtbl(howmany(x)-NC,4:6);
            end
        end
    end
    
    end
    end
    mtktconnection=kt.mt(w,2).*(kt.state(w,2)+kt.state2(w,2));%find the connections
    if mtktconnection~=0
    howmany=find(kt.mt.*(kt.state+kt.state2) == mtktconnection);%find all the connection to the mt
    if length(howmany)==1 %find how many kts are connected to the mt
    rw=((ktmtbl(w,4)-L(1)).^2+(ktmtbl(w,5)).^2+(ktmtbl(w,6)).^2).^0.5; %find the distance to the binding region
    forcew=dt.*msc.*(cpos(w,4:6)-ktmtbl(w,4:6)-ktbr*e2(w,:));%finds the force acting on the binding region
    vdrag=8.*pi.*vc.*mt.dist.L(mtktconnection)./(2.*log(mt.dist.L(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.L(mtktconnection).^3./(3.*log(mt.dist.L(mtktconnection)./dmt)));%calulates rotational drag
    angvel=cross(ktmtvect(w,4:6),rw.*forcew)./(mt.dist.L(mtktconnection).^2);(tdrag+vdrag*mt.dist.L(mtktconnection).^2);%calculates the angular velcoity
    ktmtbl(w,4:6)=cross(angvel,(ktmtbl(w,4:6)-L)).*dt+ktmtbl(w,4:6);%find the new location of the binding region
    mt.xyz.L(mtktconnection,:)=cross(angvel,mt.xyz.L(mtktconnection,:)).*dt+mt.xyz.L(mtktconnection,:);
    else
    if howmany(1)==w+NC
        forcesum=zeros(1,3);
        for x=1:length(howmany)
                rw=((ktmtbl(howmany(x)-NC,4)-L(1)).^2+(ktmtbl(howmany(x)-NC,5)).^2+(ktmtbl(howmany(x)-NC,6)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x)-NC,4:6)-ktmtbl(howmany(x)-NC,4:6)-ktbr*e2(howmany(x)-NC,:));   
                forcesum=rw.*forcew;
        end
    vdrag=8.*pi.*vc.*mt.dist.L(mtktconnection)./(2.*log(mt.dist.L(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.L(mtktconnection).^3./(3.*log(mt.dist.L(mtktconnection)./dmt)));%calulates rotational drag
        angvel=cross(ktmtvect(w,4:6),forcesum)./(tdrag+vdrag*mt.dist.L(mtktconnection).^2);%calculates the angular velcoity
        mt.xyz.L(mtktconnection,:)=cross(angvel,mt.xyz.L(mtktconnection,:)).*dt+mt.xyz.L(mtktconnection,:);%calulates the movement of the end of the mt
        for x=1:length(howmany)
               ktmtbl(howmany(x)-NC,4:6)=cross(angvel,(ktmtbl(howmany(x)-NC,4:6)-L)).*dt+ktmtbl(howmany(x)-NC,4:6);
        end
    end
    
    end
    
    mtktconnection=kt.mtR(w,1).*(kt.stateR(w,1)+kt.state2R(w,1));%find the connections
    if mtktconnection~=0
    howmany=find(kt.mtR.*(kt.stateR+kt.state2R) == mtktconnection);%find all the connection to the mt
    if length(howmany)==1 %find how many kts are connected to the mt
    rw=((ktmtblR(w,1)-L(1)).^2+(ktmtblR(w,2)).^2+(ktmtblR(w,3)).^2).^0.5; %find the distance to the binding region
    forcew=dt.*msc.*(cpos(w,1:3)-ktmtblR(w,1:3)-ktbr*e1R(w,:));%finds the force acting on the binding region
    vdrag=8.*pi.*vc.*mt.dist.R(mtktconnection)./(2.*log(mt.dist.R(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.R(mtktconnection).^3./(3.*log(mt.dist.R(mtktconnection)./dmt)));%calulates rotational drag
    angvel=cross(ktmtvectR(w,1:3),rw.*forcew)./(tdrag+vdrag*mt.dist.R(mtktconnection).^2);%calculates the angular velcoity
    ktmtblR(w,1:3)=cross(angvel,(ktmtblR(w,1:3)-R)).*dt+ktmtblR(w,1:3);%find the new location of the binding region
    ktmtblR(w,1:3);
    mt.xyz.R(mtktconnection,:)=cross(angvel,mt.xyz.R(mtktconnection,:)).*dt+mt.xyz.R(mtktconnection,:);
    else
    if howmany(1)==w
        forcesum=zeros(1,3);
        for x=1:length(howmany)
            if howmany(x)<=NC
                rw=((ktmtblR(howmany(x),1)-L(1)).^2+(ktmtblR(howmany(x),2)).^2+(ktmtblR(howmany(x),3)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x),1:3)-ktmtblR(howmany(x),1:3)-ktbr*e1R(howmany(x),:));
                forcesum=rw.*forcew;
            else
                rw=((ktmtblR(howmany(x)-NC,4)-R(1)).^2+(ktmtblR(howmany(x)-NC,5)).^2+(ktmtblR(howmany(x)-NC,6)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x)-NC,4:6)-ktmtblR(howmany(x)-NC,4:6)-ktbr*e2R(howmany(x)-NC,:));   
                forcesum=rw.*forcew;
            end
        end
    vdrag=8.*pi.*vc.*mt.dist.R(mtktconnection)./(2.*log(mt.dist.R(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.R(mtktconnection).^3./(3.*log(mt.dist.R(mtktconnection)./dmt)));%calulates rotational drag
        angvel=cross(ktmtvectR(w,1:3),forcesum)./(tdrag+vdrag*mt.dist.R(mtktconnection).^2);%calculates the angular velcoity
        mt.xyz.R(mtktconnection,:)=cross(angvel,mt.xyz.R(mtktconnection,:)).*dt+mt.xyz.R(mtktconnection,:);%calulates the movement of the end of the mt
        for x=1:length(howmany)
            if howmany(x)<=NC
               ktmtblR(howmany(x),1:3)=cross(angvel,(ktmtblR(howmany(x),1:3)-R)).*dt+ktmtblR(howmany(x),1:3);%find the new location of the binding region
            else
               ktmtblR(howmany(x)-NC,4:6)=cross(angvel,(ktmtblR(howmany(x)-NC,4:6)-R)).*dt+ktmtblR(howmany(x)-NC,4:6);
            end
        end
    end
    
    end
    end
    
    
    
    mtktconnection=kt.mtR(w,2).*(kt.stateR(w,2)+kt.state2R(w,2));%find the connections
    if mtktconnection~=0
    howmany=find(kt.mtR.*(kt.stateR+kt.state2R) == mtktconnection);%find all the connection to the mt
    if length(howmany)==1 %find how many kts are connected to the mt
    rw=((ktmtblR(w,4)-R(1)).^2+(ktmtblR(w,5)).^2+(ktmtblR(w,6)).^2).^0.5; %find the distance to the binding region
    forcew=dt.*msc.*(cpos(w,4:6)-ktmtblR(w,4:6)-ktbr*e2R(w,:));%finds the force acting on the binding region
    vdrag=8.*pi.*vc.*mt.dist.R(mtktconnection)./(2.*log(mt.dist.R(mtktconnection)./dmt)+1);%calculates per drag
    tdrag=(pi.*vc.*mt.dist.R(mtktconnection).^3./(3.*log(mt.dist.R(mtktconnection)./dmt)));%calulates rotational drag
    angvel=cross(ktmtvectR(w,4:6),rw.*forcew)./(mt.dist.R(mtktconnection).^2);(tdrag+vdrag*mt.dist.R(mtktconnection).^2);%calculates the angular velcoity
    ktmtbl(w,4:6)=cross(angvel,(ktmtblR(w,4:6)-R)).*dt+ktmtblR(w,4:6);%find the new location of the binding region
    mt.xyz.R(mtktconnection,:)=cross(angvel,mt.xyz.R(mtktconnection,:)).*dt+mt.xyz.R(mtktconnection,:);
    else
    if howmany(1)==w+NC
        forcesum=zeros(1,3);
        for x=1:length(howmany)
                rw=((ktmtblR(howmany(x)-NC,4)-R(1)).^2+(ktmtblR(howmany(x)-NC,5)).^2+(ktmtblR(howmany(x)-NC,6)).^2).^0.5;
                forcew=dt.*msc.*(cpos(howmany(x)-NC,4:6)-ktmtblR(howmany(x)-NC,4:6)-ktbr*e2R(howmany(x)-NC,:));   
                forcesum=rw.*forcew;
        end
        vdrag=8.*pi.*vc.*mt.dist.R(mtktconnection)./(2.*log(mt.dist.R(mtktconnection)./dmt)+1);%calculates per drag
        tdrag=(pi.*vc.*mt.dist.R(mtktconnection).^3./(3.*log(mt.dist.R(mtktconnection)./dmt)));%calulates rotational drag
        angvel=cross(ktmtvectR(w,4:6),forcesum)./(tdrag+vdrag*mt.dist.R(mtktconnection).^2);%calculates the angular velcoity
        mt.xyz.R(mtktconnection,:)=cross(angvel,mt.xyz.R(mtktconnection,:)).*dt+mt.xyz.R(mtktconnection,:);%calulates the movement of the end of the mt
        for x=1:length(howmany)
               ktmtblR(howmany(x)-NC,4:6)=cross(angvel,(ktmtblR(howmany(x)-NC,4:6)-R)).*dt+ktmtblR(howmany(x)-NC,4:6);
        end
    end
    
    end
    end
    
    end

bound=find(kt.state(:,1)==1);
if ~isempty(bound)
     endon=find( ((ktmtbl(bound,1)-L(1)).^2+(ktmtbl(bound,2)).^2+(ktmtbl(bound,3)).^2).^0.5>mt.dist.L(kt.mt(bound,1)));
    if ~isempty(endon)
        kt.state(bound(endon),1)=0;
        kt.state2((bound(endon)),1)=1;
    else
    ktmtspringforce=dt.*msc.*(cpos(bound,1:3)-ktmtbl(bound,1:3)-ktbr*e1(bound,:));%the force in the spring between the kinetochore and the microtuble 
    ktmtbrfpar=(ktmtvect(bound,1).*ktmtspringforce(:,1)+ktmtvect(bound,2).*ktmtspringforce(:,2)+ktmtvect(bound,3).*ktmtspringforce(:,3));% the component of the ktmt spring parellel to the mt
    %add in sping breaking here by considering perpendicular force
    notpar=find((e1(bound,1).*ktmtvect(bound,1)+e1(bound,2).*ktmtvect(bound,2)+e1(bound,3).*ktmtvect(bound,3)).^2<1-epsi);
    if(~isempty(notpar))
    m=e1(bound(notpar),1:3)-((ktmtvect(bound(notpar),1).*e1(bound(notpar),1)+ktmtvect(bound(notpar),2).*e1(bound(notpar),2)+ktmtvect(bound(notpar),3).*e1(bound(notpar),3))).*ktmtvect(bound(notpar),1:3); %calculates the vector of the perpendicular component 
    m=m./(m(:,1).^2+m(:,2).^2+m(:,3).^2); %makes m a unit vector
    ktmtfper=(m(:,1).*ktmtspringforce(notpar,1)+m(:,2).*ktmtspringforce(notpar,2)+m(:,3).*ktmtspringforce(notpar,3));% calulates perpendicular force
    %sktmtfper=ktmtfper(:,1).^2+ktmtfper(:,1).^2+ktmtfper(:,1).^2).^0.5
    %ktmtfper=(ktmtbrf.^2+ktmtbrfpar.^2).^0.5;
    breakrate=r0.*exp(ktmtfper./fperb);%calculates breaking rate
    breaking=find(rand(length(notpar),1) <= breakrate*dt); %calulates the springs that break
    
    if ~isempty(breaking)
        kt.state(bound(breaking),1)=0;
    else
    end
    detry=find(ktmtbl(bound,1)>L(1));
    if ~isempty(detry)
        speed=-(vmaxde./fstall).*(ktmtbrfpar(detry)-fstall);
        speed=min(speed,vmaxde);
        speed=max(speed,0);
        ktmtbl(bound(detry),1:3) = ktmtbl(bound(detry),1:3)+speed.*dt.*ktmtvect(bound(detry),1:3);
        %{
        check3= find((ktmtbl(bound(detry),1).^2+cktmtbl(bound(detry),2).^2+cktmtbl(bound(detry),3).^2 < r.^2)&((cktmtbl(bound(detry),1)-L(1)).^2+cktmtbl(bound(detry),2).^2+cktmtbl(bound(detry),3).^2 > pr.^2));
    if (~isempty(check3))       
        mt.state(bound(detry(check3)))=0
        ktmtbl(bound(detry(check3)),1:3)=cktmtbl(bound(detry(check3)),1:3);
    end    
        %}
        
        
        
    end
    tr=find(ktmtbl(bound,1)<=L(1));
    if ~isempty(tr)
        speed=-(vmaxtr./fstall).*(ktmtbrfpar(tr)-fstall);
        speed=min(speed,vmaxtr);
        speed=max(speed,0);
        ktmtbl(bound(tr),1:3) = ktmtbl(bound(tr),1:3)-speed.*dt.*ktmtvect(bound(tr),1:3);
        %{
        check3= find((cktmtbl(bound(tr),1).^2+cktmtbl(bound(tr),2).^2+cktmtbl(bound(tr),3).^2 < r.^2)&((cktmtbl(bound(tr),1)-L(1)).^2+cktmtbl(bound(tr),2).^2+cktmtbl(bound(tr),3).^2 > pr.^2));
    if (~isempty(check3))       
        ktmtbl(bound(tr(check3)),1:3)=cktmtbl(bound(tr(check3)),1:3);
    end    
    %}
    end 
    end
    end
end

%END ON ATTACHMENT!
bound2=find(kt.state2(:,1)==1);
if ~isempty(bound2)
    ktmtspringforce=dt.*msc.*(cpos(bound2,1:3)-ktmtbl(bound2,1:3)-ktbr*e1(bound2,:));%calculates the ktmt spring force
    dmktmtf=(e1(bound2,1).*ktmtspringforce(:,1)+e1(bound2,2).*ktmtspringforce(:,2)+e1(bound2,3).*ktmtspringforce(:,3));%calcualtes the directed maginaude of the ktmt spring force

    breakrate=r0eo.*exp(dmktmtf./fendonb)+1./(2*pi).^0.5.*exp(-(peflength(bound2,1).^2./2)); %find the breaking rate for end on mt attachments this is dependent on the force as well as the distance from the spindle pole; %calulates the breaking rate which dependednt on the force
    breaking=find( rand(length(bound2),1) <= breakrate*dt); %sees what springs break
    
    if ~isempty(breaking)
        kt.state2(bound2(breaking),1)=0; % if break removes the spring
        changetomt=kt.mt(bound2(breaking),1);
        [mt.phi.L(changetomt), mt.th.L(changetomt), mt.dist.L(changetomt)]  = cart2sph(ktmtbl(bound2(breaking),1), ktmtbl(bound2(breaking),2),ktmtbl(bound2(breaking),3));
    else
    
    eores=find(boundedktmtstate(bound2,1)==1); %finds the end on mt that are in res state
    if ~isempty(eores)
        ktmtbl(bound2(eores),1:3) = ktmtbl(bound2(eores),1:3)+eopolspeed.*dt.*ktmtvect(bound2(eores),1:3);%calculates the movement of the binding region based on its state
        mtrxyz.L(kt.mt(bound2(eores),1),1:3)=ktmtbl(bound2(eores),1:3)+L; %mt xyz relative to sphere for end on 
        mt.dist.L(kt.mt(bound2(eores),1))=(mtrxyz.L(kt.mt(bound2(eores),1),1).^2+mtrxyz.L(kt.mt(bound2(eores),1),1).^2+mtrxyz.L(kt.mt(bound2(eores),1),1).^2).^0.5;
        eocatrate(bound2(eores),1)=(r0eocatrate.*exp(dmktmtf(eores)./eocatratecon))./(30); %find the force dependednt rate of cat
        cat=find(eocatrate(bound2(eores),1)<rand(length(bound2(eores)),1));%sees what microtubles undergo cat
        if ~isempty(cat)
            boundedktmtstate(bound2(eores(cat)),1)=-1;% changes the state
        end
    end
     eocat=find(boundedktmtstate(bound2,1)==-1);%finds end on mt that are in cat state
    if ~isempty(eocat)
        ktmtbl(bound2(eocat),1:3) = ktmtbl(bound2(eocat),1:3)+eodepspeed.*dt.*ktmtvect(bound2(eocat),1:3);%calulates the movement of the binding region based on state
                mtrxyz.L(kt.mt(bound2(eocat),1),1:3)=ktmtbl(bound2(eocat),1:3)+L; %mt xyz relative to sphere for end on 
        mt.dist.L(kt.mt(bound2(eocat),1))=(mtrxyz.L(kt.mt(bound2(eocat),1),1).^2+mtrxyz.L(kt.mt(bound2(eocat),1),1).^2+mtrxyz.L(kt.mt(bound2(eocat),1),1).^2).^0.5;
        
    eoresrate(bound2(eocat),1)=(r0eoresrate.*exp(dmktmtf(eocat)./eoresratecon)./(30)); %find the force dependent rate of res
        res=find(eoresrate(bound2(eocat),1)<rand(length(bound2(eocat)),1));%sees what mt undergo res
        if ~isempty(res)
            boundedktmtstate(bound2(eocat(res)),1)=1;%changes state 
        end
    end
    end
end



bound=find(kt.state(:,2)==1);
if ~isempty(bound)
     endon=find( ((ktmtbl(bound,4)-L(1)).^2+(ktmtbl(bound,5)).^2+(ktmtbl(bound,6)).^2).^0.5>mt.dist.L(kt.mt(bound,2)));
    if ~isempty(endon)
        kt.state(bound(endon),2)=0;
        kt.state2((bound(endon)),2)=1;
    else
    ktmtspringforce=dt.*msc.*(cpos(bound,4:6)-ktmtbl(bound,4:6)-ktbr.*e2(bound,:));%the force in the spring between the kinetochore and the microtuble 
    ktmtbrfpar=(ktmtvect(bound,4).*ktmtspringforce(:,1)+ktmtvect(bound,5).*ktmtspringforce(:,2)+ktmtvect(bound,6).*ktmtspringforce(:,3));% the component of the ktmt spring parellel to the mt
    %add in sping breaking here by considering perpendicular force
    notpar=find((e2(bound,1).*ktmtvect(bound,4)+e2(bound,2).*ktmtvect(bound,5)+e2(bound,3).*ktmtvect(bound,6)).^2<1-epsi);
    if(~isempty(notpar))
    m=e2(bound(notpar),1:3)-((ktmtvect(bound(notpar),4).*e2(bound(notpar),1)+ktmtvect(bound(notpar),5).*e2(bound(notpar),2)+ktmtvect(bound(notpar),6).*e2(bound(notpar),3))).*ktmtvect(bound(notpar),4:6); %calculates the vector of the perpendicular component 
    m=m./(m(:,1).^2+m(:,2).^2+m(:,3).^2); %makes m a unit vector
    ktmtfper=(m(:,1).*ktmtspringforce(notpar,1)+m(:,2).*ktmtspringforce(notpar,2)+m(:,3).*ktmtspringforce(notpar,3));% calulates perpendicular force
    %sktmtfper=ktmtfper(:,1).^2+ktmtfper(:,1).^2+ktmtfper(:,1).^2).^0.5
    %ktmtfper=(ktmtbrf.^2+ktmtbrfpar.^2).^0.5;
    breakrate=r0.*exp(ktmtfper./fperb);%calculates breaking rate
    breaking=find(rand(length(notpar),1) <= breakrate*dt); %calulates the springs that break    

    if ~isempty(breaking)
        kt.state(bound(breaking),2)=0;
    else
    end
    detry=find(ktmtbl(bound,4)>L(1));
    if ~isempty(detry)
        speed=-(vmaxde./fstall).*(ktmtbrfpar(detry)-fstall);
        speed=min(speed,vmaxde);
        speed=max(speed,0);
        ktmtbl(bound(detry),4:6) = ktmtbl(bound(detry),4:6)+speed.*dt.*ktmtvect(bound(detry),4:6);
    end
    tr=find(ktmtbl(bound,4)<=L(1));
    if ~isempty(tr)
        speed=-(vmaxtr./fstall).*(ktmtbrfpar(tr)-fstall);
        speed=min(speed,vmaxtr);
        speed=max(speed,0);
        ktmtbl(bound(tr),4:6) = ktmtbl(bound(tr),4:6)-speed.*dt.*ktmtvect(bound(tr),4:6);
    end 
    end
    end
end

%END ON ATTACHMENT!
bound2=find(kt.state2(:,2)==1);
if ~isempty(bound2)
    ktmtspringforce=dt.*msc.*(cpos(bound2,4:6)-ktmtbl(bound2,4:6)-ktbr*e2(bound2,:));%calculates the ktmt spring force
    dmktmtf=(e2(bound2,1).*ktmtspringforce(:,1)+e2(bound2,2).*ktmtspringforce(:,2)+e2(bound2,3).*ktmtspringforce(:,3));%calcualtes the directed maginaude of the ktmt spring force
    breakrate=r0eo.*exp(dmktmtf./fendonb)+1./(2*pi).^0.5.*exp(-(peflength(bound2,2).^2./2)); %find the breaking rate for end on mt attachments this is dependent on the force as well as the distance from the spindle pole 
    breaking=find( rand(length(bound2),1) <= breakrate*dt);
    
    if ~isempty(breaking)
        kt.state2(bound2(breaking),2)=0;
        changetomt=kt.mt(bound2(breaking),2);
        [mt.phi.L(changetomt), mt.th.L(changetomt), mt.dist.L(changetomt)]  = cart2sph(ktmtbl(bound2(breaking),4), ktmtbl(bound2(breaking),5),ktmtbl(bound2(breaking),6));
    else
    
    eores=find(boundedktmtstate(bound2,2)==1);
    if ~isempty(eores)
        speed=eopolspeed;
        ktmtbl(bound2(eores),4:6) = ktmtbl(bound2(eores),4:6)+speed.*dt.*ktmtvect(bound2(eores),4:6);
        mtrxyz.L(kt.mt(bound2(eores),2),1:3)=ktmtbl(bound2(eores),4:6)+L; %mt xyz relative to sphere for end on 
        mt.dist.L(kt.mt(bound2(eores),2))=(mtrxyz.L(kt.mt(bound2(eores),2),1).^2+mtrxyz.L(kt.mt(bound2(eores),2),2).^2+mtrxyz.L(kt.mt(bound2(eores),2),3).^2).^0.5;
        eocatrate(bound2(eores),2)=(r0eocatrate.*exp(dmktmtf(eores)./eocatratecon))./(30);
        cat=find(eocatrate(bound2(eores),2)<rand(length(bound2(eores)),1));
        if ~isempty(cat)
            boundedktmtstate(bound2(eores(cat)),2)=-1;
        end
    end
     eocat=find(boundedktmtstate(bound2,2)==-1);
    if ~isempty(eocat)
        speed=eodepspeed;
        ktmtbl(bound2(eocat),4:6) = ktmtbl(bound2(eocat),4:6)+speed.*dt.*ktmtvect(bound2(eocat),4:6);
        mtrxyz.L(kt.mt(bound2(eocat),2),1:3)=ktmtbl(bound2(eocat),4:6)+L; %mt xyz relative to sphere for end on 
        mt.dist.L(kt.mt(bound2(eocat),2))=(mtrxyz.L(kt.mt(bound2(eocat),2),1).^2+mtrxyz.L(kt.mt(bound2(eocat),2),2).^2+mtrxyz.L(kt.mt(bound2(eocat),2),3).^2).^0.5;
        eoresrate(bound2(eocat),2)=(r0eoresrate.*exp(dmktmtf(eocat)./eoresratecon)./(30));
        res=find(eoresrate(bound2(eocat),2)<rand(length(bound2(eocat)),1));
        if ~isempty(res)
            boundedktmtstate(bound2(eocat(res)),2)=1;
        end
    end
    end
end






bound=find(kt.stateR(:,1)==1);
if ~isempty(bound)
     endon=find( ((ktmtblR(bound,1)-L(1)).^2+(ktmtblR(bound,2)).^2+(ktmtblR(bound,3)).^2).^0.5>mt.dist.R(kt.mtR(bound,1)));
    if ~isempty(endon)
        kt.stateR(bound(endon),1)=0;
        kt.state2R((bound(endon)),1)=1;
    else
    ktmtspringforce=dt.*msc.*(cpos(bound,1:3)-ktmtblR(bound,1:3)-ktbr*e1R(bound,:));%the force in the spring between the kinetochore and the microtuble 
    ktmtbrfpar=(ktmtvectR(bound,1).*ktmtspringforce(:,1)+ktmtvectR(bound,2).*ktmtspringforce(:,2)+ktmtvectR(bound,3).*ktmtspringforce(:,3));% the component of the ktmt spring parellel to the mt
    %add in sping breaking here by considering perpendicular force
    notpar=find((e1R(bound,1).*ktmtvectR(bound,1)+e1R(bound,2).*ktmtvectR(bound,2)+e1R(bound,3).*ktmtvectR(bound,3)).^2<1-epsi);
    if(~isempty(notpar))
    m=e1R(bound(notpar),1:3)-((ktmtvectR(bound(notpar),1).*e1R(bound(notpar),1)+ktmtvectR(bound(notpar),2).*e1R(bound(notpar),2)+ktmtvectR(bound(notpar),3).*e1R(bound(notpar),3))).*ktmtvectR(bound(notpar),1:3); %calculates the vector of the perpendicular component 
    m=m./(m(:,1).^2+m(:,2).^2+m(:,3).^2); %makes m a unit vector
    ktmtfper=(m(:,1).*ktmtspringforce(notpar,1)+m(:,2).*ktmtspringforce(notpar,2)+m(:,3).*ktmtspringforce(notpar,3));% calulates perpendicular force
    %sktmtfper=ktmtfper(:,1).^2+ktmtfper(:,1).^2+ktmtfper(:,1).^2).^0.5
    %ktmtfper=(ktmtbrf.^2+ktmtbrfpar.^2).^0.5;
    breakrate=r0.*exp(ktmtfper./fperb);%calculates breaking rate
    breaking=find(rand(length(notpar),1) <= breakrate*dt); %calulates the springs that break
    
    if ~isempty(breaking)
        kt.stateR(bound(breaking),1)=0;
    else
    end
    detry=find(ktmtblR(bound,1)>L(1));
    if ~isempty(detry)
        speed=-(vmaxde./fstall).*(ktmtbrfpar(detry)-fstall);
        speed=min(speed,vmaxde);
        speed=max(speed,0);
        ktmtblR(bound(detry),1:3) = ktmtblR(bound(detry),1:3)+speed.*dt.*ktmtvectR(bound(detry),1:3);
        %{
        check3= find((ktmtbl(bound(detry),1).^2+cktmtbl(bound(detry),2).^2+cktmtbl(bound(detry),3).^2 < r.^2)&((cktmtbl(bound(detry),1)-L(1)).^2+cktmtbl(bound(detry),2).^2+cktmtbl(bound(detry),3).^2 > pr.^2));
    if (~isempty(check3))       
        mt.state(bound(detry(check3)))=0
        ktmtbl(bound(detry(check3)),1:3)=cktmtbl(bound(detry(check3)),1:3);
    end    
        %}
        
        
        
    end
    tr=find(ktmtblR(bound,1)<=L(1));
    if ~isempty(tr)
        speed=-(vmaxtr./fstall).*(ktmtbrfpar(tr)-fstall);
        speed=min(speed,vmaxtr);
        speed=max(speed,0);
        ktmtblR(bound(tr),1:3) = ktmtblR(bound(tr),1:3)-speed.*dt.*ktmtvectR(bound(tr),1:3);
        %{
        check3= find((cktmtbl(bound(tr),1).^2+cktmtbl(bound(tr),2).^2+cktmtbl(bound(tr),3).^2 < r.^2)&((cktmtbl(bound(tr),1)-L(1)).^2+cktmtbl(bound(tr),2).^2+cktmtbl(bound(tr),3).^2 > pr.^2));
    if (~isempty(check3))       
        ktmtbl(bound(tr(check3)),1:3)=cktmtbl(bound(tr(check3)),1:3);
    end    
    %}
    end 
    end
    end
end

%END ON ATTACHMENT!
bound2=find(kt.state2R(:,1)==1);
if ~isempty(bound2)
    ktmtspringforce=dt.*msc.*(cpos(bound2,1:3)-ktmtblR(bound2,1:3)-ktbr*e1R(bound2,:));%calculates the ktmt spring force
    dmktmtf=(e1R(bound2,1).*ktmtspringforce(:,1)+e1R(bound2,2).*ktmtspringforce(:,2)+e1R(bound2,3).*ktmtspringforce(:,3));%calcualtes the directed maginaude of the ktmt spring force

    breakrate=r0eo.*exp(dmktmtf./fendonb)+1./(2*pi).^0.5.*exp(-(peflengthR(bound2,1).^2./2)); %find the breaking rate for end on mt attachments this is dependent on the force as well as the distance from the spindle pole; %calulates the breaking rate which dependednt on the force
    breaking=find( rand(length(bound2),1) <= breakrate*dt); %sees what springs break
    
    if ~isempty(breaking)
        kt.state2R(bound2(breaking),1)=0; % if break removes the spring
        changetomt=kt.mtR(bound2(breaking),1);
        [mt.phi.R(changetomt), mt.th.R(changetomt), mt.dist.R(changetomt)]  = cart2sph(ktmtblR(bound2(breaking),1), ktmtblR(bound2(breaking),2),ktmtblR(bound2(breaking),3));
    else
    
    eores=find(boundedktmtstateR(bound2,1)==1); %finds the end on mt that are in res state
    if ~isempty(eores)
        ktmtblR(bound2(eores),1:3) = ktmtblR(bound2(eores),1:3)+eopolspeed.*dt.*ktmtvectR(bound2(eores),1:3);%calculates the movement of the binding region based on its state
        mtrxyz.R(kt.mtR(bound2(eores),1),1:3)=ktmtblR(bound2(eores),1:3)+R; %mt xyz relative to sphere for end on 
        mt.dist.R(kt.mtR(bound2(eores),1))=(mtrxyz.R(kt.mtR(bound2(eores),1),1).^2+mtrxyz.R(kt.mtR(bound2(eores),1),1).^2+mtrxyz.R(kt.mtR(bound2(eores),1),1).^2).^0.5;
        eocatrate(bound2(eores),1)=(r0eocatrate.*exp(dmktmtf(eores)./eocatratecon))./(30); %find the force dependednt rate of cat
        cat=find(eocatrate(bound2(eores),1)<rand(length(bound2(eores)),1));%sees what microtubles undergo cat
        if ~isempty(cat)
            boundedktmtstateR(bound2(eores(cat)),1)=-1;% changes the state
        end
    end
     eocat=find(boundedktmtstateR(bound2,1)==-1);%finds end on mt that are in cat state
    if ~isempty(eocat)
        ktmtblR(bound2(eocat),1:3) = ktmtblR(bound2(eocat),1:3)+eodepspeed.*dt.*ktmtvectR(bound2(eocat),1:3);%calulates the movement of the binding region based on state
                mtrxyz.R(kt.mtR(bound2(eocat),1),1:3)=ktmtblR(bound2(eocat),1:3)+R; %mt xyz relative to sphere for end on 
        mt.dist.R(kt.mtR(bound2(eocat),1))=(mtrxyz.R(kt.mtR(bound2(eocat),1),1).^2+mtrxyz.R(kt.mtR(bound2(eocat),1),1).^2+mtrxyz.R(kt.mtR(bound2(eocat),1),1).^2).^0.5;
        
    eoresrate(bound2(eocat),1)=(r0eoresrate.*exp(dmktmtf(eocat)./eoresratecon)./(30)); %find the force dependent rate of res
        res=find(eoresrate(bound2(eocat),1)<rand(length(bound2(eocat)),1));%sees what mt undergo res
        if ~isempty(res)
            boundedktmtstateR(bound2(eocat(res)),1)=1;%changes state 
        end
    end
    end
end



bound=find(kt.stateR(:,2)==1);
if ~isempty(bound)
     endon=find( ((ktmtblR(bound,4)-R(1)).^2+(ktmtblR(bound,5)).^2+(ktmtblR(bound,6)).^2).^0.5>mt.dist.R(kt.mtR(bound,2)));
    if ~isempty(endon)
        kt.stateR(bound(endon),2)=0;
        kt.state2R((bound(endon)),2)=1;
    else
    ktmtspringforce=dt.*msc.*(cpos(bound,4:6)-ktmtblR(bound,4:6)-ktbr.*e2R(bound,:));%the force in the spring between the kinetochore and the microtuble 
    ktmtbrfpar=(ktmtvectR(bound,4).*ktmtspringforce(:,1)+ktmtvectR(bound,5).*ktmtspringforce(:,2)+ktmtvectR(bound,6).*ktmtspringforce(:,3));% the component of the ktmt spring parellel to the mt
    %add in sping breaking here by considering perpendicular force
    notpar=find((e2R(bound,1).*ktmtvectR(bound,4)+e2R(bound,2).*ktmtvectR(bound,5)+e2R(bound,3).*ktmtvectR(bound,6)).^2<1-epsi);
    if(~isempty(notpar))
    m=e2R(bound(notpar),1:3)-((ktmtvectR(bound(notpar),4).*e2R(bound(notpar),1)+ktmtvectR(bound(notpar),5).*e2R(bound(notpar),2)+ktmtvectR(bound(notpar),6).*e2R(bound(notpar),3))).*ktmtvectR(bound(notpar),4:6); %calculates the vector of the perpendicular component 
    m=m./(m(:,1).^2+m(:,2).^2+m(:,3).^2); %makes m a unit vector
    ktmtfper=(m(:,1).*ktmtspringforce(notpar,1)+m(:,2).*ktmtspringforce(notpar,2)+m(:,3).*ktmtspringforce(notpar,3));% calulates perpendicular force
    %sktmtfper=ktmtfper(:,1).^2+ktmtfper(:,1).^2+ktmtfper(:,1).^2).^0.5
    %ktmtfper=(ktmtbrf.^2+ktmtbrfpar.^2).^0.5;
    breakrate=r0.*exp(ktmtfper./fperb);%calculates breaking rate
    breaking=find(rand(length(notpar),1) <= breakrate*dt); %calulates the springs that break    

    if ~isempty(breaking)
        kt.stateR(bound(breaking),2)=0;
    else
    end
    detry=find(ktmtblR(bound,4)>L(1));
    if ~isempty(detry)
        speed=-(vmaxde./fstall).*(ktmtbrfpar(detry)-fstall);
        speed=min(speed,vmaxde);
        speed=max(speed,0);
        ktmtblR(bound(detry),4:6) = ktmtblR(bound(detry),4:6)+speed.*dt.*ktmtvectR(bound(detry),4:6);
    end
    tr=find(ktmtblR(bound,4)<=L(1));
    if ~isempty(tr)
        speed=-(vmaxtr./fstall).*(ktmtbrfpar(tr)-fstall);
        speed=min(speed,vmaxtr);
        speed=max(speed,0);
        ktmtblR(bound(tr),4:6) = ktmtblR(bound(tr),4:6)-speed.*dt.*ktmtvectR(bound(tr),4:6);
    end 
    end
    end
end

%END ON ATTACHMENT!
bound2=find(kt.state2R(:,2)==1);
if ~isempty(bound2)
    ktmtspringforce=dt.*msc.*(cpos(bound2,4:6)-ktmtblR(bound2,4:6)-ktbr*e2R(bound2,:));%calculates the ktmt spring force
    dmktmtf=(e2R(bound2,1).*ktmtspringforce(:,1)+e2R(bound2,2).*ktmtspringforce(:,2)+e2R(bound2,3).*ktmtspringforce(:,3));%calcualtes the directed maginaude of the ktmt spring force
    breakrate=r0eo.*exp(dmktmtf./fendonb)+1./(2*pi).^0.5.*exp(-(peflengthR(bound2,2).^2./2)); %find the breaking rate for end on mt attachments this is dependent on the force as well as the distance from the spindle pole 
    breaking=find( rand(length(bound2),1) <= breakrate*dt);
    
    if ~isempty(breaking)
        kt.state2R(bound2(breaking),2)=0;
        changetomt=kt.mtR(bound2(breaking),2);
        [mt.phi.R(changetomt), mt.th.R(changetomt), mt.dist.R(changetomt)]  = cart2sph(ktmtblR(bound2(breaking),4), ktmtblR(bound2(breaking),5),ktmtblR(bound2(breaking),6));
    else
    
    eores=find(boundedktmtstateR(bound2,2)==1);
    if ~isempty(eores)
        speed=eopolspeed;
        ktmtblR(bound2(eores),4:6) = ktmtblR(bound2(eores),4:6)+speed.*dt.*ktmtvectR(bound2(eores),4:6);
        mtrxyz.R(kt.mtR(bound2(eores),2),1:3)=ktmtblR(bound2(eores),4:6)+L; %mt xyz relative to sphere for end on 
        mt.dist.R(kt.mtR(bound2(eores),2))=(mtrxyz.R(kt.mtR(bound2(eores),2),1).^2+mtrxyz.R(kt.mtR(bound2(eores),2),2).^2+mtrxyz.R(kt.mtR(bound2(eores),2),3).^2).^0.5;
        eocatrate(bound2(eores),2)=(r0eocatrate.*exp(dmktmtf(eores)./eocatratecon))./(30);
        cat=find(eocatrate(bound2(eores),2)<rand(length(bound2(eores)),1));
        if ~isempty(cat)
            boundedktmtstateR(bound2(eores(cat)),2)=-1;
        end
    end
     eocat=find(boundedktmtstateR(bound2,2)==-1);
    if ~isempty(eocat)
        speed=eodepspeed;
        ktmtblR(bound2(eocat),4:6) = ktmtblR(bound2(eocat),4:6)+speed.*dt.*ktmtvectR(bound2(eocat),4:6);
        mtrxyz.R(kt.mtR(bound2(eocat),2),1:3)=ktmtblR(bound2(eocat),4:6)+R; %mt xyz relative to sphere for end on 
        mt.dist.R(kt.mtR(bound2(eocat),2))=(mtrxyz.R(kt.mtR(bound2(eocat),2),1).^2+mtrxyz.R(kt.mtR(bound2(eocat),2),2).^2+mtrxyz.R(kt.mtR(bound2(eocat),2),3).^2).^0.5;
        eoresrate(bound2(eocat),2)=(r0eoresrate.*exp(dmktmtf(eocat)./eoresratecon)./(30));
        res=find(eoresrate(bound2(eocat),2)<rand(length(bound2(eocat)),1));
        if ~isempty(res)
            boundedktmtstateR(bound2(eocat(res)),2)=1;
        end
    end
    end
end

    %{
for j=1:NC 
    if (kt.state(j,1) == 0 & kt.state2(j,1)==0 )
     
    posbindL =find(((cpos(j,1)-mtrxyz.L(:,1)).^2+(cpos(j,2)-mtrxyz.L(:,2)).^2+(cpos(j,3)-mtrxyz.L(:,3)).^2)<ktbr^2);
    if ~isempty(posbindL)
        mtb = posbindL(1);
        kt.mt(j,1)=mtb;
        kt.state2(j,1)=1;
        ktmtbl(j,1:3) = mtrxyz.L(mtb,:); %
        boundedktmtstate(j,1)=mt.state.L(kt.mt(j,1));
        mtvect = (mtrxyz.L(mtb,:)-L); %vector of microtuble 
        ktmtvect(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    else
    posbindL = find( ((unitmt(:,1).*disL(j,1))+(unitmt(:,2).*disL(j,2))+(unitmt(:,3).*disL(j,3))).^2-dot(disL(j,1:3),disL(j,1:3))+ktbr.^2>=0 & mt.dist.L>poldisL(j,3));
   
    if ~isempty(posbindL)
        mtb = posbindL(1);
        kt.mt(j,1)=mtb;
        kt.state(j,1)=1;
        mtvect = (mtrxyz.L(mtb,:)- L); %vector of microtuble 
        ktmtvect(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        cpoint = dot(ktmtvect(j,1:3),cpos(j,1:3)-L).*ktmtvect(j,1:3); %closest point on the mt to the kts
        ktmtbl(j,1:3) = L + cpoint; %
        mtvect = (mtrxyz.L(mtb,:)-L); %vector of microtuble 
        ktmtvect(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    end
    end
    end
    
    if (kt.state(j,2)== 0) & (kt.state2(j,2)==0)
         posbindL =find(((cpos(j,4)-mtrxyz.L(:,1)).^2+(cpos(j,5)-mtrxyz.L(:,2)).^2+(cpos(j,6)-mtrxyz.L(:,3)).^2)<ktbr^2);
    if ~isempty(posbindL)
        mtb = posbindL(1);
        kt.mt(j,2)=mtb;
        kt.state2(j,2)=1;
        ktmtbl(j,4:6) = mtrxyz.L(mtb,:);
        boundedktmtstate(j,2)=mt.state.L(kt.mt(j,2));
        mtvect = (mtrxyz.L(mtb,:)-L); %vector of microtuble 
        ktmtvect(j,4:6) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    else
    posbindL1 = find( ((unitmt(:,1).*disL(j,4))+(unitmt(:,2).*disL(j,5))+(unitmt(:,3).*disL(j,6))).^2-dot(disL(j,4:6),disL(j,4:6))+ktbr.^2>=0 & mt.dist.L>poldisL(j,6));   
    if ~isempty(posbindL1)
        mtb = posbindL1(1);
        kt.mt(j,2)=mtb;
        kt.state(j,2)=1;
        mtvect = (mtrxyz.L(mtb,:)- L); %vector of microtuble 
        ktmtvect(j,4:6)=mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        cpoint =dot(ktmtvect(j,4:6),cpos(j,4:6)-L).*ktmtvect(j,4:6); %closest point on the mt to the kts
        mtvect = (mtrxyz.L(mtb,:)-L); %vector of microtuble 
        ktmtvect(j,4:6)=mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        ktmtbl(j,4:6) = L + cpoint; %     

    end
    end
    end
end    




for j=1:NC 
    if (kt.stateR(j,1) == 0 & kt.state2R(j,1)==0 )
     
    posbindR =find(((cpos(j,1)-mtrxyz.R(:,1)).^2+(cpos(j,2)-mtrxyz.R(:,2)).^2+(cpos(j,3)-mtrxyz.R(:,3)).^2)<ktbr^2);
    if ~isempty(posbindR)
        mtb = posbindR(1);
        kt.mtR(j,1)=mtb;
        kt.state2R(j,1)=1;
        ktmtblR(j,1:3) = mtrxyz.R(mtb,:); %
        boundedktmtstateR(j,1)=mt.state.R(kt.mtR(j,1));
        mtvect = (mtrxyz.R(mtb,:)-R); %vector of microtuble 
        ktmtvectR(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    else
    posbindR = find( ((unitmtR(:,1).*disR(j,1))+(unitmtR(:,2).*disR(j,2))+(unitmtR(:,3).*disR(j,3))).^2-dot(disR(j,1:3),disR(j,1:3))+ktbr.^2>=0 & mt.dist.R>poldisR(j,3));
   
    if ~isempty(posbindR)
        mtb = posbindR(1);
        kt.mtR(j,1)=mtb;
        kt.stateR(j,1)=1;
        mtvect = (mtrxyz.R(mtb,:)- R); %vector of microtuble 
        ktmtvectR(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        cpoint = dot(ktmtvectR(j,1:3),cpos(j,1:3)-R).*ktmtvectR(j,1:3); %closest point on the mt to the kts
        ktmtblR(j,1:3) = R + cpoint; %
        mtvect = (mtrxyz.R(mtb,:)-R); %vector of microtuble 
        ktmtvectR(j,1:3) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    end
    end
    end
    
    if (kt.stateR(j,2)== 0) & (kt.state2R(j,2)==0)
         posbindR =find(((cpos(j,4)-mtrxyz.R(:,1)).^2+(cpos(j,5)-mtrxyz.R(:,2)).^2+(cpos(j,6)-mtrxyz.R(:,3)).^2)<ktbr^2);
    if ~isempty(posbindR)
        mtb = posbindR(1);
        kt.mtR(j,2)=mtb;
        kt.state2R(j,2)=1;
        ktmtblR(j,4:6) = mtrxyz.R(mtb,:);
        boundedktmtstateR(j,2)=mt.state.R(kt.mtR(j,2));
        mtvect = (mtrxyz.R(mtb,:)-r); %vector of microtuble 
        ktmtvectR(j,4:6) = mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
    else
    posbindL1 = find( ((unitmtR(:,1).*disR(j,4))+(unitmtR(:,2).*disR(j,5))+(unitmtR(:,3).*disR(j,6))).^2-dot(disR(j,4:6),disR(j,4:6))+ktbr.^2>=0 & mt.dist.R>poldisR(j,6));   
    if ~isempty(posbindL1)
        mtb = posbindL1(1);
        kt.mtR(j,2)=mtb;
        kt.stateR(j,2)=1;
        mtvect = (mtrxyz.R(mtb,:)- R); %vector of microtuble 
        ktmtvectR(j,4:6)=mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        cpoint =dot(ktmtvectR(j,4:6),cpos(j,4:6)-R).*ktmtvectR(j,4:6); %closest point on the mt to the kts
        mtvect = (mtrxyz.R(mtb,:)-R); %vector of microtuble 
        ktmtvectR(j,4:6)=mtvect./((mtvect(1).^2+mtvect(2).^2+mtvect(3).^2).^0.5);%unit microtuble vector
        ktmtblR(j,4:6) = R + cpoint; %     

    end
    end
    end
end    


    lat1=find(kt.state(:,1)==1|kt.state2(:,1)==1);
if ~isempty(lat1)
    e1(lat1,1:3)=(cpos(lat1,1:3)-ktmtbl(lat1,1:3))./((cpos(lat1,1)-ktmtbl(lat1,1)).^2+(cpos(lat1,2)-ktmtbl(lat1,2)).^2+(cpos(lat1,3)-ktmtbl(lat1,3)).^2).^0.5;
end
    lat2=find(kt.state(:,2)==1|kt.state2(:,2)==1);
if ~isempty(lat2)
    e2(lat2,1:3)=(cpos(lat2,4:6)-ktmtbl(lat2,4:6))./((cpos(lat2,4)-ktmtbl(lat2,4)).^2+(cpos(lat2,5)-ktmtbl(lat2,5)).^2+(cpos(lat2,6)-ktmtbl(lat2,6)).^2).^0.5;
end

    lat1=find(kt.stateR(:,1)==1|kt.state2R(:,1)==1);
if ~isempty(lat1)
    e1R(lat1,1:3)=(cpos(lat1,1:3)-ktmtblR(lat1,1:3))./((cpos(lat1,1)-ktmtblR(lat1,1)).^2+(cpos(lat1,2)-ktmtblR(lat1,2)).^2+(cpos(lat1,3)-ktmtblR(lat1,3)).^2).^0.5;
end
    lat2=find(kt.stateR(:,2)==1|kt.state2R(:,2)==1);
if ~isempty(lat2)
    e2R(lat2,1:3)=(cpos(lat2,4:6)-ktmtblR(lat2,4:6))./((cpos(lat2,4)-ktmtblR(lat2,4)).^2+(cpos(lat2,5)-ktmtblR(lat2,5)).^2+(cpos(lat2,6)-ktmtblR(lat2,6)).^2).^0.5;
end

    check5= find((ktmtbl(:,1).^2+ktmtbl(:,2).^2+ktmtbl(:,3).^2 > r.^2)|((ktmtbl(:,1)-L(1)).^2+ktmtbl(:,2).^2+ktmtbl(:,3).^2 < pr.^2));
    if (~isempty(check5))       
        kt.state(check5,1)=0;
        kt.state2(check5,1)=0;
    end
    check5= find((ktmtblR(:,1).^2+ktmtblR(:,2).^2+ktmtblR(:,3).^2 > r.^2)|((ktmtblR(:,1)-R(1)).^2+ktmtblR(:,2).^2+ktmtblR(:,3).^2 < pr.^2));
    if (~isempty(check5))       
        kt.stateR(check5,1)=0;
        kt.state2R(check5,1)=0;
    end
%}    
ch1 = normrnd(0,sqrt(2*D*dt),NC,3)-(dt*(sc*(cpos(:,1:3)-cpos(:,4:6)-csl*e)+kt.state(:,1).*msc.*(cpos(:,1:3)-ktmtbl(:,1:3)-ktbr*e1)+kt.state2(:,1).*msc.*(cpos(:,1:3)-ktmtbl(:,1:3)-ktbr*e1)+kt.stateR(:,1).*msc.*(cpos(:,1:3)-ktmtblR(:,1:3)-ktbr*e1R)+kt.state2R(:,1).*msc.*(cpos(:,1:3)-ktmtblR(:,1:3)-ktbr*e1R) +pefR(:,1:3)+ pef(:,1:3)))./dc;
mtktforce1=dt.*msc.*kt.state(:,1).*(cpos(:,1:3)-ktmtbl(:,1:3)-ktbr.*e1(:,:));%calculates the ktmt spring force
mtktforcetime1(i,:)=(e1(:,1).*mtktforce1(:,1)+e1(:,2).*mtktforce1(:,2)+e1(:,3).*mtktforce1(:,3));

mtktforce2=dt.*msc.*kt.state2(:,1).*(cpos(:,1:3)-ktmtbl(:,1:3)-ktbr*e1(:,:));%calculates the ktmt spring force
mtktforcetime2(i,:)=(e1(:,1).*mtktforce2(:,1)+e1(:,2).*mtktforce2(:,2)+e1(:,3).*mtktforce2(:,3));

%checks ktmtbr is in the right area 
    check5= find((ktmtbl(:,4).^2+ktmtbl(:,5).^2+ktmtbl(:,6).^2 > r.^2)|((ktmtbl(:,4)-L(1)).^2+ktmtbl(:,5).^2+ktmtbl(:,6).^2 < pr.^2)|((ktmtbl(:,4)-R(1)).^2+ktmtbl(:,5).^2+ktmtbl(:,6).^2 < pr.^2));
    if (~isempty(check5))       
        kt.state(check5,2)=0;
        kt.state2(check5,2)=0;
    end
    
        check5= find((ktmtblR(:,4).^2+ktmtblR(:,5).^2+ktmtblR(:,6).^2 > r.^2)|((ktmtblR(:,4)-L(1)).^2+ktmtblR(:,5).^2+ktmtblR(:,6).^2 < pr.^2)|((ktmtblR(:,4)-R(1)).^2+ktmtblR(:,5).^2+ktmtblR(:,6).^2 < pr.^2));
    if (~isempty(check5))       
        kt.stateR(check5,2)=0;
        kt.state2R(check5,2)=0;
    end


ch2 = normrnd(0,sqrt(2*D*dt),NC,3)+(dt*(sc*(cpos(:,1:3)-cpos(:,4:6)-csl*e)+kt.state(:,2).*msc.*(cpos(:,4:6)-ktmtbl(:,4:6)-ktbr*e2)+kt.state2(:,2).*msc.*(cpos(:,4:6)-ktmtbl(:,4:6)-ktbr*e2)+kt.stateR(:,2).*msc.*(cpos(:,4:6)-ktmtblR(:,4:6)-ktbr*e2R)+kt.state2R(:,2).*msc.*(cpos(:,4:6)-ktmtblR(:,4:6)-ktbr*e2R) + pefR(:,4:6)+ pef(:,4:6)))./dc;%+ pef(:,4:6)

ccpos(:,1:3)=cpos(:,1:3) + ch1;
ccpos(:,4:6)=cpos(:,4:6) + ch2;

check3= find((ccpos(:,1).^2+ccpos(:,2).^2+ccpos(:,3).^2 < r.^2)&((ccpos(:,1)-L(1)).^2+ccpos(:,2).^2+ccpos(:,3).^2 > pr.^2)&((ccpos(:,1)-R(1)).^2+ccpos(:,2).^2+ccpos(:,3).^2 > pr.^2));
    if (~isempty(check3))       

        cpos(check3,1:3)=ccpos(check3,1:3);
    end
check4= find(((ccpos(:,4).^2+ccpos(:,5).^2+ccpos(:,6).^2)<r.^2)&((ccpos(:,4)-L(1)).^2+ccpos(:,5).^2+ccpos(:,6).^2 > pr.^2)&((ccpos(:,4)-R(1)).^2+ccpos(:,5).^2+ccpos(:,6).^2 > pr.^2));
    if (~isempty(check4))
        
        cpos(check4,4:6)=ccpos(check4,4:6);
    end
        
chrom1(i,:)=cpos(1,:);
ktmtvecttime(i,:)=ktmtbl(1,:);
ktmtbl;
ktmtvect;
a=kt.state+2*kt.state2;
a;
kt.state2;
kt.mt.*(kt.state+kt.state2);    
   cpos 
end
  end
%}
end