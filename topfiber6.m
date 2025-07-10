%% anisotropic (4-direction) topology optimization, adapts from the 88 lines (Andreassen 2011)
% abstract model without real values of force or stiffness
% MBB-BEAM    topfiber6(60,20,0.45,.3)  topfiber6(85,24,0.4,.3)    topfiber6(120,40,0.4,.3)    %.4
% cantilever  topfiber6(120,60,0.4,.3)  topfiber6(160,80,0.4,.3)     topfiber6(200,100,0.4,.3)  

function topfiber6( nelx,nely,fibfrac, shearGxy) 
rmin =1.5; % 1.5, 2.1, 3.1
eN = nelx*nely; % number of elements
nN = (1+nelx)*(1+nely); % number of nodes

%% MATERIAL PROPERTIES  
smlth=5e-9;  % 
%      x1x2   x3x4     x1   x2   x3   x4    1
C11=[   0   shearGxy   1    0   1/4   1/4  smlth]; 
C21=[   0  -shearGxy   0    0   1/4   1/4   0];  
C31=[   0      0       0    0   1/4  -1/4   0];              
C12=[   0  -shearGxy   0    0   1/4   1/4   0];    
C22=[   0   shearGxy   0    1   1/4   1/4   smlth]; 
C32=[   0      0       0    0   1/4  -1/4   0];             
C13=[   0      0       0    0   1/4  -1/4   0];                       
C23=[   0      0       0    0   1/4  -1/4   0];     
C33=[shearGxy  0       0    0   1/4   1/4   smlth];    
CP=[C11; C21; C31;  C12; C22; C32;  C13; C23; C33]; % shape(9,7), effective elasticity matrix

zU=[4 0 3; 0 3 4; -4 0 3;  0 3 -4];
zV=[3 0 4; 0 4 3; -3 0 2;  0 2 -3];
zW=[3 0 2; 0 2 3; -3 0 4 ; 0 4 -3];
zM=[2 0 3; 0 3 2; -2 0 3;  0 3 -2];
MA=zeros(16,9);
MA(1:4,  1:3)=zU;MA(5:8,  7:9)=zU;MA(9:12, 1:3)=-zU;MA(13:16,7:9)=-zU;
MA(1:4,  7:9)=zV;MA(5:8,  4:6)=zV;MA(9:12, 7:9)=zW; MA(13:16, 4:6)=zW;
MB=zeros(16,9);
MB(1:4, 1:3)=-zM;MB(5:8, 7:9)=-zM;MB(9:12, 1:3)= zM;MB(13:16, 7:9)= zM;
MB(1:4, 7:9)=-zW;MB(5:8, 4:6)=-zW;MB(9:12, 7:9)=-zV;MB(13:16, 4:6)=-zV;
KEC = zeros(64,9);  % FEA of 2D quadrilateral element, integral of BTcB, B=LN
for i = 1:4
    for j = 1:4
         KEC(8*(i-1)+j,   :) = MA(4*(i-1)+j,:);
         KEC(8*(i-1)+j+4, :) = MB(4*(i-1)+j,:);
         KEC(8*(i-1)+j+32,:) = MB(4*(i-1)+j,:);
         KEC(8*(i-1)+j+36,:) = MA(4*(i-1)+j,:);
    end
end
KCP= KEC*CP/12;   % shape(64,7)

%PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:nN, 1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,eN,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],eN,1);
iK = reshape(kron(edofMat,ones(8,1))',64*eN,1);
jK = reshape(kron(edofMat,ones(1,8))',64*eN,1);

%% LOADS & SUPPORTS MBB-BEAM, i=1:nelx, j=1:nely, cell(i,j) maps to index 2*j + 2*(i-1)*(nely+1)
F = sparse(2,1,-1,2*nN,1);
fixeddofs = union([1:2: (2*nely+1)], [2*nN]); % HALF MBB-BEAM
U = zeros(2*nN,1);
freedofs = setdiff([1:2*nN],fixeddofs);

%% LOADS and SUPPORTS cantilever
% F = sparse( nely+2 + 2*nelx*(nely+1), 1,-1,2*nN,1); %  force = -1 on y-axis 
% fixeddofs = [1:2*(nely+1)]; 
% U = zeros(2*nN,1);
% freedofs = setdiff( [1:2*nN],fixeddofs);

%% PREPARE FILTER
iH = ones(eN*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
Ht = sparse(jH,iH,sH); % transpose of sparse(iH,jH,sH), shape(eN,eN), 
Hs = ones(4,1)*sum(Ht);  % (4,eN)

%% INITIALIZE ITERATION
rng(1998);
X= fibfrac +0.04*(rand(4, eN)-0.5);  % vol = sum(sum(phy))
phy=X; % shape(4,eN)
loop = 0;
change = 1;     
zeroE = zeros(1,eN); 
onesE = ones(1,eN);

while 80> loop &  (change > 0.01 |  10>loop )    % >0.02  >0.01
    loop = loop + 1;
    x1 = X(1,:); x2 = X(2,:); x3 = X(3,:); x4 = X(4,:);
    mc= zeros(7,eN,4); 
    matK = KCP*[ x1.*x2;  x3.*x4;   x1;   x2;    x3;   x4;   onesE]; %  KCP is(64,7), matK is (64, eN)
    mc(:,:,1)= [     x2;  zeroE; onesE; zeroE; zeroE; zeroE; zeroE]; 
    mc(:,:,2)= [     x1;  zeroE; zeroE; onesE; zeroE; zeroE; zeroE]; 
    mc(:,:,3)= [  zeroE;     x4; zeroE; zeroE; onesE; zeroE; zeroE]; 
    mc(:,:,4)= [  zeroE;     x3; zeroE; zeroE; zeroE; onesE; zeroE]; 

    K = sparse(iK,jK,matK(:));   % K is (2*nN, 2*nN), 0;
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    c = dot( K(freedofs,freedofs) * U(freedofs), U(freedofs));

    %%  Partial derivative of compliance, volume
    c_ = zeros(4,eN);
    for i = 1:4
        matK_ =KCP*mc(:,:,i);
        for ce = 1:eN
            vu=  U(edofMat(ce,:));
            colv =matK_ (:,ce); 
            c_(i,ce) = vu' * reshape(colv, 8, 8) * vu ;
        end
    end
    c_ = (c_./Hs)*Ht;
    l1 = 0; l2 = 1000; move = 0.1;   % l2 = 1e9
    newX=   zeros(4,eN);
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        tta = min(1, min(X+move,X.*sqrt(c_/lmid)) ); % see Andreassen 2011 (4) sqrt(cv/lmid)
        newX = max(0,max(X-move, tta));             %see Andreassen 2011 (3)
        phy = (newX*Ht)./Hs;  
        vol =  sum(sum(phy))/eN; 
        if vol> 4*fibfrac, l1 = lmid; else l2 = lmid; end
    end
  change = max( max(abs(newX-X)));
  X = newX;

  fprintf('It %3i  Obj %7.3f  x1 %5.3f  x2 %5.3f  x3 %5.3f  x4 %5.3f vol %5.3f lmid %5.3f\n',loop,c, mean(phy(1,:)), mean(phy(2,:)), mean(phy(3,:)), mean(phy(4,:)),vol, lmid);
  %% PLOT DENSITIES
  s1=reshape( phy(1,:) ,nely, nelx);s2=reshape( phy(2,:) ,nely, nelx);s3=reshape( phy(3,:) ,nely, nelx);s4=reshape( phy(4,:) ,nely, nelx);
  sho= [(s1+s2+s3+s4)/4; zeros(4,nelx); s1; zeros(4,nelx);  s2; zeros(4,nelx);  s3; zeros(4,nelx);  s4 ]; 
  colormap(gray); imagesc(1- sho); caxis([0 1]); axis equal; axis off; drawnow;
end
