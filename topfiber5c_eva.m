%% othotropic (rotate 45 degree) topology optimization, adapts from the 88 lines (Andreassen 2011)
%% MBB-BEAM  topfiber5c_eva(95,25,0.45,.3, true)       7.6m*1m  4000*(8.1*5.4)/0.315   
%%        %% topfiber5c_eva(120,30,0.4,.3, true)     9.6m*1.2m  4000*(10.2*6.8)/0.35 

function topfiber5c_eva(nelx,nely,fibfrac, shearGxy, runTO) 
rmin =1.5; 
eN = nelx*nely; % number of elements
nN = (1+nelx)*(1+nely); % number of nodes

%% MATERIAL PROPERTIES  
Ebam= 9.4e9;
E0= (9/20)* Ebam; % imagine depth RVE=1.0 or actual depth: 0.007*(9/20)* Ebam; 
smlth=5e-9;  % 
C11=[ shearGxy  1/4   1/4   smlth];  % x1x2, x1, x2, 1
C21=[-shearGxy  1/4   1/4   0];  
C31=[    0      1/4  -1/4   0];              
C12=[-shearGxy  1/4   1/4   0];    
C22=[ shearGxy  1/4   1/4   smlth]; 
C32=[    0      1/4  -1/4   0];             
C13=[    0      1/4  -1/4   0];                       
C23=[    0      1/4  -1/4   0];     
C33=[    0      1/4   1/4   smlth];  % x1x2, x1, x2, 1   
CP=E0*[C11; C21; C31;  C12; C22; C32;  C13; C23; C33]; % shape(9,4), effective elasticity matrix

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
KEC = zeros(64,9);
for i = 1:4
    for j = 1:4
         KEC(8*(i-1)+j,   :) = MA(4*(i-1)+j,:);
         KEC(8*(i-1)+j+4, :) = MB(4*(i-1)+j,:);
         KEC(8*(i-1)+j+32,:) = MB(4*(i-1)+j,:);
         KEC(8*(i-1)+j+36,:) = MA(4*(i-1)+j,:);
    end
end
KCP= KEC*CP/12;   % shape(64,4)

%PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:nN, 1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,eN,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],eN,1);
iK = reshape(kron(edofMat,ones(8,1))',64*eN,1);
jK = reshape(kron(edofMat,ones(1,8))',64*eN,1);

%% LOADS & SUPPORTS MBB-BEAM, i=1:nelx, j=1:nely, cell(i,j) maps to index 2*j + 2*(i-1)*(nely+1)
%centralF =  4000*(8.1*5.4)*(0.007/0.315); % floor load 4kN/m2, RVE thick(0.007), beam width = 315mm
centralF =  4000*(10.2*6.8)/0.35 ;          % floor load 4kN/m2, RVE thick(1.0), beam width = 350mm
F = sparse([2 2+2*(nely+1) 2+4*(nely+1)], [1 1 1], (-1/6)*[centralF  centralF  centralF], 2*nN,1);  %load on 3 points
fixeddofs = union([1:2:(2*nely+1)], [2*nN-4  2*nN-2  2*nN]); % fxi three points 
U = zeros(2*nN,1);
freedofs = setdiff( [1:2*nN],fixeddofs);


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
Hs = ones(2,1)*sum(Ht);  % (2,eN)

%% INITIALIZE ITERATION
if runTO
    rng(1998);
    X= fibfrac + 0.04*(rand(2, eN)-0.5);  % vol = sum(sum(phy))
else
    X= fibfrac + zeros(2, eN);
end
phy=X; % shape(2,eN)
loop = 0;
change = 1;     
zeroE = zeros(1,eN); 
onesE = ones(1,eN);

while 500> loop &  change > 0.01    % >0.02  >0.01
    loop = loop + 1;
    x1 = X(1,:); 
    x2 = X(2,:);
    matK = KCP*[ x1.*x2; x1;  x2;  onesE]; %  KCP is(64,4), matK is (64, eN)
    mc= zeros(4,eN,2);   
    mc(:,:,1)=[x2; onesE; zeroE; zeroE]; 
    mc(:,:,2)=[x1; zeroE; onesE; zeroE]; 

    K = sparse(iK,jK,matK(:));   % K is (2*nN, 2*nN), 0;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    if ~runTO
        break
    end
    c = dot( K(freedofs,freedofs) * U(freedofs), U(freedofs));

    %%  Patrial derivative of compliance, volume
    c_ = zeros(2,eN);
    for i = 1:2
        matK_ =KCP*mc(:,:,i);
        for ce = 1:eN
            vu=  U(edofMat(ce,:));
            colv =matK_ (:,ce); 
            c_(i,ce) = vu' * reshape(colv, 8, 8) * vu ;
        end
    end
    c_ = (c_./Hs)*Ht;
    l1 = 0; l2 = 10; move = 0.1;   % l2 = 1e9
    newX=   zeros(2,eN);
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        tta = min(1, min(X+move,X.*sqrt(c_/lmid)) ); % see Andreassen 2011(4) sqrt(cv/lmid)
        newX = max(0,max(X-move, tta)); %see Andreassen 2011(3)
        phy = (newX*Ht)./Hs;  
        vol =  sum(sum(phy))/eN; 
        if vol> 2*fibfrac, l1 = lmid; else l2 = lmid; end
    end
 
  change = max( max(abs(newX-X)));
  X = newX;

  fprintf('It %3i  Obj %7.3f  phy1 %5.3f  phy2 %5.3f  vol %5.3f lmid %5.3f\n',loop,c, mean(phy(1,:)), mean(phy(2,:)), vol, lmid);
  %% PLOT DENSITIES
  sa=reshape( phy(1,:) ,nely, nelx);
  sb=reshape( phy(2,:) ,nely, nelx);
  sho= [(sa+sb)/2; zeros(4,nelx); sa; zeros(4,nelx);  sb ]; 
  colormap(gray); imagesc(1- sho); caxis([0 1]); axis equal; axis off; drawnow;
end

%% post evaluation
if ~runTO
  sa=reshape( phy(1,:) , nely, nelx);
  sb=reshape( phy(2,:) , nely, nelx);
  sho= [(sa+sb)/2; zeros(4,nelx); sa; zeros(4,nelx);  sb ]; 
  colormap(gray); imagesc(1- sho); caxis([0 1]); axis equal; axis off; drawnow;
end
solColumns = X';  % could paste the data in Workspace to an excel file
maxUx = max(U(1:2:end));
minUx = min(U(1:2:end));
maxUy = max(U(2:2:end));
minUy = min(U(2:2:end));

posiLIM= 116e6/Ebam;  % limit of tensile strain 
negaLIM = -69.7e6/Ebam; % limit of compressive strain 
strainD3 = zeros(nely, nelx); 
for i = 1:nelx
  for j = 1:nely
      aID = 2*(j+1) + 2*(i-1)*(nely+1); %(i,j+1)    % 2*j + 2*(i-1)*(nely+1);
      bID = 2*j +     2*i*(nely+1);     %(i+1,j)    % 2*j + 2*(i-1)*(nely+1);
      ttx = U(bID-1) - U(aID-1);
      tty = U(bID)  -  U(aID);
      tt= (ttx+tty)/0.08;  %   ( (ttx,tty) dot (1/sqrt2, 1/sqr2t) ) / (0.04*sqrt2)
      if tt>=0 
         strainD3(j,i)= tt/posiLIM;
      else
         strainD3(j,i)= tt/negaLIM;
      end
  end
end
strainD4 = zeros(nely, nelx);
for i = 1:nelx
  for j = 1:nely
      aID =  2*j + 2*(i-1)*(nely+1);  % (i,j)  %  ij;   %  2*j + 2*(i-1)*(nely+1);
      bID =  2*(j+1) + 2*i*(nely+1);  % (i+1,j+1)
      ttx = U(bID-1) - U(aID-1);
      tty = U(bID)  -  U(aID); 
      tt = (ttx-tty)/0.08;  %   ( (ttx,tty) dot (1/sqrt2, -1/sqr2t) ) / (0.04*sqrt2)
      if tt>=0 
         strainD4(j,i)= tt/posiLIM;
      else
         strainD4(j,i)=tt/negaLIM;
      end
  end
end
X; % breakpoint this line to see data in Workspace

