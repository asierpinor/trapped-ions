% Solves classical eom of 1D Ising Hamiltonian


seed=-1;    % if negative no seed for random number generator
nosample=1; % if ==1 then initial conditions will have no noise





% System parameters and setup

Ns=10;           % number of spins
Nb=10;            % number of bosons if used for Ions experiment. Otherwise 0.
whichboson=1;   % First boson to take into account and only boson if Nb=1
iter=1;         % runs
nfile=1;
dim=1;          % Dimension of ion system (changes only type of file to be read out)

switch dim
    case 1
        F0s=[20 2 sqrt(18.55*20.33) 3 sqrt(20.84*22.58)];        % Forces
        deltas=[80 8 80 1 95];             % Detunings
        haches=[0 0.1 0.5 1 2 4 5];         % Values of h
        omeCOM=2*pi*0.477223427331887E+04;                   % COM frequency of bosons
    case 2
        F0s=[0.74159 1.75 2.5 6 10];        % Forces
        deltas=[1 5 10 50 100];             % Detunings
        haches=[0 0.1 0.5 1 2 4 5];         % Values of h
        omeCOM=2*pi*1570;                   % COM frequency of bosons
    otherwise
        disp('Problem with dimension chosen');
end

cdet=3;                             % Choose detuning and force
ch=1;                               % Choose transverse field

F0=2*pi*F0s(cdet);                  % Force: F*sqrt(hbar/(2M*omeCOM))/hbar    
M=1;                                % mass
ome=omeCOM;
detuning=deltas(cdet);              % detuning in kHz, with respect to highest frequcency!
muR=omeCOM+2*pi*detuning;           % relative freq
bi0=1/sqrt(Ns);                     % eigenvector of COM-mode 
biq=repmat(bi0,Ns,1);
F=F0*bi0*ones(Ns,Nb);                           % coupling boson-spin for COM
h=haches(ch);                              % Transverse field x

theta=0.25*pi;        % initial orientation of spin: theta
phi=0;          % initial orientation of spin: phi
ICb=0;          % choose Boson initial conditions. 0: vacuum, 1: thermal, 2: coherent state
cEOM=0;         % choose EOM. 0: BBGKYbseom, 1: BBGKYbseomCAN, 2: with transverse field, 3: with transverse field in interaction picture

n0=0.1;                 % initial mean number of bosons
T=omeCOM/log(1/n0+1);   % Temperature of COM
coh0=10.0;                 % Initial value of field for coherent state

dt=0.00002;         % step size
tmax=0.4;           % maximal time
tmeas=0.001;        % measure only in intervals of tmeas
totobs=floor(tmax/tmeas)+1;     % number of times that observables will be saved



%folder='./';
%folder='s_Nbosons_solvable/Monroe_comparison/';
%folder=sprintf('../Bosons_Spins/Nbosons_solvable/Monroe_comparison/BBGKY_nbosons%i/',Nb);
%folder='../Bosons_Spins/Nbosons_solvable/Monroe_comparison/prueba/';
%folder=sprintf('../Bosons_Spins/Nbosons_solvable/Monroe_comparison/BBGKY_nbosons%i_which%i/',Nb,whichboson);
%folder='./Monroe/';
folder=sprintf('../Bosons_Spins/Nbosons_solvable/Monroe_tilted/NoSampling_BBGKY_nbosons%i/',Nb);
%folder=sprintf('../Bosons_Spins/Nbosons_solvable/Monroe_manyMode/NoSampThirdBBGKY_nbosons%i/',Nb);


% Random number generator
if seed>0
    rng(seed,'twister');
    disp('Using manual seed.');
else
    rng('shuffle','twister');
end



Rphi=[cos(phi) -sin(phi) 0;...
      sin(phi) cos(phi) 0;...
      0 0 1];
Rtheta=[cos(theta) 0 sin(theta);...
        0 1 0;...
        -sin(theta) 0 cos(theta)];
R=Rphi*Rtheta;      % Rotation matrix for initial conditions





% Read in modes and amplitues if Nb>1, set other parameters if needed

if Nb>=1
    fr=188;     % rotation frequency in kHz
    ta=0.1844;  % Trap anisotropy
    
    formatSpec = repmat('%f ',1,Ns+1);  % Input format: each row has mode and eigenvector
    
    switch dim
        case 1
            mfname = sprintf('modes_N%i_ta%g.dat',Ns,ta);
        case 2
            mfname = sprintf('modes_N%i_fr%.2f.dat',Ns,fr);
        otherwise
            disp('Problem with dimension chosen');
    end
    Modes = fopen(mfname,'r');
    sizeInput = [ Ns+1 Ns ];
    input = fscanf(Modes,formatSpec,sizeInput); % careful: array is transposed wrt file
    
    %ome=input(1,1:Nb).*2*pi;                          % frequencies of bosons, vector
    ome=input(1,whichboson:whichboson-1+Nb).*2*pi;     % frequencies of bosons, vector
    ome=ome.';                              % -> COLUMN vector
    omeCOM=ome(1);                          % Highest frequency!!! If whichboson=1,then omeCOM is really the frequency of the COM
    muR=omeCOM+2*pi*detuning;               % relative freq, wrt highest frequency
    %biq=input(2:(Ns+1),1:Nb);           % eigenvectors of modes -> Ns x Nb
    biq=input(2:(Ns+1),whichboson:whichboson-1+Nb);           % eigenvectors of modes -> Ns x Nb
    F=F0.*sqrt(omeCOM).*(1./sqrt(repmat(ome.',Ns,1))).*biq;   % coupling boson-spin F_{j,mu}
    
    % WHAT THE HELL WERE THESE LINES USEFUL FOR???? 2D???
    %ome = ome./10;
    %omeCOM=ome(1);
    %muR=omeCOM+2*pi*detuning;           % relative freq
    
    n0=0.1.*ones(Nb,1);                 % initial mean number of bosons
    T=ome./log(1./n0+1);              % Temperature of boson modes for given mean number
    coh0=coh0.*ones(Nb,1);                 % Initial value of field for coherent state
    % All COLUMN vectors
    
    fclose(Modes);
end





    


% Output file

%discr = [];
%discr = [ 'coh=' num2str(coh0(1,1)) ];     % to discriminate simulations with different parameters
%discr = [ 'det=' num2str(detuning) 'F0=' num2str(F0/(2*pi)) ];
%discr = [ '_Ns' num2str(Ns) '_det' num2str(detuning) '_F0' num2str(F0/(2*pi)) '_h' num2str(h) '_cEOM' num2str(cEOM) ];
%discr = [ 'theta=' num2str(theta/pi) ];
%discr = [ '_Ns' num2str(Ns) '_Nb' num2str(Nb) '_dim' num2str(dim) '_det' num2str(detuning) '_F0' num2str(F0/(2*pi)) '_dt' num2str(log10(dt)) ];
discr = [ '_Ns' num2str(Ns) '_Nb' num2str(Nb) '_dim' num2str(dim) '_det' num2str(detuning) '_F0' num2str(F0/(2*pi)) ];
%discr = [ '_Ns' num2str(Ns) '_which' num2str(whichboson) '_Nb' num2str(Nb) '_dim' num2str(dim) '_det' num2str(detuning) '_F0' num2str(F0/(2*pi)) ];

parameters = ['# Parameters:\n# Nb=' num2str(Nb) '\n# Ns=' num2str(Ns) '\n# which=' num2str(whichboson)...
                '\n# iter=' num2str(iter) '\n# nfile=' num2str(nfile)...
                '\n# theta=' num2str(theta) '\n# phi=' num2str(phi)...
                '\n# ICb=' num2str(ICb) '\n# T_COM=' num2str(T(1,1)) '\n# cEOM=' num2str(cEOM) '\n# F0=' num2str(F0)...
                '\n# M=' num2str(M) '\n# detuning=' num2str(detuning) '\n# omeCOM=' num2str(omeCOM) '\n# h=' num2str(h)...
                '\n# dt=' num2str(dt) '\n# tmax=' num2str(tmax) '\n# tmeas=' num2str(tmeas)...
                '\n'];

Scoll=fopen([sprintf('%sThirdBBGKY_spincoll',folder) discr 'file' num2str(nfile) '.txt'],'w');
S1pt=fopen([sprintf('%sThirdBBGKY_spin1pt',folder) discr 'file' num2str(nfile) '.txt'],'w');
S2pt=fopen([sprintf('%sThirdBBGKY_spin2pt',folder) discr 'file' num2str(nfile) '.txt'],'w');
Bdiag=fopen([sprintf('%sThirdBBGKY_bosonsdiag',folder) discr 'file' num2str(nfile) '.txt'],'w');
B2pt=fopen([sprintf('%sThirdBBGKY_bosons2pt',folder) discr 'file' num2str(nfile) '.txt'],'w');
M2pt=fopen([sprintf('%sThirdBBGKY_mixed2pt',folder) discr 'file' num2str(nfile) '.txt'],'w');


fprintf(Scoll,['# Evolution of collective spin moments: '...
                't | Sx | Sy | Sz | Sx^2 | Sy^2 | Sz^2 | SxSy | SxSz | SySz\n']);
fprintf(S1pt,['# Evolution of spin 1 point function: '...
                't | (spin 1:) Sx | Sy | Sz | (spin 2:) ... etc. \n']);
fprintf(S2pt,['# Evolution of spin 2 point connected functions: '...
                't | (1,2) SxSx | SxSy | SxSz | SySx | SySy | SySz | SzSx | SzSy | SzSz | (1,3) ... (1,Ns) ... (2,3) ... (2,Ns) ... etc. \n']);
fprintf(Bdiag,['# Evolution of diagonal moments of bosons: '...
                't | (1st boson:) Re(a) | Im(a) | a*a^dagger | (a*a^dagger)^2 (2nd boson:) etc. \n']);
fprintf(B2pt,['# Evolution of 2pt function of bosons: '...
                't | (1,2) Re(a_1^dagger a_2) | Im(a_1^dagger a_2) | a_1^dagger a_2^dagger a_2 a_1 | '...
                '(1,3) ... (1,Nb) ... (2,3) ... (2,Nb) ... etc. \n']);
fprintf(M2pt,['# Evolution of 2pt function of mixed spin and boson products: '...
                't | (1,1) Re(a_1 Sx1) | Im(a_1 Sx1) | Re(a_1 Sy1) | Im(a_1 Sy1) | Re(a_1 Sz1) | Im(a_1 Sz1) | '...
                '(1,2) ... (1,Ns) ... (2,1) ... (2,Ns) etc. (Nb,Ns) \n']);

            
fprintf(Scoll,parameters);
fprintf(S1pt,parameters);
fprintf(S2pt,parameters);
fprintf(Bdiag,parameters);
fprintf(B2pt,parameters);
fprintf(M2pt,parameters);






% State vector: [ sx1 sy1 sz1   sx2...   S1xS2x S1xS2y S1xS2z S1yS2x... S1xS3x...    S2xS3x...   S2xS4x...   S3xS4x...
%                 a1 a2 ...   A00_11 A00_12 ...    A00_22 A00_23 ...
%                 A10_11 A10_12 ...    A10_21 A10_22 ...
%                 Mx11 My11 Mz11  Mx12 My12 Mz12 ...    Mx21 My21 Mz21 ...]

state=zeros(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2+1,iter);  % each column is state for different run

sx=1:3:3*Ns;
sy=2:3:3*Ns;
sz=3:3:3*Ns;    
bp=(3*Ns+9*Ns*(Ns-1)/2+1):(3*Ns+9*Ns*(Ns-1)/2+Nb);          % Indexing of 1pt-function
dummy=3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2+1;      % Dummy variable. Set to zero to simplify sums in EOM

SS=(3*Ns+1):(3*Ns+9*Ns*(Ns-1)/2);   % Indexing of all spin-spin 2pt-functions

XX=(3*Ns+1):9:(3*Ns+9*Ns*(Ns-1)/2);
XY=(3*Ns+2):9:(3*Ns+9*Ns*(Ns-1)/2);
XZ=(3*Ns+3):9:(3*Ns+9*Ns*(Ns-1)/2);

YX=(3*Ns+4):9:(3*Ns+9*Ns*(Ns-1)/2);
YY=(3*Ns+5):9:(3*Ns+9*Ns*(Ns-1)/2);
YZ=(3*Ns+6):9:(3*Ns+9*Ns*(Ns-1)/2);

ZX=(3*Ns+7):9:(3*Ns+9*Ns*(Ns-1)/2);
ZY=(3*Ns+8):9:(3*Ns+9*Ns*(Ns-1)/2);
ZZ=(3*Ns+9):9:(3*Ns+9*Ns*(Ns-1)/2);     % Indexing of all 2pt functions depending on {x,y,z}

sxsx=zeros(Ns,Ns);  sysx=zeros(Ns,Ns);  szsx=zeros(Ns,Ns);
sxsy=zeros(Ns,Ns);  sysy=zeros(Ns,Ns);  szsy=zeros(Ns,Ns);
sxsz=zeros(Ns,Ns);  sysz=zeros(Ns,Ns);  szsz=zeros(Ns,Ns);      % These matrices contain the index of all 2pt functions based on i,j and {x,y,z}


for ii=1:Ns
    for jj=1:Ns
        
        if ii<jj
            
                sxsx(ii,jj)=XX( (ii-1)*(Ns-ii/2)+(jj-ii) );         % Offset for a given ii is (Ns-1)+(Ns-2)+...+(Ns-(ii-1))
                sxsy(ii,jj)=XY( (ii-1)*(Ns-ii/2)+(jj-ii) );
                sxsz(ii,jj)=XZ( (ii-1)*(Ns-ii/2)+(jj-ii) );
                
                sysx(ii,jj)=YX( (ii-1)*(Ns-ii/2)+(jj-ii) );
                sysy(ii,jj)=YY( (ii-1)*(Ns-ii/2)+(jj-ii) );
                sysz(ii,jj)=YZ( (ii-1)*(Ns-ii/2)+(jj-ii) );
                
                szsx(ii,jj)=ZX( (ii-1)*(Ns-ii/2)+(jj-ii) );
                szsy(ii,jj)=ZY( (ii-1)*(Ns-ii/2)+(jj-ii) );
                szsz(ii,jj)=ZZ( (ii-1)*(Ns-ii/2)+(jj-ii) );
                
        else if ii>jj
                
                sxsx(ii,jj)=XX( (jj-1)*(Ns-jj/2)+(ii-jj) );
                sxsy(ii,jj)=YX( (jj-1)*(Ns-jj/2)+(ii-jj) );         % The position of SixSjy=sxsy(i,j) has to be the same as SjySix=sysx(j,i)
                sxsz(ii,jj)=ZX( (jj-1)*(Ns-jj/2)+(ii-jj) );
                
                sysx(ii,jj)=XY( (jj-1)*(Ns-jj/2)+(ii-jj) );
                sysy(ii,jj)=YY( (jj-1)*(Ns-jj/2)+(ii-jj) );
                sysz(ii,jj)=ZY( (jj-1)*(Ns-jj/2)+(ii-jj) );
                
                szsx(ii,jj)=XZ( (jj-1)*(Ns-jj/2)+(ii-jj) );
                szsy(ii,jj)=YZ( (jj-1)*(Ns-jj/2)+(ii-jj) );
                szsz(ii,jj)=ZZ( (jj-1)*(Ns-jj/2)+(ii-jj) );
                
            else
                
                sxsx(ii,jj)=dummy;  sysx(ii,jj)=dummy;  szsx(ii,jj)=dummy;
                sxsy(ii,jj)=dummy;  sysy(ii,jj)=dummy;  szsy(ii,jj)=dummy;         % Diagonal filled with dummy variable instead of zero -> C_ii=0
                sxsz(ii,jj)=dummy;  sysz(ii,jj)=dummy;  szsz(ii,jj)=dummy;
                
            end
        end     
        
    end
end

BB00=(3*Ns+9*Ns*(Ns-1)/2+Nb+1):(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2);         % Indexing of all (aa) functions
BB10=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+1):(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2);        % Indexing of all (a^dagger a) functions
XB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+1):3:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb);
YB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+2):3:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb);
ZB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3):3:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb);   % Indexing of all (spin a) functions

Bp00=zeros(Nb,Nb);  % (aa) function depending on sites
Bp10=zeros(Nb,Nb);  % (a^dagger a)
SxB=zeros(Ns,Nb);
SyB=zeros(Ns,Nb);
SzB=zeros(Ns,Nb);   % (spin a)

for mm=1:Nb
    for nn=1:Nb
        if mm<=nn
            Bp00(mm,nn)=BB00( (mm-1)*(Nb+1-mm/2)+(nn-mm+1) );
        else
            Bp00(mm,nn)=BB00( (nn-1)*(Nb+1-nn/2)+(mm-nn+1) );
        end
        Bp10(mm,nn)=BB10( (mm-1)*Nb+nn );
    end
end

for ii=1:Ns
    for mm=1:Nb
        SxB(ii,mm)=XB( (ii-1)*Nb+mm );
        SyB(ii,mm)=YB( (ii-1)*Nb+mm );
        SzB(ii,mm)=ZB( (ii-1)*Nb+mm );
    end
end


XXB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+1):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);
XYB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+2):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);
XZB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+3):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);

YXB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+4):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);
YYB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+5):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);  % Indexing of spin-spin-boson functions
YZB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+6):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);

ZXB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+7):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);
ZYB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+8):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);
ZZB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9):9:(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);

SSB=(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+1):(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2);  % Indexing of all spin-spin-boson functions

SxSxB = cell(1,Nb); SxSyB = cell(1,Nb); SxSzB = cell(1,Nb);
SySxB = cell(1,Nb); SySyB = cell(1,Nb); SySzB = cell(1,Nb);
SzSxB = cell(1,Nb); SzSyB = cell(1,Nb); SzSzB = cell(1,Nb);

for mm=1:Nb
    SxSxB{mm}=zeros(Ns,Ns); SxSyB{mm}=zeros(Ns,Ns); SxSzB{mm}=zeros(Ns,Ns);
    SySxB{mm}=zeros(Ns,Ns); SySyB{mm}=zeros(Ns,Ns); SySzB{mm}=zeros(Ns,Ns);
    SzSxB{mm}=zeros(Ns,Ns); SzSyB{mm}=zeros(Ns,Ns); SzSzB{mm}=zeros(Ns,Ns);
end


for mm=1:Nb

    for ii=1:Ns
        for jj=1:Ns

            if ii<jj

                    SxSxB{mm}(ii,jj)=XXB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );         % Offset for a given ii is (Ns-1)+(Ns-2)+...+(Ns-(ii-1))
                    SxSyB{mm}(ii,jj)=XYB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    SxSzB{mm}(ii,jj)=XZB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );

                    SySxB{mm}(ii,jj)=YXB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    SySyB{mm}(ii,jj)=YYB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    SySzB{mm}(ii,jj)=YZB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    
                    SzSxB{mm}(ii,jj)=ZXB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    SzSyB{mm}(ii,jj)=ZYB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );
                    SzSzB{mm}(ii,jj)=ZZB( (ii-1)*(Ns-ii/2)+(jj-ii) + (mm-1)*Ns*(Ns-1)/2 );

            else if ii>jj

                    SxSxB{mm}(ii,jj)=XXB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    SxSyB{mm}(ii,jj)=YXB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );         % The position of SixSjy=sxsy(i,j) has to be the same as SjySix=sysx(j,i)
                    SxSzB{mm}(ii,jj)=ZXB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );

                    SySxB{mm}(ii,jj)=XYB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    SySyB{mm}(ii,jj)=YYB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    SySzB{mm}(ii,jj)=ZYB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    
                    SzSxB{mm}(ii,jj)=XZB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    SzSyB{mm}(ii,jj)=YZB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );
                    SzSzB{mm}(ii,jj)=ZZB( (jj-1)*(Ns-jj/2)+(ii-jj) + (mm-1)*Ns*(Ns-1)/2 );

                else

                    SxSxB{mm}(ii,jj)=dummy;  SxSyB{mm}(ii,jj)=dummy;  SxSzB{mm}(ii,jj)=dummy;
                    SySxB{mm}(ii,jj)=dummy;  SySyB{mm}(ii,jj)=dummy;  SySzB{mm}(ii,jj)=dummy;         % Diagonal filled with dummy variable instead of zero -> C_ii=0
                    SzSxB{mm}(ii,jj)=dummy;  SzSyB{mm}(ii,jj)=dummy;  SzSzB{mm}(ii,jj)=dummy;

                end
            end     

        end
    end

end


% Grouping indices of 3-pt functions for easier implementation. For fixed
% spin index (i,j) group all mu-indices together

SxSxBfixSind = cell(Ns,Ns);  SxSyBfixSind = cell(Ns,Ns);  SxSzBfixSind = cell(Ns,Ns);
SySxBfixSind = cell(Ns,Ns);  SySyBfixSind = cell(Ns,Ns);  SySzBfixSind = cell(Ns,Ns);  % Indices of 3-point function for fixed (i,j) and all mu
SzSxBfixSind = cell(Ns,Ns);  SzSyBfixSind = cell(Ns,Ns);  SzSzBfixSind = cell(Ns,Ns);

for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        
        
        SxSxBij = zeros(1,Nb);  SxSyBij = zeros(1,Nb);  SxSzBij = zeros(1,Nb);
        SySxBij = zeros(1,Nb);  SySyBij = zeros(1,Nb);  SySzBij = zeros(1,Nb);  % Indices of 3-point function for fixed (i,j) and all mu
        SzSxBij = zeros(1,Nb);  SzSyBij = zeros(1,Nb);  SzSzBij = zeros(1,Nb);
        
        SxSxBji = zeros(1,Nb);  SxSyBji = zeros(1,Nb);  SxSzBji = zeros(1,Nb);
        SySxBji = zeros(1,Nb);  SySyBji = zeros(1,Nb);  SySzBji = zeros(1,Nb);  % Indices of 3-point function for fixed (j,i) and all mu
        SzSxBji = zeros(1,Nb);  SzSyBji = zeros(1,Nb);  SzSzBji = zeros(1,Nb);
        
        for mm=1:Nb
           SxSxBij(mm) = SxSxB{mm}(ii,jj);  SxSxBji(mm) = SxSxB{mm}(jj,ii);
           SxSyBij(mm) = SxSyB{mm}(ii,jj);  SxSyBji(mm) = SxSyB{mm}(jj,ii);
           SxSzBij(mm) = SxSzB{mm}(ii,jj);  SxSzBji(mm) = SxSzB{mm}(jj,ii);
           
           SySxBij(mm) = SySxB{mm}(ii,jj);  SySxBji(mm) = SySxB{mm}(jj,ii);
           SySyBij(mm) = SySyB{mm}(ii,jj);  SySyBji(mm) = SySyB{mm}(jj,ii);
           SySzBij(mm) = SySzB{mm}(ii,jj);  SySzBji(mm) = SySzB{mm}(jj,ii);
           
           SzSxBij(mm) = SzSxB{mm}(ii,jj);  SzSxBji(mm) = SzSxB{mm}(jj,ii);
           SzSyBij(mm) = SzSyB{mm}(ii,jj);  SzSyBji(mm) = SzSyB{mm}(jj,ii);
           SzSzBij(mm) = SzSzB{mm}(ii,jj);  SzSzBji(mm) = SzSzB{mm}(jj,ii);
        end
        
        SxSxBfixSind{ii,jj} = SxSxBij;  SxSxBfixSind{jj,ii} = SxSxBji;
        SxSyBfixSind{ii,jj} = SxSyBij;  SxSyBfixSind{jj,ii} = SxSyBji;
        SxSzBfixSind{ii,jj} = SxSzBij;  SxSzBfixSind{jj,ii} = SxSzBji;
        
        SySxBfixSind{ii,jj} = SySxBij;  SySxBfixSind{jj,ii} = SySxBji;
        SySyBfixSind{ii,jj} = SySyBij;  SySyBfixSind{jj,ii} = SySyBji;
        SySzBfixSind{ii,jj} = SySzBij;  SySzBfixSind{jj,ii} = SySzBji;
        
        SzSxBfixSind{ii,jj} = SzSxBij;  SzSxBfixSind{jj,ii} = SzSxBji;
        SzSyBfixSind{ii,jj} = SzSyBij;  SzSyBfixSind{jj,ii} = SzSyBji;
        SzSzBfixSind{ii,jj} = SzSzBij;  SzSzBfixSind{jj,ii} = SzSzBji;
        
    end
end


indices={sx,sy,sz,sxsx,sxsy,sxsz,sysx,sysy,sysz,szsx,szsy,szsz,bp,Bp00,Bp10,SxB,SyB,SzB};
indices3=[ SxSxB ; SxSyB ; SxSzB ; SySxB ; SySyB ; SySzB ; SzSxB ; SzSyB ; SzSzB ];
indicesfixSind=[ SxSxBfixSind ; SxSyBfixSind ; SxSzBfixSind ; SySxBfixSind ; SySyBfixSind ; SySzBfixSind ; SzSxBfixSind ; SzSyBfixSind ; SzSzBfixSind ];




% Initial conditions (spins,dTWA)

if nosample~=1
    state(sx,:)=2.*(randi(2,Ns,iter)-1.5);
    state(sy,:)=2.*(randi(2,Ns,iter)-1.5); 
end   
state(sz,:)=1;

for ii=1:Ns
   state([sx(ii) sy(ii) sz(ii)],:) = R*state([sx(ii) sy(ii) sz(ii)],:);
end % Rotate

state(SS,:)=0;      % Connected correlation functions initialized to 0



% Initial conditions bosons

sigma=zeros(Nb,iter);
mean=zeros(Nb,iter);

switch ICb
    case 0  % Vacuum
        sigma=0.5.*ones(Nb,iter);
    case 1  % Thermal
        nBE=1./(exp(ome./T)-1);
        sigma=repmat(sqrt((nBE+0.5)./2),1,iter);
    case 2  % Coherent
        sigma=0.5.*ones(Nb,iter);
        mean=repmat(coh0,1,iter);
    otherwise
        disp('Error with initial conditions\n');     
end

state(bp,:) = mean;
if nosample~=1
    state(bp,:) = state(bp,:) + sigma.*randn(Nb,iter) + 1i*sigma.*randn(Nb,iter);
end
%Amu0 = state(bp,:);     % save initial state of bosons

state(Bp00,:)=0;
state(Bp10,:)=0;
state(SxB,:)=0;
state(SyB,:)=0;
state(SzB,:)=0;         % connected 2-pt functions

if nosample==1 && ICb==0
    for mm=1:Nb
       state(Bp10(mm,mm),:) = 0.5;      % Quantum one-half
    end
end

state(SSB,:)=0;         % connected 3-pt functions


        




% Integrate EOM

t=0;
nextmeas=0;
tobs=1;

tic;





while t<tmax
    
    % Evolve
    
    switch cEOM
        case 0 % Original Hamiltonian
            
            % Crank-Nicholson (S0,Ns,Nb,ome,muR,F,t,Amu0,iter,indices)
            
            Sbar=state+ThirdBBGKYbseom(state,Ns,Nb,ome,muR,F,t,iter,indices,indices3,indicesfixSind).*dt.*0.5;
            Sbar=state+ThirdBBGKYbseom(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,indices3,indicesfixSind).*dt.*0.5;
            Sbar=state+ThirdBBGKYbseom(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,indices3,indicesfixSind).*dt.*0.5;
            state=2.*Sbar-state;
            
            %Sbar0=state+bseomRWA(state,Ns,Nb,ome,muR,F,t,iter).*dt.*0.5;
            %state=Sbar0+bseomRWA(Sbar0,Ns,Nb,ome,muR,F,t,iter).*dt.*0.5;
            %state=Sbar0+bseomRWA(state,Ns,Nb,ome,muR,F,t,iter).*dt.*0.5;
            %state=Sbar0+bseomRWA(state,Ns,Nb,ome,muR,F,t,iter).*dt.*0.5;
            
            % Runge-Kutta:
            
            %k1=dt.*bseomRWA(state,Ns,Nb,ome,muR,F,t,iter);
            %k2=dt.*bseomRWA(state+0.5.*k1,Ns,Nb,ome,muR,F,t+0.5*dt,iter);
            %k3=dt.*bseomRWA(state+0.5.*k2,Ns,Nb,ome,muR,F,t+0.5*dt,iter);
            %k4=dt.*bseomRWA(state+k3,Ns,Nb,ome,muR,F,t+dt,iter);
            
            %state = state + 1/6.*k1 + 1/3.*k2 + 1/3.*k3 + 1/6.*k4;
            
            
        case 1 % For IC along x, y or z directions.
            
            % Crank-Nicholson (S0,Ns,Nb,ome,muR,F,t,Amu0,iter,indices)
            
            Sbar=state+ThirdBBGKYbseomCAN(state,Ns,Nb,ome,muR,F,t,iter,indices,indices3).*dt.*0.5;
            Sbar=state+ThirdBBGKYbseomCAN(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,indices3).*dt.*0.5;
            Sbar=state+ThirdBBGKYbseomCAN(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,indices3).*dt.*0.5;
            state=2.*Sbar-state;
            
            
            
        case 2 % For transverse field. NOT IMPLEMENTED
            
            % Crank-Nicholson (S0,Ns,Nb,ome,muR,F,t,Amu0,iter,indices)
            
            Sbar=state+FullBBGKYbseomTRANSV(state,Ns,Nb,ome,muR,F,t,iter,indices,h).*dt.*0.5;
            Sbar=state+FullBBGKYbseomTRANSV(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,h).*dt.*0.5;
            Sbar=state+FullBBGKYbseomTRANSV(Sbar,Ns,Nb,ome,muR,F,t,iter,indices,h).*dt.*0.5;
            state=2.*Sbar-state;
            
            
        otherwise
            
            disp('ERROR: Problem with EOM');
            
    end
    
    
    t=t+dt;
    
 
    
    % Output observables
    
    if t>=nextmeas
        
        % Collective spin moments
        
        Stemp = [ sum(state(sx,:),1)./Ns ;...
                  sum(state(sy,:),1)./Ns ;...
                  sum(state(sz,:),1)./Ns ];     % collective spin for each iteration
              
        Si2 = [ sum(state(sx,:).^2,1)./Ns^2 ;...
              sum(state(sy,:).^2,1)./Ns^2 ;...
              sum(state(sz,:).^2,1)./Ns^2 ];     % Sum of sx_i^2, sy_i^2 etc.

        Ctemp = zeros(6,1);
              
        for ii=1:Ns
            for jj=1:Ns
                if ii~=jj
                    Ctemp = Ctemp +  [ sum(state(sxsx(ii,jj),:))./(Ns^2) ;...
                                       sum(state(sxsy(ii,jj),:))./(Ns^2) ;...
                                       sum(state(sxsz(ii,jj),:))./(Ns^2) ;...
                                       sum(state(sysy(ii,jj),:))./(Ns^2) ;...
                                       sum(state(sysz(ii,jj),:))./(Ns^2) ;...
                                       sum(state(szsz(ii,jj),:))./(Ns^2) ];
                end
            end
        end
        
        % Sum over iterations
        S1 = sum(Stemp,2)./iter;
        %S2 = sum(Stemp.^2,2)./iter;
        S2 = sum(Stemp.^2-Si2,2)./iter + 1/Ns;      % Substract self-term and add 1 instead
        SS= sum(Stemp([1 1 2],:).*Stemp([2 3 3],:),2)./iter;    % CAREFUL: need to substract self-term here too
        Css = sum(Ctemp,2)./iter;
        
        S2tot = S2 + Css([1 4 6],1);
        SStot = SS + Css([2 3 5],1);
        
        fprintf(Scoll,'%.8g\t %.8g\t %.8g\t %.8g\t',t,S1);
        fprintf(Scoll,'%.8g\t %.8g\t %.8g\t',S2tot);
        fprintf(Scoll,'%.8g\t %.8g\t %.8g\n',SStot);
        
        

        % Si expectation values
        
        fprintf(S1pt,'%.8g',t);
        %{
        for ii=1:Ns

            Si = sum(state([sx(ii) sy(ii) sz(ii)],:),2)./iter;

            fprintf(S1pt,'\t%.8g \t%.8g \t%.8g',Si); % Output of Si exp value

        end
        %}
        fprintf(S1pt,'\n');
        
        
        
        
        % SiSj correlations
        
        fprintf(S2pt,'%.8g',t);
        
        if Ns>1
        
            for ii=1:(Ns-1)
                for jj=(ii+1):Ns

                    SiSj = sum(state([sx(ii) sx(ii) sx(ii) sy(ii) sy(ii) sy(ii) sz(ii) sz(ii) sz(ii)],:)...
                              .*state([sx(jj) sy(jj) sz(jj) sx(jj) sy(jj) sz(jj) sx(jj) sy(jj) sz(jj)],:),2)./iter;
                    Cij = sum([ state(sxsx(ii,jj),:) ; state(sxsy(ii,jj),:) ; state(sxsz(ii,jj),:) ; ...
                                state(sysx(ii,jj),:) ; state(sysy(ii,jj),:) ; state(sysz(ii,jj),:) ; ...
                                state(szsx(ii,jj),:) ; state(szsy(ii,jj),:) ; state(szsz(ii,jj),:) ],2)./iter;
                    SiSjdisc = sum(state([sx(ii) sx(ii) sx(ii) sy(ii) sy(ii) sy(ii) sz(ii) sz(ii) sz(ii)],:),2)/iter...
                            .*sum(state([sx(jj) sy(jj) sz(jj) sx(jj) sy(jj) sz(jj) sx(jj) sy(jj) sz(jj)],:),2)./iter;
                    SiSjconn = SiSj + Cij - SiSjdisc;

                    fprintf(S2pt,'\t%.8g \t%.8g \t%.8g \t%.8g \t%.8g \t%.8g \t%.8g \t%.8g \t%.8g',SiSjconn); % Output of ij correlations

                end
            end
        
        end
        
        fprintf(S2pt,'\n');
        
        
        
        
        % Diagonal boson moments
        
        fprintf(Bdiag,'%.8g',t);
        
        B1 = sum(state(bp,:),2)./iter;                                    % boson 1pt function
        B2 = sum(abs(state(bp,:)).^2 + real(state(diag(Bp10),:)),2)./iter;      % boson 2pt function. diag(Bp10) should be real
        B4 = sum(abs(state(diag(Bp00),:)).^2 + 2.*abs(state(diag(Bp10),:)).^2 + 2.*real(state(diag(Bp00),:).*(conj(state(bp,:))).^2)...
                 + 4.*real(state(diag(Bp10),:)).*abs(state(bp,:)).^2 + abs(state(bp,:)).^4,2)./iter;      
                          % boson 4pt function, calculated as: C11*C00 + 2*C10^2 + 2*Re(C00*(psi^*)^2) + 3*C10|psi|^2 + |psi|^4
                     
        Ball = [ real(B1.') ; imag(B1.') ; B2.' ; B4.'];        % boson exp. values

        fprintf(Bdiag,'\t%.8g \t%.8g \t%.8g \t%.8g',Ball);
        fprintf(Bdiag,'\n');                        % output bosons
        
        
        
        
        
        % Non-diagonal mu-nu boson moments
        
        fprintf(B2pt,'%.8g',t);
        %{
        if Nb>1
        
            for mm=1:(Nb-1)
                for nn=(mm+1):Nb

                    BmuBnu = sum( conj(state(bp(mm),:)).*state(bp(nn),:) + state(Bp10(mm,nn),:) ,2)./iter;
                    BmuBnu2 = sum( abs(state(Bp00(mm,nn),:)).^2 + real(state(Bp10(mm,mm),:).*state(Bp10(nn,nn),:))...
                                   + abs(state(Bp10(mm,nn),:)).^2 + 2.*real(state(Bp00(mm,nn),:).*conj(state(bp(mm),:).*state(bp(nn),:)))...
                                   + real(state(Bp10(mm,mm),:)).*abs(state(bp(nn),:)).^2 + real(state(Bp10(nn,nn),:)).*abs(state(bp(mm),:)).^2 ...
                                   + 2.*real(state(Bp10(mm,nn),:).*conj(state(bp(nn),:)).*state(bp(mm),:)) + abs(state(bp(mm),:).*state(bp(nn),:)).^2,2)./iter;

                    fprintf(B2pt,'\t%.8g \t%.8g \t%.8g',real(BmuBnu),imag(BmuBnu),BmuBnu2); % Output of mu-nu correlations

                end
            end
        
        end
        %}
        fprintf(B2pt,'\n');
        
        
        
        
        % BmuSi correlations
        
        fprintf(M2pt,'%.8g',t);
        
        for mm=1:Nb
            for ii=1:Ns

                BmuSi = sum(state([bp(mm) bp(mm) bp(mm)],:).*state([sx(ii) sy(ii) sz(ii)],:),2)./iter;
                BmuSiconnw = sum(state([SxB(ii,mm) SyB(ii,mm) SzB(ii,mm)],:),2)./iter;
                BmuSi = BmuSi + BmuSiconnw;
                BStemp = [ real(BmuSi(1)) ; imag(BmuSi(1)) ; real(BmuSi(2)) ; imag(BmuSi(2)) ; real(BmuSi(3)) ; imag(BmuSi(3)) ];

                fprintf(M2pt,'\t%.8g \t%.8g \t%.8g \t%.8g \t%.8g \t%.8g',BStemp);   % Output of mu-j correlations

                %SxB(ii,mm)
                
            end
        end
        
        fprintf(M2pt,'\n');
        
        
        
 
        
        
        
        
        
        nextmeas = nextmeas + tmeas;
        tobs = tobs+1;
        
    end
    
    
    
end
                                                                                     

%}

fclose(Scoll);
fclose(S1pt);
fclose(S2pt);
fclose(Bdiag);
fclose(B2pt);
fclose(M2pt);
%disp(state);

disp(toc);







