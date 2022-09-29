function S = ThirdBBGKYbseomCAN(S0,Ns,Nb,ome,muR,F,t,iter,indices,indices3)
% ieom computes the rhs of EOM of Ising Hamiltonian
% 

S=zeros(3*Ns+9*Ns*(Ns-1)/2+Nb+Nb*(Nb+1)/2+Nb^2+3*Ns*Nb+9*Nb*Ns*(Ns-1)/2+1,iter);  % fill output with zeros

% Indexing

sx=indices{1};
sy=indices{2};
sz=indices{3};

sxsx=indices{4};
sxsy=indices{5};
sxsz=indices{6};

sysx=indices{7};
sysy=indices{8};
sysz=indices{9};

szsx=indices{10};
szsy=indices{11};
szsz=indices{12};

bp=indices{13};

Bp00=indices{14};  % (aa) function depending on sites
Bp10=indices{15};  % (a^dagger a)
SxB=indices{16};
SyB=indices{17};
SzB=indices{18};   % (S_z a) etc.

SxSxB=indices3(1,1:Nb);
SxSyB=indices3(2,1:Nb);
SxSzB=indices3(3,1:Nb);

SySxB=indices3(4,1:Nb);
SySyB=indices3(5,1:Nb);
SySzB=indices3(6,1:Nb);

SzSxB=indices3(7,1:Nb);
SzSyB=indices3(8,1:Nb);
SzSzB=indices3(9,1:Nb);







detuning=muR-ome;       % detuning of each mode from muR, Nb x 1
edet = exp(1i.*detuning.*t);
Fe = F.*repmat(edet.',Ns,1);    % Ns*Nb
Fec = conj(Fe);    % Ns*Nb



% EOM for 1pt functions
    
%%%%% Spins

for ii=1:Ns
    
    S(sx(ii),:)=2.*real( Fe(ii,:)*S0(SyB(ii,:),:) + S0(sy(ii),:).*( Fe(ii,:)*S0(bp,:) ) );
    S(sy(ii),:)=-2.*real( Fe(ii,:)*S0(SxB(ii,:),:) + S0(sx(ii),:).*( Fe(ii,:)*S0(bp,:) ) );
    
end

S(sz,:)=0;

%%%%% Bosons
    
S(bp,:)=1i.*0.5.*repmat(exp(-1i*detuning*t),1,iter).*((F.')*S0(sz,:));


    





% EOM for 2pt and 3pt functions



%%%%% Spin-Spin and Spin-Spin-Boson

for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        
        
        SxSxBij = zeros(1,Nb);  SxSyBij = zeros(1,Nb);
        SySxBij = zeros(1,Nb);  SySyBij = zeros(1,Nb);  % Indices of 3-point function for fixed (i,j) and all mu
        
        SxSxBji = zeros(1,Nb);  SxSyBji = zeros(1,Nb);
        SySxBji = zeros(1,Nb);  SySyBji = zeros(1,Nb);
        for mm=1:Nb
           SxSxBij(mm) = SxSxB{mm}(ii,jj);  SxSxBji(mm) = SxSxB{mm}(jj,ii);
           SxSyBij(mm) = SxSyB{mm}(ii,jj);  SxSyBji(mm) = SxSyB{mm}(jj,ii);
           SySxBij(mm) = SySxB{mm}(ii,jj);  SySxBji(mm) = SySxB{mm}(jj,ii);
           SySyBij(mm) = SySyB{mm}(ii,jj);  SySyBji(mm) = SySyB{mm}(jj,ii);
        end
        
        
        % Spin-Spin Terms
        
        % X terms
        
        S(sxsx(ii,jj),:)=2.*( real( Fe(ii,:)*S0(SySxBij,:) + S0(sysx(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(SxB(jj,:),:) ) )...
                            + real( Fe(jj,:)*S0(SySxBji,:) + S0(sysx(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sy(jj),:).*( Fe(jj,:)*S0(SxB(ii,:),:) ) ) );
                        
        S(sxsy(ii,jj),:)=2.*( real( Fe(ii,:)*S0(SySyBij,:) + S0(sysy(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(SyB(jj,:),:) ) )...
                            - real( Fe(jj,:)*S0(SxSxBji,:) + S0(sxsx(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sx(jj),:).*( Fe(jj,:)*S0(SxB(ii,:),:) ) ) );
        
        S(sxsz(ii,jj),:)=0;
                        
                        
                        
        % Y terms
                        
        S(sysx(ii,jj),:)=2.*( real( Fe(jj,:)*S0(SySyBji,:) + S0(sysy(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sy(jj),:).*( Fe(jj,:)*S0(SyB(ii,:),:) ) )...
                            - real( Fe(ii,:)*S0(SxSxBij,:) + S0(sxsx(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sx(ii),:).*( Fe(ii,:)*S0(SxB(jj,:),:) ) ) );
                        
        S(sysy(ii,jj),:)=-2.*( real( Fe(ii,:)*S0(SxSyBij,:) + S0(sxsy(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sx(ii),:).*( Fe(ii,:)*S0(SyB(jj,:),:) ) )...
                             + real( Fe(jj,:)*S0(SxSyBji,:) + S0(sxsy(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sx(jj),:).*( Fe(jj,:)*S0(SyB(ii,:),:) ) ) );
                        
        S(sysz(ii,jj),:)=0;
                        
                        
                          
        % Z terms
                        
        S(szsx(ii,jj),:)=0;
                        
        S(szsy(ii,jj),:)=0;
                        
        S(szsz(ii,jj),:)=0;
                        
        
        
        
        % Spin-Spin-Boson Terms
        
        temp = 0.5*1i*repmat(conj(edet),1,iter).*( repmat(F(ii,:).',1,iter).*repmat(S0(sz(ii),:),Nb,1) ...
                                                 + repmat(F(jj,:).',1,iter).*repmat(S0(sz(jj),:),Nb,1) );
        
        % X terms
        
        S(SxSxBij,:)=2.*( S0(SySxBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SyB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SxB(jj,:),:)),Nb,1)...
                        + S0(SySxBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SyB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SxB(ii,:),:)),Nb,1) )...
                        -temp.*repmat(S0(sxsx(ii,jj),:),Nb,1);
                        
        S(SxSyBij,:)=2.*( S0(SySyBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SyB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SyB(jj,:),:)),Nb,1)...
                        - S0(SxSxBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) - S0(SxB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SxB(ii,:),:)),Nb,1) )...
                        -temp.*repmat(S0(sxsy(ii,jj),:),Nb,1);
                        
                        
        % Y terms
                        
        S(SySxBij,:)=2.*( S0(SySyBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SyB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SyB(ii,:),:)),Nb,1)...
                        - S0(SxSxBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) - S0(SxB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SxB(jj,:),:)),Nb,1) )...
                        -temp.*repmat(S0(sxsy(jj,ii),:),Nb,1);
                        
        S(SySyBij,:)=-2.*( S0(SxSyBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SxB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SyB(jj,:),:)),Nb,1)...
                        + S0(SxSyBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SxB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SyB(ii,:),:)),Nb,1) )...
                        -temp.*repmat(S0(sysy(ii,jj),:),Nb,1);
                        
        
        
    end
end



%%%%% Boson-Boson

for mm=1:Nb
    for nn=mm:Nb
        S(Bp00(mm,nn),:)=0;
        S(Bp10(mm,nn),:)=0;
    end
end

for mm=1:Nb
    for nn=1:(mm-1)
        S(Bp10(mm,nn),:)=0;
    end
end % These variables are trivial, but are useful for the next equations.



%%%%% Spin-Boson

for ii=1:Ns
    for mm=1:Nb
        S(SxB(ii,mm),:)= S0(SyB(ii,mm),:).*( Fe(ii,:)*S0(bp,:) ) + S0(SyB(ii,mm),:).*( Fec(ii,:)*conj(S0(bp,:)) )...
                        - 1i.*0.5.*Fec(ii,mm).*S0(sz(ii),:).*S0(sx(ii),:);
                    
        S(SyB(ii,mm),:)= - S0(SxB(ii,mm),:).*( Fe(ii,:)*S0(bp,:) ) - S0(SxB(ii,mm),:).*( Fec(ii,:)*conj(S0(bp,:)) )...
                        - 1i.*0.5.*Fec(ii,mm).*S0(sz(ii),:).*S0(sy(ii),:);
                    
        S(SzB(ii,mm),:)=0;
    end
end






end















