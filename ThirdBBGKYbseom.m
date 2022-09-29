function S = ThirdBBGKYbseom(S0,Ns,Nb,ome,muR,F,t,iter,indices,indices3,indicesfixSind)
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


SxSxBfixSind = indicesfixSind(1:Ns,1:Ns);
SxSyBfixSind = indicesfixSind((Ns+1):(2*Ns),1:Ns);
SxSzBfixSind = indicesfixSind((2*Ns+1):(3*Ns),1:Ns);

SySxBfixSind = indicesfixSind((3*Ns+1):(4*Ns),1:Ns);
SySyBfixSind = indicesfixSind((4*Ns+1):(5*Ns),1:Ns);
SySzBfixSind = indicesfixSind((5*Ns+1):(6*Ns),1:Ns);  % Indices of 3-point function for fixed (i,j) and all mu

SzSxBfixSind = indicesfixSind((6*Ns+1):(7*Ns),1:Ns);
SzSyBfixSind = indicesfixSind((7*Ns+1):(8*Ns),1:Ns);
SzSzBfixSind = indicesfixSind((8*Ns+1):(9*Ns),1:Ns);





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
        
        
        SxSxBij = SxSxBfixSind{ii,jj};  SxSyBij = SxSyBfixSind{ii,jj};  SxSzBij = SxSzBfixSind{ii,jj};
        SySxBij = SySxBfixSind{ii,jj};  SySyBij = SySyBfixSind{ii,jj};  SySzBij = SySzBfixSind{ii,jj};  % Indices of 3-point function for fixed (i,j) and all mu
        SzSxBij = SzSxBfixSind{ii,jj};  SzSyBij = SzSyBfixSind{ii,jj};  SzSzBij = SzSzBfixSind{ii,jj};
        
        SxSxBji = SxSxBfixSind{jj,ii};  SxSyBji = SxSyBfixSind{jj,ii};  SxSzBji = SxSzBfixSind{jj,ii};
        SySxBji = SySxBfixSind{jj,ii};  SySyBji = SySyBfixSind{jj,ii};  SySzBji = SySzBfixSind{jj,ii};  % Indices of 3-point function for fixed (j,i) and all mu
        SzSxBji = SzSxBfixSind{jj,ii};  SzSyBji = SzSyBfixSind{jj,ii};  SzSzBji = SzSzBfixSind{jj,ii};
        
        
        % Spin-Spin Terms
        
        % X terms
        
        S(sxsx(ii,jj),:)=2.*( real( Fe(ii,:)*S0(SySxBij,:) + S0(sysx(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(SxB(jj,:),:) ) )...
                            + real( Fe(jj,:)*S0(SySxBji,:) + S0(sysx(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sy(jj),:).*( Fe(jj,:)*S0(SxB(ii,:),:) ) ) );
                        
        S(sxsy(ii,jj),:)=2.*( real( Fe(ii,:)*S0(SySyBij,:) + S0(sysy(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(SyB(jj,:),:) ) )...
                            - real( Fe(jj,:)*S0(SxSxBji,:) + S0(sxsx(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sx(jj),:).*( Fe(jj,:)*S0(SxB(ii,:),:) ) ) );
        
        S(sxsz(ii,jj),:)=2.*( real( Fe(ii,:)*S0(SySzBij,:) + S0(sysz(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(SzB(jj,:),:) ) ) );
                        
                        
                        
        % Y terms
                        
        S(sysx(ii,jj),:)=2.*( real( Fe(jj,:)*S0(SySyBji,:) + S0(sysy(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sy(jj),:).*( Fe(jj,:)*S0(SyB(ii,:),:) ) )...
                            - real( Fe(ii,:)*S0(SxSxBij,:) + S0(sxsx(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sx(ii),:).*( Fe(ii,:)*S0(SxB(jj,:),:) ) ) );
                        
        S(sysy(ii,jj),:)=-2.*( real( Fe(ii,:)*S0(SxSyBij,:) + S0(sxsy(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sx(ii),:).*( Fe(ii,:)*S0(SyB(jj,:),:) ) )...
                             + real( Fe(jj,:)*S0(SxSyBji,:) + S0(sxsy(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sx(jj),:).*( Fe(jj,:)*S0(SyB(ii,:),:) ) ) );
                        
        S(sysz(ii,jj),:)=-2.*( real( Fe(ii,:)*S0(SxSzBij,:) + S0(sxsz(ii,jj),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sx(ii),:).*( Fe(ii,:)*S0(SzB(jj,:),:) ) ) );
                        
                        
                          
        % Z terms
                        
        S(szsx(ii,jj),:)=2.*( real( Fe(jj,:)*S0(SySzBji,:) + S0(sysz(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sy(jj),:).*( Fe(jj,:)*S0(SzB(ii,:),:) ) ) );
                        
        S(szsy(ii,jj),:)=-2.*( real( Fe(jj,:)*S0(SxSzBji,:) + S0(sxsz(jj,ii),:).*( Fe(jj,:)*S0(bp,:) ) + S0(sx(jj),:).*( Fe(jj,:)*S0(SzB(ii,:),:) ) ) );
                        
        S(szsz(ii,jj),:)=0;
                        
        
        
        
        % Spin-Spin-Boson Terms
        
        temp = 0.5*1i*repmat(conj(edet),1,iter).*( repmat(F(ii,:).',1,iter).*repmat(S0(sz(ii),:),Nb,1) ...
                                                 + repmat(F(jj,:).',1,iter).*repmat(S0(sz(jj),:),Nb,1) );
        Atempi = zeros(Nb,iter);
        Atempj = zeros(Nb,iter);
        for mm=1:Nb
           Atempi(mm,:) = Fe(ii,:)*S0(Bp00(mm,:),:) + Fec(ii,:)*S0(Bp10(:,mm),:);
           Atempj(mm,:) = Fe(jj,:)*S0(Bp00(mm,:),:) + Fec(jj,:)*S0(Bp10(:,mm),:);
        end
        
        % X terms
        
        S(SxSxBij,:)=2.*( S0(SySxBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SyB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SxB(jj,:),:)),Nb,1)...
                        + S0(SySxBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SyB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SxB(ii,:),:)),Nb,1) )...
                     + Atempi.*repmat(S0(sysx(ii,jj),:),Nb,1) + Atempj.*repmat(S0(sysx(jj,ii),:),Nb,1)...
                     - temp.*repmat(S0(sxsx(ii,jj),:),Nb,1)...
                     - 0.5*1i*( ( Fec(jj,:).' )*( S0(sxsz(ii,jj),:).*S0(sx(jj),:) ) + ( Fec(ii,:).' )*( S0(sxsz(jj,ii),:).*S0(sx(ii),:) ) );
                        
        S(SxSyBij,:)=2.*( S0(SySyBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SyB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SyB(jj,:),:)),Nb,1)...
                        - S0(SxSxBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) - S0(SxB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SxB(ii,:),:)),Nb,1) )...
                     + Atempi.*repmat(S0(sysy(ii,jj),:),Nb,1) - Atempj.*repmat(S0(sxsx(jj,ii),:),Nb,1)...
                     - temp.*repmat(S0(sxsy(ii,jj),:),Nb,1)...
                     - 0.5*1i*( ( Fec(jj,:).' )*( S0(sxsz(ii,jj),:).*S0(sy(jj),:) ) + ( Fec(ii,:).' )*( S0(sysz(jj,ii),:).*S0(sx(ii),:) ) );
                 
        S(SxSzBij,:)=2.*( S0(SySzBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SyB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SzB(jj,:),:)),Nb,1) )...
                     + Atempi.*repmat(S0(sysz(ii,jj),:),Nb,1)...
                     - temp.*repmat(S0(sxsz(ii,jj),:),Nb,1)...
                     - 0.5*1i*( ( Fec(jj,:).' )*( S0(sxsz(ii,jj),:).*S0(sz(jj),:) ) );
                        
                        
        % Y terms
                        
        S(SySxBij,:)=2.*( S0(SySyBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SyB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SyB(ii,:),:)),Nb,1)...
                        - S0(SxSxBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) - S0(SxB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SxB(jj,:),:)),Nb,1) )...
                     + Atempj.*repmat(S0(sysy(jj,ii),:),Nb,1) - Atempi.*repmat(S0(sxsx(ii,jj),:),Nb,1)...
                     - temp.*repmat(S0(sxsy(jj,ii),:),Nb,1)...
                     - 0.5*1i*( ( Fec(ii,:).' )*( S0(sxsz(jj,ii),:).*S0(sy(ii),:) ) + ( Fec(jj,:).' )*( S0(sysz(ii,jj),:).*S0(sx(jj),:) ) );
                        
        S(SySyBij,:)=-2.*( S0(SxSyBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SxB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SyB(jj,:),:)),Nb,1)...
                        + S0(SxSyBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SxB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SyB(ii,:),:)),Nb,1) )...
                     - Atempi.*repmat(S0(sxsy(ii,jj),:),Nb,1) - Atempj.*repmat(S0(sxsy(jj,ii),:),Nb,1)...
                     - temp.*repmat(S0(sysy(ii,jj),:),Nb,1)...
                     - 0.5*1i*( ( Fec(jj,:).' )*( S0(sysz(ii,jj),:).*S0(sy(jj),:) ) + ( Fec(ii,:).' )*( S0(sysz(jj,ii),:).*S0(sy(ii),:) ) );
                    
        S(SySzBij,:)=-2.*( S0(SxSzBij,:).*repmat(real(Fe(ii,:)*S0(bp,:)),Nb,1) + S0(SxB(ii,:),:).*repmat(real(Fe(ii,:)*S0(SzB(jj,:),:)),Nb,1) )...
                     - Atempi.*repmat(S0(sxsz(ii,jj),:),Nb,1)...
                     - temp.*repmat(S0(sysz(ii,jj),:),Nb,1)...
                     - 0.5*1i*( ( Fec(jj,:).' )*( S0(sysz(ii,jj),:).*S0(sz(jj),:) ) );
                 
                 
        % Z terms
        
        S(SzSxBij,:)=2.*( S0(SySzBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SyB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SzB(ii,:),:)),Nb,1) )...
                     + Atempj.*repmat(S0(sysz(jj,ii),:),Nb,1)...
                     - temp.*repmat(S0(sxsz(jj,ii),:),Nb,1)...
                     - 0.5*1i*( ( Fec(ii,:).' )*( S0(sxsz(jj,ii),:).*S0(sz(ii),:) ) );
                 
        S(SzSyBij,:)=-2.*( S0(SxSzBji,:).*repmat(real(Fe(jj,:)*S0(bp,:)),Nb,1) + S0(SxB(jj,:),:).*repmat(real(Fe(jj,:)*S0(SzB(ii,:),:)),Nb,1) )...
                     - Atempj.*repmat(S0(sxsz(jj,ii),:),Nb,1)...
                     - temp.*repmat(S0(sysz(jj,ii),:),Nb,1)...
                     - 0.5*1i*( ( Fec(ii,:).' )*( S0(sysz(jj,ii),:).*S0(sz(ii),:) ) );
                 
        S(SzSzBij,:)=0;
                        
        
        
    end
end



%%%%% Boson-Boson

for mm=1:Nb
    for nn=mm:Nb
        S(Bp00(mm,nn),:)=1i.*0.5.*( (Fec(:,mm).')*S0(SzB(:,nn),:)...
                                  + (Fec(:,nn).')*S0(SzB(:,mm),:) );
        S(Bp10(mm,nn),:)=-1i.*0.5.*( (Fe(:,mm).')*S0(SzB(:,nn),:)...
                                   - (Fec(:,nn).')*conj(S0(SzB(:,mm),:)) );
    end
end

for mm=1:Nb
    for nn=1:(mm-1)
        S(Bp10(mm,nn),:)=conj(S(Bp10(nn,mm)));
    end
end % These variables are trivial, but are useful for the next equations.



%%%%% Spin-Boson

for ii=1:Ns
    for mm=1:Nb        
        S(SxB(ii,mm),:)= S0(SyB(ii,mm),:).*( Fe(ii,:)*S0(bp,:) ) + S0(sy(ii),:).*( Fe(ii,:)*S0(Bp00(:,mm),:) )...
                        + S0(SyB(ii,mm),:).*( Fec(ii,:)*conj(S0(bp,:)) ) + S0(sy(ii),:).*( Fec(ii,:)*S0(Bp10(:,mm),:) )...
                        + 1i.*0.5.*( (Fec(:,mm).')*S0(sxsz(ii,:),:) - Fec(ii,mm).*S0(sz(ii),:).*S0(sx(ii),:) );
                    
        S(SyB(ii,mm),:)= - S0(SxB(ii,mm),:).*( Fe(ii,:)*S0(bp,:) ) - S0(sx(ii),:).*( Fe(ii,:)*S0(Bp00(:,mm),:) )...
                        - S0(SxB(ii,mm),:).*( Fec(ii,:)*conj(S0(bp,:)) ) - S0(sx(ii),:).*( Fec(ii,:)*S0(Bp10(:,mm),:) )...
                        + 1i.*0.5.*( (Fec(:,mm).')*S0(sysz(ii,:),:) - Fec(ii,mm).*S0(sz(ii),:).*S0(sy(ii),:) );
                    
        S(SzB(ii,mm),:)= 1i.*0.5.*( Fec(ii,mm).*(1-S0(sz(ii),:).^2) );
    end
end






end















