function result = Mie_tetascan(m, x, nsteps)% Computation and plot of Mie Power Scattering function for given % complex refractive-index ratio m=m'+im", size parameters x=k0*a, % Inputs:%   relative complex refractive index m=k1/k=N1/N%       where N1=index of refraction for particle (N1'+iN1") and%       N = refractive index of medium (N=N'+iN'')%   size parameter x=k0*a, where k0=wave number %       in the ambient medium, a=sphere radius;%       =2*pi*N/lambda (where N = refractive index medium)%   nsteps = number of computational steps between 0..pi%%   nmax is calculated as nmax=round(2+x+4*x.^(1/3));%% according to Bohren and Huffman pg 112-113 (1983) BEWI:TDD122% Returns |S1|^2 and |S2|^2 S1 in degress 0-pi then s2 on pi-2pi% C. M�tzler, May 2002.nsteps=nsteps;m1=real(m); m2=imag(m);nx=(1:nsteps); dteta=pi/(nsteps-1);teta=(nx-1).*dteta;    for j = 1:nsteps,         u=cos(teta(j));        a(:,j)=Mie_S12(m,x,u); % Returns [S1; S2]        SL(j)= real(a(1,j)'*a(1,j)); % a' is the complex conjugate of a        SR(j)= real(a(2,j)'*a(2,j));    end;y=[teta teta+pi;SL SR(nsteps:-1:1)]'; % COMMENTED OUT DISPLAY LINES BELOW KTK 13 - 04 - 23%polar(y(:,1),y(:,2))%title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x));%xlabel('Scattering Angle')result=y; 