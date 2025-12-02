% This function solves the initial-boundary value problem for 
% the Burgers equation in a periodic domain using the
% pseudo-spectral Galerkin approach with and without dealiasing


function burgers_01
  
echo off, clc
clear all

global Nsave; global NSTART; global NADVANCE;
Nsave    = 25;
NSTART   = 0;
NADVANCE = 500;
DEALIAS  = 1;

global nu; nu = 2.5e-2;

N2 = 128;
dt = 1.0e-2;

% Generating the initial condition (a simple sine wave);
dx  = 2.0 * pi / N2;
X   = (0.0:dx:(N2-1)*dx);
uIC = sin( X );

u = Solve_Burgers( uIC, N2, dt );

figure(1); surf(X,(1:NADVANCE/Nsave)*dt,u);
xlabel('x'); ylabel('time');

figure(2); contourf(X,(1:NADVANCE/Nsave)*dt,u,10);
xlabel('x'); ylabel('time'); colorbar;


% The function solves the Burgers equation;
% Time-stepping is performed using first order explicit/implicit 
% Euler scheme
function [usave] = Solve_Burgers( u0, N2, dt )
% INPUT
% u0  - initial condition (real space)
% N2  - number of grid point (extended grid)
% dt  - time-step;
% OUTPUT
% u   - solution (in real space) saved every Nsave steps
 
global DEALIAS; global Nsave; global NSTART; global NADVANCE;
global nu;

% Determining the number of grid points
if ( DEALIAS )
  N = round(2 * N2 / 3);
else
  N = N2;
end;
if (mod(N,2) ~= 0); N = N - 1; end; 

% Assembling multiplication vectors;
Kx(1:N2/2+1)   = (0:N2/2) * j * dt;
Kx2(1:N2/2+1)  = ((0:N2/2)  .* (0:N2/2)) * nu * dt;

% Transforming the initial condition to Fourier space;
u0bar = Transform2Fourier( u0, N, N2 );
  
for i=NSTART:NADVANCE;
  
  if ( i == NSTART )
% at the initial time-step;   
    ubar(1:N/2) = u0bar(1:N/2);
    u(1:N2)     = u0(1:N2);    
  else    
% pseudospectral calculation of the nonlinear products;
    u = Transform2Real( ubar, N, N2 );
    uu(1:N2) =  0.5 * (u(1:N2) .* u(1:N2));
    uubar = Transform2Fourier( uu, N, N2 );
    
% Actual integration in time;
    ubar(1:N/2) = ( ubar(1:N/2) - Kx(1:N/2) .* uubar(1:N/2) ) ./ ...
                  ( 1.0 + Kx2(1:N/2) );
  end;

% Saving data if necessary;
  if ( mod( i, Nsave ) )
    u = Transform2Real( ubar, N, N2 );
    usave(floor(i/Nsave)+1, 1:N2 ) = u(1:N2);
  end;

end;
  
return;

% The function solves the Burgers equation;
% Time-stepping is performed using mixed Crank-Nicholson / RK
% scheme;
function [usave] = Solve_Burgers2( u0, N2, dt )
% INPUT
% u0  - initial condition (real space)
% N2  - number of grid point (extended grid)
% dt  - time-step;
% OUTPUT
% u   - solution (in real space) saved every Nsave steps
 
global DEALIAS; global Nsave; global NSTART; global NADVANCE;
global nu;

Gamma  = [8.0/15.0  5.0/12.0    3.0/4.0];
Zeta   = [0.0      -17.0/60.0  -5.0/12.0];
Beta   = [4.0/15.0   1.0/15.0   1.0/6.0];

% Determining the number of grid points
if ( DEALIAS )
  N = round(2 * N2 / 3);
else
  N = N2;
end;
if (mod(N,2) ~= 0); N = N - 1; end; 

% Assembling multiplication vectors;
Kx(1:N2/2+1)   = (0:N2/2) * j;
Kx2(1:N2/2+1)  = - ((0:N2/2)  .* (0:N2/2)) * nu;

% Transforming the initial condition to Fourier space;
u0bar = Transform2Fourier( u0, N, N2 );
  
for i=NSTART:NADVANCE;
  
  if ( i == NSTART )
% at the initial time-step;   
    ubar(1:N/2) = u0bar(1:N/2);
    u(1:N2)     = u0(1:N2);    
  else    
    for k=1:3
% Setting up the RHS vectors for the first substep;
      if ( k > 1 )
        uubar0(1:N/2) = uubar(1:N/2);
        ubar(1:N/2)   = ubar2(1:N/2);
      else
        uubar0(1:N/2) = 0.0;
      end
      
% pseudospectral calculation of the nonlinear products;
      u = Transform2Real( ubar, N, N2 );
      uu(1:N2) =  0.5 * (u(1:N2) .* u(1:N2));
      uubar = Transform2Fourier( uu, N, N2 );
    
% The explicit part;
      uubar(1:N/2) = - Kx(1:N/2) .* uubar(1:N/2);
% Advancing over the substep with the implicit terms;
      ubar2(1:N/2) = (ubar(1:N/2) + dt * (Beta(k) * Kx2(1:N/2).* ubar(1:N/2) ...
	              + Gamma(k)*uubar(1:N/2) + Zeta(k)*uubar0(1:N/2))) ./ ...
                         (1.0 - dt * Beta(k)*Kx2(1:N/2));
    end;    
    ubar(1:N/2)   = ubar2(1:N/2);
  end;

% Saving data if necessary;
  if ( mod( i, Nsave ) )
    u = Transform2Real( ubar, N, N2 );
    usave(floor(i/Nsave)+1, 1:N2 ) = u(1:N2);
  end;

end;
  
return;

%-------------------------------------------------------------
% The function performs the Fourier to real transform;
% It also adds elements with conjugate symmetry and pads the grid;  
function [ur] = Transform2Real (v, n, n2);
% INPUT
% v  - input vector (complex components);
% n  - length on the padded grid;
% n2 - length of the input vector (may be on padded grid);
% OUTPUT
% ur - the transformed real vector;
utmp(1:n/2) = v(1:n/2);
% Padding zeros;
if ( n2 ~= n )
  utmp(n/2+1:n2/2) = 0.0;
% Adding the oddball;
  utmp(n2/2+1) = 0.0;
end;
for i=2:n2/2;
  utmp(n2+2-i) = conj(utmp(i));
end;
ur = real(ifft(utmp,n2));
return;


%-------------------------------------------------------------
% The function carries out the real to Fourier transform and
% removes the aliased elements if necessary;
% The part of the spectrum with conjugate symmetry is also removed;
function [uf] = Transform2Fourier (v,n,n2);
% INPUT
% v - the input real vector;
% n2 - its length (on the padded grid - N2);
% n  - size of the output vector (may be dealiased);
% OUTPUT
% uf - the transformed vector of the Fourier coefficients;
utmp = fft(v,n2);
uf(1:n/2) = utmp(1:n/2);
% Now padding zeros;
if ( n ~= n2 )
  uf(n/2+1:n2/2) = 0.0;
end;
% In case of no dealiasing only the oddball is removed;
return;
