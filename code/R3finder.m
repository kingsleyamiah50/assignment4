function [R3] = R3finder(minV,maxV, steps)
%R3FINDER Summary of this function goes here
%   Determines the value of R3 to be used in the simulation
global C
global CuCond NoCond
global nx ny

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²

nx = 200;
ny = 100;


Curr = zeros(1,steps);
Vsweep = linspace(minV, maxV, steps);
for u = 1:steps

    Lb = 80;
    Wb = 40;
    
    CuCond = 1;
    NoCond = 10e-9;
    
    %Conductivity map
    
    cMap = zeros(nx,ny);
    
    for i = 1:nx
        for j = 1: ny
            cMap(i,j) = CuCond;
        end
    end
    
    for i = 1:nx
    for j = 1:ny
        if (j>=1 && j<=Wb && i>Lb && i<=(120))
            cMap(i,j) = NoCond;
        end
        
        if (j<=ny && j>=(ny-Wb) && i>Lb && i<=(120))
            cMap(i,j) = NoCond;
        end
    end
    end
    
    G = sparse(nx*ny, nx*ny);
    F = zeros(1, nx*ny);
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i - 1) * ny;
            
            if i == 1
                G(n,:) = 0;
                G(n,n) = 1;
                F(n) = Vsweep(u);
                
            elseif i == nx
                G(n,:) = 0;
                G(n,n) = 1;
                
            elseif j == 1
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nyp = (j+1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
                
                G(n,n) = -(rxm + rxp + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nyp) = ryp;
                
            elseif j == ny
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nym = (j-1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                rym = (cMap(i,j) + cMap(i,j-1))/2.0;
                
                G(n,n) = -(rxm + rxp + rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nym) = rym;
                
            else
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nym = (j-1) + (i-1) * ny;
                nyp = (j+1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                rym = (cMap(i,j) + cMap(i,j-1))/2.0;
                ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
                
                G(n,n) = -(rxm + rxp + rym + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nym) = rym;
                G(n, nyp) = ryp;
            end
        end
    end
    V = G\F';
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1) * ny;
            VG(i,j) = V(n);
        end
    end
    
    for i = 1:nx
        for j = 1:ny
            if i == 1
                Ex(i, j) = (VG(i + 1, j) - VG(i, j));
            elseif i == nx
                Ex(i, j) = (VG(i, j) - VG(i - 1, j));
            else
                Ex(i, j) = (VG(i + 1, j) - VG(i - 1, j)) * 0.5;
            end
            if j == 1
                Ey(i, j) = (VG(i, j + 1) - VG(i, j));
            elseif j == ny
                Ey(i, j) = (VG(i, j) - VG(i, j - 1));
            else
                Ey(i, j) = (VG(i, j + 1) - VG(i, j - 1)) * 0.5;
            end
        end
    end
    
    Ex = -Ex;
    Ey = -Ey;
    
    Jx = cMap .* Ex;
    Jy = cMap .* Ey;
    
    C0 = sum(Jx(1, :));
    Cnx = sum(Jx(nx, :));
    
    %Current
    Curr(u) = (C0 + Cnx) * 0.5;
    
end

%Current vs Voltage
fitlin = polyfit(Vsweep, Curr,1);

R3 = 1/fitlin(1);

end

