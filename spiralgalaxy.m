function [] = spiralgalaxy(dr, rmax, A0, A2, wavelength0, wavelength2, patternperiod, mu0withsign, diskscalelength, nparticles, TotalTime, dt)

graphradius = diskscalelength*3;

jsize = 300;

jjsize = 1000;
dotsizes = [2 5 8];                  nfigs = size(dotsizes,2);

stardistributiontype = 1;   
% Enter 0 for random exponential, 1 for smooth exponential, 2 for Gaussian, and 3 for uniform (all with standard deviations of 6*diskscalelength^2.

% Example command line:  spiralgalaxy(1,1000000,1,-.5,2000,1990,50000000,8.7e-13,15000,5000,250000000,20000);

% dr is the radial stepsize, and rmax is the radius at which the dark
% matter is cutoff.

% wavelength0 and wavelength2 are the wavelengths of the zeroth and second
% spherical harmonic functions f0 and f2.
% The first peak in the rotation curve due to dark matter is at 
% 4.5/2pi = 0.7 times wavelength0.

% patternperiod is the time for the rotating potential which results from
% the interference pattern to go around once.

% mu0 is the dark matter density corresponding to the standard scalar wave with amplitude one.

% The particle mass density is proportional to exp(-r/diskscalelength).
% The number of particles is nparticles.

% The stepsize in time in the simulation is dt.

% The sign of mu0 and patternperiod dictate the direction of rotation of
% the stars and dark matter, respectively.  Positive values give
% counterclockwise, and negative values give clockwise.  Generally, they
% should have the same sign (to model stars and dark matter rotating the
% same direction).

mu0 = abs(mu0withsign);

% Computer k0, k2, w0, w2, Omega (which is called Upsilon in Bray's 2010 paper).

k0 = 2*pi / wavelength0;
k2 = 2*pi / wavelength2;

w0 = patternperiod/8/pi * (k2^2 - k0^2) - 2*pi/patternperiod
w2 = patternperiod/8/pi * (k2^2 - k0^2) + 2*pi/patternperiod
dw = w2-w0;

Omega0Squared = (2*pi/patternperiod)^2 + (patternperiod/8/pi)^2 * (k2^2 - k0^2)^2 - (k2^2 + k0^2)/2;

if Omega0Squared < 0
    
    Omega0Squared
    ErrorMessage = 'Omega0Squared is negative!'
    
end

Omega0 = sqrt(Omega0Squared)
OneOverOmega0 = 1/Omega0

% Initialize f0, f2, U0, U2, U4, UU2, W0, W2, W4, and WW2 (which stands for W tilde sub 2).

nsteps = 1 + round(rmax/dr);
dr = rmax / (nsteps - 1);
rvalues = 0:dr:rmax;
graphsteps = min(1+round(graphradius/dr),size(rvalues,2));
Graphradius = dr*(graphsteps-1)

f0 = zeros(1, nsteps);
f2 = zeros(1, nsteps);

U0 = zeros(1, nsteps);
U2 = zeros(1, nsteps);
U4 = zeros(1, nsteps);
UU2 = zeros(1, nsteps);

W0 = zeros(1, nsteps);
W0p = zeros(1, nsteps);

W2 = zeros(1, nsteps);
W2p = zeros(1, nsteps);

W4 = zeros(1, nsteps);
W4p = zeros(1, nsteps);

WW2 = zeros(1, nsteps);
WW2p = zeros(1, nsteps);

% Compute f0 and f2.

f0(1) = 1;
f0p = 0;
f0pp = -k0^2 * f0(1) / 3;

f2(1) = 1;
f2p = 0;
f2pp = -k2^2 * f2(1) / 7;

for i = 2:nsteps
    
    r = (i-1)*dr;
    
    f0(i) = f0(i-1) + f0p * dr + f0pp * dr^2 / 2;
    f0p = f0p + f0pp * dr;
    f0pp = -k0^2 * f0(i) - 2 * f0p / r;
    
    f2(i) = f2(i-1) + f2p * dr + f2pp * dr^2 / 2;
    f2p = f2p + f2pp * dr;
    f2pp = -k2^2 * f2(i) - 6 * f2p / r;
    
end

% Normalize f.

f0 = -f0 / min(f0 .* rvalues);
f2 = -f2 / min(f2 .* rvalues .^ 3);

fac = max(f0);

f0 = f0 / fac;
f2 = f2 / fac;

% Plot f.

fmin = min(min(f0),min(f2));
fmax = max(max(f0),max(f2));

% Setup graphics names for initial plots.

fn = [num2str(dr) '_' num2str(rmax) '_' num2str(A0) '_' num2str(A2) '_' num2str(wavelength0) '_' num2str(wavelength2)];
fn = [fn '_' num2str(patternperiod) '_' num2str(mu0withsign) '_' num2str(diskscalelength) '_' num2str(nparticles) '_' num2str(0) '_' num2str(dt)];
fn = [fn 'P' num2str(OneOverOmega0) 'P'];
fn = nodecimals(fn);

if 1==1
   figure(1)
   plot(rvalues(1:graphsteps), f0(1:graphsteps))
   hold on
   plot(rvalues(1:graphsteps), f2(1:graphsteps) .* rvalues(1:graphsteps) .^ 2)
   hold off
   title('f_0(r) and f_2(r) * r^2')
   axis([0 Graphradius fmin fmax])
   fname = ['f1_' fn '.jpg']     
   saveas(1,fname)
end

% Compute U0, U2, U4, UU2.

for i = 1:nsteps
    
    r = (i-1)*dr;
    
    U0(i) = A0^2 * f0(i)^2 + 56/105 * A2^2 * r^4 * f2(i)^2;
    U2(i) = -40/105 * A2^2 * r^2 * f2(i)^2;
    U4(i) = 3/105 * A2^2 * f2(i)^2;
    UU2(i) = 2*A0*A2 * f0(i) * f2(i);
       
end

% Compute and plot dark matter density in xy plane.

XYDensity = zeros(jsize);
XZDensity = zeros(jsize);
YZDensity = zeros(jsize);

for jrow = 1:jsize
    for jcol = 1:jsize
        
        v1 = (jcol-1)/(jsize-1) * 2 * graphradius - graphradius;
        v2 = graphradius - 2 * graphradius * (jrow-1)/(jsize-1);
        
        XYDensity(jrow,jcol) = DMdensity(v1,v2,0,0) * (v1^2 + v2^2);
        XZDensity(jrow,jcol) = DMdensity(v1,0,v2,0) * (v1^2 + v2^2);
        YZDensity(jrow,jcol) = DMdensity(0,v1,v2,0) * (v1^2 + v2^2);
        
    end
end

if 1==1
    figure(2)
    colormap(gray);
    maximum = max(YZDensity(:));
    imagesc(YZDensity, [0  maximum]);
    title('Dark Matter Density (Times r^2) in the yz Plane at t=0')
    fname = ['f2_' fn '.jpg']     
    saveas(2,fname)
   
    figure(3)
    colormap(gray);
    maximum = max(XZDensity(:));
    imagesc(XZDensity, [0  maximum]);
    title('Dark Matter Density (Times r^2) in the xz Plane at t=0')
    fname = ['f3_' fn '.jpg']     
    saveas(3,fname)
end

figure(4)
colormap(gray);
maximum = max(XYDensity(:));
imagesc(XYDensity, [0  maximum]);
title('Dark Matter Density (Times r^2) in the xy Plane at t=0')
fname = ['f4_' fn '.jpg']     
saveas(4,fname)

% Compute W0, W2, W4, and WW2. 

W0(1) = 0;
W0p(1) = 0;
W0pp = U0(1) / 3;

W2(1) = 0;
W2p(1) = 0;
W2pp = U2(1) / 7;

W4(1) = 0;
W4p(1) = 0;
W4pp = U4(1) / 11;

WW2(1) = 0;
WW2p(1) = 0;
WW2pp = UU2(1) / 7;

for i = 2:nsteps
    
    r = (i-1)*dr;
    
    W0(i) = W0(i-1) + W0p(i-1) * dr + W0pp * dr^2 / 2;
    W0p(i) = W0p(i-1) + W0pp * dr;
    W0pp = U0(i) - 2 * W0p(i) / r;
    
    W2(i) = W2(i-1) + W2p(i-1) * dr + W2pp * dr^2 / 2;
    W2p(i) = W2p(i-1) + W2pp * dr;
    W2pp = U2(i) - 6 * W2p(i) / r;
    
    W4(i) = W4(i-1) + W4p(i-1) * dr + W4pp * dr^2 / 2;
    W4p(i) = W4p(i-1) + W4pp * dr;
    W4pp = U4(i) - 10 * W4p(i) / r;
    
    WW2(i) = WW2(i-1) + WW2p(i-1) * dr + WW2pp * dr^2 / 2;
    WW2p(i) = WW2p(i-1) + WW2pp * dr;
    WW2pp = UU2(i) - 6 * WW2p(i) / r;
    
end

% Add the correct constants to the W so that they go to zero at infinity.

W0shift = - rmax / 1 * W0p(nsteps) - W0(nsteps);
W2shift = - rmax / 5 * W2p(nsteps) - W2(nsteps);
W4shift = - rmax / 9 * W4p(nsteps) - W4(nsteps);
WW2shift = - rmax / 5 * WW2p(nsteps) - WW2(nsteps);

W0 = W0 + W0shift;
W2 = W2 + W2shift;
W4 = W4 + W4shift;
WW2 = WW2 + WW2shift;

% Plot W functions.

if 1==1
    figure(5)
    plot(rvalues(1:graphsteps), W0(1:graphsteps))
    title('W0')
    fname = ['f5_' fn '.jpg']     
    saveas(5,fname)
   
    figure(6)
    plot(rvalues(1:graphsteps), -rvalues(1:graphsteps) .^ 2 .* W2(1:graphsteps))
    title('-r^2 * W2')
    fname = ['f6_' fn '.jpg']     
    saveas(6,fname)
   
    figure(7)
    plot(rvalues(1:graphsteps), 3 * W4(1:graphsteps) .* rvalues(1:graphsteps) .^ 4)
    title('3r^4 * W4')
    fname = ['f7_' fn '.jpg']     
    saveas(7,fname)
   
    figure(8)
    plot(rvalues(1:graphsteps), WW2(1:graphsteps) .* rvalues(1:graphsteps) .^ 2)
    title('r^2 * WW2')
    fname = ['f8_' fn '.jpg']     
    saveas(8,fname)
   
end
    
% Create 3D graphs of cross sections of the potential function V.

XYPotential = zeros(jsize);
XZPotential = zeros(jsize);
YZPotential = zeros(jsize);

for jrow = 1:jsize
    for jcol = 1:jsize
        
        v1 = (jcol-1)/(jsize-1) * 2 * graphradius - graphradius;
        v2 = (jrow-1)/(jsize-1) * 2 * graphradius - graphradius;
        
        XYPotential(jrow,jcol) = Vfunct(v1,v2,0,0); 
        XZPotential(jrow,jcol) = Vfunct(v1,0,v2,0);
        YZPotential(jrow,jcol) = Vfunct(0,v1,v2,0);

    end
end

figure(9)
h = surf(XYPotential);
colormap('lines')
shading interp
title('The Potential Function in the xy Plane')
xlabel('x axis')
ylabel('y axis')
zlabel('Potential')
fname = ['f9_' fn '.jpg']     
saveas(9,fname)
   
if 1==1
    figure(10)
    h = surf(XZPotential);
    colormap('lines')
    shading interp
    title('The Potential Function in the xz Plane')
    xlabel('x axis')
    ylabel('z axis')
    zlabel('Potential')
    fname = ['f10_' fn '.jpg']     
    saveas(10,fname)

    figure(11)
    h = surf(YZPotential);
    colormap('lines')
    shading interp
    title('The Potential Function in the yz Plane')
    xlabel('y axis')
    ylabel('z axis')
    zlabel('Potential')
    fname = ['f11_' fn '.jpg']     
    saveas(11,fname)
end

% Plot the negative of the gradient vector field of the potential V.

nsize = 20;

xarray = [];
yarray = [];
vxarray = [];
vyarray = [];

for i = 1:nsize
    for j = 1:nsize
         
        x = (i-1)/(nsize-1)*2*graphradius - graphradius;
        y = (j-1)/(nsize-1)*2*graphradius - graphradius;
        z = 0;

        accelvec = -GradVfunct(x,y,z,0);

        vx = accelvec(1);
        vy = accelvec(2);
        vz = accelvec(3);
 
        xarray = [xarray x];
        yarray = [yarray y];
               
        vxarray = [vxarray vx];
        vyarray = [vyarray vy];    
      
    end
end

if 1==1
    figure(12)
    quiver(xarray, yarray, vxarray, vyarray)
    title('The Acceleration Vector Field in the xy Plane')
    xlabel('x axis')
    ylabel('y axis')
end

% Compute rotation curves relative to xaxis, yaxis, and x=y.

rrr = [];
rcxaxis = [];
rcyaxis = [];
rcxey = [];
isteps = 2000;

for i = 1:isteps
    
    r = (i-1)/(isteps-1)*graphradius;
    
    rrr = [rrr r];
    
    GradV = GradVfunct(r*sqrt(1/2), r*sqrt(1/2), 0, 0);  % GradV at radius r along x=y where potential is roughly average. 
    normGradV = sqrt(GradV(1)^2 + GradV(2)^2);
    velxey = sqrt(r * normGradV);  
    rcxey = [rcxey velxey];
    
    GradV = GradVfunct(r, 0, 0, 0);  % GradV at radius r along x axis. 
    normGradV = sqrt(GradV(1)^2 + GradV(2)^2);
    velxaxis = sqrt(r * normGradV);  
    rcxaxis = [rcxaxis velxaxis];
    
    GradV = GradVfunct(0, r, 0, 0);  % GradV at radius r along y axis.
    normGradV = sqrt(GradV(1)^2 + GradV(2)^2);
    velyaxis = sqrt(r * normGradV); 
    rcyaxis = [rcyaxis velyaxis];

end

if 1==1
    figure(13)
    plot(rrr, rcxaxis, 'b')
    hold on
    plot(rrr, rcyaxis, 'r')
    plot(rrr, rcxey, 'g')
    hold off
    title('Rotation Curves Relative to x Axis (Blue), y Axis (Red), and x=y (Green)')
    drawnow
    fname = ['f13_' fn '.jpg']     
    saveas(13,fname)
   
end
    
%----------------Begin Simulation of Galaxy-------------------------------------------

m = ones(1,nparticles);
x = zeros(1,nparticles);
y = zeros(1,nparticles);
z = zeros(1,nparticles);
vx = zeros(1,nparticles);
vy = zeros(1,nparticles);
vz = zeros(1,nparticles);

% Determine initial positions and velocities of particles.

randomphi = [0];
kount = 1;
kround = 0;

while kount < nparticles

   kround = kround + 1;
   for i = 1:kount
       randomphi(kount + i) = randomphi(i) + 1/2^kround;
   end
   kount = kount * 2;
    
end

randomphi = randomphi * (2*pi);

randomx = 0;
randomr = 0;
randomA = 0;

for i = 1:nparticles,
    
    if stardistributiontype == 0   % Random exponential decay like type 1, but not as smooth.
        
        done = 0;
        while done == 0; 
            randomr = 20*diskscalelength * rand^(1/2);  % Left unmodified, this would produce a uniform density of particles in the disk out to 20 times the disk scale length.
            density = exp(-randomr/diskscalelength);    % This is the density we are actually trying to achieve (instead of the uniform density 1).  See Binney and Tremaine 2008, page 12.
            if rand < density
                done = 1;
            end 
        end     
    
    elseif stardistributiontype == 1   % Exponential decay for star density, exp(-r/diskscalelength).  Other distributions are scaled to have the same standard deviation of 6*diskscalelength^2.
    
        randomxtarget = i/(nparticles+1);
        deltarandomx = randomxtarget - randomx;
        deltarandomA = deltarandomx * 2 * diskscalelength^2 * exp(randomr/diskscalelength);
        randomA = randomA + deltarandomA;
        randomr = sqrt(randomA);
        randomx = 1 - (1 + randomr/diskscalelength)/exp(randomr/diskscalelength); 
    
    elseif stardistributiontype == 2   % Gaussian star density, exp(-(r/sqrt(6)/diskscalelength)^2).
        
        randomr = sqrt(6)*diskscalelength*sqrt(log(1/(1 - i/(nparticles+1))));
        
    elseif stardistributiontype == 3   % Uniform star density to radius sqrt(12)*diskscalelength.
        
        randomr = sqrt(12)*diskscalelength*sqrt(i/nparticles);
        
    end
    
    phi = randomphi(i);
    
    theta = pi/2;
    
    x(i) = randomr*sin(theta)*cos(phi);
    y(i) = randomr*sin(theta)*sin(phi);
    z(i) = randomr*cos(theta);
   
    GradV = GradVfunct(randomr*sqrt(1/2), randomr*sqrt(1/2), 0, 0);  % GradV at radius randomr along x=y where potential is roughly average. 
    normGradV = sqrt(GradV(1)^2 + GradV(2)^2);
        
    vel = sqrt(randomr * normGradV);  
    % since v^2 = ra yields circular motion.  Since our potential is not spherically symmetric and the particles are not perfectly evenly distributed, this is an approximation.
    
    vx(i) = -vel*sin(phi)*sign(mu0withsign);
    vy(i) =  vel*cos(phi)*sign(mu0withsign);
    vz(i) = 0;
    
end

ax = zeros(1,nparticles);
ay = zeros(1,nparticles);
az = zeros(1,nparticles);

if 1==1
    figure(14)
    quiver3(x,y,z,vx,vy,vz)
    axis([-graphradius graphradius -graphradius graphradius -graphradius graphradius]);
    title('Initial positions and velocities')
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    drawnow
end

Maxloops = round(TotalTime/dt);

% Begin looping ---------------------------------------

for nloops = 0:Maxloops
    
    t = nloops*dt;
    
    % Define Kinetic Energy, Particle-Particle Potential Energy, 
    % Dark Matter-Particle Potential Energy, and Angular Momentum.
    
    KE = 0;
    DPE = 0;
    ANG = [0 0 0];
     
    % Note that KE + DPE and ANG would be conserved quantities if the dark
    % matter potential were static and spherically symmetric.  Since this
    % is not the case, their lack of conservation represents transfers of
    % energy and angular momentum to the test particles.
    
    for j = 1:nparticles
        
        KE = KE + (1/2)*m(j)*(vx(j)^2 + vy(j)^2 + vz(j)^2);        
        DPE = DPE + m(j)*Vfunct(x(j),y(j),z(j),t);
        ANG = ANG + m(j)*cross([x(j) y(j) z(j)], [vx(j) vy(j) vz(j)]);
        
        accel = -GradVfunct(x(j),y(j),z(j),t);
        ax(j) = accel(1);
        ay(j) = accel(2);
        az(j) = accel(3);
   
    end
    
    if nloops == 0;
        ANG0 = ANG;
        KE0 = KE;
        DPE0 = DPE;
        TE0 = KE+DPE;
    end 
    
    x = x + vx*dt + ax*dt*dt/2;
    y = y + vy*dt + ay*dt*dt/2;
    z = z + vz*dt + az*dt*dt/2;
    
    vx = vx + ax*dt;
    vy = vy + ay*dt;   
    vz = vz + az*dt;
    
    if nloops == round(nloops/10)*10
        
        % graphradius = max([graphradius, max(abs(x)), max(abs(y)), max(abs(z))]);
        
        figure(15)
        quiver3(x,y,z,vx,vy,vz)
        %scatter3(x,y,z,1,'filled')
        angle = dw*t/2;
        xbar = graphradius*cos(angle);
        ybar = graphradius*sin(angle);
        hold on
        quiver3(0,0,0,xbar, ybar, 0)
        hold off
        axis([-graphradius graphradius -graphradius graphradius -graphradius graphradius]);
        title(sprintf('%g %g %g %g %g %g', [dr rmax A0 A2 wavelength0 wavelength2 ]))
        xlabel(sprintf('%g %g %g %g %g %g', [patternperiod mu0withsign diskscalelength nparticles t dt]))
        ylabel(sprintf('%g',[OneOverOmega0]))
        view([0,0,1])
        set(gca,'CameraUpVector',[0 1 0])       
        drawnow
        
        fn = [num2str(dr) '_' num2str(rmax) '_' num2str(A0) '_' num2str(A2) '_' num2str(wavelength0) '_' num2str(wavelength2)];
        fn = [fn '_' num2str(patternperiod) '_' num2str(mu0withsign) '_' num2str(diskscalelength) '_' num2str(nparticles) '_' num2str(t) '_' num2str(dt)];
        fn = [fn 'P' num2str(OneOverOmega0) 'P'];
        fn = nodecimals(fn);

        fname = ['f15_' fn '.jpg']     
        saveas(15,fname)
        
        for k = 1:nfigs
            
            [dd,maxdd] = spiraldensity(jjsize, dotsizes(k)); 
            dd = dd/maxdd;
            
            figure(20 + k)
            colormap(gray);
            imagesc(dd, [0  1]);
            title(sprintf('%g %g %g %g %g %g', [dr rmax A0 A2 wavelength0 wavelength2 ]))
            xlabel(sprintf('%g %g %g %g %g %g', [patternperiod mu0withsign diskscalelength nparticles t dt]))
            ylabel(sprintf('%g %g %g',[OneOverOmega0 jjsize dotsizes(k)]))
            drawnow
                
            if k==2        
                    
                fn = [num2str(dr) '_' num2str(rmax) '_' num2str(A0) '_' num2str(A2) '_' num2str(wavelength0) '_' num2str(wavelength2)];
                fn = [fn '_' num2str(patternperiod) '_' num2str(mu0withsign) '_' num2str(diskscalelength) '_' num2str(nparticles) '_' num2str(t) '_' num2str(dt)];
                fn = [fn 'P' num2str(OneOverOmega0) '_' num2str(jjsize) '_' num2str(dotsizes(k)) 'P'];
                fn = nodecimals(fn);

                fname = [fn '.jpg']
                imwrite(dd, fname, 'jpg')
                
                fnameinvert = [fn 'inv.jpg']   
                imwrite(dd', fnameinvert, 'jpg')
              
            end
           
        end
        
        TE = KE + DPE;
        
        nloops
        ElapseTime = nloops*dt
        ANG0
        ANG
        TE0
        TE
        
    end
    
end

Omega0
OneOverOmega0





%-----------------------END OF PROGRAM-------------------------------






% Functions

    function [density, maxdensity] = spiraldensity(jjsize, dotradius)
        
       
        density = zeros(jjsize);
        maxdensity = 0;

        for i = 1:nparticles
    
            jcol0 = round((x(i)+graphradius)/2/graphradius*(jjsize-1) + 1);
            jrow0 = (jjsize+1) - round((y(i)+graphradius)/2/graphradius*(jjsize-1) + 1);

            for djcol = -dotradius:dotradius
            for djrow = -dotradius:dotradius
        
                if djcol^2 + djrow^2 <= dotradius^2
            
                    jcol = jcol0 + djcol;
                    jrow = jrow0 + djrow;
                    if ((jcol > 0) & (jcol <= jjsize) & (jrow > 0) & (jrow <= jjsize))
                
                    density(jrow,jcol) = density(jrow,jcol)+1;
                    maxdensity = max(maxdensity,density(jrow,jcol));
                
                    end
            
                end
        
            end
            end
    
        end

    end

    function output = DMdensity(x,y,z,t)
        
        r = sqrt(x^2 + y^2 + z^2);
        
        U0value = 0;
        U2value = 0;
        U4value = 0;
        UU2value = 0;

        if r <= rmax
            ival = 1 + round(r/dr);
            U0value = U0(ival);
            U2value = U2(ival);
            U4value = U4(ival);
            UU2value = UU2(ival);
        end
       
        temp = U0value + U2value * (3*z^2 - r^2) + U4value * (35*z^4 - 30*r^2*z^2 + 3*r^4);
        temp2 = temp + UU2value * (cos(dw*t) * (x^2 - y^2) + sin(dw*t) * (2*x*y));
        output = mu0 * temp2;
        
    end

    function output = Vfunct(x,y,z,t)
        
        r = sqrt(x^2 + y^2 + z^2);
        
        if r <= rmax
            ival = 1 + round(r/dr);
            W0value = W0(ival);
            W2value = W2(ival);
            W4value = W4(ival);
            WW2value = WW2(ival);
        else
            W0value = W0(nsteps) * (rmax/r) ^ 1;
            W2value = W2(nsteps) * (rmax/r) ^ 5;
            W4value = W4(nsteps) * (rmax/r) ^ 9;
            WW2value = WW2(nsteps) * (rmax/r) ^ 5;
        end
        
        temp = W0value + W2value*(3*z^2 - r^2) + W4value*(35*z^4 - 30*r^2*z^2 + 3*r^4);
        temp2 = temp + WW2value * (cos(dw*t) * (x^2 - y^2) + sin(dw*t) * (2*x*y));
        output = 4*pi*mu0 * temp2;
        
    end

    function output = GradVfunct(x,y,z,t)
        
        r = sqrt(x^2 + y^2 + z^2);
        if r == 0
            output = [0 0 0];
        else
        
            if r <= rmax
                ival = 1 + round(r/dr);
            
                W0value = W0(ival);
                W2value = W2(ival);
                W4value = W4(ival);
                WW2value = WW2(ival);
            
                W0pvalue = W0p(ival);
                W2pvalue = W2p(ival);
                W4pvalue = W4p(ival);
                WW2pvalue = WW2p(ival);
            
            else
                W0value = W0(nsteps) * (rmax/r) ^ 1;
                W2value = W2(nsteps) * (rmax/r) ^ 5;
                W4value = W4(nsteps) * (rmax/r) ^ 9;
                WW2value = WW2(nsteps) * (rmax/r) ^ 5;
            
                W0pvalue = -W0value / r * 1;
                W2pvalue = -W2value / r * 5;
                W4pvalue = -W4value / r * 9;
                WW2pvalue = -WW2value / r * 5;   
            
            end
        
            term1 = W0pvalue + W2pvalue * (3*z^2 - r^2) + W4pvalue * (35*z^4 - 30*r^2*z^2 + 3*r^4);
            term2 = WW2pvalue * (cos(dw*t) * (x^2 - y^2) + sin(dw*t) * (2*x*y));
            vector1 = (term1 + term2) * [x/r, y/r, z/r];
        
            vector2 = W2value * [-2*x, -2*y, 4*z];
        
            vector3 = W4value * [12*x*(r^2 - 5*z^2), 12*y*(r^2 - 5*z^2), 16*z*(5*z^2 - 3*r^2)];
            
            vector4 = cos(dw*t) * WW2value * [2*x, -2*y, 0];
            
            vector5 = sin(dw*t) * WW2value * [2*y, 2*x, 0];
        
            output = 4*pi*mu0*(vector1 + vector2 + vector3 + vector4 + vector5);
            
        end
        
    end 

    function [goodfilename] = nodecimals(badfilename)
        
        stringsize = size(badfilename);
        stringlength = stringsize(2);
        
        goodfilename = [];
        for k = 1:stringlength
            badchar = badfilename(k);
            if badchar == '.'
                goodchar = 'p';
            else
                goodchar = badchar;
            end
            goodfilename = [goodfilename goodchar];
        end
           
    end

end