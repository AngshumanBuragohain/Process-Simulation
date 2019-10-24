%% Authored by:
%  1. Angshuman Buragohain :    221552
%  2. Devidas Khatri       :    221549

%% Problem 5.6: Colliding particles
%  Two particles at positions x1=(0,0,0)m and x2=(0,2,0)m have diameters=1m. 
%  Particle 1 is stationary and fixed at the origin and particle 2 with a mass 
%  of 0.01kg moves at a velocity of (0,-0.5,0)m/s. Objectives:
%  1. Designing a function which takes positions of two particles and their diameters
%     and returns the overlap (scalar) and the normal (vector)
%  2. Using initial condition of the particles, move the particle 2 with a
%     fixed time step=1e-4s using Leapfrog algorithm (Using predictor -corrector)
%  3. Moving the particle with timestep and evaluating the overlap. kn=1000
%     Computation of the normal force based on the linear model F=-kn*del*n
%  4. Verification of the working of the program using 40,000 time-steps
%  5. Plotting Kinetic energy of particle 2 as a function of time
%  6. Perform the objectives 3 and 4 using coefficient of restitution e=0.8

clear all;
close all;

%% Initializing the variables for the given data about the mass, position and velocity of particles.

   % Position of particles
   x1  = [0, 0, 0];
   x2_i= [0, 2, 0];
   
   % Initial velocity of the particles
   v1   = [0, 0, 0];
   v2_i = [0, -0.5, 0];
   
   % Diameter of the particles
   d1 = 1;
   d2 = 1;
   
   % Initial mass of the particles
   m1   = [];
   m2   = 0.01;
   
   % Elastic constant
   kn_loading = 1000;
   
   F=[0, 0, 0];
   
   % Time step
   dt = 1e-4;
   
   % Number of timesteps
   N = 4e+4;
   
   % Time vector
   t=0:dt:(dt*N);
   
   % Positions of particle 2
   x2 = zeros(N,length(x2_i));
   x2(1,:) = x2_i;
   
   % Velocities of particle 2
   v2 = zeros(N,length(v2_i));
   v2(1,:) = v2_i;
   
   % Kinetic energy matrix initialization to store K.E for e=1
   KE = zeros(1,N);
   KE(1) = 0.5*m2*dot(v2_i,v2_i);
   
   % Kinetic energy matrix initialization to store K.E for e=0.8
   KEe = zeros(1,N);
   KEe(1) = 0.5*m2*dot(v2_i,v2_i);
   
   % Coefficient of restitution
   e_rest=0.8;
   
   % Iteratios for varying time step
   num=4;
   
%% Objective 3 & 4: Moving particle 2 with the given fixed time step dt using the initial conditions given 
%  and calculate normal force and if there is an overlap

   Force{1}=[0 0 0];
   for j=1:N
       [del, n] = Overlap(x1, x2(j,:), d1, d2);  
       if (del<0)
           Force{j+1}=[0 0 0];
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);
       
       elseif (del>0)
           Force{j+1}=-kn_loading*del*n;
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);
       end
       KE(j+1) = 0.5*m2*dot(v2(j+1,:),v2(j+1,:));   
   end

%% Objective 5: Plotting the Kinetic energy of particle 2 with time and varying time step from 10^(-4) to 10^(-5)
   
   figure (1);
   plot(t,KE);
   title(['Kinetic Energy of particle 2 as a function of time with e=1 with time step=',num2str(dt),'s'])
   xlabel("Time [s]")
   ylabel("Kinetic Energy [J]")   
   
    % Time step
   dt = 1e-5;
   
   % Number of timesteps
   N = 4e+5;
   
   % Time vector
   t=0:dt:(dt*N);
   
   Force{1}=[0 0 0];
   for j=1:N
       [del, n] = Overlap(x1, x2(j,:), d1, d2);  
       if (del<0)
           Force{j+1}=[0 0 0];
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);
       
       elseif (del>0)
           Force{j+1}=-kn_loading*del*n;
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);
       end
       KE(j+1) = 0.5*m2*dot(v2(j+1,:),v2(j+1,:));   
   end
   
   figure (2);
   plot(t,KE);
   title(['Kinetic Energy of particle 2 as a function of time with e=1 with time step=',num2str(dt),'s'])
   xlabel("Time [s]")
   ylabel("Kinetic Energy [J]")   

%% Objective 6: Coefficient of resitution, loading and unloading
   Force{1} = [0 0 0];
   for j=1:N
       [del, n] = Overlap(x1, x2(j,:), d1, d2);
       
       % Before and after overlapping
       if (del<0)
           Force{j+1}=[0 0 0];
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);
       
       % During overlap
       elseif (del>0)
           v21 = v2(j+1,:)-v1;
           
           % Loading condition
           if (dot(v21,n)>0)        
                kn = kn_loading;
           
           % Unloading condition (e=sqrt(k_unloading/k_loading))
           elseif (dot(v21,n)<0)    
                kn = kn_loading*(e_rest^2);
           end
           Force{j+1}=-kn*del*n;
           F_new=Force{j+1};
           F_old=Force{j};
           [x2(j+1,:) , v2(j+1,:)] = Particle_move(x2(j,:),v2(j,:), F_new, F_old, m2, dt);

       end
       KEe(j+1) = 0.5*m2*dot(v2(j+1,:),v2(j+1,:));            
   end
   
   figure (3);
   plot(t,KEe);
   title("Kinetic Energy of particle 2 as a function of time with e=0.8")
   xlabel("Time [s]")
   ylabel("Kinetic Energy [J]")
   
%% Objective 1: Function to take arbitary position and diameter of two particles and return the overlap and the normal 

    function [del, n] = Overlap(x1, x2, d1, d2)
        
        % Radius of the particles
        r1=d1/2;
        r2=d2/2;
        
        % Relative displacement between the two particles
        x12 = x1-x2;
        
        % Distance between the two partilces
        dist_part=norm(x12);
        
        % Overlap 
        del = r1+r2-dist_part;
        
        % Normal vector
        n = x12/norm(x12);
    
    end
    
%% Objective 2: Function to move the particle 2 with a fixed time step using Leapfrog mechanism
   
    function [x_new,v_new] = Particle_move(x_old,v_old,F_new, F_old,m,dt)
            
    % *********************** Leap frog mechanism ***********************************
    % For leap frog mechanism the velocities for calculation of new position
    % are calculated at the mid-points of the time intervals (dt/2) which are
    % assumed as vstar values. Using vstar values x_new is calculated from
    % x_old (further explanation is provided in the report). 
    % The correction of velocity is done by calculating the velocities
    % v_new (using v_old values) at the exact points of the iterations i.e at 
    % (dt) intervals using vstar values and the average value of acceration before
    % and after interval dt.
        
        % Accelaration from forces before and after time interval dt 
        a_new = F_new./m;
        a_old = F_old./m;
        
        % Prediction of velocity at mid interval using old velocity and old acceleration
        vstar = v_old + 0.5*dt.*a_old;
        
        % Position shifting calculation 
        x_new= x_old + vstar*dt;
        
        % Correction of velocity after interval dt
        v_new = vstar + 0.5*dt.*(a_old+a_new);
    
    end
    