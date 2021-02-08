   % PART 2: COLLISIONS WITH MEAN FREE PATH

    num_of_particles = 1000;
    
    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons

    %% Plot Canvas Setup ----------------------------------------------
    
    x_plane = 20e-8;
    y_plane = 10e-8;
    
    step_increment = 1e-9;
    timestep_increase = 1000;
    
    Temp = 300;
    Thermal_V = sqrt(2*C.kb*Temp/C.m_n);

    ChangeIn_V = step_increment/Thermal_V;
    
    % set the mean time collisions between 0.2e-12
    Collision_mean_time = 0.2e-12;
    
    %% Initializing Random Particles -------------------------------------
    
    % generates two rows containing random numbers for X coordinates of each particle.
    % Using the two rows, the position and angles are acquired
    X = rand(2,num_of_particles);
    Y = rand(1,num_of_particles);
    
    % Giving postitional coordinates for x and y points and multiply by
    % plane to ensure restrictions
    X_posn(1,:) = X(1,:)*x_plane;
    Y_posn(1,:) = Y(1,:)*y_plane;
    
    % Apply an angle for each particle (rads unit)
    angle_value(1,:) = X(2,:)*2*pi;
    
    % Setting up the Maxwell_Blotzmann distribution using the global
    % components. This is for the components which contain the thermal
    % velocity distribution
    % Generates a histogram based off calculated values
    sig = sqrt(C.kb*Temp/C.m_n)/4;
    m_b_d = makedist('Normal',Thermal_V,sig);
    v = random(m_b_d,1,num_of_particles);
    figure(1)
    hist(v)
    title('Histogram for Avg Velocity')
    V_x = ChangeIn_V*v(1,:).*cos(angle_value(1,:));
    V_y = ChangeIn_V*v(1,:).*sin(angle_value(1,:));
    
    % Setting up scattering percent
    Perc_scat = 1 - exp(-ChangeIn_V/Collision_mean_time);
    vector_mfp = zeros(1,num_of_particles);
    
    %% Canvas design ------------------------------------------------
    

    for k = 1:timestep_increase
        
        % sets up the Scattering electrons
        scattered_particles = rand(1,num_of_particles);
        scattered_P = scattered_particles < Perc_scat;
        
        % If the scattered percentage oversees the scatter then the
        % particle scatters. In such circumstance the logical set array
        % becomes all 1's which will also scatter.
        
        % generate a new random angle
        angle_value(scattered_P) = rand*2*pi;
        
        % generate a new random velocity:
        v = random(m_b_d,1,num_of_particles);
        V_x(scattered_P) = ChangeIn_V*v(scattered_P) ...
            .*cos(angle_value(scattered_P));
        V_y(scattered_P) = ChangeIn_V*v(scattered_P) ... 
            .*sin(angle_value(scattered_P));
        
        % The follow code finds the mean free path by keeping track of the
        % particles and their movements
        vector_mfp(~scattered_P) = vector_mfp(~scattered_P)+step_increment;
        
        % Any scattered electron will be moved to 0.
        vector_mfp(scattered_P) = 0;
        
        % Using initial position, addition velocity, and the logical array
        % particles were moved around.
        
        % The particles are kept within their canvas regions and the code
        % checks if the particle attemps to exit such. If this is true then
        % the particle will be reset to the opposite side of such boundary
        X_R = X_posn + V_x > x_plane;
        X_posn(X_R) = X_posn(X_R) ... 
            + V_x(X_R) - x_plane;

        X_L = X_posn + V_x < 0;
        X_posn(X_L) = X_posn(X_L) ... 
            + V_x(X_L) + x_plane;
        
       % When the logical array sets to zero, we are informed that nothing
       % has passed the field boundaries. To do this the position and
       % velocity are added.
       
       % bound particles are the displayed particles which aren't in x_left or x_right
       bounds_particles = ~(X_L | X_R);
       X_posn(bounds_particles) = X_posn(bounds_particles) + V_x(bounds_particles);
       
       % The following code checks the bounds for the Y coordinate to make
       % sure the particle does not go over its specified limits. If this
       % is true then the particles will "bounce" and turn around and go in
       % the opposite direction
       Y_bounds = (Y_posn + V_y > y_plane | Y_posn + V_y < 0);
       
       % this checks the bounded particles and changes them to send them
       % back into the canvas
       V_y(Y_bounds) = -1*V_y(Y_bounds);
       Y_posn(1,:) = Y_posn(1,:) + V_y(1,:);
       
       % The following cache code is a feature to save the timestepping for
       % looping reference
       X_cache(k,:) = X_posn(1,:);
       Y_cache(k,:) = Y_posn(1,:);
       

       % Make sure the temperature stays at 300 by calculation
       V_avg = sum(v)/num_of_particles;
       present_temp(1,k) =  V_avg^2*C.m_n/(2*C.kb);
       
       % Calculate the mean free path given values that were obtained
       MFP = sum(vector_mfp)/num_of_particles;
       TBC_avg = MFP/V_avg;
       
    end
    
    %% Plot ----------------------------------------------------
    % Everything is calculated! Now let's start plotting...
    figure(2)
    for paths = 1:num_of_particles
        subplot(3,1,1);
        plot(X_cache(:,paths),Y_cache(:,paths),'-')
        xlim([0 x_plane])
        ylim([0 y_plane])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Final Electron Paths')
        hold on
    end
    for row = 1:timestep_increase

        subplot(3,1,2);
        plot(row,present_temp(row),'.k');
        title('Temperature_Caused_Electrons')
        ylabel('Temperature (K)')
        xlabel('Time-step')

        legend(['Current Temperature:' num2str(present_temp(row))], ...
            ['Avg Time Between Collisions:' num2str(TBC_avg)], ...
            ['Mean Free Path:' num2str(MFP)])
        hold on
        xlim([0 timestep_increase])
        ylim([250 350])
        pause(0.0001)
    end