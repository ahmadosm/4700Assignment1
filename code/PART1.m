% PART 1 ELECTRON MODELLING

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
    
    % The given planary size for the particle region is 200nm × 100nm
    x_plane = 20e-8;
    y_plane = 10e-8;
    
    % The chosen time step increment is given to be smaller than 1/100 the region
    step_increment = 1e-9;
    timestep_increase = 1000;
    
    % Assume the Temperature is 300 Kelvins units
    Temp = 300;
    
    % Thermal velocity equation calculation
    Thermal_V = sqrt(2*C.kb*Temp/C.m_n);
    
    % calculating the change in velocity for each given timestep
    ChangeIn_V = step_increment/Thermal_V;
    
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
    
    % Set up X and Y using thermal velocity
    V_x = Thermal_V*ChangeIn_V*cos(angle_value(1,:));
    V_y = Thermal_V*ChangeIn_V*sin(angle_value(1,:));
    
    
    %% Canvas design ------------------------------------------------
    
    for k = 1:timestep_increase
        % Using a logical array the electrons move starting at an initial
        % position and continue to speed up when velocity is applied to the
        % speed of the electron
        
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
       X_posn(bounds_particles) = X_posn(bounds_particles) ...
           + V_x(bounds_particles);
       
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
       X_Cache(k,:) = X_posn(1,:);
       Y_Cache(k,:) = Y_posn(1,:);
       
       % Make sure the temperature stays at 300 by calculation
       total_X = sum( (V_x/ChangeIn_V).^2);
       total_Y = sum( (V_y/ChangeIn_V).^2);
       temp = (total_X+total_Y)*C.m_n/(2*C.kb);
       present_temp(k) = temp/num_of_particles;
    end
    
    %% Plots ----------------------------------------------------

    figure(1)
    for paths = 1:num_of_particles
        subplot(3,1,1);
        plot(X_Cache(:,paths),Y_Cache(:,paths),'-')
        xlim([0 x_plane])
        ylim([0 y_plane])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Final Electron Paths')
        hold on
    end
    
    for row = 1:timestep_increase
        subplot(3,1,3);
        plot(row,present_temp(row),'.k');
        title('The Current Temperature Caused by Those Electrons')
        ylabel('Temperature (K)')
        xlabel('Time-step')

        legend(['Present Temperature:' num2str(present_temp(row))], ...
            ['Thermal Velocity:' num2str(Thermal_V)], ...
            ['The free path mean:' num2str(Thermal_V*0.2e-12)])
        hold on
        xlim([0 timestep_increase])
        ylim([Temp-1 Temp+1])
        pause(0.01)
    end
