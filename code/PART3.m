% PART 3: ENHANCEMENTS

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
    
    % The following code checks the status of a particle within the X
    % boundary
    check_X_L = X_posn > 0.8e-7;
    check_X_R = X_posn < 1.2e-7;
    check_X = check_X_L & check_X_R;
    
    % The following code checks the status of the Y position within that X
    % boundary's location
    check_T = Y_posn > 0.6e-7;
    check_B = Y_posn < 0.4e-7;
    
    % Box top, bottom and In determine which box the particles is located
    % in
    Box_Top = check_T & check_X;
    Box_Bottom = check_B & check_X;
    InBox = Box_Top | Box_Bottom;
    
    % Random generation selection for particles
    while(sum(InBox) > 0)
        
        X_temp = rand(1,sum(InBox));
        Y_temp = rand(1,sum(InBox));
        X_posn(InBox) = X_temp*x_plane;
        Y_posn(InBox) = Y_temp*y_plane;
        
        % Refer to code outside of loop to determine positions
        check_X_L = X_posn > 0.8e-7;
        check_X_R = X_posn < 1.2e-7;
        check_X = check_X_L & check_X_R;
        check_T = Y_posn > 0.6e-7;
        check_B = Y_posn < 0.4e-7;
        Box_Top = check_T & check_X;
        Box_Bottom = check_B & check_X;
        InBox = Box_Top | Box_Bottom;
    end
    
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
    title('histogram for particle velocity')
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
        
        % Anything that scattered will be set to 0
        vector_mfp(scattered_P) = 0;
        
        % Using initial position, addition velocity, and the logical array
        % particles were moved around.
        
        % The particles are kept within their canvas regions and the code
        % checks if the particle attemps to exit such. If this is true then
        % the particle will be reset to the opposite side of such boundary.
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
       % Check if the Y coordinate isn't cool:
       bounds_particles = (Y_posn + V_y > y_plane | Y_posn + V_y < 0);
       
       % The following code checks the bounds for the Y coordinate to make
       % sure the particle does not go over its specified limits. If this
       % is true then the particles will "bounce" and turn around and go in
       % the opposite direction
       V_y(bounds_particles) = -1*V_y(bounds_particles);
       
       % this checks the bounded particles and changes them to send them
       % back into the canvas       
       Y_posn(1,:) = Y_posn(1,:) + V_y(1,:);
       
       
       
       % BOX REFLECTIONS -------------------------------------------------
        check_X_L = (X_posn + V_x) > 0.8e-7;
        check_X_R = (X_posn + V_x) < 1.2e-7;
        check_X = check_X_L & check_X_R;
        
        % Check the boundary conditions
        check_Box = (Y_posn + V_y) < 0.4e-7;
        Box_Bottom = check_Box & check_X;
        V_x(Box_Bottom) = -1*V_x(Box_Bottom);
        
        check_X_L = (X_posn + V_x) > 0.8e-7 + step_increment;
        check_X_R = (X_posn + V_x) < 1.2e-7 - step_increment;
        check_X = check_X_L & check_X_R;
        check_Y_bottombox_top = Y_posn < 0.4e-7 - step_increment;
        Ybox_Bottom = check_X & check_Y_bottombox_top;
        V_y(Ybox_Bottom) = -1*V_y(Ybox_Bottom);

        check_X_L = (X_posn + V_x) > 0.8e-7;
        check_X_R = (X_posn + V_x) < 1.2e-7;
        check_X = check_X_L & check_X_R;

        check_T = (Y_posn + V_y) > 0.6e-7;
        Box_Top = check_T & check_X;
        V_x(Box_Top) = -1*V_x(Box_Top);
        
        check_X_L = (X_posn + V_x) > 0.8e-7 + step_increment;
        check_X_R = (X_posn + V_x) < 1.2e-7 - step_increment;
        check_X = check_X_L & check_X_R;
        check_Y_Box_Top_bottom = Y_posn > 0.6e-7 + step_increment;
        Y_BOX_TOP = check_X & check_Y_Box_Top_bottom;
        V_y(Y_BOX_TOP) = -1*V_y(Y_BOX_TOP);
        
       % The following code checks the bounds for the Y coordinate to make
       % sure the particle does not go over its specified limits. If this
       % is true then the particles will "bounce" and turn around and go in
       % the opposite direction
       
       % this checks the bounded particles and changes them to send them
       % back into the canvas
       Particles_P = ~(X_L | X_R | Box_Top | Box_Bottom);
       X_posn(Particles_P) = X_posn(Particles_P) ...
           + V_x(Particles_P);
       
       % The following cache code is a feature to save the timestepping for
       % looping reference
       X_Cache(k,:) = X_posn(1,:);
       Y_Cache(k,:) = Y_posn(1,:);
       

       % Make sure the temperature stays at 300 by calculation
       V_avg = sum(v)/num_of_particles;
       present_temp(1,k) =  V_avg^2*C.m_n/(2*C.kb);
       
       % Calculate the mean free path given values that were obtained
       MFP = sum(vector_mfp)/num_of_particles;
       TBC_avg = MFP/V_avg;
       
    end
    
    %% Plot----------------------------------------------------
    
    figure(2)
    electron_density = [X_posn',Y_posn'];
    hist3(electron_density,'CdataMode','auto')
    xlabel('X (m)')
    ylabel('Y (m)')
    title('Electron Density Map')
    colorbar
    view(2)
    
    X_field = linspace(0,x_plane,10);
    Y_field = linspace(0,y_plane,10);

    x_bin = discretize(X_posn,X_field);
    y_bin = discretize(Y_posn,Y_field);

    temperatureBin = zeros(10,10);

    for x = 1:10
        for y=1:10
            logicX = x_bin == x;
            logicY = y_bin == y;
            logic = logicX & logicY;

            sum_x = sum(V_x(logic)) / ChangeIn_V;
            sum_y = sum(V_y(logic)) / ChangeIn_V;
            
            Vavg = sqrt((sum_x)^2 + (sum_y)^2);

            temperatureBin(x,y) = Vavg^2*C.m_n/(2*C.kb);
        end
    end
    
    figure(3)
    surf(temperatureBin)
    title('Electron Temperature Map')
    colorbar
    
    boxXpic = [0.8e-7 0.8e-7  1.2e-7 1.2e-7];
    boxTopY = [1e-7 0.6e-7 0.6e-7 1e-7];
    boxBottomy = [0 0.4e-7 0.4e-7 0];
    for Row = 1:timestep_increase
        figure(4)
        subplot(2,1,1);
        plot(Row,present_temp(Row),'.k');
        title('caused temperature')
        ylabel('Temperature (K)')
        xlabel('Time-step')
        legend(['Current Temperature:' num2str(present_temp(Row))], ...
            ['Avg Time Between Collisions:' num2str(TBC_avg)], ...
            ['Mean Free Path:' num2str(MFP)])
        hold on
        xlim([0 timestep_increase])
        ylim([250 350])
        pause(0.00000001)
        
        figure(5)
        plot(X_Cache(Row,:),Y_Cache(Row,:),'o')
        xlim([0 x_plane])
        ylim([0 y_plane])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('2-D Plot of Particle Trajectories')
        legend(['Number of Particles:' num2str(num_of_particles)])
        hold on
        plot(boxXpic,boxTopY,'k')
        plot(boxXpic,boxBottomy,'k')
    end
    
    