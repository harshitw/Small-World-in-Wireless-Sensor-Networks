Min_Diam = 12.0;
Max_Diam = 11.0;
Intv = 1.00;

rad_omni = 1.0;
omni_area = pi*(rad_omni^2.0);

base_diam = 8.0;
base_N = 300;
base_area = base_diam^2;

Num_Sim = 1 + (Max_Diam - Min_Diam)/Intv; % number of simulations 

N = 300;        
SimCount = 0;

% Antenna Model
lambda = (3*10^8)/(2.4*10^9);
Delta = lambda/2;

for Diam = Min_Diam:Intv:Max_Diam
    
    x_low = 0.0;
    Y_low = 0.0;
    
    x_high = Diam;
    y_high = Diam;
    
    sim_area = x_high * y_high; % area of rectangular simulation region 
    
    N = ceil(base_N * (sim_area/base_area));
    
    density_omni = N * (omni_area/sim_area); % no of nodes in omnidirectional area
    
    SimCount = SimCount + 1;
    
    % computation on different values of theta on 0 to pi/2
    % compute theta for whole number values of r
    
    r_theta = 10.0 * rad_omni; % -- Assuming range to be 10 hops
    
    theta = (2.0 * pi) * ((rad_omni^2.0)/(r_theta^2.0));

    dpmax = 0.0;
    theta_opt = 0.0;
    
    while theta <= pi/2.0
        
        ar_pl = (theta/2.0) * ((ceil(r_theta) - 1.0)^2.0) * (rad_omni^2.0); % area of sectoral region upto the second last one 
        Ar_last = omni_area - ar_pl; % maintains the condition that omni directional area is same as the directional area
        
        % Obtain the expected distance of a node in the last region
        r_prev = ceil(r_theta) - 1.0;
        k = r_theta - r_prev;
        
        x = r_prev;
        delx = k/1000.0;
        
        pdf = 0.0;
        cdf = 0.0;
        
        while x <= r_theta
            pdf = delx * ((x^2.0) * theta)/Ar_last;
            cdf = cdf + pdf;

            x = x + delx;
        end
        
        E_dlast = cdf;
        
        E_nlast = density_omni * (Ar_last/omni_area); % Compute number of nodes in the last sector 
        
        pnl = 1.0 - (1.0 - (Ar_last/omni_area))^density_omni; % probability of finding a node in the last region 
         
        Ar_first = theta/2.0;
        pnf = 1.0 - (1.0 - (Ar_first/omni_area))^density_omni;% probability of finding a node in the first sector 
        
        dir_prod = E_dlast * pnf * pnl;
        
        if dpmax < dir_prod 
            dpmax = dir_prod;
            theta_opt = theta;
        end
        
        % fprintf('theta\t%.2f\trtheta\t%.2f\tArL\t%.4f\tEdl\t%.4f\tEnl\t%.2f\tpnl\t%.4f\tpnf\t%.2f\tdpr\t%.4f\n', theta, r_theta, Ar_last, E_dlast, E_nlast, pnl, pnf, dir_prod/dpmax);
        
        r_theta = r_theta - 1.0;
        theta = (2.0 * pi) * ((rad_omni^2.0)/(r_theta^2.0));
    
    end
    
    Adj_Omni = zeros(N,N);    % -- Adjacency Matrix that will store neighbourhood information for the initial omnidirectional case
    Adj_Dir  = zeros(N,N);    % -- Adjacency Matrix after nodes orient their beams
    node_loc = zeros(N,4);         % -- Store information of nodes - location & beam info: angle and width
    
    % Intermediate variables to be used for coordinates
    x_temp = 0.0;
    y_temp = 0.0;
    
    % Generate node locations
    for i=1:N
        % Generate random relative location of the node i
        x_temp = rand(1) * (x_high - x_low);
        y_temp = rand(1) * (y_high - y_low);

        node_loc(i,1) = x_low + x_temp;         % -- Absolute X-coordinate of node
        node_loc(i,2) = y_low + y_temp;         % -- Absolute Y-coordinate of node
        node_loc(i,3) = 0;                      % -- Beam direction - initially along the X-axis
        node_loc(i,4) = 2 * pi;                 % -- Beam width - omnidirectional initially
    end
    
    distij = 0.0;
    
    % Generate adjacency matrix - iterate through all nodes twice O(N^2)
        for i=1:N
            for j=1:N
                if( i == j )            % -- 0 hops to myself
                    Adj_Omni(i,j) = 0;
                else
                    distij = get_distance(node_loc(i,1), node_loc(i,2), node_loc(j,1), node_loc(j,2));

                    if( distij < rad_omni )
                        Adj_Omni(i,j) = 1;
                    end
                end
            end
        end
        
        PathList_Omni = [];             % -- List of shortest path lengths between each node pair - omnidirectional
        NextHops_Omni = [];
        
        %sh_paths = simple_dijkstra(Adj_Omni, i, N);
        for i=1:N
            [sh_paths, nxt_hops] = dijkstra_paths(Adj_Omni, i, N);
            PathList_Omni = [PathList_Omni; sh_paths];
            NextHops_Omni = [NextHops_Omni; nxt_hops];
        end
        
        rth_sim = r_thopt ;
        theta_sim = (2.0 * pi) * ((rad_omni^2.0)/(rth_sim^2.0));
        
        order = 0;
        A = [];
        zeta = 0;

        for i = 1:N
            for j = 1:N
                if Adj_Omni(i,j) == 1
                    order = order + 1;
                    A = [A j];
                end
            end    
            [e,F] = size(A);
            for k = 1:F
                for t = (k+1):F
                    if Adj_Omni(A(1,k),A(1,t)) == 1
                        zeta = zeta + 1;
                    end    
                end
            end
            if order > 1
                cc = [cc ; (zeta*2)/(order*(order-1))];
            elseif order == 1
                cc = [cc ; 0];
            end     
            order = 0;
            zeta = 0;
            A = [];
        end
        
        [cl_cf,index] = sort(cc,'descend');
        
        % dil_list = zeros(1,N);
        
        deg = sum(Adj_Omni,1) + sum(Adj_Omni,2) - diag(Adj_Omni);
        
        dil_list = betwn(Adj_Omni,deg,N);
        
        p = 0.1 ;
        
        for i = 1 : N
            if (dil_list(i) == i)
                nbd_dir = index(1,i);
                beam_dir = get_line_angle(node_loc(i,1), node_loc(i,2), node_loc(nbr_dir,1), node_loc(nbr_dir,2));
                node_loc(i,3) = beam_dir;       % -- Beam direction
                node_loc(i,4) = theta_sim;      % -- Beam width
            end
            for j = 1 : N
                if (dil_list(i) == 1)
                    if( i == j )
                        Adj_Dir(i,j) = 0;               % -- 0 hops to myself
                    else
                        beam_dir = node_loc(i,3);
                        beam_width = node_loc(i,4);
                        line_angle = get_line_angle(node_loc(i,1), node_loc(i,2), node_loc(j,1), node_loc(j,2));    % -- Angle of the line joining i & j to the X-axis

                        chi  = ((pi * Delta)/lambda) * (cos(line_angle) - cos(beam_dir));
                        Gain1 = (1/rth_sim) * (sin(rth_sim*chi)./sin(chi)).^2;
                        distij = get_distance(node_loc(i,1), node_loc(i,2), node_loc(j,1), node_loc(j,2));
                        
                        if( distij < Gain1 )
                            Adj_Dir(i,j) = 1;
                        end
                        
                    end
                else
                    if( i == j )            % -- 0 hops to myself
                        Adj_Dir(i,j) = 0;
                    else
                        distij = get_distance(node_loc(i,1), node_loc(i,2), node_loc(j,1), node_loc(j,2));

                        if( distij < rad_omni )
                            Adj_Dir(i,j) = 1;
                        end
                    end
                end
            end
        end
                
                
                