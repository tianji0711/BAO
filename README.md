# AW-BAO
%% Parameter Setting
Dim=3;      % Dimension
SwarmSize=10; % Swarm size (number of particles)
ObjFun=@fmin; % Objective function handle
MaxIter=N;     % Maximum number of iterations
MinFit=0.001;  % Minimum fitness value
Vmax=1;
Vmin=-1;
Ub=[10 10 10]; % Upper bound
Lb=[0 0 0];    % Lower bound

c=5;           % ratio between step and d0
step=1;        % Step size
chaos_init = 0.7; % Initial value for chaotic map

%% Particle Swarm Initialization
Range = ones(SwarmSize,1)*(Ub-Lb);
Swarm = rand(SwarmSize,Dim).*Range + ones(SwarmSize,1)*Lb;  
VStep = rand(SwarmSize,Dim)*(Vmax-Vmin) + Vmin;             
fSwarm = zeros(SwarmSize,1);
for i=1:SwarmSize
    fSwarm(i,:) = feval(ObjFun,Swarm(i,:));                 
end

%% Individual and Global Extremes
[bestf, bestindex]=min(fSwarm);
zbest=Swarm(bestindex,:);   
gbest=Swarm;                
fgbest=fSwarm;              
fzbest=bestf;               

%% Iterative Optimization
iter = 0;

while( (iter < MaxIter) && (fzbest > MinFit) )
    
    % Adaptive weights and constants
    w = 0.9 - iter * 0.5 / MaxIter;
    c1 = 2 + sin(iter * pi / MaxIter);
    c2 = 2 - sin(iter * pi / MaxIter);
    c3 = 2 - cos(iter * pi / MaxIter);
    
    for j=1:SwarmSize
        d0=step/c;
        dir(j,:)=rand(1,Dim); 
        dir(j,:)=dir(j,:)/(eps+norm(dir(j,:))); 
        xleft(j,:)=Swarm(j,:)+dir(j,:)*d0/2;    
        fleft(j,:)=fmin(xleft(j,:));            
        xright(j,:)=Swarm(j,:)-dir(j,:)*d0/2;   
        fright(j,:)=fmin(xright(j,:));          
        force_bas(j,:)=-dir(j,:)*sign(fleft(j,:)-fright(j,:));
   
        % Chaotic mapping
        chaos_init = 4 * chaos_init * (1 - chaos_init);
        chaos_factor = chaos_init;
        
        % Update velocity with adaptive weights and chaotic mapping
        VStep(j,:) = w*VStep(j,:) + c1*rand*(gbest(j,:) - Swarm(j,:)) + c2*rand*(zbest - Swarm(j,:)) + c3*rand*force_bas(j,:) + chaos_factor;
        
        % Boundary checks for velocity
        VStep(j,:) = max(min(VStep(j,:), Vmax), Vmin);
        
        % Update position
        Swarm(j,:)=Swarm(j,:)+VStep(j,:);
        
        % Boundary checks for position
        Swarm(j,:) = max(min(Swarm(j,:), Ub), Lb);
        
        % Fitness value
        fSwarm(j,:) = feval(ObjFun,Swarm(j,:));
        
        % Update individual best
        if fSwarm(j) < fgbest(j)
            gbest(j,:) = Swarm(j,:);
            fgbest(j) = fSwarm(j);
        end
        
        % Update global best
        if fSwarm(j) < fzbest
            zbest = Swarm(j,:);
            fzbest = fSwarm(j);
        end
    end 
    
    iter = iter+1;           % Update iteration count
end


