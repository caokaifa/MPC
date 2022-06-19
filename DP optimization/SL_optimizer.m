function controller_output = SL_optimizer(p,x0,controller_output_previous)

controller_output = struct;
timer = tic;
obstacle=[57;0];
if isempty(controller_output_previous) % First time step (start of the race)
    % Initial guessed trajectory: Standing still at x0.
    x = repmat(x0,1,p.Hp);
else
    % Shift previous trajectory by one time step.
    x = controller_output_previous.x;    
    x(:,1:end-1) = x(:,2:end);
end

for i = 1:1
    % For each trajectory point, find the closest track checkpoint.i 
    checkpoint_indices = nan(1, p.Hp);
    for k=1:p.Hp
        ff=sum(([p.checkpoints.center]-repmat(x(1:2,k),1,length(p.checkpoints))).^2);
        [bb,index]=min(sum(([p.checkpoints.center]-repmat(x(1:2,k),1,length(p.checkpoints))).^2));
        [~,checkpoint_indices(k)] = min(sum(([p.checkpoints.center]-repmat(x(1:2,k),1,length(p.checkpoints))).^2));
    end

    left_points = [p.checkpoints(:).left];
    right_points = [p.checkpoints(:).right];
    center_points =[p.checkpoints(:).center];
    forward_vectors = [p.checkpoints.forward_vector];
%     plot(left_points(1,checkpoint_indices) ,left_points(2,checkpoint_indices),'.r','LineWidth',5);
%     hold on
%     plot(right_points(1,checkpoint_indices),right_points(2,checkpoint_indices),'.b','LineWidth',5);
    Dx=right_points(1,checkpoint_indices)-left_points(1,checkpoint_indices);
    Dy=right_points(2,checkpoint_indices)-left_points(2,checkpoint_indices);
    grid_width=28;
    grid_x=zeros(p.Hp,grid_width);
    grid_y=zeros(p.Hp,grid_width);
    for i=1:p.Hp
        for j=1:grid_width
            % equally spaced grid
            grid_x(i,j) = left_points(1,checkpoint_indices(i)) + j/(grid_width+1)*Dx(i);
            grid_y(i,j) = left_points(2,checkpoint_indices(i)) + j/(grid_width+1)*Dy(i);
%             hold on
%             plot(grid_x ,grid_y,'*r','LineWidth',6);
        end
    end  
%     hold on
%     plot(grid_x ,grid_y,'*r','LineWidth',5);
    grid_isOccupiedBy = zeros(p.Hp,grid_width); %initialize
    % for c, type cast this to integers
    % see documents: mxGetData, mxClassID and cast
    grid_isOccupiedBy = cast(grid_isOccupiedBy, 'int32');
    obstacle = repmat(obstacle,1,p.Hp);
%     if CarIsClose(left_points(:,checkpoint_indices),obstacle)
   [grid_isOccupiedBy,vertices] = getGridOccupancy(grid_isOccupiedBy,grid_x,grid_y,obstacle,p);
   
%    hold on
%    plot(vertices(:,1),vertices(:,2),'g','LineWidth',10)
%     end
    if max(max(abs(grid_isOccupiedBy)))==0 %in C -1, first car is car0
    % if all elements of grid_isOccupiedBy are zeros (ie. no cars nearby),
    % then use default borders
        new_b_array= [p.checkpoints(:,checkpoint_indices).left;p.checkpoints(:,checkpoint_indices).right]';
    else
%         try
        new_b_array = getNewBorders_helper(center_points(:,checkpoint_indices),grid_x,grid_y,grid_isOccupiedBy,p,checkpoint_indices,obstacle);
    end
%替换更新的边界
   for m=1:p.Hp
         p.checkpoints(:,checkpoint_indices(m)).left=new_b_array(m,1:2)';
         if m<p.Hp
            if(checkpoint_indices(m+1)-checkpoint_indices(m)>1)
                obstacle_flag=grid_isOccupiedBy(m,:);
                obstacle_flag_1=grid_isOccupiedBy(m+1,:);
                if(max(abs(obstacle_flag))>0 && max(abs(obstacle_flag_1))>0)
                    for n=checkpoint_indices(m):1:checkpoint_indices(m+1)
                        p.checkpoints(:,n).left=new_b_array(m,1:2)';
                    end
                end
            end
         end
                   
         p.checkpoints(:,checkpoint_indices(m)).right=new_b_array(m,3:4)';
   end
    % Formulate and solve the linearized trajectory optimization problem.
    [x, U, optimization_log] = SL_QP(p,x0,checkpoint_indices,x);
end

opt_time = toc(timer);
fprintf('optT %6.0fms slack %6.2f fval %6.1f\n', opt_time*1000, optimization_log.slack_lateral, optimization_log.fval);

controller_output.checkpoint_indices = checkpoint_indices;
controller_output.x = x;
controller_output.U = U;
controller_output.optimization_log = optimization_log;
controller_output.opt_time = opt_time;
controller_output.vertices=vertices;
controller_output.left=[p.checkpoints(:).left];
controller_output.right=[p.checkpoints(:).right];
end


function bool = CarIsClose(X_curr,X_oth) 
% This implements a quick and cheap method to check if two cars are close.
% Simply compare their theta values at time 0. If they are within a
% threshold of each other, returns true. Otherwise, returns false. Note
% that wrap-around is implemented with respect to tracklength.
%
% Note: assumed theta index is nx

threshold=5.0;


bool=false;

if   sqrt(X_curr(1,end) - X_oth(1,end)^2+X_curr(2,end) - X_oth(2,end)^2) < threshold
   
    % 2nd line of above if-condition checks for track wrap-around
    bool=true;
    return;
end

end


function [grid_isOccupiedBy,vertices] = getGridOccupancy(grid_isOccupiedBy,grid_x,grid_y,X_oth,p)
% check if (grid_x,grid_y)(i,j) is inside polygon of {othercaridx}'th car at
% corresponding time.
%
% if grid(i,j) is occupied by the car, grid_isOccupiedBy(i,j)
% is assigned the index of the occupying car (othercaridx). Otherwise it is
% untouched.

carwidth=1;
carlength=2;

for i=1:p.Hp
    % position of car
    posn = [X_oth(1,i) X_oth(2,i)];
    
    % orientation of car (zero is positive x axis)
    phi = 0;
    cosphi=cos(phi);
    sinphi=sin(phi);
    
    % tangent and normal vectors of bounding box (normal vector is 90 deg
    % clockwise of tangent vector)
    t = [cosphi,sinphi]*carlength/2;
    n = [sinphi,-cosphi]*carwidth/2;
    
    % corners of bounding box of the car with respect to the point of
    % rotation used in the model
    FL = posn + t + n;
    FR = posn + t - n;
    BL = posn - t + n;
    BR = posn - t - n;
    
    % check gridpoints of i'th row is inside bounding box of car
    vertices = [FR;BR;BL;FL;FR]; % [number of corners]by[2] array
 
    points = [grid_x(i,:);grid_y(i,:)]'; % [number of points to test]by[2] array

    in = inpoly(points,vertices);
    
    grid_isOccupiedBy(i,in) = 1;

end

end

function [in] = inpoly(points,vertices)

% determinse whether a point is inside the polytop given by the vertices
n = length(points);
CorrectSum = sum(convhull(vertices(:,1),vertices(:,2)));
in = false(n,1);
for i = 1:n
    
	tmp_vertices = [vertices;points(i,:)];
    if sum(convhull(tmp_vertices(:,1),tmp_vertices(:,2))) == CorrectSum
        in(i) = 1;
    else
        in(i) = 0;
    end
end

end


function new_b_array = getNewBorders_helper(X,grid_x,grid_y,grid_isOccupiedBy,p,checkpoint_indices,obstacle)
% to do


priorPath = Lock2Grid(X,grid_x,grid_y,p);
[Optm_Path_Mtx Optm_Cost_Mtx]=getOptimalPathMatrix(X,grid_x,grid_y,grid_isOccupiedBy,priorPath,p);
% trace optimal policy matrix to get optimal path
OptmPath=TraceOptimalPath(Optm_Path_Mtx,p);
new_b_array= Find_Nearest_Borders(X,grid_x,grid_y,grid_isOccupiedBy,OptmPath,p,checkpoint_indices,obstacle);
end

function priorPath = Lock2Grid(X,grid_x,grid_y,p)
% priorPath - [N]by[1] (dimensionless)


grid_width = 28;


priorPath = zeros(p.Hp,1);
% for c, type cast this to integers
priorPath = cast(priorPath, 'int32');
    
for i=1:p.Hp
    min=100; % choose some big number that distance^2 can never reach
    ind=1;
    for j=1:grid_width
        dx = grid_x(i,j) - X(1,i);
        dy = grid_y(i,j) - X(2,i);
        d2=dx*dx+dy*dy;
        if d2<min
            min=d2;
            ind=j;
        end
    end
    priorPath(i) = ind;
end
end


function [Optm_Path_Mtx Optm_Cost_Mtx]=getOptimalPathMatrix(X,grid_x,grid_y,grid_isOccupiedBy,priorPath,p)
% for i>1: jnxt=Optm_Path_Mtx(i+1,j,j1) is the optimal next position
%          (i+1,jnxt) if you are at (i,j) and you came from (i-1,j1)
%          (goes up to i=N-1)
% for i=1: Optm_Path_Mtx(2,j,1) is the optimal next position if you are at
%          (1,j)
% for i=0: Optm_Path_Mtx(1,1,1) is the optimal next position from the
%          starting position
%
% Optm_Cost_Mtx(i+1,j,j1) is the optimal cost if you are at (i,j) and came
% from (i-1,j1) and choose to go to (i+1,jnxt)
%
% note: if (i+1,jnxt) is invalid for all jnxt, then Optm_Path_Mtx(i+1,j,j1)
% is assigned value of DEAD_END_PATH. Note that the position (i,j) itself
% is a valid position.


N=p.Hp;
grid_width=28;
Cost_dead_end=5000;
time_step_cost_discount_factor=0.95;
UNOCCUPIED=0;
DEAD_END_PATH=-2;

Optm_Path_Mtx=zeros(N,grid_width,grid_width);
Optm_Path_Mtx = cast(Optm_Path_Mtx, 'int32');


% figure(4); title('getNewBorders.m\getOptimalPathMatrix()\grid_isOccupiedBy');
% spy(grid_isOccupiedBy)
% asdf=5;

for i=N-1:-1:1
    for j1=1:grid_width
        % no need to look for a valid path if previous position (i-1,j1) is
        % invalid. Also no need to check for i=1, since i=0 will always
        % have a valid j
        if i>1 && grid_isOccupiedBy(i-1,j1)~=UNOCCUPIED
            continue;
        end
        
        % change to true if *any* valid jnxt is found for at least 1 j
        valid_path_is_found = false; 
        
        for j=1:grid_width
            % no need to look for a valid path if current position (i,j) is
            % invalid
            if grid_isOccupiedBy(i,j)~=UNOCCUPIED
                continue;
            end

            
            % change to true if a valid jnxt is found for *current* j
            valid_jnxt_is_found = false;
        

            for jnxt=1:grid_width
                %if (i<=N-2 && Optm_Path_Mtx(i+2,jnxt,j)==DEAD_END_PATH) || ...
                if    PathIsObstructed(grid_isOccupiedBy,i,j,jnxt,UNOCCUPIED)
                    
                    if jnxt==grid_width && ~valid_jnxt_is_found
                        %uh oh: no valid jnxt found for this j
                        %assign a really big cost to this state
                        Optm_Cost_Mtx(i+1,j,j1) = Cost_dead_end;
                        Optm_Path_Mtx(i+1,j,j1) = DEAD_END_PATH;

                    end
                    
                    % if not yet at last point in row, continue checking.
                    % if at last point, move onto the next j.
                    continue;
                end
                
                % if code reaches here, both j and jnxt are valid

                % calculate costs
                Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt);
%                 Cost_PathLength=getCost_PathLength(grid_x,grid_y,X,i,j,jnxt,stateindex_x,stateindex_y);
%                 Cost_AngleChange=getCost_AngleChange(grid_x,grid_y,X,i,j,j1,jnxt,DesignParameters,ModelParams);
                
                % add costs
                if i==N-1
                    Cost_Total = 0.5*Cost_PathDeviation;
                else
                    % also add on (discounted) cost of (i+2,jnxt,j)
                    Cost_Total = 0.5*Cost_PathDeviation ...
                        + time_step_cost_discount_factor * Optm_Cost_Mtx(i+2,jnxt,j);
                end
                
%                 if(i==12 && j==4 && j1==1)
%                     tmp=[i j j1 jnxt;
%                         Cost_PathDeviation Cost_PathLength Cost_AngleChange Optm_Cost_Mtx(i+2,jnxt,j);
%                         Cost_Total zeros(1,3)];
%                     asdf=5;
%                 end
                
                if ~valid_jnxt_is_found
                    % this is first valid path found, so set values for
                    % later comparison
                    minCost = Cost_Total;
                    jnxt_minCost = jnxt;
                    valid_jnxt_is_found = true;
                else
                    % compare with previous minimum cost
                    if Cost_Total<minCost
                        minCost = Cost_Total;
                        jnxt_minCost=jnxt;
                    end
                end
            end %end of jnxt loop

            
            

            if valid_jnxt_is_found
                % minCost is now the minimum cost for (i,j,j1)
                Optm_Cost_Mtx(i+1,j,j1) = minCost;
                Optm_Path_Mtx(i+1,j,j1) = jnxt_minCost;
                valid_path_is_found = true;
            end
                       
            
        end % end of j loop
        
        if ~valid_path_is_found
            %uh oh: no valid path is found!
            % return an error code and exit
            error('no_valid_path');
        end
        
        if i==1
            % for i=1, does not use j1, so arbitrarily choose j1=1 to
            % store data in. Hence, can break out of the j1 loop.
            break;
        end
    end % end of j1 loop
    
end % end of i loop


% i=0 must be dealt with separately
i=0;
valid_jnxt_is_found=false; %since only the single j, find a valid jnxt IFF find a valid path

for jnxt=1:grid_width
    % check if grid(1,jnxt) is occupied
    if grid_isOccupiedBy(1,jnxt)~=UNOCCUPIED
        if jnxt==grid_width && ~valid_jnxt_is_found
            % uh oh: no valid path for any jnxt!
            % return an error code and exit
            error('no_valid_path');
        end
        % if not yet reached last point in row, continue checking
        continue;
    end
    
    % if code reaches here, valid jnxt has been found
    
    % costs
    Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt);
%     Cost_PathLength=getCost_PathLength(grid_x,grid_y,X,i,j,jnxt,stateindex_x,stateindex_y);
%     Cost_AngleChange=getCost_AngleChange(grid_x,grid_y,X,i,j,j1,jnxt,DesignParameters,ModelParams); %j1 NOT SET...
    
    % total cost
    Cost_Total = 0.5*Cost_PathDeviation ...
        + time_step_cost_discount_factor ...
        * Optm_Cost_Mtx(i+2,jnxt,1); % 3rd dim=1 here because j1 for i=1 is always 1
    
    if ~valid_jnxt_is_found
        % first valid path found, so set values for later comparison
        minCost=Cost_Total;
        jnxt_minCost=jnxt;
        valid_jnxt_is_found=true;
    else
        % compare with previous minimum cost
        if Cost_Total<minCost
            minCost = Cost_Total;
            jnxt_minCost=jnxt;
        end
    end
    
end % end of jnxt loop

% minCost is now the minimum cost for Optm_Cost_Mtx(i+1=1,:,:)
% jnxt_minCost is the choice of path associated with this
Optm_Path_Mtx(1,1,1) = jnxt_minCost; % choose to store this in first element
Optm_Cost_Mtx(1,1,1) = minCost; 

% clean up and exit (no error)


% note: Optm_Path_Mtx_mex not implemented (due to 1-based vs 0-based
% counting), so just rely on Optm_Cost_Mtx_mex comparison

end



function bool=PathIsObstructed(grid_isOccupiedBy,i,j,jnxt,UNOCCUPIED)
% returns true if path from (i,j) to (i+1,jnxt) is obstructed.
% True if any of the points in the box defined by its diagonal corners
% at (i,j) and (i+1,jnxt) is occupied.


if j<=jnxt
    for jj=j:jnxt
        if grid_isOccupiedBy(i,jj) ~= UNOCCUPIED || ...
           grid_isOccupiedBy(i+1,jj) ~= UNOCCUPIED
            bool = true; return;
        end
    end
else
    for jj=j:-1:jnxt
        if grid_isOccupiedBy(i,jj) ~= UNOCCUPIED || ...
           grid_isOccupiedBy(i+1,jj) ~= UNOCCUPIED
            bool = true; return;
        end
    end
end

bool = false;

end

function Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt)
% cost of deviation from previously planned path when at (i+1,jnxt)
% i goes from 0 to N-1

dj = abs( priorPath(i+1) - jnxt );

% cast to double
dj = cast(dj, 'double');

Cost_PathDeviation = dj*100;


end


function OptmPath=TraceOptimalPath(Optm_Path_Mtx,p)
% Optm_Path_Mtx is [N]by[grid_width]by[grid_width]
% OptmPath is [N]by[1]

N=p.Hp;
DEAD_END_PATH = -2;

OptmPath = zeros(N,1); % initialize.
% OptmPath = cast(OptmPath, 'int32'); % cast to int

OptmPath(1) = Optm_Path_Mtx(1,1,1);
OptmPath(2) = Optm_Path_Mtx(2,OptmPath(1),1);

for i1=3:N
    % i1 represents i+1

    OptmPath(i1) = Optm_Path_Mtx(i1,OptmPath(i1-1),OptmPath(i1-2));
    if OptmPath(i1)==DEAD_END_PATH
        % exit function. Find_Nearest_Borders() will take care of the rest.
        disp('BLAAAAAAAAAAAAAAAAAAAAAAAAAAAAH');
        return;
    end

end

end

function new_b=Find_Nearest_Borders(X_all,grid_x,grid_y,grid_isOccupiedBy,OptmPath,p,checkpoint_indices,obstacle)

% new_b_default is [Lx, Ly, Rx, Ry]
% note: have made assumption that state indices for x,y,psi are 1,2,3



grid_width=28;
N=40;
carlength=2;
carwidth=1;
UNOCCUPIED=0;
DEAD_END_PATH = -2;
LAST_BORDER_XMAX=10;

new_b_L = zeros(N,2);
new_b_R = zeros(N,2);

% for each row
for i=1:N
    j_optm = OptmPath(i);
    
    if j_optm==DEAD_END_PATH
        
%         % (i,j_optm) is invalid for any choice of j_optm
%         % note that (i-1,j_prev) is still a valid position
%         if i==1
%             % cannot use previous borders, what to do?
%             error('no valid first point');
%         end
%         
%         j_prev=OptmPath(i-1);
%         
%         % if not first point, use borders from previous point for rest of
%         % horizon
%         for ii=i:N
%             new_b_L(ii,:)=new_b_L(i-1,:);
%             new_b_R(ii,:)=new_b_R(i-1,:);
%         end
%         new_b = [new_b_L new_b_R];
%         
%         % set the last border
%         P1=[grid_x(i-1,j_prev) grid_y(i-1,j_prev)];
%         P2=[grid_x(  i,j_prev) grid_y(  i,j_prev)];
%         last_border=calculate_last_border(P1,P2);
%         
% %         figure(5);hold on;
% %         plot(grid_x,grid_y,'go');
% %         plot([P1(1) P2(1)],[P1(2) P2(2)],'r*');
% %         xx=[-2,2];
% %         yy=(last_border(3)-last_border(1)*xx)/last_border(2);
% %         plot(xx,yy,'b');
% %         axis([-2 2 -2 2]);
        
        return; % exit from function
    end
    left_points = [p.checkpoints(:).left];
    right_points = [p.checkpoints(:).right];
    center_points =[p.checkpoints(:).center];
    % x,y coordinates of left and right track borders
    Lx=left_points(1,checkpoint_indices(i));
    Ly=left_points(2,checkpoint_indices(i));
    Rx=right_points(1,checkpoint_indices(i));
    Ry=right_points(2,checkpoint_indices(i));
    
    % identify which is first point to left that is occupied
    for j=j_optm:-1:1
        k = grid_isOccupiedBy(i,j); % car index
        if k~=UNOCCUPIED
            % this is the first occupied point to left of car  
            % get states of car
            xypsi = obstacle(:,k);
            % find closest intersect point to right track border
            new_b_L(i,:) = FindClosestIntersect(Rx,Ry,Lx,Ly,xypsi,carlength,carwidth);
            break;
        elseif j==1
            % has hit track border
            new_b_L(i,:) = [Lx Ly];
        end
    end
    
    % identify which is first point to right that is occupied
    for j=j_optm+1:grid_width
        k = grid_isOccupiedBy(i,j); % car index
        if k~=UNOCCUPIED
            % this is the first occupied point to right of car
            
            % get states of car
            xypsi =  obstacle(:,k);
                        
            % find closest intersect point to left track border
            new_b_R(i,:) = FindClosestIntersect(Lx,Ly,Rx,Ry,xypsi,carlength,carwidth);       
            break;
        elseif j==grid_width
            % has hit track border
            new_b_R(i,:) = [Rx Ry];
        end
    end
end

new_b = [new_b_L new_b_R];

% % set the last border as a "trivial" last border
% last_border=[1;0;LAST_BORDER_XMAX]; % this constraint means: x<=LAST_BORDER_XMAX
% 
% if run_Find_Nearest_Borders_mex
% if((norm(new_b - new_b_mex)>1e-12) || norm(last_border-last_border_mex)>1e-12)
%    asdf=5;
%    error('mex function Find_Nearest_Borders_mex not working.');
% end
% end

end

function last_border=calculate_last_border(P1,P2)
% Calculates the half-plane constraint: [a1 b2][x;y]<=b
% last_border=[a1;a2;b]
%
% The half-plane is defined by the line passing through the midpoints
% between P1 and P2 and perpendicular to the line through P1 and P2.
%
% choose your ordering of P1 and P2 s.t. P1 will fulfill the constraint, 
% while P2 does not.

x1=P1(1);
x2=P2(1);
y1=P1(2);
y2=P2(2);

a1=x2-x1;
a2=y2-y1;
b=0.5*(-x1^2-y1^2+x2^2+y2^2);

% check if direction of half-plane is correct
if a1*x1+a2*y1<=b
    last_border=[a1;a2;b];
else
    % change the direction of the half-plane
    last_border=[-a1;-a2;-b];
end

end


function intersect = FindClosestIntersect(xA,yA,xB,yB,xypsi,carlength,carwidth)
% finds the intersection between the 4 sides of the car and the line AB.
% in general there are up to two: this finds the intersection closest to A.
%
% note: assumed other car could be backwards, hence need to check all 4 sides
% 
% http://paulbourke.net/geometry/lineline2d/

center=xypsi(1:2);
psi=pi/2;

cospsi=cos(psi);
sinpsi=sin(psi);

R = [cospsi -sinpsi;sinpsi cospsi]; % rotation matrix

% first calculate position of 4 corners
% corners is [2]by[4]
w=carlength/2; l=carwidth/2;
corners(:,1) = center + R*[ l; w]; %FL
corners(:,2) = center + R*[ l;-w]; %FR
corners(:,3) = center + R*[-l;-w]; %BR
corners(:,4) = center + R*[-l; w]; %BL
corners(:,5) = center + R*[ l; w]; %FL
% hold on
% plot(corners(1,:),corners(2,:))

% fAB is the normalized fraction from point A to point B that the intersect
% lies on.
fABinit=2.0;
fABmin=fABinit; % initialize (need only be larger than 1)


% for each adjacent pair of corners, find intersect with line AB
for i=1:4
    if i==4
        j=1;
    else
        j=i+1;
    end
    xi = corners(1,i);
    yi = corners(2,i);
    xj = corners(1,j);
    yj = corners(2,j);
%     hold on
%     plot([xi,xj],[yi,yj],'LineWidth',10)
%     hold on 
%     plot([xA,xB],[yA,yB],'LineWidth',10)
    % http://paulbourke.net/geometry/lineline2d/
%     numerij = (xB-xA)*(yA-yi) - (yB-yA)*(xA-xi);
%     denom = (yj-yi)*(xB-xA) - (xj-xi)*(yB-yA);
    
    numerij = (yA-yB)*(xi-xA) + (xB-xA)*(yi-yA);
    denom = (xB-xA)*(yi-yj) - (xi-xj)*(yB-yA);
    
    if denom==0
        %two lines are parallel
        continue;
    end
    
    % normalized distance between point i and point j where intersect lies
    fij = numerij/denom; 
    
    if ~( fij<=1 && fij>=0 )
        % intersect does not lie between point i and j
        continue;
    end
    
    % intersect lies between i,j so see how close it is to A
%     numerAB = (xj-xi)*(yA-yi) - (yj-yi)*(xA-xi);

    numerAB = (yi-yj)*(xi-xA) + (xj-xi)*(yi-yA);
    
    fAB = numerAB/denom;
    
    if ~( fAB<=1 && fAB>=0 )
        % intersect does not lie between point A and B
        continue;
    end
    
    if fAB<fABmin
        fABmin=fAB;
    end
end

if fABmin==fABinit,error('no intersect between line and polygon'),end
intersect = [xA;yA] + fABmin*[xB-xA;yB-yA];
intersect=intersect'; % return a row vector

% plotaid_FindClosestIntersect(xA,yA,xB,yB,corners,intersect)

end
