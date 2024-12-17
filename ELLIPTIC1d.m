function [xc,nsol] = ELLIPTIC1d (nxdx,bcond,ifexact,ifplot,ifdemo) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELLIPTIC1d (nxdx,bcond,ifexact) 
%% 1d FD cell-centered (FV) solution to linear diffusion/flow equation
%% - (k u_x)_x  = f(x)  on (a,b)
%% default domain is (a,b)=(0,1)     %% HARDCODED {a}, {b}
%% k is conductivity                 %% HARDCODED {kcof}  
%% f(.)                              %% HARDCODED {rhsfun}
%% with boundary conditions given in <bcond> or {exactfun}
%% user must code/change {hardcoded data} 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT
%% <nxdx> is either nx or x position of grid nodes 
%% <bcond>: vector of 4 numbers: 
%%   <flag at x=a,value at x=a,flag at x=b,value at x=b> or {bflag_left,bval_left,bflag_right,bval_right}
%%   bflag* == 0: Dirichlet b.c.; bflag*~=0: Neumann b.c. 
%%   bval* is the Dirichlet value or the value -ku'. 
%% <ifexact>=0: no exact sol known. 
%% <ifexact>~=- ignore boundary values; plot numerical and exact; calculate err
%% <ifplot>: plot or no [default=1]
%% <ifdemo>: output some values or not [default=0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT/solution
%% returns [xc, nsol]: 
%% xc: location of cell centers; solution at cell centers
%% if <ifplot", solution [nsol] is plotted at [xc]
%% if <ifexact> sol is known, L2 grid norm of the error is calculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%% nonuniform grid, no exact solution known: Dirichlet conditions 0,1 
%ELLIPTIC1d([0 1/3 1/2 2/3 0.9 1],[0 0 0 1],0,1,0);
%% uniform grid, exact solution known; Dirichlet bc, boundary values from exact
%ELLIPTIC1d(5,[0 0 0 0],1,1,0);
%% uniform grid, exact solution known; mixed bc, boundary values from exact
%ELLIPTIC1d(5,[0 1 0 0],1,1,0);
%Elliptic1d(5, [1 0 1 0], 1, 1, 0); Newman
%% same but no plot, and info on matrices etc is provided
%ELLIPTIC1d(5,[0 1 0 0],1,1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 5, ifdemo = 0;end
    if nargin < 4, ifplot = 1;end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DATA on PDE is here
    %% Domain and grid: 
    % x are the nodes of the grid, xc are cell centers
    if length(nxdx)>1
        x = nxdx'; nx=length(x)-1; dx=diff(x); 
    else
        nx = nxdx; dx =1/nx*ones(nx,1); a=0; b=1; %%default domain is (a,b)=(0,1)
        x = 0*dx; x(1)=a; 
        for j=2:nx+1 
            x(j)=x(j-1)+dx(j-1);
        end 
    end 
    % setup the cell centers in xc
    xc = x(1:nx)+dx/2; h = max(dx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% kcof = permeability/diffusivity/conductivity 
    kcof = ones(nx,1); 
    % for j = 1:length(x)-1 (for -(ku_x)_x = f, u piecewise)
    %     if x(j)<1/3
    %         kcof(j) = 1;
    %     else
    %         kcof(j) = 10;
    %     end
    % end    
    %% Boundary data: unpack <bcond>
    bflag_left=bcond(1);
    bval_left=bcond(2);
    bflag_right=bcond(3);
    bval_right=bcond(4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CCFD = FV algorithm starts here 
    %% calculate transmissibilities k/dx on interior edges by harmonic weighting
    tx = zeros(nx+1,1); 
    for j=2:nx,  tx(j)=2/(dx(j-1)/kcof(j-1)+dx(j)/kcof(j));end
    %% account for Dirichlet boundary conditions in {tx}. For Neumann, do nothing 
    %0 is 'D', otherwise is 'N'
    if bflag_left == 0,  j = 1;         tx(j)=2/(dx(j)/kcof(j));end
    if bflag_right == 0, j = nx + 1;    tx(j)=2/(dx(j-1)/kcof(j-1));end
    %
    if ifdemo, fprintf('Transmissibilities\n'); tx',end
    %% set up the discrete diffusion matrix
    diffmat = sparse (nx,nx); 
    % for -(ku_x)_x +cu = f
    for j=2:nx % loop on interior edges                         
      gl = j-1; gr = j;
      diffmat(gl,gl) = diffmat(gl,gl) + tx(j);                                    
      diffmat(gl,gr) = diffmat(gl,gr) - tx(j);                                                              
      diffmat(gr,gl) = diffmat(gr,gl) - tx(j);                                    
      diffmat(gr,gr) = diffmat(gr,gr) + tx(j);     
    end
  
    %% contributions of bconditions to the matrix (if Dirichlet)
    if bflag_left == 0, j = 1; gr = 1; 
        diffmat(gr,gr) = diffmat(gr,gr) + tx(j);
    end
    % if bflag_left == 0, j = 1; gr = 1; 
    %   diffmat(gr,gr) = diffmat(gr,gr) + tx(j); 
    % end
    if bflag_right == 0, j = nx+1; gl = nx; 
        diffmat(gl,gl) = diffmat(gl,gl) + tx(j);  
    end
    if ifdemo, full(diffmat), pause; end
    %%%MODIFY  set up the other matrices (if reaction terms, time dependent etc.
    %for -(ku_x)_x +cu = f. 
    cmat = sparse(nx,nx); 
    ccof = ones(nx,1);
    for j=1:nx
        cmat(j,j) = ccof(j)*dx(j);
    end
    %%%MODIFY calculate the matrix of linear system to be solved in time loop
    %allmat = diffmat; % for -(ku_x)_x = f.
    allmat = diffmat + cmat; % for -(ku_x)_x +cu = f.
    %%
    if ifdemo,full(allmat), end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%MODIFY e.,g. if TIME LOOP solver is needed 
    %t = 0; n = 0; nsol = init; maxl2err = 0;
    % set up time stepping e.g.  while 1 .... loop ... end 
    clf;
    rhs = dx.* rhsfun (xc); 
        
    %% contributions of bdary conditions to qleft, qright 
    %% to be included in the rhs or in the residual
    % evaluate current qleft,qright at current time
    %LEFT
    if bflag_left == 0 %% Dirichlet contribution from dval_left
        if ifexact %% use exfun and ignore input values
            [dval_left,~]= exfun(x(1));
        else, dval_left = bval_left; end
        qleft = tx(1)*dval_left;
    else %% Neumann condition specifying -k du/dx
        if ifexact
            [~,nval_left]= exfun(x(1));
            qleft = -kcof(1)*nval_left;
        else, qleft = bval_left; %% qleft is given
        end
    end
    %RIGHT
    if bflag_right == 0 %% Dirichlet contribution from dval_right
        if ifexact %% use exfun and ignore input values
            [dval_right,~]= exfun(x(nx+1));
        else, dval_right = bval_right; end
        qright = - tx(nx+1)*dval_right;
    else %% Neumann condition specifying -k du/dx
        if ifexact
            [~,nval_right]= exfun(x(nx+1));
            qright = -kcof(nx)*nval_right;
        else, qright = bval_right; %% qleft is given
        end
    end
    if ifdemo, qleft, qright, end
    %% UPDATE rhs or residual (if residual, use different signs)        
    rhs(1) = rhs(1) + qleft;    
    rhs(nx) = rhs(nx) - qright; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE the linear system
    nsol = allmat \ rhs;
    %% Post-process: plot, compute the error if true solution available
    if ifplot
        plot(xc,nsol,'b*-'); title(sprintf('Solution')); 
    end
    %% ... plot exact solution and check the error in l2 norm
    if ifexact 
        exsol = exfun(xc);
        if ifplot
            plot(xc,exsol,'k',xc,nsol,'r*-'); 
            legend('exact','numerical','location','best');
            title(sprintf('Solution '));
        end
        %
        l2err = sqrt(sum((nsol-exsol).^2.*dx));
        fprintf('Error is err_l2=%g\n',l2err); 
    end
    %ylim([ 0 3]); choose some bounds to see the evolution nicely
    pause(0.05);

 end

%%%%%%%%%%%%%%%%%%%%%%%%% data for the problem 
function [v,dv] = exfun(x)
    v = (exp(10*x.^2))/exp(10);
    dv = (20*x.*exp(10*x.^2))/exp(10);
end

function v = rhsfun(x) 
       v = -((400*x.^2.*exp(10*x.^2))/exp(10) + (20*exp(10.*x.^2))/exp(10)) + 2*exp(10*x.^2)/exp(10);

end

%%%%%%%%%%%%%%%%%%%%%%%%% data for the problem 
% function [v, dv] = exfun(x)
%     % Piecewise function for u(x)
%     v = 0*x; 
%     dv = 0*x;
%     a = 1/3;
% 
%     % Compute v and dv for each piece
%     for i = 1:length(x)
%         if x(i) <= a
%             v(i) = -5/2 * x(i) + 1;        
%             dv(i) = -5/2;                   
%         else
%             v(i) = -5/20 * x(i) + 5/20;     
%             dv(i) = -1/4;                   
%         end
%     end
% end
% 
% function v = rhsfun(x)
%     v = 0*x;  % f(x) = 0 for all x
% end





