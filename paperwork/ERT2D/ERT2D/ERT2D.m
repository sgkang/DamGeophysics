function [ sens_data Jacobian ] = ERT2D (cond_disb, evaluation, ERTParams)

% [sens_data Jacobian] = ERT2D(cond_disb, evaluation, ERTParams)

% ERT2D is a Matlan function that simulates electrical resistance
% tomography in a 2D rectangular domain. The setup is rather similar to
% Geophysics application where the top side of the domain corresponds to
% the air interface (Neumann B.C.) and the sides undergo a mixed boundary
% condition to simulate infinite half-space.

% Aside from solving the Poisson's equation over a set of predefined sensor
% configurations, ERT2D also calculates the sensitivity (Jacobian) of the
% model outputs to the domain conductivity. The program uses an Adjoint
% technique to calculate the sensitivities.

% =========================================================================
% For mode details about the Adjoint technique please refer to the
% following papers:

% * Aghasi, A. and Miller, E., ?Sensitivity Calculations for the Poisson?s
% Equation via the Ad- joint Field Method?, IEEE Geoscience and Remote
% Sensing Letters, Vol. 9, no. 2, pp. 237?241, 2012

% * Polydorides, N., Aghasi, A. and Miller, E., ?High-Order Regularized
% Regression in Electrical Impedance Tomography?, SIAM Journal on Imaging
% Sciences, Vol. 5, no. 3, pp. 912?943, 2012
% =========================================================================

% The boundary conditions are applied through the gradient matrices
% calculated in paramPackGenerator program. You can easily modify the code
% to handle other type of boundary conditions and setups

% INPUTS:
% =========================================================================
% cond_disb: is a matrix that contains the conductivity values over the
% rectangular domain. Clearly the size of cond_disb will be used to
% determine the number of grids in solving the Poisson's equation
% numerically over the 2D region. Note that this code is merely designed
% for resistance tomography and therefore cond_disb needs to be real,
% however, with slight modifications this code is applicable to the complex
% case, i.e., impedance tomography

% evaluation: only takes values of 1 and 2. If 1, the program only
% calculates the voltage values at the designated sensors. If 2, the
% program calculates the potential values, as well as the sensitivity
% (derivative) of each sensor potential to a given pixel within the domain

% ERTParams: To make the program faster, all variables and matrices that
% need to be calculated only once are encapsulated inside this single
% structure

% OUTPUTS:
% =========================================================================
% sens_data: is a vector that contains the potential values at the
% designated sensors. If there are ne experiments where at each experiment
% nm sensors measure the potential at different points in the domain,
% sens_data would be of length nm*ne

% Jacobian: is a matrix the columns of which correspond to the sensitivity
% (derivative) of the sensor potentials to a certain pixel in the domain.
% This is an important peice of information in performing the inversion.
% If there are ne experiments where at each experiment nm sensors measure
% the potential at different points in the domain, and the domain contains
% np pixels, Jacobian would be of size nm*ne by np

% Please see the user's guide for more technical details


%% Feel free to redistribute and/or modify this software as long as you
% make reference to the original source (or cite the corresponding papers):
% Created by: Alireza Aghasi, Georgia Tech
% Email: aghasi@gatech.edu
% Created: Summer 2013
%%
%%
% This program is based upon work supported by the National Science 
% Foundation under grant no. EAR 0838313. Any opinions, findings, and
% conclusions or recommendations expressed in this material are those of 
% the authors and do not necessarily reflect the views of theNational 
% Science Foundation
%%

if nargin ~= 3
    error('ERT2D:argChk', ['Wrong number of inputs: make sure to use ERT2D in the following format:\n'...
        '[sens_data Jacobian] = ERT2D(cond_disb, evaluation, ERTParams)']);
end

% decapsulating the variables from ERTParams
disbe = ERTParams.disbe;
G1x = ERTParams.G1x;
G2x = ERTParams.G2x;
G1y = ERTParams.G1y;
G2y = ERTParams.G2y;

boundary_conductivity = ERTParams.boundary_conductivity;
ind_b_nodes = ERTParams.ind_b_nodes;
sensInd = ERTParams.sensInd;
sources_index = ERTParams.sources_index;
f = ERTParams.f;
Q = ERTParams.Q;
Ix = ERTParams.Ix;
Iy = ERTParams.Iy;

ExNo_sides = ERTParams.ExNo_sides;
m = ERTParams.m;
n = ERTParams.n;

m0 = ERTParams.m0;
n0 = ERTParams.n0;

nsrc = ERTParams.nsrc;
nsens = ERTParams.nsens;


% Cehcking if the size of cond_disb is m0 by n0
siz = size(cond_disb);

if m0 ~= siz(1) || n0 ~= siz(2)
    error('The size of conductivity matrix does not match m0 and n0 in ERTParams.');
end


switch evaluation
    case 1
        
        % placing the given conductivity distribution in the middle of
        % the extended domain
        disbe(2 : m0 + 1,(ExNo_sides + 2) : (n0 + ExNo_sides + 1)) = cond_disb;
        
        temp_x = speye(m * (n + 1));
        temp_y = speye((m + 1) * n);
        
        data = zeros(nsens,nsrc);
        
        
        % Checking if Matlabpool is available for parallel loops
        
        if matlabpool('size') ~= 0
            
            parfor i = 1 : nsrc
                disbee = disbe;
                disbee(ind_b_nodes) = boundary_conductivity(:,i);
                
                disbfx=[disbee,disbee(:,end)];
                disbfy=[disbee;disbee(end,:)];
                
                % The finite difference matrix, basically the r.h.s of the
                % Poisson's equation in discrete form
                K = G2x*(spdiags(disbfx(:),0,temp_x) * G1x)...
                    + G2y * (spdiags(disbfy(:),0,temp_y) * G1y);
                
                % and here we determine the voltages
                v = K \ f(:,sources_index(i));
                
                data(:,i) = v(sensInd);
                
                % In case you like to see how the Green's functions over the
                % extended domain look, uncomment the following, but I warn you
                % it would be tons of figures!
                % figure;surf(reshape(full(v),[m n]));pause(.1);
            end
            
        else
            
            for i = 1 : nsrc
                disbee = disbe;
                disbee(ind_b_nodes) = boundary_conductivity(:,i);
                
                disbfx=[disbee,disbee(:,end)];
                disbfy=[disbee;disbee(end,:)];
                
                % The finite difference matrix, basically the r.h.s of the
                % Poisson's equation in discrete form
                K = G2x*(spdiags(disbfx(:),0,temp_x) * G1x)...
                    + G2y * (spdiags(disbfy(:),0,temp_y) * G1y);
                
                % and here we determine the voltages
                v = K \ f(:,sources_index(i));
                
                data(:,i) = v(sensInd);
                
                % In case you like to see how the Green's functions over the
                % extended domain look, uncomment the following, but I warn you
                % it would be tons of figures!
                % figure;surf(reshape(full(v),[m n]));pause(.1);
            end
        end
        
        
        sens_data = Q*data(:);
        Jacobian = [];
        
    case 2
        
        % placing the given conductivity distribution in the middle of
        % the extended domain
        disbe(2 : m0 + 1,(ExNo_sides + 2) : (n0 + ExNo_sides + 1)) = cond_disb;
        
        temp_x = speye(m * (n + 1));
        temp_y = speye((m + 1) * n);
        
        s = speye(m0 * n0);
        sensitivity = zeros(nsens,m0 * n0,nsrc);
        data = zeros(nsens,nsrc);
        
        % Checking if Matlabpool is available for parallel loops
        
        if matlabpool('size') ~= 0
            
            parfor i = 1 : nsrc
                disbee = disbe;
                disbee(ind_b_nodes) = boundary_conductivity(:,i);
                
                disbfx = [disbee,disbee(:,end)];
                disbfy = [disbee;disbee(end,:)];
                
                % The left side finite difference matrix
                K = G2x*(spdiags(disbfx(:),0,temp_x) * G1x)...
                    + G2y * (spdiags(disbfy(:),0,temp_y) * G1y);
                
                % Finally!
                V = K \ f;
                
                data(:, i) = V(sensInd,sources_index(i));
                
                % The gradient of the voltage later to be used in calculating
                % the Jacobian through the Adjoint technique
                
                sens_x = G1x * V;
                sens_y = G1y * V;
                
                sensitivity(:,:,i)=full((spdiags(sens_x(Ix,sources_index(i)),0,s) * sens_x(Ix,:) +...
                    spdiags(sens_y(Iy,sources_index(i)),0,s) * sens_y(Iy,:) ))';
            end
            
        else
            
            for i = 1 : nsrc
                disbee = disbe;
                disbee(ind_b_nodes) = boundary_conductivity(:,i);
                
                disbfx = [disbee,disbee(:,end)];
                disbfy = [disbee;disbee(end,:)];
                
                % The left side finite difference matrix
                K = G2x*(spdiags(disbfx(:),0,temp_x) * G1x)...
                    + G2y * (spdiags(disbfy(:),0,temp_y) * G1y);
                
                % Finally!
                V = K \ f;
                
                data(:, i) = V(sensInd,sources_index(i));
                
                % The gradient of the voltage later to be used in calculating
                % the Jacobian through the Adjoint technique
                
                sens_x = G1x * V;
                sens_y = G1y * V;
                
                sensitivity(:,:,i)=full((spdiags(sens_x(Ix,sources_index(i)),0,s) * sens_x(Ix,:) +...
                    spdiags(sens_y(Iy,sources_index(i)),0,s) * sens_y(Iy,:) ))';
            end
        end
        
        
        sens_flat = zeros(nsens * nsrc,m0 * n0);
        
        for i = 1 : nsrc
            sens_flat((i-1) * nsens + 1 : nsens * i,:) = sensitivity(:,:,i);
        end
        
        sens_data = Q*data(:);
        Jacobian = Q*sens_flat;
        
    otherwise
        error('ERT2D:argChk', ['Wrong choice for the evaluation parameter, make sure to\n'...
            'use either 1 or 2.']);
end




