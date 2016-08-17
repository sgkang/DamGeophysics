function ERTParams = paramPackGenerator( domainGridSize, sigmaBackground, domainDimX, domainDimY, ExNo_sides, ExNo_bottom)

% ERTParams = paramPackGenerator( domainGridSize, sigmaBackground, ...
%                          domainDimX, domainDimY, ExNo_sides, ExNo_bottom)

% paramPackGenerator generates a structure of vaiables (ERTParams), which
% is used in the electrical resistance tomography simulator, ERT2D. 
% The pack contains the parameter settings and matrices that need to be
% calculated only once. Among the most important ones are the gradient 
% matrices.

% INPUTS:
% =========================================================================
% domainGridSize: The size of the rectangular domain where the Poisson
% equation is solved, in pixels, e.g.: [50 50] 

% sigmaBackground: The conductivity of the domain far from the sensors. You
% will specifically provide the domain conductivity map inside the ERT2D
% program. However this variable (which is only a constant) is the 
% conductivity of the homogenous background and merely used in extending 
% the imaging domain for more accurate solutions and realistic boundary 
% conditions

% domainDimX: The horizontal length of the domain in meters

% domainDimY: The vertical (depth) of the domain in meters

% ExNo_sides: The number of pixels in extending the domain horizontally to 
% apply an appropriate mixed boundary condition that models infinite half-
% space. The domain will be extended as much as ExNo_sides pixels from each
% side (left and right)

% ExNo_bottom: The number of pixels in extending the domain vertically to 
% apply an appropriate mixed boundary condition that models infinite half-
% space. The domain will be extended as much as ExNo_bottom only at the
% bottom of the domain. Note that the top portion of the domain undergoes a
% Neumann B.C. and no extension is required.

% OUTPUT:
% =========================================================================
% ERTParams is a structure that contains several variables and matrices
% needed to be calculated only once for the ERT2D simulator

% Please see the user's guide for more details

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

if nargin ~= 6
    error('ERT2D:argChk', ['Wrong number of inputs: make sure to use paramPackGenerator in the following format:\n'...
            'ERTParams = paramPackGenerator( domainGridSize, sigmaBackground, domainDimX, domainDimY, ExNo_sides, ExNo_bottom)']);
end


% loading the sensor configuration files
load nodes.mat;
load dipole_configuration.mat;
load dipole_voltages.mat;


% determining the number of grids
m0 = domainGridSize(1);
n0 = domainGridSize(2);

% determining the grid spacing
dx = domainDimX / (n0 + 1);
dy = domainDimY / (m0 + .5);

% determing the coordinates of the grid points
hx0 = linspace(- domainDimX/2 - dx / 2, domainDimX/2 + dx / 2, n0 + 3);
hy0 = linspace(- dy / 2, domainDimY, m0 + 2);

% extending the domain to a larger domain to later apply the
% absorbing boundary conditions (Robin B.C.)
nodes_x = [(hx0(1) - ExNo_sides * dx) : dx : (hx0(1) - dx), hx0, ...
                 (-hx0(1) + dx) : dx : (- hx0(1) + ExNo_sides * dx) ]';
             
nodes_y = [hy0, (hy0(end) + dy): dy : (hy0(end) + ExNo_bottom * dy) ]';

% hx and hy contain the sizes of the grid blocks in each
% direction x and y
hx = dx * ones(size(nodes_x));
hy = dy * ones(size(nodes_y));

% n and m are the number of grid points in the extended domain
n = n0 + 2 + 2 * ExNo_sides;
m = m0 + 1 + ExNo_bottom;

% building a uniform conductivity distribution in the new extended domain
disbe = sigmaBackground * ones(m,n);

% Generating the dimension matrices
Del_hx = repmat(hx',m,1);
Del_hy = repmat(hy,1,n);

% Generating the gradient and divergence matrix operators
[G1x G2x] = grad_x(m,n,Del_hx);
[G1y G2y] = grad_y(m,n,Del_hy);


% The following lines are mainly in charge of setting the parameters for an
% absorbing (mixed) boundary condition on the sides and the bottom portion
% of the domain. For more details please refer to the manual.
%
%
% ************* absorbing B.C. *************
%
%
% Extracting the boundary nodes (only boundaries with mixed B.C.)
I = [1 : m - 1, m * ones(1,n), m-1 : -1 : 1 ];
J = [ones(1,m ),2 : n ,n * ones(1,m-1)];


% Extracting the single-index representation of the boundary nodes
ind_b_nodes = sub2ind([m n],I,J);

% The normal vector corresponding to such nodes
boundaries_normal_vector = [[-ones(m - 1,1);zeros(n,1);ones(m - 1,1)] , ...
    [zeros(m - 1,1);ones(n,1);zeros(m - 1,1)]];

% obaining the boundary nodes (the center points)
cnx = (nodes_x(1 : end - 1) + nodes_x(2 : end)) / 2;
cny = (nodes_y(1 : end - 1) + nodes_y(2 : end)) / 2;


% the coordinates of the boundary nodes
boundaries_coords = [[cnx(1) * ones(m - 1,1); cnx; cnx(end) * ones(m - 1,1)],...
    [cny; cny(end) * ones(n - 2,1); flipud(cny)]];

% determining the actual position of the sensors obtained from the data
% file nodes.mat (note that the positions in the file nodes.mat are in
% normalized form)
sens_nodes = [ abs(domainDimX) * nodes(:,1), - abs(domainDimY) * nodes(:,2)];

nsens = size(sens_nodes,1);
boundary_scale_factor = zeros(length(ind_b_nodes),nsens);

% determing the scale factor for every point on the extended domain that
% needs to undergo an absorbing boundary condition
for i = 1 : nsens
    ri = [boundaries_coords(:,1) - sens_nodes(i,1),boundaries_coords(:,2) ...
                                                        - sens_nodes(i,2)];
                                                    
    boundary_scale_factor(:,i) = sum((ri .* boundaries_normal_vector)')' ...
                                      ./ ((ri(:,1)) .^ 2 + (ri(:,2)) .^ 2);
end

% The scaled conductivity values are saved in a matrix where each column
% corresponds to the proper boundary conductivity values for a given sensor
boundary_conductivity = sigmaBackground ./ boundary_scale_factor;


% ************* Somehow done with the absorbing B.C. *************


% determining the sensor positions in the extended domain grid
sensor_index_x = 1 + floor((sens_nodes(:,1) - nodes_x(1)) / dx)';
sensor_index_y = 1 + floor((sens_nodes(:,2) - nodes_y(1)) / dy)';

sensInd = sub2ind([m n],sensor_index_y,sensor_index_x);

% number of experiments
ne = size(dipole_configuration,1);

% finding the sensors that are used as current the source in at least one
% experiment
sources_index = sort( union(dipole_configuration(:,1),dipole_configuration(:,2)) );
nsrc = size(sources_index,1);

% building the current source vector, f is in fact the delta function
% finding the response to which leads to the Green's function.
f = sparse(m*n,nsens);
for i = 1 : nsens
    f(sensInd(i),i) = 1;
end

% Q is the matrix which picks the data and takes care of the dipole
% combinations. In every experiment, two sensors work as diople nodes
% and the remainder serve as voltage measuring sensors. In the main program
% once the measurement vector corresponding to all sensors is obtained, an
% operation as data=Q*sensor_data(:) would give us the measurement vector
% resulted from the designated experiment and dipole combinations.

Q = sparse((nsens-2) * ne,nsrc * nsens);

k = 1;
for i = 1 : ne
    for j = 1 : nsens
        if((j == dipole_configuration(i,1)) || (j == dipole_configuration(i,2)))
            continue;
        else
            ind1 = find(dipole_configuration(i,1) == sources_index);
            ind2 = find(dipole_configuration(i,2) == sources_index);
            Q(k,(ind1 - 1) * nsens + j) = dipole_voltages(i,1);
            Q(k,(ind2 - 1) * nsens + j) = dipole_voltages(i,2);
        end
        k = k + 1;
    end
end


% Since Gx and Gy extend the disb matrix by one pixel and also the
% sensitivity values are only required for the pixels on the original disb
% and not the extended one, here we find those indices to be extracted.
Ix = [1 : m * (n + 1)]';
Iy = [1 : (m + 1) * n]';
Ix = reshape(Ix,[m,n + 1]);
Iy = reshape(Iy,[m + 1,n]);


Ix = Ix(2 : m0 + 1,(ExNo_sides + 2) : (n0 + ExNo_sides + 1) );
Iy = Iy(2 : m0 + 1,(ExNo_sides + 2) : (n0 + ExNo_sides + 1) );

Ix = Ix(:);
Iy = Iy(:);


% encapsulating all the variables inside a single structure
ERTParams.disbe = disbe;
ERTParams.G1x = G1x;
ERTParams.G2x = G2x;
ERTParams.G1y = G1y;
ERTParams.G2y = G2y;

ERTParams.boundary_conductivity = boundary_conductivity;
ERTParams.ind_b_nodes = ind_b_nodes;
ERTParams.sensInd = sensInd;
ERTParams.sources_index = sources_index;
ERTParams.f = f;
ERTParams.Q = Q;
ERTParams.Ix = Ix;
ERTParams.Iy = Iy;

ERTParams.ExNo_sides = ExNo_sides;
ERTParams.m = m;
ERTParams.n = n;

ERTParams.m0 = m0;
ERTParams.n0 = n0;

ERTParams.nsrc = nsrc;
ERTParams.nsens = nsens;





function [Gx_1 Gx_2] = grad_x(m,n,Del_hx)
% grad_x generates x-derivative matrices
% m and n are the image size and Del_hx of size m.(n+1) is the actual
% spacing between the grids
% Gx_1 and Gx_2 are the x-derivative sparse matrices. For a given image I,
% reshape(Gx_1*I(:), [m,n+1]) returns the derivative of the image in the x
% (horizontal) direction, where a Dirichlet boundary condition is imposed
% on both sides of the image. Gx_2 is again an x-derivative matrix which
% applies to the result of Gx_1*I(:). More specifically Gx_2*Gx_1*I(:)
% returns the second derivative in the x direction where the result is the
% same size as I(:).

%% Feel free to redistribute and/or modify this software as long as you
% make reference to the original source (or cite the corresponding papers):
% Created by: Alireza Aghasi, Georgia Tech
% Email: aghasi@gatech.edu
% Created: Summer 2013
%%
hx = Del_hx(:);
b = zeros(m * n,2);
Gx_1 = sparse(m * (n + 1),m * n);
b(:,1) = 1 ./ hx(1 : n * m);
b(:,2) = - 1 ./ hx(m + 1 : end);
Gx_1 = spdiags(b,[0 -m],Gx_1);

Gx_2 = sparse(m * n,m * (n + 1));
b(:,1) = 1 ./ hx(1 : n * m);
b(:,2) = -1 ./ hx(1 : n * m);
Gx_2 = spdiags(b,[m 0],Gx_2);



function [Gy_1 Gy_2] = grad_y(m,n,Del_hy)
% grad_y generates derivative matrices
% m and n are the image size and Del_hx of size (m+1).n is the actual
% spacing between the grids
% Gy_1 and Gy_2 are the y-derivative sparse matrices. For a given image I,
% reshape(Gy_1*I(:), [m+1,n]) returns the derivative of the image in the y
% (vertical) direction, where a Dirichlet boundary condition is imposed
% on the bottom and an Neumann B.C is applied on the top. Gy_2 is again a
% y-derivative matrix which applies to the result of Gy_1*I(:).
% More specifically Gy_2*Gy_1*I(:) returns the second derivative in the y
% direction where the result is the same size as I(:).

%% Feel free to redistribute and/or modify this software as long as you
% make reference to the original source (or cite the corresponding papers):
% Created by: Alireza Aghasi, Georgia Tech
% Email: aghasi@gatech.edu
% Created: Summer 2013
%%
Gy_1 = sparse((m + 1) * n,(m + 1) * n);
hy = Del_hy(:);
b1 = 1./hy(1 : (m + 1) * n);
b2 = [-1 ./ hy(2 : (m + 1) * n);0];
b1(1 : m + 1 : end) = 0;
Gy_1 = spdiags([b1 b2],[0 -1],Gy_1);
Gy_1(:,(m + 1) : (m + 1) : end) = [];


Gy_2 = sparse((m + 1) * n,(m + 1) * n);
b1 = -1 ./ hy(1 : (m + 1) * n);
b2 = [0;1 ./ hy(1 : (m + 1) * n - 1)];
Gy_2 = spdiags([b1 b2],[0 1],Gy_2);
Gy_2((m + 1) : (m + 1) : end,:) = [];

