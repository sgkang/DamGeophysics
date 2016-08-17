% A short Demo to test the 2D electrical resistance tomography simulator
% Created by: Alireza Aghasi, Georgia Tech
% Email: aghasi@gatech.edu
% Created: Summer 2013

%%
% This program is based upon work supported by the National Science 
% Foundation under grant no. EAR 0838313. Any opinions, findings, and
% conclusions or recommendations expressed in this material are those of 
% the authors and do not necessarily reflect the views of theNational 
% Science Foundation
%%

%==========================================================================
clc;
clear all;
% We consider a rectangular domain of length 10m and depth 5m, where the
% top portion corresponds to the air interface (Neumann B.C.) and the
% boundary conditions on the remaining sides is somehow that models an
% infinte half space

domainDimX = 10;
domainDimY = 5;

% The domain is discretized into 40 by 50 pixels and the conductivity of
% the domain is defined by a 40x50 matrix. The conductivity is 0.1
% throughout the domain except a rectangular region which has conductivity
% 0.3

domainGridSize = [40 50];
cond_disb = 0.1 * ones(40,50);
cond_disb(5:30, 15:35) = 0.3;
% adding some inhomogenuity
cond_disb = cond_disb + .02*randn(40, 50);

sigmaBackground = 0.1;

% Although we apply absorbing boundary conditions on the the left, right
% and bottom portion of the domain, for more accurate results we can
% consider few pixels of domain extension. Increasing the number of pixel
% extensions increases the simulator accuracy, but of course at the expense
% of computational load

ExNo_sides = 3;
ExNo_bottom = 3;

% The sensor configurations and the experiment setup are passed to the
% program through the mat files:
% - dipole_configuration
% - dipole_voltages
% - nodes
% NOTE THAT THE COORDINATES IN THE MAT FILE NODES ARE AND SHOULD BE
% NORMALIZED TO THE DOMAIN SIZE. THE SENSORS ARE CONSIDERED TO BE AROUND A
% RECTANGULAR DOMAIN WHICH EXTENDS FROM -0.5 TO 0.5 IN THE X-DIRECTION AND
% EXTENDS FROM 0 TO -1 IN THE Y-DIRECTION

% Now we are all set to generate the ERTParams pack:

ERTParams = paramPackGenerator( domainGridSize, sigmaBackground, domainDimX,...
    domainDimY, ExNo_sides, ExNo_bottom);

% And we can proceed with generating the ERT data
[ data, ~ ] = ERT2D (cond_disb, 1 , ERTParams);

% plotting the problem setup:
% =========================================================================
load nodes;
load dipole_configuration;
figure;

x = linspace(-domainDimX / 2, domainDimX / 2, domainGridSize(2));
y = linspace(0, -domainDimY, domainGridSize(1));

[X Y] = meshgrid(x,y);
surf(X,Y, cond_disb); view([0,90]);
axis([-domainDimX/2-1, domainDimX/2+1, -domainDimY-.5, .5]);
hold on;
box on;
xlabel('x (meters)');
ylabel('y (meters)');
for i = 1 : size(nodes,1)
    scatter3(domainDimX * nodes(i,1),domainDimY * nodes(i,2),1);
    hold on;
    text(domainDimX * nodes(i,1)+.3 * sign(nodes(i,1)), domainDimY * ...
                      nodes(i,2) + .1,1, num2str(i),'color',[.3 .3 .3]);
    pause(.05);
end
for i = 1 : size(dipole_configuration,1)
    line([domainDimX * nodes(dipole_configuration(i,1),1) domainDimX * ...
                                  nodes(dipole_configuration(i,2),1) ],...
        [domainDimY * nodes(dipole_configuration(i,1),2) domainDimY * ...
            nodes(dipole_configuration(i,2),2) ] , [1 1], 'color','yellow');
    hold on;
    pause(.05);
end
title('Problem Setup: Sensors indicated by dots & index, Yellow lines represent the dipoles per experiment');



% plotting the data
% =========================================================================
figure;
plot(data,'MarkerSize',4,'Marker','*',...
    'Color',[0.0392156876623631 0.141176477074623 0.415686279535294]);
title('Potential Measurements Over the Sensors at Different Experiments');

% A Jacobian test plot
% =========================================================================
% A Jacobian accuracy check: we calculate the derivative of the 
% measurements to a given pixel, once by perturbing the pixel value and
% approximating the derivative and once by extracting the corresponding
% column of the Jacobian matrix and plot them for comparison

% pixel index, here can be any number between 1 to 40*50 
i = 100;


del = zeros(size(cond_disb));
del(i) = .00001;
[datap J] = ERT2D (cond_disb + del, 2, ERTParams);

figure;

subplot(2,1,1);plot((datap - data) / .00001);
title('Sensitivity of measurements to a specific pixel using fine perturbation')
subplot (2,1,2); plot(J(:,i),'r');
title('Sensitivity vector extracted from the Jacobian matrix')

