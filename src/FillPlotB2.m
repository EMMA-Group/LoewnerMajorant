function h = FillPlotB2(BETA, c)
%% FILLPLOTB2 plots the critical beta_1/beta_2 curve
% and fills the admissible region with color c
% into a new figure and returns the figure handle
%
% INPUTS
% BETA          2 x N or N x 2 matrix containing beta_1/beta_2 pairs
% c             [OPTIONAL] color of the curve (default:  light blue)
%
% OUTPUTS
% h             handle to the new figure
%
%% related to the article
% Sharp Loewner majorants for a set of symmetric matrices
% by Felix Fritzen and Mauricio Fernandez
% [ submitted to AMCS (Aug 2018) ]
%
% CONTENT OF THIS FILE
% - example the majorization of a set of random symmetric matrices
% - illustrates how to access the different tools including
%   a discrete set majorization
%
% LICENSE
% This software is distributed under the GNU GPL v3.
% An active contribution to this software is welcome (see CONTACT below).
%
% Further licensing information is contained in the file 'license'
% distributed with this code.
%
% CONTACT
% Felix Fritzen        felix.fritzen@mechbau.uni-stuttgart.de
% Mauricio Fernandez   mauricio.fernandez@mechbau.uni-stuttgart.de
% http://www.mechbau.uni-stuttgart.de/EMMA
%
% GITHUB
% This software is hosted on GITHUB. It is accessible via
% XXXXXXXXXXXXXXXXXXXX
%
% CHANGELOG
% Aug 09, 2018      initial creation
%

if( size(BETA,1) > 2 )
    BETA=BETA';
end

if( ~exist('c','var') )
    c = [0,190,255]/255;
end

h = figure;
hold on;
fill(   [BETA(1,:), BETA(1,end), BETA(1,1)], ...
        [BETA(2,:), BETA(2,1),   BETA(2,1)], ...
        c, 'EdgeColor', 'none' );
plot( BETA(1,:), BETA(2,:), 'LineWidth', 3, 'Color', 'black' );    

end