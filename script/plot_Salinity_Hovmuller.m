clear all; close all; clc;
 
%% Load data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input/'))
load yearly_salinity.mat

%% Best fit graph (95% confidence level)
depths = Z(1:19);
FramStrait = FramStrait(1:size(depths,1),1:288); 
BarentsSea = BarentsSea(1:size(depths,1),1:288);
BeaufortGyre = BeaufortGyre(1:size(depths,1),1:288);

clr = flipud(brewermap(15,'BrBG'));
clri= interp1(1:1:15,clr,1:0.25:15,'linear');

% Assuming FramStrait is your 20x288 matrix

% Create the tiled layout
t = tiledlayout(1,8);
t.TileSpacing = 'tight';
t.Padding = 'compact';

% Loop over each 3-year period
for i = 1:8
    % Calculate the starting and ending indices for each 3-year chunk
    startIdx = 1 + (i-1)*36;
    endIdx = i*36;

    % Select the data for the current 3-year period
    currentData = BeaufortGyre(:, startIdx:endIdx);

    % Create a new tile
    nexttile
    
    % contour
    contourf(startIdx:endIdx, depths, currentData,15,'edgecolor','k'); 
    shading interp; 
    colormap(gca, clri); hold on; 
    box on; set(gca, 'Layer', 'top');
%     [C,L] = contour(startIdx:endIdx, depths, currentData, 15,'linewidth',1.5,'LineColor','k');
%     clabel(C, L, 'FontSize', 12); L.LevelList = round(L.LevelList,1);

    if i~=1
        set(gca,'YTickLabel',[])
    end
    if i==1
        ylabel('Depth (m)');
    end
    % Set the title for each tile
    startYear = 1992 + (i-1)*3;
    endYear = startYear + 2;
    
    
    title(sprintf('%d-%d', startYear, endYear));
    set(gcf, 'Position', [-5856,941,2419,342]); 
    set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold on
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel');
    %caxis([33.5 35.5])
    caxis([29 34]) % for Beaufort
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
end

c = colorbar;
c.FontSize = 30;
c.FontWeight = 'bold';
c.Layout.Tile = 'east';