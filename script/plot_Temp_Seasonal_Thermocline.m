clear all; close all; clc;

%% Data
% HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
% % sourceDir = strcat(HOME,'input/winter/');
% % sourceDir = strcat(HOME,'input/spring/');
% % sourceDir = strcat(HOME,'input/summer/');
% % sourceDir = strcat(HOME,'input/fall/');
% 
% cd(strcat(sourceDir));
% 
% 
% files = dir(fullfile(sourceDir, '*.nc'));
% FramStrait = zeros(33, length(files)); % Pre-allocating matrix
% BarentsSea = zeros(33, length(files)); 
% BeaufortGyre = zeros(33, length(files)); 
% 
% for ii = 1:length(files)
%     fname = files(ii).name;
%     readncfile;
%     tile = 7; % arctic slice
%     depth = 1:33; Z = Z(depth); % upto 2 km
%     THETA = squeeze(THETA(:,:,tile,depth)); % considering Arctic slice upto 2 km
%     THETA = flip(THETA);
%     clear A i i_g j j_g k k_l k_p1 k_u ll XC_bnds ...
%         SALT XC YC nYG YC_bnds Z_bnds Zl Zp1 Zu XG YG
%     FramStrait(:,ii) = squeeze(THETA(80,53,:));
%     BarentsSea(:,ii) = squeeze(THETA(53,27,:));
%     BeaufortGyre(:,ii) = squeeze(THETA(3,66,:));
% end

%% load thermocline data directly
clear all; clc;
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
% Ask the user for the season
seasonInput = input('Enter the season (w for winter, sp for spring, sm for summer, f for fall): ', 's');
% Determine the source directory and filename based on the user input
switch seasonInput
    case 'w'
        sourceDir = strcat(HOME, 'input/winter/');
        fileName = 'winter_thermocline.mat';
    case 'sp'
        sourceDir = strcat(HOME, 'input/spring/');
        fileName = 'spring_thermocline.mat';
    case 'sm'
        sourceDir = strcat(HOME, 'input/summer/');
        fileName = 'summer_thermocline.mat';
    case 'f'
        sourceDir = strcat(HOME, 'input/fall/');
        fileName = 'fall_thermocline.mat';
    otherwise
        error('Invalid season input. Please enter w, sp, sm, or f.');
end
% Change to the appropriate directory
cd(sourceDir);
% Load the data file for the selected season
load(fileName);


%% Plot
figure(1); clf; clc;
t = tiledlayout(1,3);
t.TileSpacing = 'tight';
t.Padding = 'compact';
L = 33; % L=23 --> upto 500 m

% Fram Strait
ax1 = nexttile;
cmap1 = slanCM('dense'); 
colormap(ax1, cmap1);
numMonths = size(FramStrait, 2); 
for i = 1:numMonths 
    colorIndex = floor((i-1) * size(cmap1, 1) / numMonths) + 1;
    plotColor = cmap1(colorIndex, :); ylim([-1400 -5])
    plot(ax1, FramStrait(1:L,i), Z(1:L), '-', 'linewidth', 1, 'Color', plotColor); hold on;
end
PLOT(ax1, 'Fram Strait', [-3 2], 'Temperature (°C)', 'Depth (m)', numMonths, cmap1);

% Barents Sea
ax2 = nexttile;
cmap2 = slanCM('matter'); 
colormap(ax2, cmap2);
for i = 1:numMonths 
    colorIndex = floor((i-1) * size(cmap2, 1) / numMonths) + 1;
    plotColor = cmap2(colorIndex, :); ylim([-1400 -5])
    plot(ax2, BarentsSea(1:L,i), Z(1:L), '-', 'linewidth', 1, 'Color', plotColor); hold on;
end
PLOT(ax2, 'Barents Sea', [-3 2], 'Temperature (°C)', 'Depth (m)', numMonths, cmap2);

% Beaufort Gyre
ax3 = nexttile;
cmap3 = flipud(slanCM('nuclear')); 
colormap(ax3, cmap3);
for i = 1:numMonths 
    colorIndex = floor((i-1) * size(cmap3, 1) / numMonths) + 1;
    plotColor = cmap3(colorIndex, :); ylim([-1400 -5])
    plot(ax3, BeaufortGyre(1:L,i), Z(1:L), '-', 'linewidth', 1, 'Color', plotColor); hold on;
end
PLOT(ax3, 'Beaufort Gyre', [-3 2], 'Temperature (°C)', 'Depth (m)', numMonths, cmap3);

% PLOT function
function PLOT(ax, titleText, xLimits, xLabel, yLabel, numMonths, cmap)
    set(ax, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3, 'FontSize', 15);
    box on; title(ax, titleText);
    xlim(ax, xLimits); xlabel(ax, xLabel); ylabel(ax, yLabel);
    colormap(ax, cmap); caxis(ax, [1 numMonths]); 
    cb = colorbar(ax,'horizontal'); ylabel(cb, 'Number of months');
    cb.Ticks = linspace(1, numMonths, 10);
    cb.TickLabels = round(linspace(1, numMonths, 10));
    set(gcf, 'Position', [-2543,831,1630,506]);
end

