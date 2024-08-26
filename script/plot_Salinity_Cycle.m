clear all; close all; clc;

%% Data
% HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
% cd(strcat(HOME,'input/data/'));
% sourceDir = (strcat(HOME,'input/data/'));
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
%     SALT = squeeze(SALT(:,:,tile,depth)); % considering Arctic slice upto 2 km
%     SALT = flip(SALT);
%     clear A i i_g j j_g k k_l k_p1 k_u ll XC_bnds ...
%         THETA XC YC nYG YC_bnds Z_bnds Zl Zp1 Zu XG YG
%     FramStrait(:,ii) = squeeze(SALT(80,53,:));
%     BarentsSea(:,ii) = squeeze(SALT(53,27,:));
%     BeaufortGyre(:,ii) = squeeze(SALT(3,66,:));
% end
%  
%% Load data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input/'))
load yearly_salinity.mat
FramStrait = FramStrait_year;
BarentsSea = BarentsSea_year;
BeaufortGyre = BeaufortGyre_year;

%% Plot salinity cycle
figure(1), close all;
% Considering first three layers
L = 6;

% Create a datetime vector for each month from 1992 to 2017
startDate = datetime(1992, 1, 1);  % Start in January 1992
endDate = datetime(2017, 12, 1);   % End in December 2017
months = dateshift(startDate, 'start', 'month', 0:311);

FramStrait_m = mean(FramStrait(1:L,:),1);
FramStrait_std = std(FramStrait(1:L,:),1);
BarentsSea_m = mean(BarentsSea(1:L,:),1);
BarentsSea_std = std(BarentsSea(1:L,:),1);
BeaufortGyre_m = mean(BeaufortGyre(1:L,:),1);
BeaufortGyre_std = std(BeaufortGyre(1:L,:),1);

% Define x-axis and the upper and lower bounds for the patch
x = 1:size(FramStrait, 2);
FramStrait_up = FramStrait_m + FramStrait_std;
FramStrait_low = FramStrait_m - FramStrait_std;
BarentsSea_up = BarentsSea_m + BarentsSea_std;
BarentsSea_low = BarentsSea_m - BarentsSea_std;
BeaufortGyre_up = BeaufortGyre_m + BeaufortGyre_std;
BeaufortGyre_low = BeaufortGyre_m - BeaufortGyre_std;

% Create the plot
plot(months, FramStrait_m, 'b', 'LineWidth', 3); hold on;
plot(months, BarentsSea_m, 'k', 'LineWidth', 3); hold on;
plot(months, BeaufortGyre_m, 'r', 'LineWidth', 3); hold on;

% Create the patch
fill([months fliplr(months)], [FramStrait_up fliplr(FramStrait_low)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([months fliplr(months)], [BarentsSea_up fliplr(BarentsSea_low)], ...
    'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([months fliplr(months)], [BeaufortGyre_up fliplr(BeaufortGyre_low)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
       

% Format the x-axis
ax = gca;
ax.XTick = months(1:6:end);  % Adjust this to change the frequency of labels
ax.XTickLabelRotation = 45;  % Rotate labels for better readability
set(gcf, 'Position', [-2257,775,1630,506]); 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold on
set(gca,'XTickLabel');set(gca,'YTickLabel');
set(gca,'XTick');set(gca,'YTick');
title(sprintf('Mean salinity upto  %.2f m deep',Z(L))); 
ylabel('Salinity (pss)');ylim([28 36])
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

%% Spectrum analysis
yg = zeros(3,size(FramStrait_m,2));
yg(1,:) = FramStrait_m';
yg(2,:) = BarentsSea_m';
yg(3,:) = BeaufortGyre_m';
obsv = datetime(months);


[pxx1,f] = plomb(yg(1,:),obsv,[],10,'power'); hold on
[pxx2,f] = plomb(yg(2,:),obsv,[],10,'power'); hold on
[pxx3,f] = plomb(yg(3,:),obsv,[],10,'power'); hold on
f = f*86400*365;


figure(2), clf;
loglog(f,pxx1,'b-','linewidth',2); hold on
loglog(f,pxx2,'k-','linewidth',2); hold on
loglog(f,pxx3,'r-','linewidth',2); hold on
set(gcf, 'Position', [-2257,439,1630,256]); 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold on
set(gca,'XTickLabel');set(gca,'YTickLabel');
set(gca,'XTick');set(gca,'YTick');
ylabel('Magnitude'); 
xlabel('Frequency (Year^{-1})');
set(gca, 'FontSize', 16, 'FontWeight', 'bold');
xlim([10e-3 10]);
ylim([10e-12 1000]); 
title('Power Spectrum and Prominent Peak')


