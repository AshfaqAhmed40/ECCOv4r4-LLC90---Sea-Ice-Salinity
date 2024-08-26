clear all; close all; clc;

%% Data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input'));
fname = "ExampleFile-Salinity-Temperature-2017-12-01.nc";
readncfile;
tile = 7; % arctic slice
depth = 1:33; Z = Z(depth); % upto 2 km
SALT = squeeze(SALT(:,:,tile,depth)); % considering Arctic slice upto 2 km
SALT = flip(SALT);
n = 1;

clear A fname i i_g j j_g k k_l k_p1 k_u ll ...
    XC_bnds YG YC_bnds Z_bnds Zl Zp1 Zu XG


%% Thgermohaline profile
% Fram Strait (80,53)
FramStrait_SSS = squeeze(SALT(80,53,:));

% Barents Sea (53,27)
BarentsSea_SSS = squeeze(SALT(53,27,:));

% Beaufort Gyre (3,66)
BeaufortGyre_SSS = squeeze(SALT(3,66,:));


%% Plot thermohaline profile
% Sample data
salinity = {FramStrait_SSS;BarentsSea_SSS;BeaufortGyre_SSS};

figure(1), clf;
t = tiledlayout(1,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile

for i = 1:3
    plot(salinity{i}, Z, '-','linewidth',5); hold on
end

set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3);
box on; set(gca, 'GridLineStyle', '-','LineWidth', 3);
hold on; set(gca, 'Layer', 'top');
title('Halocline'); xlim([30 36])
set(gca,'XTickLabel'); set(gca,'YTickLabel');
xlabel('Salinity (pss)'); ylabel('Depth (m)');
set(gca,'FontSize', 30)
legend('Fram Strait','Barents Sea','Beaufort Sea',...
    'fontsize',30,'location','southwest')
set(gcf, 'Position', [-2433,355,844,982]);
