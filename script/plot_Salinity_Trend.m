clear all; close all; clc;
 
%% Load data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input/'))
load yearly_salinity.mat

%% Best fit graph (95% confidence level)
figure(1), close all;
L = 6;
% Create a datetime vector for each month from 1992 to 2017
startDate = datetime(1992, 1, 1);  % Start in January 1992
endDate = datetime(2017, 12, 1);   % End in December 2017
months = dateshift(startDate, 'start', 'month', 0:311);

% Calculate the mean temperature up to layer L for each variable
FramStrait_m = mean(FramStrait(1:L,:),1);
FramStrait_std = std(FramStrait(1:L,:),1);
BarentsSea_m = mean(BarentsSea(1:L,:),1);
BarentsSea_std = std(BarentsSea(1:L,:),1);
BeaufortGyre_m = mean(BeaufortGyre(1:L,:),1);
BeaufortGyre_std = std(BeaufortGyre(1:L,:),1);

% Define x-axis and the upper and lower bounds for the patch
FramStrait_up = FramStrait_m + FramStrait_std;
FramStrait_low = FramStrait_m - FramStrait_std;
BarentsSea_up = BarentsSea_m + BarentsSea_std;
BarentsSea_low = BarentsSea_m - BarentsSea_std;
BeaufortGyre_up = BeaufortGyre_m + BeaufortGyre_std;
BeaufortGyre_low = BeaufortGyre_m - BeaufortGyre_std;


% Plotting for each variable
t = tiledlayout(3,1);
t.TileSpacing = 'tight';
t.Padding = 'compact';

nexttile
CalculateTREND(months, FramStrait_m, 'b'); hold on;
fill([months fliplr(months)], [FramStrait_up fliplr(FramStrait_low)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
ylabel('Salinity (pss)'); ylim([33.5 35])
title('Annual trend (95% Confidence)');
subtitle(sprintf('Mean salinity upto %.1f m deep',Z(L))); 
box on; set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold off
set(gca,'XTickLabel',[]); 
legend('Fram Strait','location','southwest','fontsize',20);


nexttile
CalculateTREND(months, BarentsSea_m,'k'); hold on;
fill([months fliplr(months)], [BarentsSea_up fliplr(BarentsSea_low)], ...
    'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
ylabel('Salinity (pss)'); ylim([32.5 34.5])
box on; set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold off
set(gca,'XTickLabel',[]);
legend('Barents Sea','location','southwest','fontsize',20);


nexttile
CalculateTREND(months, BeaufortGyre_m, 'r'); hold on;
fill([months fliplr(months)], [BeaufortGyre_up fliplr(BeaufortGyre_low)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
ylabel('Salinity (pss)'); ylim([29.5 32.5])
box on; set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold off
ax = gca; ax.XTick = months(1:12:end); ax.XTickLabelRotation = 45;
legend('Beaufort Sea','location','southwest','fontsize',20);


% Function to plot data, best-fit line, and confidence interval
function CalculateTREND(months, data, color)
    x = (1:length(data))';
    [fitResult, ~] = polyfit(x, data, 1);
    yFit = polyval(fitResult, x);
    stdResiduals = std(data - yFit);
    ci = 1.96 * stdResiduals;
    lowerBound = yFit - ci;
    upperBound = yFit + ci;

    plot(datetime(months), data, 'Color', color, 'LineStyle', '-', 'LineWidth', 3, 'MarkerSize', 8); hold on;
    plot(datetime(months), yFit, [color '-'], 'LineWidth', 4); hold on;
    plot(datetime(months), lowerBound, [color '--'], 'LineWidth', 2); hold on;
    plot(datetime(months), upperBound, [color '--'], 'LineWidth', 2); hold on;
    set(gca, 'FontSize', 20, 'FontWeight', 'bold');
    set(gcf, 'Position', [-4677,198,924,1072]); 
end


