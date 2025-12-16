function porkchop_plot(MISSION)
% PORKCHOPPLOT: Delta-V contour plot vs departure and arrival dates (MJD2000)
%
% INPUT: MISSION.Solution.t1_v  : vector MJD2000 of departures (N elements)
%        MISSION.Solution.t2_v  : vector MJD2000 of arrivals   (M elements)
%        MISSION.Solution.DELTAV_TOT : matrice N x M con deltaVtot [km/s]
%
% OUTPUT : Contour plot : X-axis = departure date
%                         Y-axis = arrival date ; 
%                         Colorbar = Δv [km/s]

% Checking correctness 
if ~isfield(MISSION,'Solution') || ...
                ~all(isfield(MISSION.Solution, {'t1_v','t2_v','DELTAV_TOT'}))
    error('MISSION.Solution does not have t1_v, t2_v and DELTAV_TOT.');
end

t1_v = MISSION.Solution.t1_v(:);         % N x 1
t2_v = MISSION.Solution.t2_v(:);         % M x 1
Z = MISSION.Solution.DELTAV_TOT;         % N x M (dep x arr)

% Check dimensions
if size(Z,1) ~= numel(t1_v) || size(Z,2) ~= numel(t2_v)
    error('Dimensions of DELTAV_TOT non coherent with t1_v (rows) e t2_v (columns).');
end

% Convert from MJD2000 -> datetime
dep_dt = datetime(2000,1,1,12,0,0) + days(t1_v);    % N x 1
arr_dt = datetime(2000,1,1,12,0,0) + days(t2_v);    % M x 1

% Convert from datetime to serial number
dep_num = datenum(dep_dt);                   % N x 1
arr_num = datenum(arr_dt);                   % M x 1

% Mesh for X-axis and Y-axis
[X, Y] = meshgrid(dep_num(:).', arr_num(:).');

% Transpose in order to have X = departures, Y = arrival
Zplot = Z'; 

% Checking all NaN values
Zvalid = Zplot(~isnan(Zplot));
if isempty(Zvalid)
    error('DELTAV_TOT is NaN: no point is valid.');
end

levels = 1:1:30;

% Contour plot
figure
[C, h] = contour(X, Y, Zplot, levels, 'LineWidth', 2);
% clabel(C, h, 'Color', [0.1 0.1 0.1], 'FontSize', 14, 'LabelSpacing', 2000);
hold on;

% Colorbar
cb = colorbar;
title(cb, '\DeltaV [km/s]');

% Axis and title
xlabel('Departure date');
ylabel('Arrival date');
title('Porkchop plot');
grid on

% Dates on axis
ax = gca;
xticks_dt = dep_dt(1):calmonths(3):dep_dt(end);
yticks_dt = arr_dt(1):calmonths(3):arr_dt(end);
ax.XTick = datenum(xticks_dt);
ax.YTick = datenum(yticks_dt);
ax.XTickLabel = cellstr(datestr(xticks_dt,'yyyy mmm dd'));
ax.YTickLabel = cellstr(datestr(yticks_dt,'yyyy mmm dd'));
ax.XTickLabelRotation = 45;
ax.YTickLabelRotation = 45;

% MINIMUM DELTAV_TOT
[min_deltaV, ~] = min(Zvalid);
[row, col] = find(Zplot == min_deltaV, 1, 'first');
if ~isempty(row) && ~isempty(col)
    plot(X(row, col), Y(row, col), 'o', 'MarkerSize', 7, ...
         'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', 'w');
end

% % Constant TOF lines (in days)
% TOF_const = [60, 120, 180, 240, 300]; % in days 
% ymin = datenum(arr_dt(1));
% ymax = datenum(arr_dt(end));
% for tof = TOF_const
%     t2_line = dep_num + tof;                  % t2 = t1 + TOF
%     mask = t2_line >= ymin & t2_line <= ymax;
%     x_plot = dep_num(mask);
%     y_plot = t2_line(mask);
%     plot(x_plot, y_plot, 'LineWidth', 0.8, 'Color', [0.1 0.1 0.1]);
%     if ~isempty(x_plot)
%         mid_idx = round(numel(x_plot)/2);
%         text(x_plot(mid_idx), y_plot(mid_idx), [num2str(tof) 'D'], ...
%              'Color', [0.1 0.1 0.1], 'FontSize', 16);
%     end
% end

hold off
end