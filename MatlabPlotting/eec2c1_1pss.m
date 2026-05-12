% Define network name pattern
networkPattern = 'eec2c1_1pss';

% Options
startYatZero = true;
manualTsRange = false;   % use full data range
manualTsmin = 18.5;
manualTsmax = 20;
plotLineWidth = 3;

% -----------------------------------------------------------------------
% Color and style convention:
%   Positive steady states:  solid (stable) / dashed (unstable), BLUE
%   Living boundary:         solid (stable) / dashed (unstable), ORANGE
%   Dead boundary:           solid (stable) / dashed (unstable), RED
%
% Plot order: dashed (unstable) first, solid (stable) on top,
% so solid always wins visually when lines overlap.
% -----------------------------------------------------------------------
colorPositive = [0.00 0.45 0.74];
colorLiving   = [0.93 0.69 0.13];
colorDead     = [0.85 0.10 0.10];

% Find all matching files
dataDir = '../BifDataForPlotting/';
files = dir([dataDir networkPattern '*_bifdata.csv']);
if isempty(files), error('No files found.'); end

disp('Found files:');
for i = 1:length(files), disp(['  ' files(i).name]); end

% Read and combine data
allData = [];
for i = 1:length(files)
    currentData = readtable(fullfile(files(i).folder, files(i).name));
    if contains(files(i).name, 'living_boundary')
        dataType = 'living_boundary';
    elseif contains(files(i).name, 'dead_boundary')
        dataType = 'dead_boundary';
    else
        dataType = 'positive';
    end
    if ~ismember('u', currentData.Properties.VariableNames)
        currentData.u = nan(height(currentData), 1);
    end
    currentData.type = repmat({dataType}, height(currentData), 1);
    requiredCols = {'Ts','u','c1','c2','c3','c4','e','s0','s1','s2',...
                    'stable','stability_flag','type'};
    for col = requiredCols(1:end-1)
        if ~ismember(col{1}, currentData.Properties.VariableNames)
            currentData.(col{1}) = nan(height(currentData), 1);
        end
    end
    allData = [allData; currentData(:, requiredCols)];
end
disp(['Total rows: ' num2str(height(allData))]);

% Separate categories
stableLiving     = allData(allData.stability_flag==1 & strcmp(allData.type,'living_boundary'),:);
unstableLiving   = allData(allData.stability_flag==0 & strcmp(allData.type,'living_boundary'),:);
stablePositive   = allData(allData.stability_flag==1 & strcmp(allData.type,'positive'),:);
unstablePositive = allData(allData.stability_flag==0 & strcmp(allData.type,'positive'),:);
stableDead       = allData(allData.stability_flag==1 & strcmp(allData.type,'dead_boundary'),:);
unstableDead     = allData(allData.stability_flag==0 & strcmp(allData.type,'dead_boundary'),:);

% LaTeX defaults
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% Ts range
if manualTsRange
    Tsmin = manualTsmin;
    Tsmax = manualTsmax;
else
    Tsmin = min(allData.Ts);
    Tsmax = max(allData.Ts);
end
fprintf('Ts range: [%.4f, %.4f]\n', Tsmin, Tsmax);
fprintf('>> If plot looks too small, set manualTsRange=true and use the values above.\n');

% Substrate variables only
variables = {'s0','s1','s2'};
labels    = {'$s_0$','$s_1$','$s_2$'};

% Datasets in plot order: DASHED first, SOLID on top
datasetData   = {unstableDead, unstableLiving, unstablePositive, ...
                 stableDead,   stableLiving,   stablePositive};
datasetColors = {colorDead,    colorLiving,    colorPositive, ...
                 colorDead,    colorLiving,    colorPositive};
datasetStyles = {'--',         '--',           '--', ...
                 '-',          '-',            '-'};
datasetLabels = {'Unstable dead boundary',  'Unstable living boundary', 'Unstable positive', ...
                 'Stable dead boundary',    'Stable living boundary',   'Stable positive'};

% -----------------------------------------------------------------------
% Create figure: 1 x 3
% -----------------------------------------------------------------------
figure('Position',[100,100,1400,460]);

t = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

allH    = [];
allLbls = {};
seen    = {};

for i = 1:3
    nexttile(t, i);
    hold on;

    for k = 1:length(datasetData)
        d = datasetData{k};
        d = d(~isnan(d.(variables{i})), :);
        if isempty(d), continue; end
        d = sortrows(d, 'Ts');

        h = plot(d.Ts, d.(variables{i}), datasetStyles{k}, ...
            'Color', datasetColors{k}, 'LineWidth', plotLineWidth);

        lbl = datasetLabels{k};
        if ~ismember(lbl, seen)
            seen{end+1} = lbl;
            allH(end+1)    = h;
            allLbls{end+1} = lbl;
        end
    end

    hold off;

    xlabel('$T_s$','Interpreter','latex','FontSize',14,'FontWeight','bold');
    ylabel(labels{i},'Interpreter','latex','FontSize',14,'FontWeight','bold');
    title(labels{i},'Interpreter','latex','FontSize',16,'FontWeight','bold');
    xlim([Tsmin Tsmax]);

    % Y limits from all data in Ts range
    allVals = [];
    cats = {stableLiving,unstableLiving,stablePositive,unstablePositive,...
            stableDead,unstableDead};
    for c = 1:length(cats)
        d = cats{c};
        d = d(d.Ts>=Tsmin & d.Ts<=Tsmax & ~isnan(d.(variables{i})), :);
        allVals = [allVals; d.(variables{i})(:)];
    end

    if ~isempty(allVals)
        ylo = min(allVals);
        yhi = max(allVals);
    else
        ylo = 0; yhi = 1;
    end

    if yhi > ylo
        % For s1 (i==2): use slightly negative lower bound so living
        % boundary steady state at s1=0 is clearly visible
        if i == 2 && ylo >= 0
            ylim([-0.05*yhi, yhi*1.15]);
        elseif startYatZero && ylo >= 0
            ylim([0, yhi*1.15]);
        else
            ylim([ylo - 0.05*(yhi-ylo), yhi*1.15]);
        end
    else
        ylim([max(0, ylo*0.85), ylo*1.15]);
    end

    grid on;
    set(gca,'FontSize',13,'LineWidth',1.5,'TickLabelInterpreter','latex');
end

% Single shared legend placed inside the s2 tile (top left)
if ~isempty(allH)
    nexttile(t, 3);
    lgd = legend(allH, allLbls, ...
        'Orientation','vertical', ...
        'Interpreter','latex', ...
        'FontSize',11, ...
        'Location','northwest');
    lgd.Box = 'on';
end

% Save
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3),pos(4)]);

savePath = 'MatlabPlotting';
print([savePath networkPattern '_substrates.png'],'-dpng','-r300');
print([savePath networkPattern '_substrates.pdf'],'-dpdf','-r300','-painters');

disp(['Saved: ' savePath]);
disp(['Network: ' networkPattern]);
disp(['Ts range: [' num2str(Tsmin) ', ' num2str(Tsmax) ']']);
