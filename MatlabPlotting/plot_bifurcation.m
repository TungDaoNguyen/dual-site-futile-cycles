% Define network name pattern
networkPattern = 'eec1c2_4pss';

% Options
startYatZero = true;
manualTsRange = false;
manualTsmin = 18.5;
manualTsmax = 20;
plotLineWidth = 3;

colorPositive = [0.00 0.45 0.74];
colorLiving   = [0.93 0.69 0.13];
colorDead     = [0.85 0.10 0.10];

dataDir = '../BifDataForPlotting/';
files = dir([dataDir networkPattern '*_bifdata.csv']);
if isempty(files), error('No files found.'); end

disp('Found files:');
for i = 1:length(files), disp(['  ' files(i).name]); end

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

    allData = [allData; currentData(:, requiredCols)]; %#ok<AGROW>
end

disp(['Total rows: ' num2str(height(allData))]);

stableLiving     = allData(allData.stability_flag==1 & strcmp(allData.type,'living_boundary'),:);
unstableLiving   = allData(allData.stability_flag==0 & strcmp(allData.type,'living_boundary'),:);
stablePositive   = allData(allData.stability_flag==1 & strcmp(allData.type,'positive'),:);
unstablePositive = allData(allData.stability_flag==0 & strcmp(allData.type,'positive'),:);
stableDead       = allData(allData.stability_flag==1 & strcmp(allData.type,'dead_boundary'),:);
unstableDead     = allData(allData.stability_flag==0 & strcmp(allData.type,'dead_boundary'),:);

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

if manualTsRange
    Tsmin = manualTsmin;
    Tsmax = manualTsmax;
else
    Tsmin = min(allData.Ts);
    Tsmax = max(allData.Ts);
end

variables = {'s0','s1','s2'};
labels    = {'$s_0$','$s_1$','$s_2$'};

% Plot order: dashed first, solid on top
datasetData   = {unstableDead, unstableLiving, unstablePositive, ...
                 stableDead,   stableLiving,   stablePositive};

datasetColors = {colorDead,    colorLiving,    colorPositive, ...
                 colorDead,    colorLiving,    colorPositive};

datasetStyles = {'--',         '--',           '--', ...
                 '-',          '-',            '-'};

datasetLabels = {'Unstable dead boundary',  'Unstable living boundary', 'Unstable positive', ...
                 'Stable dead boundary',    'Stable living boundary',   'Stable positive'};

% Variables used to identify branches. This is the important part: branch
% tracing is done in the full state space, not from a single plotted y-value.
branchVars = {'u','c1','c2','c3','c4','e','s0','s1','s2'};

% Branch colors: each traced branch/component gets its own color.
% Stability is encoded only by line style: dashed = unstable, solid = stable.
% Cool palette only: blue/teal/green hues, avoiding orange/red because those
% are reserved for boundary steady states in related figures.
branchPalette = [
    0.00 0.45 0.74;   % branch 1: MATLAB blue
    0.00 0.40 0.00;   % branch 2: dark green (clearly distinguishable from blue)
    0.00 0.15 0.85;   % branch 3: vivid deep blue
    0.13 0.55 0.13;   % branch 4: forest green
    0.00 0.28 0.55;   % branch 5: dark slate blue
    0.10 0.65 0.40;   % branch 6: medium green
    0.05 0.05 0.45;   % branch 7: indigo blue
    0.00 0.25 0.08;   % branch 8: very dark green
    0.18 0.50 0.95;   % branch 9: light saturated blue
    0.18 0.72 0.18;   % branch 10: bright green
    0.00 0.30 0.12;   % branch 11: dark forest green
    0.20 0.60 0.90;   % branch 12: sky blue
];

% -----------------------------------------------------------------------
% Create figure using tiledlayout
% -----------------------------------------------------------------------
fig = figure('Position',[100,100,1400,520]);

t = tiledlayout(fig, 1, 3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

axArr   = gobjects(1,3);

% Precompute traced branches for each dataset so all panels use the same
% branch decomposition.
datasetBranches = cell(size(datasetData));
for k = 1:length(datasetData)
    d = datasetData{k};
    if isempty(d)
        datasetBranches{k} = {};
    else
        datasetBranches{k} = traceBranchesByNearestState(d, branchVars);
    end
end

for i = 1:3
    ax = nexttile(t, i);
    axArr(i) = ax;
    hold(ax, 'on');

    for k = 1:length(datasetData)
        d = datasetData{k};
        branches = datasetBranches{k};

        if isempty(d) || isempty(branches)
            continue;
        end

        for b = 1:numel(branches)
            idx = branches{b};
            idx = idx(~isnan(d.(variables{i})(idx)));

            if numel(idx) < 2
                continue;
            end

            x = d.Ts(idx);
            y = d.(variables{i})(idx);

            branchColor = branchPalette(mod(b-1, size(branchPalette,1)) + 1, :);

            plot(ax, x, y, datasetStyles{k}, ...
                'Color', branchColor, ...
                'LineWidth', plotLineWidth);
        end
    end

    hold(ax, 'off');

    xlabel(ax, '$T_s$', ...
        'Interpreter','latex', ...
        'FontSize',14, ...
        'FontWeight','bold');

    ylabel(ax, labels{i}, ...
        'Interpreter','latex', ...
        'FontSize',14, ...
        'FontWeight','bold');

    title(ax, labels{i}, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'FontWeight','bold');

    xlim(ax, [Tsmin Tsmax]);

    allVals = [];
    cats = {stableLiving,unstableLiving,stablePositive,unstablePositive,...
            stableDead,unstableDead};

    for c = 1:length(cats)
        d = cats{c};
        d = d(d.Ts>=Tsmin & d.Ts<=Tsmax & ~isnan(d.(variables{i})), :);
        allVals = [allVals; d.(variables{i})]; %#ok<AGROW>
    end

    if ~isempty(allVals)
        ylo = min(allVals);
        yhi = max(allVals);
    else
        ylo = 0;
        yhi = 1;
    end

    if yhi > ylo
        if startYatZero && ylo >= 0
            ylim(ax, [0, yhi*1.15]);
        else
            ylim(ax, [ylo - 0.05*(yhi-ylo), yhi*1.15]);
        end
    else
        ylim(ax, [max(0, ylo*0.85), ylo*1.15]);
    end

    grid(ax, 'on');

    set(ax, ...
        'FontSize',13, ...
        'LineWidth',1.5, ...
        'TickLabelInterpreter','latex');
end

% -----------------------------------------------------------------------
% Minimal legend inside middle tile: stability only
% -----------------------------------------------------------------------
axMid = axArr(2);
hold(axMid, 'on');

hStable = plot(axMid, nan, nan, '-', ...
    'Color', 'k', ...
    'LineWidth', plotLineWidth);

hUnstable = plot(axMid, nan, nan, '--', ...
    'Color', 'k', ...
    'LineWidth', plotLineWidth);

lgd = legend(axMid, [hStable, hUnstable], ...
    {'Stable','Unstable'}, ...
    'Orientation','vertical', ...
    'Interpreter','latex', ...
    'FontSize',11, ...
    'Location','north');

lgd.Box = 'on';

% -----------------------------------------------------------------------
% Save
% -----------------------------------------------------------------------
set(fig,'Units','Inches');
pos = get(fig,'Position');

set(fig, ...
    'PaperPositionMode','Auto', ...
    'PaperUnits','Inches', ...
    'PaperSize',[pos(3),pos(4)]);

savePath = '/Users/bjoshi/Library/CloudStorage/Dropbox/BifuncMSS/MatlabPlotting/';

print(fig, [savePath networkPattern '_substrates.png'], '-dpng', '-r300');
print(fig, [savePath networkPattern '_substrates.pdf'], '-dpdf', '-r300', '-painters');

disp(['Saved: ' savePath]);
disp(['Network: ' networkPattern]);
disp(['Ts range: [' num2str(Tsmin) ', ' num2str(Tsmax) ']']);

% -----------------------------------------------------------------------
% Local helper functions
% -----------------------------------------------------------------------
function branches = traceBranchesByNearestState(d, branchVars)
    % Build continuous branches by matching points at neighboring Ts values
    % using nearest neighbor distance in normalized state space.

    n = height(d);
    branches = {};
    if n == 0
        return;
    end

    % Keep only variables that exist and are not entirely NaN.
    keep = false(size(branchVars));
    for j = 1:numel(branchVars)
        keep(j) = ismember(branchVars{j}, d.Properties.VariableNames) && any(~isnan(d.(branchVars{j})));
    end
    branchVars = branchVars(keep);

    if isempty(branchVars)
        branches = {transpose(1:n)};
        return;
    end

    X = zeros(n, numel(branchVars));
    for j = 1:numel(branchVars)
        v = d.(branchVars{j});
        vmin = min(v, [], 'omitnan');
        vmax = max(v, [], 'omitnan');
        vrng = vmax - vmin;
        if isnan(vrng) || vrng == 0
            vrng = 1;
        end
        X(:,j) = (v - vmin) / vrng;
    end

    % Replace remaining NaNs after normalization with zero so distances do
    % not become NaN. These should be rare because missing columns were
    % created as NaN, but branchVars filters all-NaN columns above.
    X(isnan(X)) = 0;

    TsVals = unique(d.Ts);
    TsVals = sort(TsVals);

    active = struct('branchIndex', {}, 'lastRow', {}, 'lastTs', {});
    maxStepNorm = 0.30;
    maxTsGapFactor = 3.5;

    if numel(TsVals) > 1
        dTs = diff(TsVals);
        typicalTsGap = median(dTs(dTs > 0));
        if isempty(typicalTsGap) || isnan(typicalTsGap) || typicalTsGap == 0
            typicalTsGap = inf;
        end
    else
        typicalTsGap = inf;
    end

    for tIdx = 1:numel(TsVals)
        ts = TsVals(tIdx);
        currentRows = find(d.Ts == ts);

        % Drop stale active branches if there is a real Ts gap.
        if isfinite(typicalTsGap)
            fresh = true(size(active));
            for a = 1:numel(active)
                fresh(a) = (ts - active(a).lastTs) <= maxTsGapFactor * typicalTsGap;
            end
            active = active(fresh);
        end

        if isempty(active)
            for r = 1:numel(currentRows)
                branches{end+1} = currentRows(r); %#ok<AGROW>
                active(end+1).branchIndex = numel(branches); %#ok<AGROW>
                active(end).lastRow = currentRows(r);
                active(end).lastTs = ts;
            end
            continue;
        end

        activeRows = [active.lastRow];
        D = zeros(numel(activeRows), numel(currentRows));
        for a = 1:numel(activeRows)
            for r = 1:numel(currentRows)
                D(a,r) = norm(X(activeRows(a),:) - X(currentRows(r),:));
            end
        end

        assignedActive = false(1, numel(activeRows));
        assignedCurrent = false(1, numel(currentRows));

        % Greedy global nearest-neighbor assignment. Counts are small, so a
        % simple loop is adequate and avoids optimization toolbox functions.
        while true
            Dwork = D;
            Dwork(assignedActive,:) = inf;
            Dwork(:,assignedCurrent) = inf;
            [bestVal, linearIdx] = min(Dwork(:));

            if isinf(bestVal) || bestVal > maxStepNorm
                break;
            end

            [aBest, rBest] = ind2sub(size(Dwork), linearIdx);
            bIdx = active(aBest).branchIndex;
            row = currentRows(rBest);

            branches{bIdx}(end+1) = row; %#ok<AGROW>
            active(aBest).lastRow = row;
            active(aBest).lastTs = ts;

            assignedActive(aBest) = true;
            assignedCurrent(rBest) = true;
        end

        % Any unassigned current point starts a new branch.
        for r = 1:numel(currentRows)
            if ~assignedCurrent(r)
                branches{end+1} = currentRows(r); %#ok<AGROW>
                active(end+1).branchIndex = numel(branches); %#ok<AGROW>
                active(end).lastRow = currentRows(r);
                active(end).lastTs = ts;
            end
        end
    end

    % Remove one-point branches because they cannot form a line.
    keepBranch = false(size(branches));
    for b = 1:numel(branches)
        keepBranch(b) = numel(branches{b}) >= 2;
    end
    branches = branches(keepBranch);
end
