% -----------------------------------------------------------------------
% Generate bifurcation data and plot for the (E,C1) single-site network
%
% ODE system (after eliminating e = Te-c1-c2, s0 = Ts-s1-c1-2*c2):
%   dc1/dt = k1p*s0*e - (k1m+k1)*c1 - k2p*s1*c1 + (k2m+k2)*c2
%   dc2/dt = k2p*s1*c1 - (k2m+k2)*c2
%   ds1/dt = k1*c1 - k2p*s1*c1 + k2m*c2
%
% ACR value: s1* = k1/k2*  where k2* = k2*k2p/(k2+k2m)
% Transcritical bifurcation at Ts = k1/k2*
% -----------------------------------------------------------------------

% --- Rate constants ---
k1p = 2.0;    % k1+
k1m = 1.0;    % k1-
k1  = 1.0;    % k1 (catalytic)
k2p = 3.0;    % k2+
k2m = 0.5;    % k2-
k2  = 1.5;    % k2 (catalytic)
Te  = 1.0;    % total enzyme

% Derived quantities
k2star = k2*k2p / (k2 + k2m);    % k2*
acrValue = k1 / k2star;           % ACR value = k1/k2*
fprintf('k2* = %.4f\n', k2star);
fprintf('ACR value (Ts threshold) = k1/k2* = %.4f\n', acrValue);

% --- Ts range: show below and above transcritical bifurcation ---
TsMin  = 0;
TsMax  = acrValue * 4;
TsStep = acrValue / 100;
TsVals = TsMin : TsStep : TsMax;

% --- Options ---
plotLineWidth = 3;
colorPositive = [0.00 0.45 0.74];   % blue
colorDead     = [0.85 0.10 0.10];   % red

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% --- Storage ---
% Dead boundary SS: e=Te, s1=Ts, c1=c2=0, s0=0
% Positive SS: solve 3-variable system numerically

deadTs  = []; deadS0 = []; deadS1 = []; deadStab = [];
posTs   = []; posS0  = []; posS1  = []; posStab  = [];

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

for Ts = TsVals

    % ---- Dead boundary steady state ----
    % s1 = Ts, c1 = c2 = 0, s0 = 0, e = Te
    deadTs(end+1)  = Ts;
    deadS0(end+1)  = 0;
    deadS1(end+1)  = Ts;
    % Stability: linearize reduced 3D system (c1,c2,s1) at BSS
    % At BSS: c1=0, c2=0, s1=Ts, s0=0, e=Te
    J_bss = jacobianEC1(0, 0, Ts, Ts, Te, k1p,k1m,k1,k2p,k2m,k2);
    ev_bss = eig(J_bss);
    % BSS is stable if all eigenvalues have negative real part
    if all(real(ev_bss) < 0)
        deadStab(end+1) = 1;   % stable
    else
        deadStab(end+1) = 0;   % unstable
    end

    % ---- Positive steady state ----
    if Ts > acrValue
        % Solve for steady state: F(c1,c2,s1) = 0
        % Use parameterization: at PSS, s1 = k1/k2* (ACR)
        s1ss = acrValue;
        % From dc2/dt=0: c2 = k2p*s1*c1/(k2m+k2)
        % From dc1/dt=0 + dc2/dt=0: k1p*s0*e = (k1m+k1)*c1 + k2p*s1*c1 - (k2m+k2)*c2
        %   => k1p*s0*e = (k1m+k1)*c1
        % Conservation: e = Te-c1-c2, s0 = Ts-s1ss-c1-2*c2
        % Solve for c1 numerically
        try
            c1ss = fzero(@(c1) steadyStateResidual(c1, s1ss, Ts, Te, k1p,k1m,k1,k2p,k2m,k2), ...
                         [1e-10, min(Te, Ts)/2]);
            c2ss = k2p*s1ss*c1ss / (k2m+k2);
            s0ss = Ts - s1ss - c1ss - 2*c2ss;
            ess  = Te - c1ss - c2ss;

            if s0ss > 0 && ess > 0 && c1ss > 0 && c2ss > 0
                posTs(end+1)  = Ts;
                posS0(end+1)  = s0ss;
                posS1(end+1)  = s1ss;

                % Stability via Jacobian eigenvalues (3x3 reduced system)
                J = jacobianEC1(c1ss, c2ss, s1ss, Ts, Te, k1p,k1m,k1,k2p,k2m,k2);
                ev = eig(J);
                if all(real(ev) < 0)
                    posStab(end+1) = 1;
                else
                    posStab(end+1) = 0;
                end
            end
        catch
        end
    end
end

fprintf('Dead BSS points: %d\n', length(deadTs));
fprintf('Positive SS points: %d\n', length(posTs));

% --- Separate stable/unstable ---
stableDeadTs  = deadTs(deadStab==1);  stableDeadS0  = deadS0(deadStab==1);  stableDeadS1  = deadS1(deadStab==1);
unstableDeadTs= deadTs(deadStab==0);  unstableDeadS0= deadS0(deadStab==0);  unstableDeadS1= deadS1(deadStab==0);
stablePosTs   = posTs(posStab==1);    stablePosS0   = posS0(posStab==1);    stablePosS1   = posS1(posStab==1);
unstablePosTs = posTs(posStab==0);    unstablePosS0 = posS0(posStab==0);    unstablePosS1 = posS1(posStab==0);

% --- Plot: 1x2 ---
figure('Position',[100,100,950,460]);

vars   = {'s0','s1'};
labels = {'$s_0$','$s_1$'};

deadData = {stableDeadS0, unstableDeadS0; stableDeadS1, unstableDeadS1};
posData  = {stablePosS0,  unstablePosS0;  stablePosS1,  unstablePosS1};
deadTsData = {stableDeadTs, unstableDeadTs};
posTsData  = {stablePosTs,  unstablePosTs};

allH = []; allLbls = {}; seen = {};

for i = 1:2
    subplot(1,2,i);
    hold on;

    % Dashed (unstable) first, solid (stable) on top
    % Unstable dead boundary
    if ~isempty(unstableDeadTs)
        h = plot(unstableDeadTs, deadData{i,2}, '--', 'Color', colorDead, 'LineWidth', plotLineWidth);
        lbl = 'Unstable boundary';
        if ~ismember(lbl,seen), seen{end+1}=lbl; allH(end+1)=h; allLbls{end+1}=lbl; end
    end
    % Unstable positive
    if ~isempty(unstablePosTs)
        h = plot(unstablePosTs, posData{i,2}, '--', 'Color', colorPositive, 'LineWidth', plotLineWidth);
        lbl = 'Unstable positive';
        if ~ismember(lbl,seen), seen{end+1}=lbl; allH(end+1)=h; allLbls{end+1}=lbl; end
    end
    % Stable dead boundary
    if ~isempty(stableDeadTs)
        h = plot(stableDeadTs, deadData{i,1}, '-', 'Color', colorDead, 'LineWidth', plotLineWidth);
        lbl = 'Stable boundary';
        if ~ismember(lbl,seen), seen{end+1}=lbl; allH(end+1)=h; allLbls{end+1}=lbl; end
    end
    % Stable positive
    if ~isempty(stablePosTs)
        h = plot(stablePosTs, posData{i,1}, '-', 'Color', colorPositive, 'LineWidth', plotLineWidth);
        lbl = 'Stable positive';
        if ~ismember(lbl,seen), seen{end+1}=lbl; allH(end+1)=h; allLbls{end+1}=lbl; end
    end

    hold off;
    xlabel('$T_s$','Interpreter','latex','FontSize',14,'FontWeight','bold');
    ylabel(labels{i},'Interpreter','latex','FontSize',14,'FontWeight','bold');
    title(labels{i},'Interpreter','latex','FontSize',16,'FontWeight','bold');
    xlim([TsMin TsMax]);

    % For s0, set slightly negative lower bound so boundary SS (s0=0) is visible
    ydata = [deadData{i,1}(:); deadData{i,2}(:); posData{i,1}(:); posData{i,2}(:)];
    ydata = ydata(~isnan(ydata));
    if ~isempty(ydata)
        yhi = max(ydata) * 1.15;
        if i == 1  % s0 panel
            ylo = -0.05 * max(ydata);
        else
            ylo = 0;
        end
        ylim([ylo, yhi]);
    end

    grid on;
    set(gca,'FontSize',13,'LineWidth',1.5,'TickLabelInterpreter','latex');
end

% Legend inside s1 panel
subplot(1,2,2);
lgd = legend(allH, allLbls, 'Orientation','vertical','Interpreter','latex',...
    'FontSize',11,'Location','northwest');
lgd.Box = 'on';

% Save
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);

savePath = '/Users/bjoshi/Library/CloudStorage/Dropbox/BifuncMSS/MatlabPlotting/';
print([savePath 'ec1_substrates.png'],'-dpng','-r300');
print([savePath 'ec1_substrates.pdf'],'-dpdf','-r300','-painters');

fprintf('Saved to %s\n', savePath);
fprintf('k1/k2* = %.4f  (transcritical bifurcation)\n', acrValue);
fprintf('Ts range: [%.2f, %.2f]\n', TsMin, TsMax);

% -----------------------------------------------------------------------
% Local helper: residual for finding c1 at positive SS
% -----------------------------------------------------------------------
function res = steadyStateResidual(c1, s1, Ts, Te, k1p,k1m,k1,k2p,k2m,k2)
    c2  = k2p*s1*c1 / (k2m+k2);
    s0  = Ts - s1 - c1 - 2*c2;
    e   = Te - c1 - c2;
    % From dc1/dt = 0: k1p*s0*e - (k1m+k1)*c1 - k2p*s1*c1 + (k2m+k2)*c2 = 0
    % But we also used dc2/dt=0, so just check dc1/dt = 0 gives:
    % k1p*s0*e = (k1m+k1)*c1  (after substituting dc2/dt=0)
    res = k1p*s0*e - (k1m+k1)*c1;
end

% -----------------------------------------------------------------------
% Local helper: 3x3 Jacobian of reduced system (c1,c2,s1)
% -----------------------------------------------------------------------
function J = jacobianEC1(c1, c2, s1, Ts, Te, k1p,k1m,k1,k2p,k2m,k2)
    s0 = Ts - s1 - c1 - 2*c2;
    e  = Te - c1 - c2;
    % dc1/dt = k1p*s0*e - (k1m+k1)*c1 - k2p*s1*c1 + (k2m+k2)*c2
    % dc2/dt = k2p*s1*c1 - (k2m+k2)*c2
    % ds1/dt = k1*c1 - k2p*s1*c1 + k2m*c2
    % Partial derivatives w.r.t. (c1, c2, s1)
    % ds0/dc1 = -1, ds0/dc2 = -2, ds0/ds1 = -1
    % de/dc1  = -1, de/dc2  = -1, de/ds1  =  0
    J = zeros(3,3);
    % d(dc1/dt)/dc1
    J(1,1) = k1p*((-1)*e + s0*(-1)) - (k1m+k1) - k2p*s1;
    % d(dc1/dt)/dc2
    J(1,2) = k1p*((-2)*e + s0*(-1)) + (k2m+k2);
    % d(dc1/dt)/ds1
    J(1,3) = k1p*((-1)*e) - k2p*c1;
    % d(dc2/dt)/dc1
    J(2,1) = k2p*s1;
    % d(dc2/dt)/dc2
    J(2,2) = -(k2m+k2);
    % d(dc2/dt)/ds1
    J(2,3) = k2p*c1;
    % d(ds1/dt)/dc1
    J(3,1) = k1 - k2p*s1;
    % d(ds1/dt)/dc2
    J(3,2) = k2m;
    % d(ds1/dt)/ds1
    J(3,3) = -k2p*c1;
end