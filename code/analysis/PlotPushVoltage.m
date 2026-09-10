function PlotPushVoltage()
%PLOTPUSHVOLTAGE  Hydrophone peak-to-peak & peak-negative voltage vs TX voltage
%   for the S5-1 79-element push, with the saturation zone shaded.
%
%   Peak-to-peak = max(y)-min(y); peak-negative = -min(y) (baseline-removed).
%   The saturation onset is detected where the peak-negative voltage first
%   falls >5% below the low-voltage linear trend (the AH-2010 preamp clips at
%   ~2 V peak / ~4 V peak-to-peak). Data path is external (edit `d` if moved).

    addpath(fileparts(mfilename('fullpath')));                 % readHWS
    d = miData('S5-1','2026-08-13_sessionA_preamp_50ohm','push_79el_sweep');
    V = [2 3 4 5 6 7 8 9 10 11 12]';

    Vpp = nan(numel(V),1); Vneg = nan(numel(V),1); Vpos = nan(numel(V),1);
    for i = 1:numel(V)
        w = readHWS(fullfile(d, sprintf('S5-1_%dV_50ohm_fund_push.hws', V(i))));
        yb = w(1).y - median(w(1).y);
        Vneg(i) = -min(yb); Vpos(i) = max(yb); Vpp(i) = Vpos(i) + Vneg(i);
    end

    % --- saturation onset: deviation of Vneg from the low-V linear trend ---
    anc = V<=6;  pl = polyfit(V(anc), Vneg(anc), 1);           % anchor on 2-6 V
    dev = Vneg ./ polyval(pl, V);
    Vsat = V(end)+1;
    hi = find(V>6);
    for j = 1:numel(hi)
        k = hi(j);
        if dev(k) < 0.95 && (k==numel(V) || dev(k+1) < 0.95), Vsat = V(k); break; end
    end
    m  = V < Vsat;                                             % unsaturated
    ppf = polyfit(V(m), Vpp(m), 1); ngf = polyfit(V(m), Vneg(m), 1);

    % --- plot ---
    fig = figure('Color','w','Position',[100 100 760 500]); hold on; grid on; box on
    yl = [0 max(Vpp)*1.12]; xl = [0 13];
    patch([Vsat xl(2) xl(2) Vsat], [yl(1) yl(1) yl(2) yl(2)], [0.55 0.55 0.55], ...
        'FaceAlpha',0.15, 'EdgeColor','none', 'DisplayName','saturation zone');
    Vx = xl(1):0.25:xl(2);
    plot(Vx, polyval(ppf,Vx), '--', 'Color',[0 0.45 0.74],  'LineWidth',1.1, 'HandleVisibility','off');
    plot(Vx, polyval(ngf,Vx), '--', 'Color',[0.85 0.33 0.10],'LineWidth',1.1, 'HandleVisibility','off');
    plot(V, Vpp,  'o-', 'Color',[0 0.45 0.74],   'MarkerFaceColor',[0 0.45 0.74],   'MarkerSize',7, 'LineWidth',1.4, 'DisplayName','peak-to-peak');
    plot(V, Vneg, 's-', 'Color',[0.85 0.33 0.10],'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',7, 'LineWidth',1.4, 'DisplayName','peak-negative');
    xline(Vsat, ':', sprintf('saturation onset ~%g V', Vsat), 'Color',[0.4 0.4 0.4], ...
        'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    xlabel('TX voltage (V)'); ylabel('hydrophone voltage (V)');
    title('S5-1 push, 79 elements: hydrophone V_{pp} / V_{neg} vs TX voltage');
    legend('Location','northwest'); xlim(xl); ylim(yl);
    text(0.98,0.05, 'dashed = linear fit of unsaturated region', 'Units','normalized', ...
        'HorizontalAlignment','right', 'FontSize',9, 'Color',[0.35 0.35 0.35]);

    outdir = miRoot('results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    out = fullfile(outdir,'PushVoltage_79el_Vpp_Vneg.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('saturation onset ~%g V | Vpp(max)=%.2f V, Vneg(max)=%.2f V\n', Vsat, max(Vpp), max(Vneg));
    fprintf('saved: %s\n', out);
end
