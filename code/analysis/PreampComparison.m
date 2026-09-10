function PreampComparison()
%PREAMPCOMPARISON  S5-1 push: with vs without pre-amplifier analysis + figures.
%
%   Answers:
%     1. Do MI with/without preamp agree, and does the load capacitance vary?
%     2. Axial sweep MI vs depth, with and without preamp.
%     3. Chosen location, TX voltage vs MI, with and without preamp.
%     4. Chosen location (no preamp), TX voltage vs MI, for 41/61/79 elements.
%
%   Method: depth = distance from Home (transducer centre); f_awf = measured
%   -6 dB fundamental; 0.3 dB/cm/MHz derating; Onda parallel-circuit
%   sensitivity (HydrophoneTables.mat).
%     - WITH preamp: 50 ohm (the notes say 1e6 but the ~2.2 V clip ceiling and
%       the ~73x ratio prove 50 ohm) -> gain + 6.1 pF divider, no /2.
%     - WITHOUT preamp: hydrophone straight into 1 MOhm scope -> no gain, and
%       the capacitive divider by the cable+scope load C_load (estimated from
%       the preamp/no-preamp ratio over the unsaturated 3-12 V range).
%
%   Data folder is external; edit `d` if it moves.  Run:  PreampComparison

    addpath(fileparts(mfilename('fullpath')));
    d = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp');
    Home = [55.8 105.84 -17.3];  rho=1000; c=1500;

    % --- calibration at f_awf ---
    fawf = measureFawf(fullfile(d,'S5-1_8V_Preamp_peak.hws'));
    [Moc,Ch,gain,Ca] = cal(fawf);                 % V/Pa, pF, -, pF
    sens_with = gain*Moc*Ch/(Ch+Ca)*1e6;          % V/MPa (50-ohm, preamp)

    % ================= Q1: C_load + MI agreement =================
    Vp = [2 3 4 5 6 7 8 9 10 12 15 20];
    Vn = [2 3 4 5 6 7 8 9 10 12 15 20 25 30 35 40 45 50];
    VnegP = arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',v))), Vp);
    VnegN = arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt.hws',v))), Vn);
    [tf,loc] = ismember(Vn,Vp); ratio = nan(size(Vn));
    ratio(tf) = VnegP(loc(tf)) ./ VnegN(tf);
    unsat = Vn>=3 & Vn<=12;                        % preamp unsaturated
    Cload = mean(ratio(unsat & tf),'omitnan')*(Ch+Ca)/gain - Ch;
    sens_no = Moc*Ch/(Ch+Cload)*1e6;              % V/MPa (no preamp)
    depth = 35;  der = db2mag(0.3*depth/10*fawf); % chosen-location depth
    fprintf('f_awf=%.2f MHz | sens_with=%.3f V/MPa | C_load=%.1f pF | sens_no=%.3f mV/MPa\n', ...
        fawf, sens_with, Cload, sens_no*1e3);
    fprintf('C_load: this (S5-1, %.0f pF) vs L11-5 (121 pF) -> %+.0f%%\n', Cload, 100*(Cload-121)/121);
    fprintf('preamp/no-preamp ratio 3-12V: mean %.1f, std %.1f (constant => same beam point)\n', ...
        mean(ratio(unsat&tf),'omitnan'), std(ratio(unsat&tf),'omitnan'));

    MIp = VnegP(:)/sens_with/der/sqrt(fawf);
    MIn = VnegN(:)/sens_no  /der/sqrt(fawf);

    % ---------- Figure 3 (Q3): MI vs TX voltage, preamp vs no preamp ----------
    Vsat = Vp(find(diff(VnegP)<=0,1));  if isempty(Vsat), Vsat=Vp(end); end
    f3=figure('Color','w','Position',[80 80 720 500]); hold on; grid on; box on
    yl=[0 max(MIn)*1.12];
    patch([Vsat 52 52 Vsat],[0 0 yl(2) yl(2)],[.55 .55 .55],'FaceAlpha',.13,'EdgeColor','none','DisplayName','preamp saturation');
    plot(Vn,MIn,'s-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74],'MarkerSize',6,'LineWidth',1.4,'DisplayName','no preamp (1 MOhm)');
    plot(Vp,MIp,'o-','Color',[.85 .33 .10],'MarkerFaceColor',[.85 .33 .10],'MarkerSize',6,'LineWidth',1.4,'DisplayName','preamp (50 Ohm)');
    xlabel('TX voltage (V)'); ylabel('MI (derated)'); xlim([0 52]); ylim(yl);
    title(sprintf('S5-1 push, 61 el, chosen location: MI vs TX voltage\npreamp vs no preamp (C_{load}=%.0f pF)',Cload));
    legend('Location','northwest');
    save_fig(f3,d,'MI_vs_TXvoltage_preamp_vs_nopreamp.png');

    % ================= Q2: axial MI vs depth =================
    [depP,MIaxP] = axial(d,'S5-1_4V_Preamp_%d.hws',1:10,sens_with,fawf);
    [depN,MIaxN] = axial(d,'S5-1_4V_NoPreamp_%d.hws',1:8,sens_no,fawf);
    f1=figure('Color','w','Position',[80 80 720 500]); hold on; grid on; box on
    plot(depP,MIaxP,'o-','Color',[.85 .33 .10],'MarkerFaceColor',[.85 .33 .10],'MarkerSize',6,'LineWidth',1.4,'DisplayName','preamp (50 Ohm)');
    plot(depN,MIaxN,'s-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74],'MarkerSize',6,'LineWidth',1.4,'DisplayName','no preamp (1 MOhm)');
    [~,ip]=max(MIaxP); plot(depP(ip),MIaxP(ip),'p','MarkerSize',14,'MarkerFaceColor',[.85 .33 .10],'MarkerEdgeColor','k','HandleVisibility','off');
    [~,in]=max(MIaxN); plot(depN(in),MIaxN(in),'p','MarkerSize',14,'MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor','k','HandleVisibility','off');
    xlabel('axial distance from transducer (mm)'); ylabel('MI (derated), 4 V');
    title('S5-1 push, 61 el: axial MI vs depth (local lateral+elev. max)');
    legend('Location','south');
    save_fig(f1,d,'MI_vs_axial_preamp_vs_nopreamp.png');

    % ================= Q4: element comparison (no preamp) =================
    Vel=[15 20 25 30 35 40 45 50];
    MI41=arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt_41el.hws',v))),Vel)/sens_no/der/sqrt(fawf);
    MI79=arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt_79el.hws',v))),Vel)/sens_no/der/sqrt(fawf);
    MI61=MIn;  % 61 el over full range
    f4=figure('Color','w','Position',[80 80 720 500]); hold on; grid on; box on
    plot(Vn, MI61,'s-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74],'MarkerSize',6,'LineWidth',1.4,'DisplayName','61 elements');
    plot(Vel,MI41,'^-','Color',[.47 .67 .19],'MarkerFaceColor',[.47 .67 .19],'MarkerSize',6,'LineWidth',1.4,'DisplayName','41 elements');
    plot(Vel,MI79,'d-','Color',[.85 .33 .10],'MarkerFaceColor',[.85 .33 .10],'MarkerSize',6,'LineWidth',1.4,'DisplayName','79 elements');
    yline(1.9,'k--','MI limit 1.9','HandleVisibility','off','LabelHorizontalAlignment','left');
    xlabel('TX voltage (V)'); ylabel('MI (derated)');
    title('S5-1 push, no preamp: MI vs TX voltage at the chosen location');
    legend('Location','northwest'); xlim([0 52]);
    save_fig(f4,d,'MI_vs_TXvoltage_41_61_79el.png');

    % ================= INTENSITY: Isppa.3 & Ispta.3 =================
    PRFeff = 24/(1.2+30);            % 20 Hz for 1.2 s, then >=30 s off -> 0.77 Hz
    Ilim = 190; Tlim = 720;         % FDA limits: W/cm2, mW/cm2
    cO=[.85 .33 .10]; cB=[0 .45 .74]; cG=[.47 .67 .19];

    % --- voltage sweep (preamp vs no preamp) ---
    [IsP,ItP]=arrayfun(@(v) intens(fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',v)),sens_with,depth,fawf,PRFeff),Vp);
    [IsN,ItN]=arrayfun(@(v) intens(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt.hws',v)),sens_no,depth,fawf,PRFeff),Vn);
    f5=figure('Color','w','Position',[60 60 1000 440]);
    subplot(1,2,1); hold on; grid on; box on
    patch([Vsat 52 52 Vsat],[0 0 Ilim*1.3 Ilim*1.3],[.55 .55 .55],'FaceAlpha',.13,'EdgeColor','none','DisplayName','preamp saturation');
    yline(Ilim,'k--','limit 190','HandleVisibility','off');
    plot(Vn,IsN,'s-','Color',cB,'MarkerFaceColor',cB,'DisplayName','no preamp'); plot(Vp,IsP,'o-','Color',cO,'MarkerFaceColor',cO,'DisplayName','preamp');
    xlabel('TX voltage (V)'); ylabel('I_{sppa.3} (W/cm^2)'); title('I_{sppa.3}'); legend('Location','northwest'); xlim([0 52]); ylim([0 Ilim*1.3]);
    subplot(1,2,2); hold on; grid on; box on
    patch([Vsat 52 52 Vsat],[0 0 max(ItN)*1.2 max(ItN)*1.2],[.55 .55 .55],'FaceAlpha',.13,'EdgeColor','none');
    plot(Vn,ItN,'s-','Color',cB,'MarkerFaceColor',cB); plot(Vp,ItP,'o-','Color',cO,'MarkerFaceColor',cO);
    xlabel('TX voltage (V)'); ylabel('I_{spta.3} (mW/cm^2)'); title(sprintf('I_{spta.3}  (burst; limit 720, far above)')); xlim([0 52]);
    sgtitle('S5-1 push, 61 el, chosen location: intensity vs TX voltage');
    save_fig(f5,d,'Intensity_vs_TXvoltage_preamp_vs_nopreamp.png');

    % --- axial (preamp vs no preamp) ---
    [daP,IsaxP,ItaxP]=axialInt(d,'S5-1_4V_Preamp_%d.hws',1:10,sens_with,fawf,PRFeff);
    [daN,IsaxN,ItaxN]=axialInt(d,'S5-1_4V_NoPreamp_%d.hws',1:8,sens_no,fawf,PRFeff);
    f6=figure('Color','w','Position',[60 60 1000 440]);
    subplot(1,2,1); hold on; grid on; box on
    plot(daP,IsaxP,'o-','Color',cO,'MarkerFaceColor',cO,'DisplayName','preamp'); plot(daN,IsaxN,'s-','Color',cB,'MarkerFaceColor',cB,'DisplayName','no preamp');
    xlabel('axial distance (mm)'); ylabel('I_{sppa.3} (W/cm^2)'); title('I_{sppa.3}, 4 V'); legend('Location','south');
    subplot(1,2,2); hold on; grid on; box on
    plot(daP,ItaxP,'o-','Color',cO,'MarkerFaceColor',cO); plot(daN,ItaxN,'s-','Color',cB,'MarkerFaceColor',cB);
    xlabel('axial distance (mm)'); ylabel('I_{spta.3} (mW/cm^2)'); title('I_{spta.3}, 4 V');
    sgtitle('S5-1 push, 61 el: intensity vs axial depth');
    save_fig(f6,d,'Intensity_vs_axial_preamp_vs_nopreamp.png');

    % --- element comparison (no preamp) ---
    [Is41,It41]=arrayfun(@(v) intens(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt_41el.hws',v)),sens_no,depth,fawf,PRFeff),Vel);
    [Is79,It79]=arrayfun(@(v) intens(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt_79el.hws',v)),sens_no,depth,fawf,PRFeff),Vel);
    f7=figure('Color','w','Position',[60 60 1000 440]);
    subplot(1,2,1); hold on; grid on; box on
    yline(Ilim,'k--','limit 190','HandleVisibility','off');
    plot(Vn,IsN,'s-','Color',cB,'MarkerFaceColor',cB,'DisplayName','61 el'); plot(Vel,Is41,'^-','Color',cG,'MarkerFaceColor',cG,'DisplayName','41 el'); plot(Vel,Is79,'d-','Color',cO,'MarkerFaceColor',cO,'DisplayName','79 el');
    xlabel('TX voltage (V)'); ylabel('I_{sppa.3} (W/cm^2)'); title('I_{sppa.3}'); legend('Location','northwest'); xlim([0 52]);
    subplot(1,2,2); hold on; grid on; box on
    plot(Vn,ItN,'s-','Color',cB,'MarkerFaceColor',cB); plot(Vel,It41,'^-','Color',cG,'MarkerFaceColor',cG); plot(Vel,It79,'d-','Color',cO,'MarkerFaceColor',cO);
    xlabel('TX voltage (V)'); ylabel('I_{spta.3} (mW/cm^2)'); title('I_{spta.3} (burst; limit 720, far above)'); xlim([0 52]);
    sgtitle('S5-1 push, no preamp: intensity vs TX voltage, 41/61/79 el');
    save_fig(f7,d,'Intensity_vs_TXvoltage_41_61_79el.png');

    % ================= numeric summary =================
    fprintf('\nMI at chosen location (derated, depth %.0f mm):\n V | preamp | noPreamp | 41el | 79el\n',depth);
    for i=1:numel(Vn)
        s1=sprintf('%.3f',MIn(i)); s2='   -  '; j=find(Vp==Vn(i)); if ~isempty(j), s2=sprintf('%.3f',MIp(j)); end
        s3='  -  '; s4='  -  '; k=find(Vel==Vn(i));
        if ~isempty(k), s3=sprintf('%.3f',MI41(k)); s4=sprintf('%.3f',MI79(k)); end
        fprintf('%2d | %6s | %8s | %s | %s\n',Vn(i),s2,s1,s3,s4);
    end
end

% ===================================================================
function v = vneg(file)
    addpath(fileparts(mfilename('fullpath')));
    w = readHWS(file); v = -min(w(1).y - median(w(1).y));
end
function [dep,MI] = axial(d,pat,ks,sens,fawf)
    Hs=[]; n=numel(ks); dep=nan(n,1); MI=nan(n,1);
    for i=1:n
        [w,c]=readHWS(fullfile(d,sprintf(pat,ks(i)))); yb=w(1).y-median(w(1).y);
        h=regexp(c.notes,'Home\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        p=regexp(c.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        H=str2double(h); P=str2double(p); dep(i)=norm(P-H);
        MI(i)=(-min(yb))/sens/db2mag(0.3*dep(i)/10*fawf)/sqrt(fawf);
    end
    [dep,o]=sort(dep); MI=MI(o);
end
function [Isppa,Ispta] = intens(file,sens,depth,fawf,PRFeff)
%INTENS  Derated Isppa.3 [W/cm2] and Ispta.3 [mW/cm2] from one capture.
    addpath(fileparts(mfilename('fullpath'))); rho=1000; c=1500;
    w=readHWS(file); yb=w(1).y-median(w(1).y); dt=w(1).dt;
    p=yb/sens*1e6; I=p.^2/(rho*c); cR=rescale(cumsum(I),0,1);
    i10=find(cR>0.1,1); i90=find(cR>0.9,1); PD=1.25*(i90-i10)*dt; PII=trapz(I(i10:i90))*dt;
    Ider=(1/db2mag(0.3*depth/10*fawf))^2; PII3=PII*Ider;
    Isppa=PII3/PD/1e4; Ispta=Isppa*PD*PRFeff*1000;
end
function [dep,Isppa,Ispta] = axialInt(d,pat,ks,sens,fawf,PRFeff)
    n=numel(ks); dep=nan(n,1); Isppa=nan(n,1); Ispta=nan(n,1);
    for i=1:n
        [~,c]=readHWS(fullfile(d,sprintf(pat,ks(i))));
        h=regexp(c.notes,'Home\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        pp=regexp(c.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        dep(i)=norm(str2double(pp)-str2double(h));
        [Isppa(i),Ispta(i)]=intens(fullfile(d,sprintf(pat,ks(i))),sens,dep(i),fawf,PRFeff);
    end
    [dep,o]=sort(dep); Isppa=Isppa(o); Ispta=Ispta(o);
end

function fawf = measureFawf(file)
    addpath(fileparts(mfilename('fullpath')));
    w=readHWS(file); y=w(1).y-median(w(1).y); dt=w(1).dt;
    Y=abs(fft(y.*hann(numel(y)))); f=(0:numel(Y)-1)'/(numel(Y)*dt); hb=1:floor(numel(Y)/2);
    bf=f(hb)>1.3e6&f(hb)<3.5e6; pk=max(Y(hb).*bf); idx=find(Y(hb)>=pk/2&bf); fawf=(f(idx(1))+f(idx(end)))/2/1e6;
end
function [Moc,Ch,gain,Ca] = cal(ft)
    S=load(miRoot('calibration','HydrophoneTables.mat'),'AmplifierTable','HydrophoneTable');
    HT=S.HydrophoneTable; AT=S.AmplifierTable;
    Moc=db2mag(HT{find(HT.FREQ_MHz>=ft,1),"SENS_DB"})*1e6; Ch=HT{find(HT.FREQ_MHz>=ft,1),"CAP_PF"};
    gain=db2mag(AT{find(AT.FREQ_MHZ>=ft,1),"GAIN_DB"}); Ca=AT{find(AT.FREQ_MHZ>=ft,1),"CAP_PF"};
end
function save_fig(f,~,name)
    d = miRoot('results'); if ~isfolder(d), mkdir(d); end
    out=fullfile(d,name); exportgraphics(f,out,'Resolution',150); fprintf('saved: %s\n',out);
end
