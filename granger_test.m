cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';

cfg.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];
                      
cfg.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];
                        
cfg.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];

data              = ft_connectivitysimulation(cfg);
%%
%% Fit autoregressive model

cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);
%%
%%
cfg           = [];
cfg.method    = 'granger';
granger1       = ft_connectivityanalysis(cfg, mfreq);


%% Nonparametric freq analysis (MTMFFT)
cfg           = [];
cfg.method    = 'mtmfft';
%cfg.method    = 'mtmconvol';
%cfg.pad = 'nextpow2';
cfg.pad = 10;

cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi=[0:1:100];
cfg.tapsmofrq = 2;
freq_mtmfft          = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq_mtmfft);
%%
%% Nonparametric freq analysis (MTMconvol)
cfg           = [];
%cfg.method    = 'mtmfft';
cfg.method    = 'mtmconvol';
%cfg.pad = 'nextpow2';
cfg.pad = 10;

cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi=[0:2:100];
 cfg.tapsmofrq = 2;
cfg.t_ftimwin=ones(length(cfg.foi),1).*(0.5);
%cfg.t_ftimwin=1000./cfg.foi;
%cfg.tapsmofrq = 0.4*cfg.foi;


%cfg.t_ftimwin=7./cfg.foi;

% cfg.toi='50%';
cfg.toi=linspace(min(data.time{1}),max(data.time{1}),5);
freq_mtmfft          = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq_mtmfft);
%%
%% Nonparametric freq analysis (Wavelet)
cfg           = [];
%cfg.method    = 'mtmfft';
cfg.method    = 'wavelet';
%cfg.pad = 'nextpow2';
cfg.pad = 10;

cfg.width=2;

cfg.taper     = 'dpss';
cfg.output    = 'powandcsd';
cfg.foi=[0:5:100];
cfg.tapsmofrq = 2;
cfg.t_ftimwin=ones(length(cfg.foi),1).*(0.5);
%cfg.t_ftimwin=7./cfg.foi;

% cfg.toi='50%';
cfg.toi=linspace(min(data.time{1}),max(data.time{1}),10);
freq_mtmfft          = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq_mtmfft);
%%
%% Nonparametric freq analysis (Wavelet OTHER)
cfg           = [];
%cfg.method    = 'mtmfft';
cfg.method    = 'tfr';
%cfg.pad = 'nextpow2';
cfg.pad = 10;

%cfg.width=2;

cfg.taper     = 'dpss';
cfg.output    = 'powandcsd';
cfg.foi=[0:5:100];
cfg.tapsmofrq = 2;
cfg.t_ftimwin=ones(length(cfg.foi),1).*(0.5);
%cfg.t_ftimwin=7./cfg.foi;

% cfg.toi='50%';
cfg.toi=linspace(min(data.time{1}),max(data.time{1}),10);
freq_mtmfft          = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq_mtmfft);



%%
% % % % % freq = 
% % % % % 
% % % % %   struct with fields:
% % % % % 
% % % % %             label: {3×1 cell}
% % % % %            dimord: 'rpttap_chan_freq'
% % % % %              freq: [1×101 double]
% % % % %     fourierspctrm: [1500×3×101 double]
% % % % %         cumsumcnt: [500×1 double]
% % % % %         cumtapcnt: [500×1 double]
% % % % %               cfg: [1×1 struct]
%% Nonparametric freq analysis 
cfg           = [];
cfg.method    = 'mtmconvol';
cfg.foi          = 1:1:100;
cfg.taper     = 'dpss';
cfg.output    = 'fourier';

%cfg.taper     = 'hanning';
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window = 0.5 sec
%cfg.toi=0:0.0050:0.9950;
%cfg.toi=0.50:0.0050:0.9000;
%cfg.toi='50%';
cfg.toi=0:0.05:0.9950;
% % 
% cfg.output    = 'pow';

cfg.tapsmofrq = 10;
% cfg.toi='50%'; 
%cfg.toi          = linspace(min(data.time{1}),max(data.time{1}),100);
%cfg.t_ftimwin    =10./cfg.foi;

freq_mtmconvol          = ft_freqanalysis(cfg, data);
 %%
% freq_mtmconvol.powspctrm
chinga=freq_mtmconvol.powspctrm;

%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq_mtmconvol);
%%
chingada=granger.grangerspctrm
%% Nonparametric freq analysis (WAVELET)
cfg           = [];
cfg.method    = 'wavelet';
cfg.width=7;
cfg.pad=2;
cfg.output='powandcsd';
cfg.foi          = 0:1:100;
% cfg.taper     = 'dpss';
% cfg.output    = 'fourier';
cfg.toi=0:0.5:0.9950;
%cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window = 0.5 sec
freq          = ft_freqanalysis(cfg, data);
%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);
%%
%% Nonparametric freq analysis (OTHER WAVELET)
cfg           = [];
cfg.method    = 'tfr';
% cfg.width=15;
% cfg.gwidth=5;

cfg.foi          = 0:5:100;
%cfg.foilim=[10 100];
% cfg.taper     = 'dpss';
cfg.output    = 'powandcsd';
%andcsd
cfg.toi=0:0.1:0.9950;
%cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window = 0.5 sec
freq          = ft_freqanalysis(cfg, data);


%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);


%%
%cfg.pad         =  'nextpow2';
 %cfg.toi=(min(data.time{1})::max(data.time{1}),0.005);
cfg.foi          = 20:20:100;
%cfg.toi=7./cfg.foi;
% cfg.tapsmofrq = 10;

cfg.toi='all';
%cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window = 0.5 sec

% cfg.taper     = 'dpss';
cfg.output    = 'freq';
% cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);
%%
% 
% T=(length(data1)*(1/fs)); % Signal duration
% W=1; %Bandwidth

movingwin=[.25 0.01]; %Moving time window 

% nw=T*W; %Time-bandwidth product
% k=2*nw-1; %Number of tapers to use

%Function parameters
params.tapers=[3 5]; %Taper specs
params.fpass=[0 100]; %Visualization specs
params.Fs=data.fsample; %Sampling freq.
%params.err=[1 0.05]; %Chi-squared error with p=0.05


params.pad=3; %Zero-padding  
params.trialave=0; %No averaging of trials.

% ese=cell(1,3,1);
% for h=1:3
ch=data.trial{1};
[S,t,f] = mtspecgramc(ch(h,:) , movingwin, params ); %Main function computing spectrogram
% ese{h}=S;
% end
ese=nan(length(data.label),length(f),length(t));

for h=1:length(data.label)
 [S,t,f] = mtspecgramc(ch(h,:) , movingwin, params ); %Main function computing spectrogram
S=S.';
 ese(h,:,:)=S;   
end
%%
freqq=freq;
freqq.time=t;
freqq.freq=f;
freqq.powspctrm=ese;
%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freqq);

%%
% plot_matrix(S,t,f,'n') %Visualization of Spectrogra

%%
% cfg           = [];
% cfg.method    = 'coh';
% coh           = ft_connectivityanalysis(cfg, freq);
% cohm          = ft_connectivityanalysis(cfg, mfreq);
% 
% cfg           = [];
% cfg.parameter = 'cohspctrm';
% cfg.zlim      = [0 1];
% ft_connectivityplot(cfg, coh, cohm);

%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);
%%
granger1       = ft_connectivityanalysis(cfg, mfreq);




%%
for j=1:39
for row=1:3
for col=1:3
  subplot(3,3,(row-1)*3+col);
  %plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:,j)),'LineWidth',1,'Color',[0 0 0])
   hold on
%   %plot(granger.freq, squeeze(chocho(row,col,:)))
%  %hold on
    %plot(granger1.freq, squeeze(granger1.grangerspctrm(row,col,:)))
    ylim([0 1])
end
end
end
%%
ssn=cell(21,1);
% sumame=0;
for j=1:21;
   ssn{j}= squeeze(granger.grangerspctrm(:,:,:,j));
% sumame=sumame+squeeze(granger.grangerspctrm(:,:,:,j));
end


%%
chocho=zeros(3,3,100);
conta=0;
for i=1:3
    for j=1:3
        
        for l=1:21
            cha=ssn{l};
               Q=cha(i,j,:);
               Q(isnan(Q))=0;
            if l==1            
            cho=(Q*0);
            else
            cho=cho+(Q);
            end
            if l==21
                conta=conta+1;
               chocho(i,j,:)=(cho)/21; 
            end
        end
    end
end

%%
chale=cell2mat(ssn);
%%
ssnm=max(ssn);
plot(granger.freq,ssnm)
%%
averno=granger.grangerspctrm;
averno=squeeze(granger.grangerspctrm(:,:,:,j));
avern=squeeze(granger.grangerspctrm(row,col,:));
%%
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, granger, granger1);
%%
% % % %%
% % % figure
% % % plot(data.time{1}, data.trial{1}) 
% % % legend(data.label)
% % % xlabel('time (s)')
% % % %%
% % % cfg = [];
% % % cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
% % % ft_databrowser(cfg, data);
% % % %% Fit autoregressive model
% % % 
% % % cfg         = [];
% % % cfg.order   = 5;
% % % cfg.toolbox = 'bsmart';
% % % mdata       = ft_mvaranalysis(cfg, data);
% % % 
% % % cfg        = [];
% % % cfg.method = 'mvar';
% % % mfreq      = ft_freqanalysis(cfg, mdata);
% % % %% Calculate the cross-spectral density . THe spectral factorization of this gives the spectral transfer function. 
% % % 
% % % cfg           = [];
% % % cfg.method    = 'mtmfft';
% % % cfg.taper     = 'dpss';
% % % cfg.output    = 'fourier';
% % % cfg.tapsmofrq = 2;
% % % freq          = ft_freqanalysis(cfg, data);
% % % %%
% % % cfg           = [];
% % % cfg.method    = 'coh';
% % % coh           = ft_connectivityanalysis(cfg, freq);
% % % cohm          = ft_connectivityanalysis(cfg, mfreq);
% % % %%
% % % cfg           = [];
% % % cfg.parameter = 'cohspctrm';
% % % cfg.zlim      = [0 1];
% % % ft_connectivityplot(cfg, coh, cohm);
% % % %% 
% % % cfg           = [];
% % % cfg.method    = 'granger';
% % % granger       = ft_connectivityanalysis(cfg, mfreq);
% % % 
% % % cfg           = [];
% % % cfg.parameter = 'grangerspctrm';
% % % cfg.zlim      = [0 1];
% % % ft_connectivityplot(cfg, granger);
% % % %%
% % % figure
% % % for row=1:3
% % % for col=1:3
% % %   subplot(3,3,(row-1)*3+col);
% % %   plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
% % %   ylim([0 1])
% % % end
% % % end
% % % 
% % % %%
% % % sigcell=cell(1,63);
% % % timecell=cell(1,63);
% % % 
% % % for i=1:63
% % % sig9=Mono9{i};
% % % sig12=Mono12{i};
% % % sig17=Mono17{i};
% % % tt=0:length(sig9)-1;
% % % tt=tt*(1/fn);
% % % timecell{1,i}=tt;
% % % 
% % % sig=[sig9 sig12 sig17];
% % % sig=sig.';
% % % sigcell{1,i}=sig;
% % % end
% % % 
% % % 
% % % %%
% % % % sig9=Mono9{1};
% % % % sig12=Mono12{1};
% % % % sig17=Mono17{1};
% % % % sig=[sig9 sig12 sig17];
% % % % sig=sig.';
% % % data1.trial=sigcell;
% % % data1.time=timecell;
% % % data1.fsample=fn;
% % % data1.label=cell(3,1);
% % % data1.label{1}='PFC';
% % % data1.label{2}='Parietal';
% % % data1.label{3}='Hippocampus';
% % % %%
% % % epoch=32;
% % % figure
% % % plot(data1.time{epoch}, data1.trial{epoch}) 
% % % legend(data1.label)
% % % xlabel('time (s)')
% % % %%
% % % 
% % % cfg = [];
% % % cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
% % % ft_databrowser(cfg, data1);
% % % %% Fit autoregressive model
% % % 
% % % cfg         = [];
% % % cfg.order   = 5;
% % % cfg.toolbox = 'bsmart';
% % % mdata1       = ft_mvaranalysis(cfg, data1);
% % % 
% % % cfg        = [];
% % % cfg.method = 'mvar';
% % % mfreq1      = ft_freqanalysis(cfg, mdata1);
% % % %%
% % % cfg           = [];
% % % cfg.method    = 'mtmfft';
% % % cfg.taper     = 'dpss';
% % % cfg.output    = 'fourier';
% % % cfg.tapsmofrq = 2;
% % % freq1          = ft_freqanalysis(cfg, data1);
%% RESULTADO DE GRANGER 
%% 
ro=1000;
[p,q,timecell]=getwin(carajo,veamos,sig1,sig2,label1,label2,ro);
%%
gc(q,timecell,'Bandpassed',ro)

gc(p,timecell,'Widepass',ro)
%%
ro=1000;
[p,q,timecell]=getwinbip(carajo,veamos,sig1,sig2,label1,label2,ro);

gcbip(q,timecell,'Bandpassed',ro)
gcbip(p,timecell,'Widepass',ro)
%%
Ro=[200 500 1000];
for i=1:3
  ro=Ro(i);  
  [p,q,timecell]=getwin(carajo,veamos,sig1,sig2,label1,label2,ro);
  close all
  gc(q,timecell,'Bandpassed',ro)
  
 
end


%%
ro=1000;

if ro==200;
i=1;
[p1, p4, z1, z4]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});


i=3;
[p2, p5, z2, z5]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=7;
[p3, p6, z3, z6]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=11;
[p7, p8, z7, z8]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
end
    
if ro==500;
i=1;
[p1, p4, z1, z4]=generate500(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});


i=3;
[p2, p5, z2, z5]=generate500(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=7;
[p3, p6, z3, z6]=generate500(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=11;
[p7, p8, z7, z8]=generate500(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
end

if ro==1000;
i=1;
[p1, p4, z1, z4]=generate1000(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});


i=3;
[p2, p5, z2, z5]=generate1000(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=7;
[p3, p6, z3, z6]=generate1000(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=11;
[p7, p8, z7, z8]=generate1000(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
end
    
    
    

%% NUEVO. AQUI EMPIEZA LO CHIDO 
p=cell(length(z1),1);
q=cell(length(z1),1);
timecell=cell(length(z1),1);


for i=1:length(z1)
%    p{i,1}=[z1{i}.';z2{i}.';z3{i}.'];
%    q{i,1}=[z4{i}.';z5{i}.';z6{i}.'];
   p{i,1}=[z1{i}.';z2{i}.';z3{i}.';z7{i}.']; %Widepass
   q{i,1}=[z4{i}.';z5{i}.';z6{i}.';z8{i}.']; %Bandpassed
   timecell{i,1}=[0:400]*(1/fn)-0.2;

end

p=p.';
q=q.';
timecell=timecell.';

%%
gc(q,timecell,'Bandpassed')

gc(p,timecell,'Widepass')
%%
sigcell1=cell(3,1);
sigcell2=cell(3,1);

vecon=[1 3 7]';
i=2;
[p1, p4, z1, z4]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=5;
[p2, p5, z2, z5]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

i=9;
[p3, p6, z3, z6]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});

p=cell(length(z1),1);
q=cell(length(z1),1);
timecell=cell(length(z1),1);


for i=1:length(z1)
   p{i,1}=[z1{i}.';z2{i}.';z3{i}.'];
   q{i,1}=[z4{i}.';z5{i}.';z6{i}.'];
%    p{i,1}=[z1{i}.';z2{i}.';z3{i}.';z7{i}.']; %Widepass
%    q{i,1}=[z4{i}.';z5{i}.';z6{i}.';z8{i}.']; %Bandpassed
   timecell{i,1}=[0:400]*(1/fn)-0.2;

end

p=p.';
q=q.';
timecell=timecell.';

%%
gcbip(q,timecell,'Bandpassed')

gcbip(p,timecell,'Widepass')


% % % %%
% % % data1.trial=q;
% % % data1.time= timecell; %Might have to change this one 
% % % data1.fsample=fn;
% % % data1.label=cell(3,1);
% % % data1.label{1}='Hippocampus';
% % % data1.label{2}='Parietal';
% % % data1.label{3}='PFC';
% % % data1.label{4}='Reference';
% % % 
% % % 
% % % %%
% % % 
% % % % %%
% % % % pe1=[p1;p2; p3];
% % % % pe2=[p4;p5; p6];
% % % % 
% % % % cellpe=cell(1,1);
% % % % cellpe{1}=pe1;
% % % % %data1.trial=cellpe;
% % % % data1.trial=pe1;
% % % % 
% % % % data1.time=[0:length(p3)-1]*(1/fn); %Might have to change this one 
% % % % data1.fsample=fn;
% % % % data1.label=cell(3,1);
% % % % data1.label{1}='Hippocampus';
% % % % data1.label{2}='Parietal';
% % % % data1.label{3}='PFC';
% % % %% VISUALIZATION
% % % 
% % % %%
% % % figure
% % % plot(data1.time{1}, data1.trial{1}) 
% % % legend(data1.label)
% % % xlabel('time (s)')
% % % %%
% % % 
% % % cfg = [];
% % % cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
% % % ft_databrowser(cfg, data1);
% % % %% Fit autoregressive model
% % % 
% % % cfg         = [];
% % % cfg.order   = 5;
% % % cfg.toolbox = 'bsmart';
% % % mdata       = ft_mvaranalysis(cfg, data1);
% % % 
% % % cfg        = [];
% % % cfg.method = 'mvar';
% % % mfreq      = ft_freqanalysis(cfg, mdata);
% % % %% Non parametric 
% % % cfg           = [];
% % % cfg.method    = 'mtmfft';
% % % cfg.taper     = 'dpss';
% % % cfg.output    = 'fourier';
% % % cfg.tapsmofrq = 2;
% % % freq          = ft_freqanalysis(cfg, data1);
% % % %% Coherence
% % % cfg           = [];
% % % cfg.method    = 'coh';
% % % coh           = ft_connectivityanalysis(cfg, freq);
% % % cohm          = ft_connectivityanalysis(cfg, mfreq);
% % % %%
% % % cfg           = [];
% % % cfg.parameter = 'cohspctrm';
% % % cfg.zlim      = [0 1];
% % % ft_connectivityplot(cfg, coh, cohm);
% % % %%
% % % cfg           = [];
% % % cfg.method    = 'granger';
% % % granger       = ft_connectivityanalysis(cfg, mfreq);
% % % 
% % % cfg           = [];
% % % cfg.parameter = 'grangerspctrm';
% % % cfg.zlim      = [0 1];
% % % ft_connectivityplot(cfg, granger);
% % % %%
% % % lab=cell(16,1);
% % % lab{1}='Hippo -> Hippo';
% % % lab{2}='Hippo -> Parietal';
% % % lab{3}='Hippo -> PFC';
% % % lab{4}='Hippo -> Reference';
% % % 
% % % lab{5}='Parietal -> Hippo';
% % % lab{6}='Parietal -> Parietal';
% % % lab{7}='Parietal -> PFC';
% % % lab{8}='Parietal -> Reference';
% % % 
% % % lab{9}='PFC -> Hippo';
% % % lab{10}='PFC -> Parietal';
% % % lab{11}='PFC -> PFC';
% % % lab{12}='PFC -> Reference';
% % % 
% % % lab{13}='Reference -> Hippo';
% % % lab{14}='Reference -> Parietal';
% % % lab{15}='Reference -> PFC';
% % % lab{16}='Reference -> Reference';
% % % 
% % % 
% % % figure
% % % conta=0;
% % % for row=1:length(data1.label)
% % % for col=1:length(data1.label)
% % %   subplot(length(data1.label),length(data1.label),(row-1)*length(data1.label)+col);
% % %   plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
% % %   ylim([0 1])
% % %   xlim([0 300])
% % %   xlabel('Frequency (Hz)')
% % %   grid minor
% % %   conta=conta+1;
% % %   
% % %  if conta==1 || conta==6 || conta==11 || conta==16 
% % %  set(gca,'Color','k')
% % %  end
% % % %   if conta==4
% % % %       error('stop')
% % % %   end
% % %   title(lab{conta})
% % % end
% % % end
% % % mtit('Monopolar','fontsize',14,'color',[1 0 0],'position',[.5 1 ])
% % % mtit('Bandpassed','fontsize',14,'color',[1 0 0],'position',[.5 0.75 ])
% % % mtit('(+/-200ms)','fontsize',14,'color',[1 0 0],'position',[.5 0.5 ])
% % % 
%% REAL DATA 
label='Bandpassed';
% (q,timecell,label,ro)
fn=1000;
data1.trial=q;
data1.time= timecell; %Might have to change this one 
data1.fsample=fn;
data1.label=cell(3,1);
data1.label{1}='Hippocampus';
data1.label{2}='Parietal';
data1.label{3}='PFC';
data1.label{4}='Reference';
%%
%Parametric model
cfg         = [];
cfg.order   = 10;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data1);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);
%%
%Non parametric
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data1);
%%
%Non parametric- Multitaper
cfg           = [];
cfg.method    = 'mtmconvol';
%cfg.foi          = 1:10:500;
cfg.foi          = 10:10:50;

cfg.taper     = 'hanning';
cfg.output    = 'fourier';
%cfg.t_ftimwin    = ones(length(cfg.foi),1).*.5;   % length of time window = 0.5 sec
cfg.t_ftimwin    = [0.3 0.15 0.1 1 1 ];
%cfg.t_ftimwin    = [0.6 0.30 0.2];

%cfg.t_ftimwin    =[.2 .4];
cfg.toi='50%'; 
%cfg.toi=-0.2:.01:0.2;

cfg.tapsmofrq = 2;
cfg.pad='nextpow2',30;

freq1          = ft_freqanalysis(cfg, data1);
%freq1.freq=round(freq1.freq);
% freq1.freq=10:10:300;
freq1.freq
%%
cfg = [];
cfg.baseline     = [-0.2 0.2]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-3e-27 3e-27];	        
cfg.showlabels   = 'yes';	
 cfg.layout       = 'CTF151.lay';
figure 
ft_multiplotTFR(cfg, freq1);

%%
% % % % %  label: {4×1 cell}
% % % % %            dimord: 'rpttap_chan_freq_time'
% % % % %              freq: [1×40 double]
% % % % %              time: [-0.2000 -0.1500 -0.1000 -0.0500 0 0.0500 0.1000 0.1500 0.2000]
% % % % %     fourierspctrm: [77×4×40×9 double]
% % % % %         cumtapcnt: [77×40 double]
% % % % %               cfg: [1×1 struct]

%%


cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);
granger1       = ft_connectivityanalysis(cfg, mfreq);
granger2       = ft_connectivityanalysis(cfg, freq1);
%%

lab=cell(16,1);
lab{1}='Hippo -> Hippo';
lab{2}='Hippo -> Parietal';
lab{3}='Hippo -> PFC';
lab{4}='Hippo -> Reference';

lab{5}='Parietal -> Hippo';
lab{6}='Parietal -> Parietal';
lab{7}='Parietal -> PFC';
lab{8}='Parietal -> Reference';

lab{9}='PFC -> Hippo';
lab{10}='PFC -> Parietal';
lab{11}='PFC -> PFC';
lab{12}='PFC -> Reference';

lab{13}='Reference -> Hippo';
lab{14}='Reference -> Parietal';
lab{15}='Reference -> PFC';
lab{16}='Reference -> Reference';

%%
% figure
conta=0;
figure('units','normalized','outerposition',[0 0 1 1])
% for j=1:length(freq1.time)
%     conta=0;
for row=1:length(data1.label)
for col=1:length(data1.label)
    
  subplot(length(data1.label),length(data1.label),(row-1)*length(data1.label)+col);
  
  plot(granger2.freq, squeeze(granger2.grangerspctrm(row,col,:,2)),'LineWidth',.1,'Color',[0 0 0])
  hold on
  plot(granger1.freq, squeeze(granger1.grangerspctrm(row,col,:)))
  %hold on
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
 
  %plot(granger2.freq, squeeze(granger2.grangerspctrm(row,col,:)))
  
  ylim([0 1])
  xlim([0 300])
  xlabel('Frequency (Hz)')
  grid minor
  conta=conta+1;
  
 if conta==1 || conta==6 || conta==11 || conta==16
   legend('NP:Multitaper','Parametric: AR(10)','NP:MTMFFT')  
%  legend('Parametric: AR(10)','NP:MTMFFT')    
 set(gca,'Color','k')
 end
%   if conta==4
%       error('stop')
%   end
  title(lab{conta})
end
end
%end
mtit('Monopolar','fontsize',14,'color',[1 0 0],'position',[.5 1 ])
mtit(label,'fontsize',14,'color',[1 0 0],'position',[.5 0.75 ])
mtit(strcat('(+/-',num2str(ro),'ms)'),'fontsize',14,'color',[1 0 0],'position',[.5 0.5 ])

