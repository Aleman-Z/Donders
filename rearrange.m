%% Rearrange

%Make labels
label1=cell(11,1);
label1{1}='Hippocampus';
label1{2}='Hippocampus';
label1{3}='Parietal';
label1{4}='Parietal';
label1{5}='Parietal';
label1{6}='Parietal';
label1{7}='PFC';
label1{8}='PFC';
label1{9}='PFC';
label1{10}='PFC';
label1{11}='Reference';

label2=cell(11,1);
label2{1}='Monopolar';
label2{2}='Bipolar';
label2{3}='Monopolar';
label2{4}='Mono+ICA';
label2{5}='Bipolar';
label2{6}='Bip+ICA';
label2{7}='Monopolar';
label2{8}='Mono+ICA';
label2{9}='Bipolar';
label2{10}='Bip+ICA';
label2{11}='Monopolar';

sig1=cell(11,1);

sig1{1}=Mono17;
sig1{2}=Bip17;
sig1{3}=Mono12;
sig1{4}=MonoR12;
sig1{5}=Bip12;
sig1{6}=BipSSS12;
sig1{7}=Mono9;
sig1{8}=MonoR9;
sig1{9}=Bip9;
sig1{10}=BipSSS9;
sig1{11}=Mono6;


sig2=cell(11,1);

sig2{1}=V17;
sig2{2}=S17;
sig2{3}=V12;
sig2{4}=R12;
sig2{5}=S12;
sig2{6}=SSS12;
sig2{7}=V9;
sig2{8}=R9;
sig2{9}=S9;
sig2{10}=SSS9;
sig2{11}=V6;
%% 200ms SWR from Bipolar Hippocampus
for i=1:11
[p3, p5]=generate(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
cd NewFolderSinFiltro
string=strcat('200_WAV_',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..     
% [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% plotwave(D1, D2, D3, D4,D5)
% mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
% cd WaveDec
% string=strcat('200_WD_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

end
close all


%% 500ms SWR from Bipolar Hippocampus
for i=1:11
[p3,p5]=generate500(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
cd NewFolderSinFiltro
string=strcat('500_WAV_',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..  

% [D1, D2, D3, D4,D5 ]=deco500(p3,p5);
% plotwave500(D1, D2, D3, D4,D5)
% mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
% cd WaveDec
% string=strcat('500_WD_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

end
close all


%% 1000ms SWR from Bipolar Hippocampus 
for i=1:11
[p3,p5]=generate1000(carajo,veamos, sig1{i},sig2{i},label1{i},label2{i});
%cd probando
cd NewFolderSinFiltro
string=strcat('1000_WAV_',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..     
% 
% [D1, D2, D3, D4,D5 ]=deco1000(p3,p5);
% plotwave1000(D1, D2, D3, D4,D5)
% mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
% cd WaveDec
% string=strcat('1000_WD_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

end
close all



%% 200ms Common 
for i=1:11
[p3, p5]=generate(car,veamos3, sig1{i},sig2{i},label1{i},label2{i});
% cd NewFolder
% string=strcat('200_CN_WAV_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

[D1, D2, D3, D4,D5 ]=deco(p3,p5);
plotwave(D1, D2, D3, D4,D5)
mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
cd WaveDec
string=strcat('200_WD_CD_',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..     

end
close all
%% 500ms Common 
for i=1:11
[p3,p5]=generate500(car,veamos3, sig1{i},sig2{i},label1{i},label2{i});
% cd NewFolder
% string=strcat('500_CN_WAV_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

[D1, D2, D3, D4,D5 ]=deco500(p3,p5);
plotwave500(D1, D2, D3, D4,D5)
mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
cd WaveDec
string=strcat('500_WD_CD_',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..     


end
close all
%% 1000ms Common  
for i=1:11
[p3,p5]=generate1000(car,veamos3, sig1{i},sig2{i},label1{i},label2{i})
% cd probando
% string=strcat('1000_CN_WAV_',label1{i},label2{i},'.png');
% saveas(gcf,string)
% cd ..     

[D1, D2, D3, D4,D5 ]=deco1000(p3,p5);
plotwave1000(D1, D2, D3, D4,D5)
mtit(strcat(label1{i},' (',label2{i},')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
cd WaveDec
string=strcat('1000_WD_CD',label1{i},label2{i},'.png');
saveas(gcf,string)
cd ..     

end
close all

%% Wavelet decomposition 
% % 
% % label1='Hippocampus';
% % label2='Bipolar';
% % [p3, p5]=generate(carajo,veamos, Bip17,S17,label1,label2);
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % 
% % %%
% % plotwave(D1, D2, D3, D4,D5)
% % %%
% % fn=1000;
% % data=D5;
% % 
% % [cfs,f] = cwt(data,fn,'TimeBandwidth',35);
% % 
% % t=linspace(-0.2,0.2, length(cfs));
% % helperCWTTimeFreqPlot(cfs,t,f,'surf','CWT of Quadratic Chirp','Seconds','Hz')
% % title('Spectrogram Wideband (High frequency resolution)')
% % 
% % %%
% % plot(D5)
% % 
% % %%
% % label1='Hippocampus';
% % label2='Monopolar';
% % [p3, p5]=generate(carajo,veamos, Mono17,V17,label1,label2)
% % 
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % 
% % %%
% % 
% % label1='Hippocampus';
% % label2='Bipolar';
% % generate(carajo,veamos, Bip17,S17,label1,label2)
% % cd NewFolder
% % string=strcat('WAV_',label1,label2,'_200');
% % saveas(gcf,string)
% % cd ..
% % %%
% % label1='Parietal';
% % label2='Monopolar';
% % [p3, p5]=generate(carajo,veamos, Mono12,V12,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='Parietal';
% % label2='Mon+ICA';
% % [p3, p5]=generate(carajo,veamos, MonoR12,R12,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='Parietal';
% % label2='Bipolar';
% % [p3, p5]=generate(carajo,veamos, Bip12,S12,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='Parietal';
% % label2='Bip+ICA';
% % [p3, p5]=generate(carajo,veamos, BipSSS12,SSS12,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='PFC';
% % label2='Monopolar';
% % [p3, p5]=generate(carajo,veamos, Mono9,V9,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='PFC';
% % label2='Mon+ICA';
% % [p3, p5]=generate(carajo,veamos, MonoR9,R9,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % 
% % %%
% % label1='PFC';
% % label2='Bipolar';
% % [p3, p5]=generate(carajo,veamos, Bip9,S9,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)
% % mtit(strcat(label1,' (',label2,')'),'fontsize',14,'color',[1 0 0],'position',[.5 1 ])
% % 
% % %%
% % label1='PFC';
% % label2='Bip+ICA';
% % [p3, p5]=generate(carajo,veamos, BipSSS9,SSS9,label1,label2)
% % [D1, D2, D3, D4,D5 ]=deco(p3,p5);
% % plotwave(D1, D2, D3, D4,D5)

%%
label1='Reference';
label2='Monopolar';
[p3, p5]=generate(carajo,veamos, Mono6 ,V6,label1,label2)
[D1, D2, D3, D4,D5 ]=deco(p3,p5);
plotwave(D1, D2, D3, D4,D5)
%% GRANGER CONNECTIVITY 

[GRANGER, V, N] = ft_connectivity_granger(H, Z, S, key1, value1, ...);

%   H is the spectral transfer matrix, Nrpt x Nchan x Nchan x Nfreq (x Ntime),
%        or Nrpt x Nchancmb x Nfreq (x Ntime). Nrpt can be 1.
%     Z is the covariance matrix of the noise, Nrpt x Nchan x Nchan (x Ntime),
%        or Nrpt x Nchancmb (x Ntime).
%     S is the cross-spectral density matrix, same dimensionality as H
 %%   
% Monopolar 
 Ro=[200 500 1000];
for i=1:3
  ro=Ro(i);  
  [p,q,timecell,cfs,f]=getwin(carajo,veamos,sig1,sig2,label1,label2,ro);
  close all
   q=cut(q);
   p=cut(p);
  timecell=cut(timecell);
   gc(q,timecell,'Bandpassed',ro)
   
 string=strcat(num2str(ro),'_GC_','Monopolar','Bandpassed','.png');
cd Nuevo
fig=gcf;
fig.InvertHardcopy='off';
saveas(gcf,string)
cd ..  
 close all
 
 gc(p,timecell,'Widepass',ro)

string=strcat(num2str(ro),'_GC_','Monopolar','Widepass','.png');
cd Nuevo
fig=gcf;
fig.InvertHardcopy='off';

saveas(gcf,string)
cd ..  
 close all
 
end

%%
%Bipolar
 Ro=[200 500 1000];
for i=1:3
  ro=Ro(i);  
  [p,q,timecell]=getwinbip(carajo,veamos,sig1,sig2,label1,label2,ro);
  close all
  q=cut(q);
  p=cut(p);
  timecell=cut(timecell);
  
%   gcbip(q,timecell,'Bandpassed',ro)
%   
%   string=strcat(num2str(ro),'_GC_','Bipolar','Bandpassed','.png');
% cd G
% fig=gcf;
% fig.InvertHardcopy='off';
% 
% saveas(gcf,string)
% cd ..  
%  close all
 
 gcbip(p,timecell,'Widepass',ro)
error('Stop here')
 string=strcat(num2str(ro),'_GC_','Bipolar','Widepass','.png');
cd G
fig=gcf;
fig.InvertHardcopy='off';

saveas(gcf,string)
cd ..  
 close all
 
end
