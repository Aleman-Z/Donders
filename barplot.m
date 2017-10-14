function barplot(p4,p3)
fn=1000;
data=p4;

clear params
params.tapers=[3 9];
params.Fs=fn;
params.err=0;
params.fpass=[0 30];
params.pad=5; %-1 means no zero-padding  
params.trialave=0; %No averaging of trials.
%movingwin=[.015 .002];
movingwin=[.2 .002];


[S,t,f] = mtspecgramc( data, movingwin, params );
%subplot(1,2,1)
subplot(2,2,3)

%t=linspace(-0.2,0.2, length(t));

t=linspace(-0.5,0.5, length(t));
plotter1(S,t,f)
title('Spectrogram Wideband')
xticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5])
xticklabels({'-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0', '0.1', '0.2', '0.3', '0.4', '0.5'})


data=p3;

clear params
params.tapers=[3 5];
params.Fs=fn;
params.err=0;
params.fpass=[100 300];
params.pad=8; %-1 means no zero-padding  
params.trialave=0; %No averaging of trials.
%movingwin=[.03 .001];
movingwin=[.2 .002];

[S,t,f] = mtspecgramc( data, movingwin, params );
%subplot(1,2,2)
subplot(2,2,4)
%t=linspace(-0.2,0.2, length(t));

t=linspace(-0.5,0.5, length(t));
plotter1(S,t,f)
xticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5])
xticklabels({'-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0', '0.1', '0.2', '0.3', '0.4', '0.5'})

title('Spectrogram Bandpass')

end