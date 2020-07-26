%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   test_mrGEDIfb
%   Test code for mr-GEDIfb
%
%   Irino, T.,  Katsuihiko Yamamoto 
%   Created : 16 May 2020   IT, based on test_mrGEDI (KY)
%   Modified: 16 May 2020   IT, introduction of Frame-base cGCFB -->  mrGEDI_OutdcGCfb (IT)
%   Modified: 21 May 2020   IT, using GCFBv221pack
%   Modified: 24 May 2020   IT, adding note
%   Modified: 26 May 2020   IT, debug
%
%   Note (24 & 26 Jul 2020 by IT): 
%     This mrGEDIfb is an initial version for fast processing. 
%     The results derived by the original mrGEDI and the mrGEDIfb were
%     very similar but not exactly the same in the preliminary tests
%     at a moderate SPL.  The parallel shift of percent correct line can be
%     easily compensated by the parameters of the ideal observer.
%     Please use this carefully for your purpose.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Environment settings
% Root
DirRoot = [getenv('HOME') '/Desktop/GEDI/'];
try
    chdir(DirRoot);
catch
    DirRoot = [pwd '/'];
end;
% Sounds
DirData = [DirRoot 'wav_sample/'];

% Package of dynamic compressive gammachirp filterbank
% Please download the "GCFBv221pack" and put into the "package" directory 
DirGCFB = [DirRoot 'package/GCFBv221pack/']; % essential
addpath(DirGCFB)  % modified by IT

%StrPath = path;
% if ~contains(StrPath,'GCFBv221pack/') == 1
%    addpath(genpath(DirRoot));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEDI and materials
% Parameters of dcGC filterbank
GCparam.Ctrl    = 'dynamic';
GCparam.NumCh   = 100;
GCparam.fs      = 48000;
% introduction of frame-base processing
% GCparam.DynHPAF.StrPrc = 'sample';  % sample-by-sample
GCparam.DynHPAF.StrPrc = 'frame';    % frame-base introduced by IT

% Parameter settings for materials
ParamSNR = [-6 -3 0 3]; %SNR between clean speech and noise
SPL = 65; % sound pressure level

% Parameter settings for materials
fs = NaN;  % setting it later. IT
Conditions = [1.50 0.5 20000 1.64 fs]; % [k q m sigma_s fs]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start simulation
Pcorrects = zeros(1,length(ParamSNR));

for i = 1:length(ParamSNR)
    
    % Index number of a sample speech file
    idxSp = i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test signal (enhanced/Unprocessed speech); the speech intelligiblity is calculated
    % Name of wav-file (example: '*/GEDI_Software/wav_sample/sample_sp1')
    strTest = [DirData 'sample_sp' num2str(idxSp)];
    % Read wav-file of test speech
    [SndTest, fs] = audioread([strTest '.wav']);
    disp(strTest);
    
    %% Reference signal (Clean speech)
    % Name of wav-file
    strClean = [DirData 'sample_sp_clean'];
    % Read wav-file of clean speech
    [SndClean, fs2] = audioread([strClean '.wav']);
    disp(strClean);
    if fs ~= fs2  %  IT 
        error('Inconsistency of sampling rate.');
    end;
    Conditions(5)  = fs; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-processing with a half-cosine function
    TimeHalfCos     = 0.02;    % [sec]    
    LenSndClean     = length(SndClean);
    LenHalfCos      = round(TimeHalfCos * fs);
    SigHalfCos      = hann(LenHalfCos*2);
    SigHalfCosUp    = SigHalfCos(1:LenHalfCos)';
    CompHalfCosFunc = [SigHalfCosUp ones(1,LenSndClean-LenHalfCos*2) fliplr(SigHalfCosUp)];
    SndClean    = SndClean .* CompHalfCosFunc';

    % Extract a speech segment for sample data
    TimeSndBefore   = 0.35;
    SndTest     = SndTest(fs*TimeSndBefore+1:fs*TimeSndBefore+LenSndClean) .* CompHalfCosFunc';    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Speech intelligibility prediction by mr-GEDI
    Result = mrGEDIfb(SndTest, SndClean, GCparam, Conditions, SPL);
    Pcorrect = Result.Pcorrect;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Pcorrects(idxSp) = Pcorrect;
    disp('==========================================');
    disp(['Percent correct:' num2str(Pcorrect) '(%)']);
    disp('==========================================');
    
end

%% Plot results
figure
plot(ParamSNR,Pcorrects);
xlim([-0.5+min(ParamSNR) max(ParamSNR)+0.5]);
ylim([-0.5+0 100+0.5]);
xlabel('SNR (dB)');
ylabel('Percent correct (%)')
legend('Unprocessed')
title('Examples of mr-GEDI');
grid on;

% Keep results for comparison    20 May 20
if 0
    NameRslt = 'Rslt1';
    save([NameRslt '_Val'])
    print([NameRslt '_Fig'],'-depsc')
end;
