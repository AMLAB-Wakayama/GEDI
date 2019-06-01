%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A multi-resolution gammachirp envelope distortion index (mr-GEDI) for open software
%   Objective mesure index for speech intelligibility
%
%   Yamamoto
%   Created : 26 Dec. 2018; based on GEDI.m
%   Modified: 01 June 2019; added the GCresp for inputs of mrGEDI_OutdcGC.m
%
%   Inputs:
%       SndTest:  input signal of enhanced/unprocessed noisy speech
%       SndClean: input signal of clean speech reference
%       GCparam:  parameters of dynamic compressive gammachirp filterbank
%       (dcGC-FB)
%       Conditions: parameters for a material, [k q m sigma_s fs]
%       SPL: sound pressure level of input signals in dB
%
%   Outputs:
%       Result.
%           Pcorrect: percent correct of speech intelligibility
%           SDRenv: SDRenv in dB
%           SDRenvsModFB: SDRenvs in all modulation filter channels
%           ParamModFB: parameters of modulation filterbank
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Result = mrGEDI(SndTest, SndClean, GCparam, Conditions, SPL)

%% Default parameters of dcGC filterbank
if isfield(GCparam,'fs')  == 0, GCparam.fs = 48000; end
if isfield(GCparam,'NumCh')  == 0, GCparam.NumCh = 100; end
if isfield(GCparam,'FRange')  == 0, GCparam.FRange = [100, 6000]; end
if isfield(GCparam,'OutMidCrct')  == 0,  GCparam.OutMidCrct = 'ELC'; end
if isfield(GCparam, 'Ctrl') == 0,GCparam.Ctrl = 'dynamic'; end     % Cf: GCparam.Ctrl = 'static'; % or 'fixed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processsing for input signals
% Extract material information
fs = Conditions(5); % sampling frequency of original materials

% Time alignment
[SndTest, SndClean] = GEDI_TimeAlign(SndTest, SndClean);

% Upsampling to 48 kHz
SndTest = interp(SndTest,GCparam.fs/fs);
SndClean = interp(SndClean,GCparam.fs/fs);

% Calibrate input level of SndTest and SndClean by using Meddis hair cell level
[SndTest, ~]   = Eqlz2MeddisHCLevel(SndTest, SPL);
[SndClean, ~]   = Eqlz2MeddisHCLevel(SndClean, SPL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyzing by dynamic compressive gammachirp filterbank
% Test signal (enhanced/unprocessed speech)
[OutdcGCTest, ~, ~, GCresp] = GCFBv211(SndTest',GCparam);

% Reference signal (Clean speech; S)
[OutdcGCClean, ~, ~, ~] = GCFBv211(SndClean',GCparam);

%% Main processing part of mr-GEDI
Result = mrGEDI_OutdcGC(OutdcGCTest, OutdcGCClean, GCparam, GCresp, Conditions);

end

function [OutSndTest, OutSndClean] = GEDI_TimeAlign(InSndTest, InSndClean)
%% Align signals using cross-correlation

[CorrSnd, Lag] = xcorr(InSndTest, InSndClean);
[~,IdxMaxLag] = max(CorrSnd);
TimeLag = Lag(IdxMaxLag);

if TimeLag > 0
    OutSndTest = InSndTest(TimeLag:end);
    LenOutSndTest = length(OutSndTest);
    OutSndClean = InSndClean(1:LenOutSndTest);
elseif TimeLag < 0
    OutSndClean = InSndClean(-TimeLag:end);
    LenOutSndClean = length(OutSndClean);
    OutSndTest = InSndTest(1:LenOutSndClean);
else
    OutSndTest = InSndTest;
    OutSndClean = InSndClean;
end

end