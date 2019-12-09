%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mrGEDI_OutdcGC
%   Katsuhiko Yamamoto,  Irino, T.
%   Created:  12 Dec. 2017; based on GEDI_OutdcGC_v3d
%   Modified: 07 Feb. 2018; renamed mrGEDI_OutdcGC_v1h -> mrGEDI_OutdcGC
%   Modified: 01 Jul. 2018; limitations for modulation filter outputs
%             01 June 2019; added an inputs and corrected a weighting function
%   Modified: 9 Dec 2019; norminv_erf and normcdf_erf are included to avoid
%                Statistics Toolbox, (Irino, T.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = mrGEDI_OutdcGC(OutdcGCMix, OutdcGCClean, GCparam, GCresp, Conditions)

%%%%%%%%%%%%%%%%
% Define parameters of dcGC filterbank
if isfield(GCparam,'fs')  == 0, GCparam.fs = 48000; end
if isfield(GCparam,'NumCh')  == 0, GCparam.NumCh = 100; end
if isfield(GCparam,'FRange')  == 0, GCparam.FRange = [100, 6000]; end
if isfield(GCparam,'OutMidCrct')  == 0,  GCparam.OutMidCrct = 'ELC'; end
if isfield(GCparam, 'Ctrl') == 0,GCparam.Ctrl = 'dynamic'; end    % Cf: GCparam.Ctrl = 'static'; % or 'fixed
[Fr1, ~]  = EqualFreqScale('ERB',GCparam.NumCh,GCparam.FRange);
GCparam.Fr1 = Fr1;

% Frequency definition of each domaion
fsOrg   = Conditions(5);    % Sampling rate of original source (16000 Hz)
fsdcGC  = GCparam.fs;       % Sampling rate of the dcGC outputs (48000 Hz)
fsEnv   = fsOrg/10;         % Sampling rate of the temporal envelopes (1600 Hz)
fcEnv   = 150;              % Cutoff frequency of envelope LPF
%fcEnvDg = 30;               % delta-gamma—p@ cutoff freq. for delta-gamma filter ‚±‚ê‚Í‚È‚¾‚ç‚©‚È•Ï‰»

% Channel numbers of auditory filterbank
numChAud = 1:GCparam.NumCh; % 100ch

% Cutoff frequency of LP and center frequencies of BP mod filterbank
numFcMod = [1 2 4 8 16 32 64 128 256];  % 1 LPF + 8 BPFs
%%%%%%%%%%%%%%%%

% Lowpass filtering at fcMod
[bzLPF, apLPF] = butter(1, fcEnv/(fsEnv/2));

% Flame processings
lenAudModEnvs = ceil(length(OutdcGCClean)*(fsEnv/fsdcGC));
WinDurs      = 1./numFcMod; % (sec) :The window duration is the inverse of the centerfrequency of the modulation channel

lenWin = floor(WinDurs * fsEnv); % (samples)
numSegWin = ceil(lenAudModEnvs./lenWin); % The total number of segments is Nframes plus any additional "leftover"

PenvMix     = zeros(numSegWin(end),length(numFcMod),length(numChAud));
PenvClean   = PenvMix;
PenvDist    = PenvMix;

if find(lenWin == lenAudModEnvs) % If the duration of the stimulus is exactly equalt to the window duration
    segIdx = find(lenWin == lenAudModEnvs);
    numSegWin(segIdx) =  numSegWin(segIdx)-1;
end



%% -------------------------------------------------------------------%
% Envelope extraction and modulation filterbank
% All temporal envelopes filtered by auditory and modulation filters
AudModEnvsMix   = zeros(length(numFcMod), lenAudModEnvs, length(numChAud));
AudModEnvsClean = AudModEnvsMix;
AudModEnvsDist  = AudModEnvsMix;

cnt = 0;
for nAud = numChAud % 1:GCparam.NumCh
    cnt = cnt +1;
    
    %% Calculate temporal envelope
    % Mix/enhanced (target)
    tmp1    = abs(hilbert(OutdcGCMix(nAud,:)));      % Envelope mix
    tmp2    = decimate(tmp1,fsdcGC/fsOrg);          % Resampling to the original fs (48k -> 16k)
    tmp3    = decimate(tmp2,fsOrg/fsEnv);           % Resampling to the envlope fs (16k -> 1.6k)
    EnvMix  = filter(bzLPF,apLPF,tmp3);             % LPF
    
    % Clean (reference)
    tmp1        = abs(hilbert(OutdcGCClean(nAud,:)));   % Envelope Clean
    tmp2        = decimate(tmp1,fsdcGC/fsOrg);          % Resampling to the original fs (48k -> 16k)
    tmp3        = decimate(tmp2,fsOrg/fsEnv);           % Resampling to the envelope fs (16k -> 1.6k)
    EnvClean    = filter(bzLPF,apLPF,tmp3);             % LPF
    
    % Distortion (Yamamoto et al., 2017)
    % Temporal envelope distiortion calculated by sample-by-sample
    EnvDist = sqrt(abs(EnvClean.^2 - EnvMix.^2));
    
    
    %% Modulation filterbank
    %% Extract 9 tmporal envelopes filtered by 9 IIR-modulation filters
    ModEnvsMix   = modFbank_YK_v2(EnvMix,fsEnv,numFcMod);   % Mix/enhanced (target)
    ModEnvsClean = modFbank_YK_v2(EnvClean,fsEnv,numFcMod); % Clean (reference)
    ModEnvsDist  = modFbank_YK_v2(EnvDist,fsEnv,numFcMod); % Distortion (Yamamoto et al., 2017)
    
    % Save all outputs from modulation filters (length(numFcMod), length(EnvMix), length(numChAud))
    AudModEnvsMix(:,:,cnt)      = ModEnvsMix;
    AudModEnvsClean(:,:,cnt)    = ModEnvsClean;
    AudModEnvsDist(:,:,cnt)     = ModEnvsDist;
    
    % Caluculte DC powers of temporal envelopes
    DCPowModEnvMix   = (mean(EnvMix).^2)/2;
    
    for nMod = 1:length(numFcMod) %For each modulation channel
        %  Initialize temporary variables:
        tmpEnvMix   = zeros(lenWin(nMod),numSegWin(nMod));
        tmpEnvClean = tmpEnvMix;
        tmpEnvDist  = tmpEnvMix;
        segLengths  = zeros(1,numSegWin(nMod));
        
        for iSeg = 1:numSegWin(nMod) % For each temoral segment of the signal
            % find the start and end index of the frame
            if iSeg > (numSegWin(nMod)-1)
                startIdx = 1 + (iSeg-1)*lenWin(nMod);
                endIdx = lenAudModEnvs;
            else
                startIdx = 1 + (iSeg-1)*lenWin(nMod);
                endIdx = startIdx + lenWin(nMod)-1;
            end
            
            idxSeg = startIdx:endIdx;
            segLengths(iSeg) = length(idxSeg);
            
            tmpEnvMix(1:segLengths(iSeg),iSeg)   = AudModEnvsMix(nMod,idxSeg,nAud)-(sum(AudModEnvsMix(nMod,idxSeg,nAud))/segLengths(iSeg));
            tmpEnvClean(1:segLengths(iSeg),iSeg) = AudModEnvsClean(nMod,idxSeg,nAud)-(sum(AudModEnvsClean(nMod,idxSeg,nAud))/segLengths(iSeg));
            tmpEnvDist(1:segLengths(iSeg),iSeg)  = AudModEnvsDist(nMod,idxSeg,nAud)-(sum(AudModEnvsDist(nMod,idxSeg,nAud))/segLengths(iSeg));
            
        end % iSeg
        
        %% Normalize envelope powers based on the DC component
        PenvMix(1:numSegWin(nMod),nMod,nAud)    = sum(tmpEnvMix.*tmpEnvMix,1)./segLengths ./ DCPowModEnvMix; % computing the envelope power
        PenvClean(1:numSegWin(nMod),nMod,nAud)  = sum(tmpEnvClean.*tmpEnvClean,1)./segLengths ./ DCPowModEnvMix;
        PenvDist(1:numSegWin(nMod),nMod,nAud)   = sum(tmpEnvDist.*tmpEnvDist,1)./segLengths ./ DCPowModEnvMix;
        
        % Change NaN data to zero value
        if sum(sum(isnan(PenvMix(:,nMod,nAud))))
            PenvMix(isnan(PenvMix(:,nMod,nAud)),nMod,nAud) = 0;
        end
        
        if sum(sum(isnan(PenvClean(:,nMod,nAud))))
            PenvClean(isnan(PenvClean(:,nMod,nAud)),nMod,nAud) = 0;
        end
        
        if sum(sum(isnan(PenvDist(:,nMod,nAud))))
            PenvDist(isnan(PenvDist(:,nMod,nAud)),nMod,nAud) = 0;
        end
        
        % The envelope power of the noise is the minimum of Penv of the mixture or the noise.
        idxNonZeroPenvMix     =  PenvMix(:,nMod,nAud)~=0;
        idxNonZeroPenvClean   =  PenvClean(:,nMod,nAud)~=0;
        idxNonZeroPenvDist    =  PenvDist(:,nMod,nAud)~=0;
        
        % The envelope power cannot go below 0.001 (-30 dB) reflecting our minimum threshold of sensitivity to modulation detection
        threshold = 0.001;
        PenvClean(idxNonZeroPenvClean,nMod,nAud) = max(PenvClean(idxNonZeroPenvClean,nMod,nAud),threshold);
        PenvMix(idxNonZeroPenvMix,nMod,nAud)     = max(PenvMix(idxNonZeroPenvMix,nMod,nAud),threshold);
        PenvDist(idxNonZeroPenvDist,nMod,nAud)   = max(PenvDist(idxNonZeroPenvDist,nMod,nAud),threshold);
        
        % The modulation powers exceeding a quaprter of the center
        % frequency of the corresponding auditory fillter are not
        % considered in the model (v3, 01 July 18)
        if numFcMod(nMod) > GCparam.Fr1(nAud)/4
            PenvClean(:,nMod,nAud)  = 0;
            PenvMix(:,nMod,nAud)    = 0;
            PenvDist(:,nMod,nAud)   = 0;
        end
        
    end % nMod
    
end % nAud = numChAud

%% Calculate signal-to-distortion ratios of modulation domain, SDRenvs
% Weighting for auditory filter channels
Weight = zeros(numSegWin(end),length(numFcMod),length(numChAud));
tmpWeight = makeCoefERBwidth_v2(GCparam, GCresp, 0);
for nAud = numChAud
    Weight(:,:,nAud) = tmpWeight(nAud);
end

% average between auditory filters
SumPenvClean = sum(PenvClean.*Weight,3); % average between modulation filters
SumPenvDist  = sum(PenvDist.*Weight,3);
SDRenvSegMod = SumPenvClean ./ SumPenvDist;
SDRenvSegMod = max(0.001,SDRenvSegMod); % Truncated at -30 dB for numerical reasons.

SDRenvMod = zeros(1,length(numFcMod));
for nMod = 1:length(numFcMod)
    tmpSDRenvSegMod = SDRenvSegMod(:,nMod);
    idxSegWin = 1:numSegWin(nMod);
    SDRenvMod(1,nMod) = mean(tmpSDRenvSegMod(idxSegWin));
end

% Integrating the SNRenv across audio and modulation bands
SDRenv = sqrt(sum(SDRenvMod.^2));   % Combine across modulation filters
SDRenvdB = 10*log10(SDRenv);

% Convert the SDRenv to percent correct through an ideal observer model
Pcorrect = IdealObserver_v1(SDRenv,Conditions);

%% Save information
Output.Pcorrect = Pcorrect;
Output.SDRenvdB = SDRenvdB;

% Parameter of modulation filterbank (ModFB)
Output.mod_fcs = numFcMod;
Output.aud_fcs = numChAud;
Output.GCparam = GCparam;

end

function x_filt = modFbank_YK_v2(Env,fsEnv,cf_mod)
%
% Inputs:
%   Env:  The envelope to be filtered
%   fsMod: sampling frequency of the envelope
%   cf_mod:  centerfrequencies of the modulation filters
%
% Outputs:
%   x_filt:  Temporal outputs for each of the modulation filters
%
% Katsuhiko Yamamoto
% Created:  05 Sep 2017; based on modFbank_v3
% Modified: 07 Sep 2017; plot the transfer function of modFbank_YK
%           12 Sep 2017; add the version number of modFbank_YK (v1)
%           15 Feb 2018; 2nd IIR LPF -> 3rd IIR LPF (v2)

if nargin<3
    %band center frequencies
    cf_mod=[1 2 4 8 16 32 64 128 256];
end

if size(Env,1) > 1
    Env = Env';
end

x_filt = zeros(length(cf_mod),length(Env));

IIR_b = zeros(length(cf_mod),4);
IIR_a = zeros(length(cf_mod),4);

for i = 1:length(cf_mod)
    
    if cf_mod(i) == 1
        
        % Third order lowpass filter
        [b, a] = butter(3, cf_mod(i)/(fsEnv/2));
        b4 = b/a(1);
        a4 = a/a(1);
        
        IIR_b(i,:) = b4;
        IIR_a(i,:) = a4;
        
        % Filtering
        x_filt(i,:) = filter(b4, a4, Env);
        
    else
        
        % Pre-warping
        w0 = 2*pi*cf_mod(i)/fsEnv;
        
        % Bilinear z-transform
        W0 = tan(w0/2);
        
        % Second order bandpass filter
        Q = 1;
        B0 = W0/Q;
        b = [B0; 0; -B0];
        a = [1 + B0 + W0^2; 2*W0^2 - 2; 1 - B0 + W0^2];
        b3 = b/a(1);
        a3 = a/a(1);
        
        IIR_b(i,1:3) = b3;
        IIR_a(i,1:3) = a3;
        
        % Filtering
        x_filt(i,:) = filter(b3, a3, Env);
        
    end
    
end


% Plot frequency response of the digital filter
Sw = 0;
%Sw = 1;
if Sw == 1
    
    hold on
    for i = 1:length(cf_mod)
        
        % Frequency response of the digital filter
        w = 0:1/fsEnv:pi;
        
        if cf_mod(i) == 1
            IIR_TF = freqz(IIR_b(i,:),IIR_a(i,:),w);
        else
            IIR_TF = freqz(IIR_b(i,1:3),IIR_a(i,1:3),w);
        end
        
        wf = w*fsEnv/(2*pi);
        plot(wf,20*log10(abs(IIR_TF))); % Filter attenuation (dB)
        
    end
    
    hold off
    box on
    axis([0.25 max(cf_mod)*2 -40 5]);
    grid;
    set(gca,'xscale','log');
    set(gca,'xtick',cf_mod);
    xlabel('Frequency (Hz)');
    ylabel('Filter attenuation (dB)');
    Str_cf_mod = num2str(cf_mod');
    legend(Str_cf_mod,'location','southwest');
    title('Modulation filterbank');
    
end

end


function weight = makeCoefERBwidth_v2(GCparam, GCresp, SwPlot)
%
%   Katsuhiko YAMAMOTO
%   Created:  21 Dec. 2017
%   Modified: 21 Dec. 2017 (v1b); Added arguments
%             31 Dec. 2017 (v1c); Corrected a numerator value
%             01 June 2019 (v2);  Corrected a weighting function
%

if isfield(GCresp,'Fr1')  == 0
    [~, GCresp] = GCFBv211_SetParam(GCparam);
end

% Convert linear frequency to ERB
[~, ERBwidth] = Freq2ERB(GCresp.Fr1);
[~, ERBwidth1kHz] = Freq2ERB(1000);

% Weighting
weight = ERBwidth1kHz./ERBwidth;

if SwPlot == 1
    figure
    %plot(1:100,weight);
    plot(1:100,weight./max(weight));
    % In this plot, weight coefficients are normalized by the maximum value
    % but SumPenv***s don't change
    xlabel('Channel of dcGC filterbank');
    ylabel('Coefficient')
    legend({...
        '$\frac{\textrm{ERBwidth(1000 [Hz])}}{\textrm{ERBwidth(f [Hz])}}$'},...
        'interpreter','latex','fontsize',20);
end

end


function [ERBrate, ERBwidth] = Freq2ERB(cf)
%
%	Frequency -> ERB_N-rate and ERB_N-Bandwidth (Glasberg and Moore, 1990)
%	Toshio IRINO
%	11 Mar. 1998
%       Nodified: 26 Jul 2004 (no warning)
%       Nodified: 17 Nov 2006 (modified the comments only. ERB-> ERB_N)
%
%	function [ERBrate, ERBwidth] = Freq2ERB(cf),
%	INPUT	cf:       Center frequency
%	OUTPUT  ERBrate:  ERB_N rate
%		ERBwidth: ERB_N Bandwidth
%
%	Ref: Glasberg and Moore: Hearing Research, 47 (1990), 103-138
%            For different formulae (years), see Freq2ERBYear.m
%

if nargin < 1,  help Freq2ERB; end

ERBrate		= 21.4.*log10(4.37*cf/1000+1);
ERBwidth	= 24.7.*(4.37*cf/1000 + 1);

return % no warning

%%% Warning for Freq. Range %%%
cfmin = 50;
cfmax = 12000;
if (min(cf) < cfmin | max(cf) > cfmax)
    disp(['Warning : Min or max frequency exceeds the proper ERB range:']);
    disp(['          ' int2str(cfmin) '(Hz) <= Fc <=  ' int2str(cfmax) '(Hz).']);
end

end


function Pcorrect  = IdealObserver_v1(SDRenv_lin,parameters)
%%
% IdealObserver: Converts the overall SDRenv to percent correct.
%
% Usage: Pcorrect = IdealObserver(SDRenv_lin,parameters)
% Parameters :  vector with the parameters for the ideal Observer formatted as [k q m sigma_s]
%
% Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
% In Signal Detection and Recognition by Human Observers,
% edited by John A. Swets (John Wiley & Sons, New York)
%%

if nargin < 2
    error('You have to specify the k,q,m,sigma_s parameters for the IdealObserver')
end
k = parameters(1);
q = parameters(2);
m = parameters(3);
sigma_s = parameters(4);


% ---------- Converting from SNRenv to d_prime  --------------
d_prime = k*(SDRenv_lin).^q;

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
% Un = 1*norminv(1-(1/m));   % 8 Dec 2019
Un = 1*norminv_erf(1-(1/m));
mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
% Pcorrect = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;  % 8 Dec 2019
Pcorrect = normcdf_erf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;

end



%
%    norminv_erf  
%    Equivalent to norminv in Statistics Toolbox.
%    Irino, T
%    Created: 2 Dec 19
%    Modified: 2 Dec 19
%
function ni = norminv_erf(p)

ni = sqrt(2)*erfinv(2*p-1);

end

%    normcdf_erf   
%    Equivalent to normcdf in Statistics Toolbox.
%    Irino, T
%    Created: 2 Dec 19
%    Modified: 2 Dec 19
%    Modified: 9 Dec 19  % debug
%
%
function nc =normcdf_erf(x,mu,sigma)

if nargin==1
    mu=0;
    sigma=1;
elseif nargin==2
    sigma=1;
end
nc = (1+erf((x-mu)./(sigma*sqrt(2))))./2;

end
