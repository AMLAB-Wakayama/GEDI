%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GEDI_OutdcGC for open software
%   Main processing part of gammachirp envelope distortion index (GEDI)
%   using outputs of dynamic compressive gammachirp filterbank (dcGC-FB)
%
%   Katsuhiko Yamamoto,  Irino, T.
%   Created:    26 Jan 2018 based on GEDI_OutdcGC_v1a
%   Modified:    9 Dec 2019  norminv_erf and normcdf_erf are included 
%                                     to avoid Statistics Toolbox, (Irino, T.)
%               12 Jul 2020 add makeCoefERBwidth to calculate with
%               weighting function of auditory filter channels
%
%   Inputs:
%       OutdcGCTest:  output of dcGC-FB (enhanced/unprocessed noisy speech)
%       OutdcGCClean: output of dcGC-FB (clean speech reference)
%       GCparam: parameters of dcGC-FB
%       Conditions: parameters for a material, [k q m sigma_s fs]
%
%   Outputs:
%       Pcorrect: percent correct of speech intelligibility
%       SDRenv: SDRenv in dB
%       SDRenvsModFB: SDRenvs in all modulation filter channels
%       ParamModFB: parameters of modulation filterbank
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pcorrect, SDRenvdB, SDRenvsModFB, ParamModFB] = GEDI_OutdcGC(OutdcGCTest, OutdcGCClean, GCparam, Conditions)

%%%%%%%%%%%%%%%%
fsdcGC = GCparam.fs;    % sampling frequency of OutdcGCTest/OutdcGCClean: 48 kHz
fs = Conditions(5);     % sampling frequency of original material
fcMod = 150;            % cutoff frequency of envelope filter
%%%%%%%%%%%%%%%%

nChCal = 1:GCparam.NumCh; % dcGC-100ch

% Center frequencies of modulation filterbank
mod_fcs=[ 2 4 8 16 32 64 ]; % BPFs
mod_Chs = [1 mod_fcs];  % LPF + BPFs

[bzLPF, apLPF] = butter(1, fcMod/(fs/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Envelope extraction and modulation filterbank
% Calculating envelopes of temporal outputs
PenvsTest    = -20*ones(length(mod_Chs),length(nChCal));
PenvsClean  = -20*ones(length(mod_Chs),length(nChCal));
PenvsDist   = -20*ones(length(mod_Chs),length(nChCal));

for nCh = nChCal % 1:GCparam.NumCh
    
    % Calculate envelope
    tmp1        = abs(hilbert(OutdcGCTest(nCh,:)));   %  envelope Test
    tmp2        = decimate(tmp1,fsdcGC/fs);
    EnvTest      = filter(bzLPF,apLPF,tmp2);
    
    tmp1        = abs(hilbert(OutdcGCClean(nCh,:))); %  envelope Clean
    tmp2        = decimate(tmp1,fsdcGC/fs);
    EnvClean    = filter(bzLPF,apLPF,tmp2);
    
    %% Calculate envelope sample-by-sample
    EnvDist = sqrt(abs(EnvClean.^2 - EnvTest.^2));
    
    % Modulation filterbank
    % Calculation of envelope power in 7 modulation filter
    [~ , ExcPtnTest, DCpower]  =  modfilterbankepsmYK_v1(EnvTest, 0, 0,fs);
    % Normalized by the DCpower of the EnvTest
    [~ , ExcPtnClean, ~]    =  modfilterbankepsmYK_v1(EnvClean, 0, DCpower,fs);
    [~ , ExcPtnDist, ~]     =  modfilterbankepsmYK_v1(EnvDist, 0, DCpower,fs);
    
    % The noisefloor restricted to minimum 0.01 (-20 dB) reflecting and internal noise threshold
    ExcPtnClean    = max(ExcPtnClean,0.01);
    ExcPtnTest      = max(ExcPtnTest,0.01);
    ExcPtnDist     = max(ExcPtnDist,0.01);
    
    PenvsClean(:,nCh)   = ExcPtnClean';
    PenvsTest(:,nCh)     = ExcPtnTest';
    PenvsDist(:,nCh)    = ExcPtnDist';
    
end

%% Calculate signal-to-distortion ratios of modulation domain, SDRenvs
% Weighting for auditory filters
Weight = zeros(length(mod_Chs),length(nChCal));
tmpWeight = makeCoefERBwidth_v2(GCparam, GCresp, 0);
for nMod = 1:length(mod_Chs) 
    Weight(nMod,:) = tmpWeight;
end

%% Calculate signal-to-distortion ratios of modulation domain, SDRenvs
% Initialization
SumPenvsClean = zeros(7,1);
SumPenvsTest   = zeros(7,1);
SumPenvsDist  = zeros(7,1);

% Sum all outputs from each auditory filter in each modulation filters
for iModCh = 1:7
    SumPenvsClean(iModCh,1) = sum(PenvsClean(iModCh,:).*Weight(iModCh,:));
    SumPenvsTest(iModCh,1)   = sum(PenvsTest(iModCh,:).*Weight(iModCh,:));
    SumPenvsDist(iModCh,1)  = sum(PenvsDist(iModCh,:).*Weight(iModCh,:));
end

% SDRenvs in all modulation filter channels
SDRenvs      = SumPenvsClean ./ SumPenvsDist;

%% Combine the SDRenvs into a total of SDRenv
tmpSDRenvs                = sqrt(sum(SDRenvs.^2, 1));
LinearCombinedSDRenv      = sqrt(sum(tmpSDRenvs.^2, 2));

%% Ideal observer stage
% Converting SNRenv to percent correct
[Pcorrect, SDRenvdB] = speechpercentcorrect_YK_v3(LinearCombinedSDRenv,Conditions);

%% Save information
% Parameter of modulation filterbank (ModFB)
ParamModFB.mod_fcs = mod_fcs;
ParamModFB.aud_fcs = nChCal;
ParamModFB.GCparam = GCparam;

% Output of SDRenvs
SDRenvsModFB.SDRenvs = SDRenvs;
SDRenvsModFB.PenvsTest = PenvsTest;
SDRenvsModFB.PenvsClean = PenvsClean;
SDRenvsModFB.PenvsDist  = PenvsDist;

end

function [FcsModBPF, TotalOutModFilterPosPowerSpecSig, DC_power] = modfilterbankepsmYK_v1(EnvSig,fig,RefEnvDC,fs)

if mod(length(EnvSig),2) == 0
    %number is even
    EnvSig = EnvSig(1:end-1);
else
    %number is odd
end


Qs = 1;  % Quality factor of resonance circuit

LenEnv = length(EnvSig);
SpecSig = fft(EnvSig);
MagSpecSig  = abs(SpecSig) ;
PowerSpecSig = MagSpecSig.^2/LenEnv ;% power spectrum.
PosPowerSpecSig = PowerSpecSig(1:fix(LenEnv/2)+1) ;
%take positive frequencies only and mulitply by two-squared to get the same total energy
PosPowerSpecSig(2:end) = PosPowerSpecSig(2:end).* (2)  ;

FreqPosPowerSpecSig = linspace(0,fs/2,length(PosPowerSpecSig));
AllFreq = [FreqPosPowerSpecSig -1*fliplr(FreqPosPowerSpecSig(2:end))];


% Bandpass filter
% ELECTRICAL ENGINEERING: PRINCIPLES AND APPLICATIONS, Fourth Edition,
% by Allan R. Hambley, 2008 Pearson Education, Inc.
% band center frequencies of BPF
FcsModBPF = [2 4 8 16 32 64];

% Initialize transfer function
TransFuncs = zeros(length(FcsModBPF),length(AllFreq));

% Calculating frequency-dmoain transferfunction for each center frequency:
for k = 1:length(FcsModBPF)
    TransFuncs(k+1,2:end) = 1./(1+ (1j*Qs*(AllFreq(2:end)./FcsModBPF(k) - FcsModBPF(k)./AllFreq(2:end)))); % p287 Hambley.
end

% squared BPF filter magnitude transfer functions
WeightModFilter = (abs(TransFuncs)).^2;

% Lowpass filter squared transfer function: third order butterworth filter
% TF from: http://en.wikipedia.org/wiki/Butterworth_filter
% cutoff frequency of lowpassfilter:
FreqCutoffLPF = 1;
% order:
OrderLPF = 3;
WeightModFilter(1,:) =  1./(1+((2*pi*AllFreq/(2*pi*FreqCutoffLPF)).^(2*OrderLPF)));

%TransFuncs(1,:) = sqrt(WeightModFilter(1,:));

% initialize output product:
OutModFilterPosPowerSpecSig = zeros(length(FcsModBPF),length(FreqPosPowerSpecSig));
TotalOutModFilterPosPowerSpecSig = zeros(1,7);

% ------------ DC-power, --------------------------
% %here devided by two such that a fully modulated tone has an AC-power of 1.
if RefEnvDC == 0
    % None
    DC_power =  PosPowerSpecSig(1)/ LenEnv /2;
else
    % Use other reference DCpower
    DC_power = RefEnvDC;
end

% ------------------------------------------------

for k = 1:size(WeightModFilter,1)
    OutModFilterPosPowerSpecSig(k,:) = PosPowerSpecSig.*WeightModFilter(k,1:floor(end/2)+1);
    % Integration estimated as a sum from f > 0
    % integrate envelope power in the
    % passband of the filter. Index goes from 2:end since
    TotalOutModFilterPosPowerSpecSig(k) =  sum(OutModFilterPosPowerSpecSig(k,2:end))/ LenEnv / DC_power;
    % integration is for f>0
end

FcsModBPF = [1 FcsModBPF];

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

function [Pc_est, SDRenvdB] = speechpercentcorrect_YK_v3(SDRenv_lin,parameters)

k = parameters(1);
q = parameters(2);
m = parameters(3);
sigma_s = parameters(4);


%--------  Only used for output:
SDRenvdB = 10*log10(SDRenv_lin);

% ---------- Converting from SNRenv to d_prime  --------------
d_prime = k*(SDRenv_lin).^q;

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
% Un = 1*norminv(1-(1/m));  % 8 Dec 2019
Un = 1*norminv_erf(1-(1/m));
mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
% Pc_est = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100; % 8 Dec 2019
Pc_est = normcdf_erf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;
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



