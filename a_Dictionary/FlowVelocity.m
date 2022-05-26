function velocity = FlowVelocity(ParaSet)
%Generate 1 period of biphasic flow velocity for simulation
% using 1ms step size
% ---------------------------
% INPUTS:
% ---------------------------
% * ParaSet:    * peakIn: mm/s
%               * peakOut:mm/s
%               * dcFlow: background constant flow mm/s
%               * timeRatio: of out and in flow, >1          
%               * period, default 1000 (ms)    
% ---------------------------
% OUTPUTS:
% ---------------------------
% * velocity in mm/s. vector. 1second with 1ms resolution
%
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%

peakIn = 6; peakOut = 3; dcFlow = 1; timeRatio = 3;

if isfield(ParaSet,'peakIn')
    peakIn = ParaSet.peakIn;
end
if isfield(ParaSet,'peakOut')
    peakOut = ParaSet.peakOut;
end  
if isfield(ParaSet,'dcFlow')
    dcFlow = ParaSet.dcFlow;
end      
if isfield(ParaSet,'timeRatio')
    timeRatio = ParaSet.timeRatio;
end    

if isfield(ParaSet,'period')
    period = ParaSet.period;
else
    period = 1000;
end  

time_in=round(period/(timeRatio+1));
time_out=period-time_in;

velocity = [peakIn*sin((1:time_in)/time_in*pi)...
    -peakOut*sin((1:time_out)/time_out*pi)]+dcFlow;


end

