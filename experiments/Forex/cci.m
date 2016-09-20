function out = cci(data,period)
% Function to calculate the Commodity Channel Index of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the number of periods over which to calculate the CCI
%
% Example:
% out = cci(data,period)
%

% Error check
if nargin ~= 2
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end
if (numel(period) ~= 1) || (mod(period,1)~=0)
    error('The period must be a scalar integer.');
end
if length(data) < (2*period-1)
    error('The length of the data must be at least (2*period-1).');
end

% calculate the CCI
smavg = sma(data,period);
ld = length(data);
out = nan*ones(ld,1);
for idx = 2*period-1:ld
    out(idx) = sum(abs(smavg(idx)-data(idx-period+1:idx)));
end
idx = 2*period-1:ld;
out(idx) = (data(idx)-smavg(idx))./(0.015*out(idx)/period);
