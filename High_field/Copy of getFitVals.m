function [vals, unc] = getFitVals(params, conf, parameter_name)
% After 'lorentzian_fit' get a value from the fit and its uncertainty.
switch parameter_name
    case 'Peak'
        idx=1;
    case 'Width'
        idx=2;
    case 'Amplitude'
        idx=3;
    case 'Offset'
        vals = params(end);
        unc = abs(conf(end,2)-conf(end,1));
        return
    otherwise
        error('Unknown Paramter of Lorentzian.')
end
if isempty(conf)
    unc=[];
else
    unc = 0.5*abs(conf(idx:3:(end-1),2)-conf(idx:3:(end-1),1));
end

vals = params(idx:3:(end-1));
end