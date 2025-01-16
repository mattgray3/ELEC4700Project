%ELEC  4700 Project Milestone 1 Source Function Matt Gray 101183570

function E = SourceFct(t, InputParas) %taking inputs for source function

if isfield(InputParas, 'rep') %Checking to see if InputParas are fields or structs
    n = floor(t/InputParas.rep);
    t = t - n*InputParas.rep;
end

if ~isstruct(InputParas) %Checking to see if InputParas are fields or structs
    E = InputParas;
else
    E = InputParas.E0 * exp(-(t-InputParas.t0)^2 / InputParas.wg^2)* ... 
        exp(1i*(InputParas.we*t + InputParas.phi));
end
end