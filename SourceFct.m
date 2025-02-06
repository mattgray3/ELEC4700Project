%ELEC  4700 Project Milestone 1 Source Function Matt Gray 101183570

function E = SourceFct(t, InputParas) %taking inputs for source function, time and parameters for source

if isfield(InputParas, 'rep') %Checking to see if InputParas has a field named rep
    n = floor(t/InputParas.rep);% Calculates how many full repetitions of the period have occurred in time 't'
    t = t - n*InputParas.rep; % Adjusts 't' to be within the current repetition period
end

if ~isstruct(InputParas) %Checking to see if InputParas are fields or NOT structs
    E = InputParas; % If InputParas is not a structure, assigns its value directly to E
else % If InputParas is a struct it does this loop instead of the one above
    E = InputParas.E0 * exp(-(t-InputParas.t0)^2 / InputParas.wg^2)* ... 
        exp(1i*(InputParas.we*t + InputParas.phi));  % Multiplies the amplitude e0 by a gaussian envelope centered at t0 with width wg ...
                                                     %Applies a complex exponential modulation with angular frequency 'we' and phase 'phi'
end
end