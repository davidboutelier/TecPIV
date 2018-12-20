% Stand-Alone_tracers
function StandAloneLagTracers
close all
clear all
clc

% Folow a number of tracers
NTracers = 0;

% initialise a cell array containing teh initial positions of the tracers
TracerInit = {};

TracerInit = AddLagrangianTracer(NTracers, TracerInit)



end


