function TecPIV_TrackLagrangianTracers(TracerInit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% find list of unique parameters to track

TableParam = TracerInit(:, 5:end)  % Tracer num, Frame num, then x and y coord mean that we start list of param names at column 5
ListParam=[TableParam{:}];
ListUniqueParam = unique(ListParam);

for i = 1:length(ListUniqueParam)
    Param = ListUniqueParam(i)
    
    % find the tracers that include this param
    S = size(TableParam);
    ListTracers = [];
    ListStartFrame = [];
    for j=1:S(2)
        for k=1:S(1)
            if TableParam{k,j} == Param
                ListTracers = [ListTracers k];
                ListStartFrame = [ListStartFrame TracerInit{2,j}]
            end
        end
    end
    
    ListTracers
    ListStartFrame
    
    
end

end

