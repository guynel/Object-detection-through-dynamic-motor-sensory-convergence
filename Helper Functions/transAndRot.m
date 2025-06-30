function [varargout] = transAndRot(transformMat,varargin)
% Used in Fig. 1; For video frames that were rotated, this also rotates associated traces
% Given a 3x3 transformation matrix 'transformMat', for each variable in the input,
% transform it then output it. If input is a cell matrix, perform
% transformation for each cell. All input should be [x,y] coordinates in vector
% of size Nx2. 

for i = 1:length(varargin)
    if ~iscell(varargin{i}) % If cell-matrix, do for each cell separately
        varargout{i} = transAndRotSingleArg(varargin{i},transformMat);
    else % Otherwise, just do for this argument
        varargout{i} = cellfun(@(x) transAndRotSingleArg(x,transformMat),varargin{i},'Uni',0);
    end
end
end

function [output] = transAndRotSingleArg(input,transformMat)

output = input'; % Shift dimensions for multiplication
output(3,:) = 1; % Pad third dimension with ones, as transformMat is 3x3
output = transformMat*output;
output = output(1:2,:); % Get back [x,y] observations
output = output'; % Go back to original orientation

end