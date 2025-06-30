function [ islands, output_marked ] = getOneIslands( binaryVec, connComp, padLeft, padRight)
% Get islands of one observations, helper indexing function

% Minimal island size
if nargin<2
    connComp = 0;
end

% We're finding flips from zeros to one and vice versa, so pad to faciliate
% indexing
binDiff = diff([zeros(size(binaryVec,1),1), binaryVec, zeros(size(binaryVec,1),1)],[],2);
islands = cell(1,size(binDiff,1));
% Generate for each row a cell with indices of islands
for i = 1:size(binDiff,1) % For each row
    islands{i}(:,1) = find( binDiff(i,:)>0 );
    islands{i}(:,2) = find( binDiff(i,:)<0 )-1;
end

% For each island found, only keep it if it exceeds minimal length
for i = 1:length(islands)
    counter = 1;
    while true
        if counter+1>size(islands{i},1) || counter>size(islands{i},1)
            break
        elseif abs(islands{i}(counter,2)-islands{i}(counter+1,1))<(connComp+1)
            islands{i}(counter,2) = islands{i}(counter+1,2);
            islands{i}(counter+1,:) = [];
        else
            counter = counter+1;
        end
        
    end
    
end

% Padding adjusment, if requested by user
if nargin>2
    for i = 1:length(islands)
        islands{i}(:,1) = islands{i}(:,1)-padLeft;
        islands{i}(:,2) = islands{i}(:,2)+padRight;
        islands{i}(:,1) = max(islands{i}(:,1),1);
        islands{i}(:,2) = min(islands{i}(:,2),size(binaryVec,2));
    end
end

% Return not just the islands as sets of incides, but also a marked output,
% with numbered islands in relevant places
output_marked = zeros(size(binaryVec));
for i = 1:length(islands)
    for j = 1:size(islands{i},1)
        output_marked(i,islands{i}(j,1):islands{i}(j,2)) = 1;
    end
end


end

