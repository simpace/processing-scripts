function X = deal_names(str,dirname)

if nargin==1;
    d = dir(str);
else
    if ~isempty(strmatch(dirname(end),'/'));%==1;
        d = dir([dirname str]);
    else %strmatch(dirname(end),'/')==0;
        d = dir([dirname '/' str]);
    end
end

if ~isempty(d);
    [X{1:length(d)}] = deal(d.name);
    X = X';
    if length(X)==1;
        X = X{1};
    end
else
    warning('No desired files in directory')
    X = {};
end
    