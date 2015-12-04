function X = deal_names_full(str,dirname)

if nargin==1;
    d = dir(str);
    dirname = [pwd '/'];
else
    if isempty(strmatch(dirname(end),'/'));%==1;
        dirname = [dirname '/'];
    end
    
    d = dir([dirname str]);
    
end

if ~isempty(d);
    [X{1:length(d)}] = deal(d.name);
    X = X';
    if length(X)==1;
        X = [dirname X{1}];
    else
        for irow = 1:length(X)
            X{irow} = [dirname X{irow}];
        end
    end
else
    warning('No desired files in directory')
    X = {};
end
    