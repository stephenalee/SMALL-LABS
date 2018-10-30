function [out]=MakeSkipFrames(mainfold)
%Makes skipframes and the movie .mat files for SMALL labs


%% find all the movie files
[datalist,dataloc,~]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin;*.mat'],...
    'Select the movie files','multiselect','on');

if ~dataloc
    display('no data selected')
    return
end

if ~iscell(datalist); datalist={datalist}; end
datalist = cellfun(@(x)[dataloc, x],datalist,'uniformoutput',false);
[dlocs,dnames,exts]=cellfun(@fileparts,datalist,'uniformoutput',false);

%%
for ii=1:numel(dlocs)
    try
        load([dlocs{ii},filesep,dnames{ii},'.mat'],'mov');
    catch
        Movie2mat([dlocs{ii},filesep,dnames{ii},exts{ii}]);
        load([dlocs{ii},filesep,dnames{ii},'.mat'],'mov');
    end
    ints=squeeze(sum(sum(mov)));
    figure(2)
    plot(ints); 
    title('click above the base value to remove activation frames');axis tight;
    clickY=[];
    [~,clickY]=ginput(1);
    if ~isempty(clickY)
        SkipFrames=find(ints>clickY);
    end
    goodframe=true(1,size(mov,3));
    goodframe(SkipFrames)=0;
    
    save([dlocs{ii},filesep,dnames{ii},'.mat'],'goodframe','-append')
end
close all force
end