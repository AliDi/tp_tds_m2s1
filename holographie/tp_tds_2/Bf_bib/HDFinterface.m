function varargout = HDFinterface(varargin)
%HDFinterface
%
% PURPOSE
%    Interface to control HDF file for AS applications
%
% SYNOPSIS
%      [pathname,filename,datagroupname] = HDFinterface('select',HDFFilePath)
%      ischecked = HDFinterface('checkfile',pathfilename,AStype,datasetlist) ...,color) 
%      HDFinterface('showinfo',pathfilename) ...,color) 
%      parstruct = HDFinterface('getappparam',pathfilename)
%      [status,errtype] = HDFinterface('savefreqlist',pathfilename,cellfreqlist) ...,questbox)
%      cellfreqlist = HDFinterface('loadfreqlist',pathfilename)
%      [status,errtype] = HDFinterface('savedrawlist',pathfilename,drawstruct) 
%      drawstruct = HDFinterface('loaddrawlist',pathfilename)
%      AppStruct = HDFinterface('measdata',pathfilename) or ...,AppStruct)
%      [micropnts,coordsys,arraysys] = HDFinterface('micropnts',pathfilename)
%      [spectrum/quadsignal,frequency/time,axename] = HDFinterface('averageddata',pathfilename)
%      [specto,frequency,time] = HDFinterface('specgram',pathfilename) ...,prbvect) ...,seqvect) ...,nfft) ...,window) ...,numoverlap)
%
%
% DESCRIPTION
%     action -> 'select'       : select, open and control an HDF file
%               'showinfo'     : show file information
%               'getappparam'  : get application parameters need to the application process
%               'savefreqlist' : save the frequency list in (desc) annotation
%               'loadfreqlist' : load the frequency list from (desc) annotation
%               'savedrawlist' : save the drawstruct in (desc) annotation
%               'loadfreqlist' : load the drawstruct from (desc) annotation
%               'measdata'     : create or actualize the measdata field of an HDF data field
%               'micropnts'    : return coordinates of probes and informations on theirs arangements
%               'averageddata' : return the averaged spectrum/quadsignal data on the array and the corresponding frequency/time vector
%               'specgram'     : return the specto matrix and the corresponding frequency and time vectors (only with 'MeasTime' AStype)
%
%   for description see for each 'case'
%
% SEE ALSO
%        

% EXAMPLES   
%
% ALGORITHM
%
% REFERENCES
%    Matlab, Using Matlab, Version 6, 2002 (R13).
%
% Comments:   
%
% Matlab Toolbox functions used: makexyvector (measdata)
%
% External librairies and languages: none
%
% Test procedure: none
%
% Copyright 1997-2008 VisualVibroAcoustics
% jcp 20/12/02                                                                                       version 02/04/08
%--------------------------------------------------------------------------------------------------------------------
%  name  |   date   |  modification description
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 17/11/04 | add 'checkfile' action
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 12/11/05 | use uigetfile(searchpath,... instead of cd(searchpath), 'getappparam' : add parstruct.Fs
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 27/12/05 | add 'savefreqlist', 'loadfreqlist' (frequency list)
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 01/01/06 | add 'savedrawlist', 'loaddrawlist' (drawing)
%--------------------------------------------------------------------------------------------------------------------
%  jfl   | 17/01/06 | (1) case 'savefreqlist': in quesdlg, delete 'Cancel' button and correct error: status == 0 
%        |          | (2) case 'savedrawlist': correct error: status == 0 
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 11/05/06 | add 'measdata' (AppStruct)
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 25/07/06 | turn around the problem when the user enter '25°C' in ANstruct.ApplTemp instead '25' (getparam)
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 17/08/07 | change ANStationaryNAH to ASannotations, add 'micropnts', 'averagespect', 'spectgram' 
%--------------------------------------------------------------------------------------------------------------------
%  jcp   | 02/04/08 | add the field CreateStruct.Interpreter in internal function badfilewarn need for run R14 SP2
%--------------------------------------------------------------------------------------------------------------------


action = varargin{1};

switch lower(action)

%--------------------------------------------------------------------------------------------------------------------
case 'select'
%
%  [pathname,filename,datagroupname] = HDFinterface('select',HDFFilePath)
%  select, open and control an HDF file
%  HDFFilePath -> path to search the HDF file
HDFFilePath = varargin{2};
varargout{1} = []; varargout{2} = []; varargout{3} = [];

if ~isempty(HDFFilePath), searchpath = [HDFFilePath '*.hdf']; else, searchpath = '*.hdf'; end 
[filename, pathname] = uigetfile(searchpath,'Select an HDF-file');

if ~(isequal(filename,0)|isequal(pathname,0))
    HDFfileAScontrol('closeall');
    HDF_struct = HDFfileAScontrol('getinfo',[pathname,filename],'short');
    if ~isempty(HDF_struct)
        varargout{1} = pathname;
        varargout{2} = filename;
        varargout{3} = HDF_struct.attributes.name;
    else
        badfilewarn;    
    end
end
%--------------------------------------------------------------------------------------------------------------------
case 'measdata'
%
%      AppStruct = HDFinterface('measdata',pathfilename)   create AppStruct with field 'measdata'
%      AppStruct = HDFinterface('measdata',AppStruct)      actualize AppStruct.measdata
%
%    measured data             AppStruct.measdata.HDFFilePath    -> file path [set in ASSelectBarapp]
%                              AppStruct.measdata.HDFFileName    -> file name [set in ASSelectBarapp]
%                              AppStruct.measdata.DataGroupName  -> data group name [set in ASSelectBarapp] 
%                              AppStruct.measdata.Astype         -> AS type of HDF file
%                              AppStruct.measdata.SampleFrequency-> sampling frequency [set in 'reset')
%                              AppStruct.measdata.Frequency/Time -> frequency or time vector (depend on Astype)
%                              AppStruct.measdata.CoordSys       -> coordinates system
%                              AppStruct.measdata.ArraySys       -> array system ('rectangular' or 'cluster')
%                              AppStruct.measdata.Xvect          -> original x vector or cells of vectors
%                              AppStruct.measdata.Yvect          -> original y vector or cells of vectors
%                              AppStruct.measdata.Zvect          -> original z vector or cells of vectors
%
if isstruct(varargin{2})
    AppStruct = varargin{2}; 
    [logic,pathname] = isallfields(AppStruct,'measdata','HDFFilePath');
    [logic,filename] = isallfields(AppStruct,'measdata','HDFFileName');
    pathfilename = [pathname filename];
else
    pathfilename = varargin{2};
end
if ~isempty(pathfilename)
    HDF_struct = HDFfileAScontrol('getinfo',pathfilename);
    if ~isempty(HDF_struct)
        Nchild = length(HDF_struct.children);
        % determine file identifiers
        AppStruct.bidon = []; AppStruct = rmfield(AppStruct,'bidon');   % need when AppStruct don't exist
        if ~isallfields(AppStruct,'measdata'), AppStruct.measdata = []; end
        if ~isfield(AppStruct.measdata,'HDFFileName')|isempty(AppStruct.measdata.HDFFileName)
            AppStruct.measdata.HDFFileName = HDF_struct.properties.filename;
        end
        ind = findstr(HDF_struct.properties.filename,'.hdf');
        if isempty(ind)
            AppStruct.measdata.DataGroupName = HDF_struct.properties.filename(1:end);
        else
            AppStruct.measdata.DataGroupName = HDF_struct.properties.filename(1:ind(1)-1);
        end
        if ~isfield(AppStruct.measdata,'HDFFilePath')|isempty(AppStruct.measdata.HDFFilePath)
            AppStruct.measdata.HDFFilePath = [cd '\'];
        end
        
        %  determine AS type of HDF file
        AppStruct.measdata.Astype = HDF_struct.attributes.AStype;
        
        %  read sampling frequency
        if isallfields(HDF_struct,'attributes','sampling_frequency')&~isempty(HDF_struct.attributes.sampling_frequency)
            AppStruct.measdata.SampleFrequency = HDF_struct.attributes.sampling_frequency;
        else
            AppStruct.measdata.SampleFrequency = [];
        end
        
        % read the frequency or time vector
        if strcmp(AppStruct.measdata.Astype,'MeasSpct')
            axename = 'Frequency'; 
        elseif strcmp(AppStruct.measdata.Astype,'MeasTime')
            axename = 'Time';
        else
            error('incorrect AS type of HDF file'),
        end
        localflag = 0;
        for ii=1:Nchild
            if localflag==0&strcmp(HDF_struct.children(ii).children(1).properties.name,axename)
                AppStruct = setfield(AppStruct,'measdata',axename,HDF_struct.children(ii).children(1).properties.axevector);
                localflag = 1;
            end
        end
        
        % read coordsys attribute
        AppStruct.measdata.CoordSys = [];
        for ii=1:Nchild
            if strcmp(HDF_struct.children(ii).properties.name,'Coordinates')
                AppStruct.measdata.CoordSys = HDF_struct.children(ii).attributes.coordsys;
            end
        end
        
        % read arraysys attribute
        AppStruct.measdata.ArraySys = [];
        for ii=1:Nchild
            if strcmp(HDF_struct.children(ii).properties.name,'IndexMap')
                AppStruct.measdata.ArraySys = HDF_struct.children(ii).attributes.arraysys;  
            end
        end
        
        % determine Xvect Yvect Zvect
        Coordinates = HDFfileAScontrol('getdata',pathfilename,'Coordinates',':',':',':');
        IndexMap = HDFfileAScontrol('getdata',pathfilename,'IndexMap',':',':',':');
        [x,y,z] = makexyvector(Coordinates,IndexMap,AppStruct.measdata.CoordSys,AppStruct.measdata.ArraySys);
        AppStruct.measdata.Xvect = x;
        AppStruct.measdata.Yvect = y;
        AppStruct.measdata.Zvect = z;
    else
        badfilewarn;    
    end
end    
varargout{1} = AppStruct;

%--------------------------------------------------------------------------------------------------------------------
case 'checkfile'
%
%  ischecked = HDFinterface('checkfile',pathfilename,AStype,datasetlist) ...,color) 
%  check file in regard to application
%  pathfilename -> file name (in optional with path)
%  AStype       -> to compare with HDF_struct.attributes.AStype
%  datasetlist  -> cell of the names of needed datasets to compare with HDF_struct.properties.datasetnames
%  color        -> color of the message box, if [], the box is not shown

pathfilename = varargin{2};
AStype = varargin{3};
datasetlist = varargin{4};
if nargin>4, color = varargin{5}; else, color = [0.8 0.8 0.8]; end

ischecked = 0;
if ~isempty(pathfilename)
    HDF_struct = HDFfileAScontrol('getinfo',pathfilename);   % datasetnames is not set with 'short' (to modify)
    if ~isempty(HDF_struct)
        V = []; 
        for ii=1:length(datasetlist)
            ind = find(strcmp(HDF_struct.properties.datasetnames,datasetlist{ii})); 
            V = [V ind]; 
        end
        if strcmp(HDF_struct.attributes.AStype,AStype)&length(V)==length(datasetlist), ischecked = 1; end
    else
        file_AStype = 'unknow type';
    end
    if ischecked==0&~isempty(color) % open message box
        cmsgbox = cell(2,1); 
        cmsgbox{1} = ['This file (' file_AStype ') is not compatible'];
        cmsgbox{2} = 'with the current application.';
        CreateStruct.WindowStyle = 'modal';
        CreateStruct.Interpreter = 'tex';
        
        h = warndlg(cmsgbox,'Selected HDF file',CreateStruct);
        set(h,'Color',color);
    end
end
varargout{1} = ischecked;

%--------------------------------------------------------------------------------------------------------------------
case 'showinfo'
%
%  HDFinterface('showinfo',pathfilename) ...,color)
%  show file information
%  pathfilename -> file name (in optional with path)
%  color        -> color of the message box

pathfilename = varargin{2};
if nargin>2, color = varargin{3}; else, color = [0.75 0.87 1]; end

HDF_struct = HDFfileAScontrol('getinfo',pathfilename);   % datasetnames is not set with 'short' (to modify)
if ~isempty(HDF_struct)
    Name = HDF_struct.attributes.name;
    ind = findstr(Name,'\'); if ~isempty(ind), Name = Name(ind(end)+1:end); end
        cmsgbox = cell(9,1);      
        cmsgbox{1} = ['\rmName   : \bf',Name];
        cmsgbox{2} = ['\rmAStype : ',HDF_struct.attributes.AStype];
        cmsgbox{3} = ' ';
        cmsgbox{4} = ['frequency raws         : ',num2str(HDF_struct.attributes.frequency_raws)];
        cmsgbox{5} = ['probe number           : ',num2str(HDF_struct.attributes.probe_num)];
        cmsgbox{6} = ['auxiliary signals       : ',num2str(HDF_struct.attributes.aux_num)];
        cmsgbox{7} = ['sequence number   : ',num2str(HDF_struct.attributes.seq_num)];
        cmsgbox{8} = ['reference number    : ',num2str(HDF_struct.attributes.ref_num)];
        ind = find(strcmp(HDF_struct.properties.datasetnames,'spSpct'));
        if isempty(ind), workstr = 'none'; else, workstr = 'spSpct'; end
        cmsgbox{9} = ['work Data Set           : ',workstr];
        CreateStruct.WindowStyle = 'modal';
        CreateStruct.Interpreter = 'tex';
        
    h = msgbox(cmsgbox,'Selected HDF Data Set','help',CreateStruct);
    set(h,'Color',color);
else
    badfilewarn;    
end

%--------------------------------------------------------------------------------------------------------------------
case 'getappparam'
%
%  parstruct = HDFinterface('getappparam',pathfilename)
%  get application parameters need to the application process
%  parstruct    -> structure of parameters (returns parameters fields in all cases)
%                    ProjDistance : projection distance in cm
%                            Temp : temperature (°C)
%                              Fs : sampling frequency in Hz
%  pathfilename -> file name (in optional with path)

pathfilename = varargin{2};
parstruct = [];

%-- read annotation
HDF_struct = HDFfileAScontrol('getinfo',pathfilename); 

%-- read annotation
ANstruct = ASannotations('get',pathfilename);

%  determine the source-hologram distance
if ~isempty(ANstruct.ApplSrcp)
    Srcp = str2num(ANstruct.ApplSrcp{:}); 
    CoordMat = HDFfileAScontrol('getdata',pathfilename,'Coordinates',':',':',':');
    [Srcp,dist] = obtainSrcp(CoordMat,Srcp);
    parstruct.ProjDistance = dist;  
else
    parstruct.ProjDistance = [];
end

%  read the temperature
if ~isempty(ANstruct.ApplTemp)
    ApplTempStr = ANstruct.ApplTemp{:};
    n = 0;  while ~isempty(str2num(ApplTempStr(n+1)))&n<length(ApplTempStr), n = n + 1; end
    if n>1, parstruct.Temp = str2num(ApplTempStr(1:n)); else, parstruct.Temp = 20; end
else 
    parstruct.Temp = 20; 
end

%  read sampling frequency
if isallfields(HDF_struct,'attributes','sampling_frequency')&~isempty(HDF_struct.attributes.sampling_frequency)
    parstruct.Fs = HDF_struct.attributes.sampling_frequency;
else
    parstruct.Fs = [];
end

varargout{1} = parstruct;

%--------------------------------------------------------------------------------------------------------------------
case 'savefreqlist'
%
%  [status,errtype] = HDFinterface('savefreqlist',pathfilename,cellfreqlist,questbox)
%  [status,errtype] = HDFinterface('savefreqlist',pathfilename,cellfreqlist) equivalent to questbox = 0
%
%  store selected frequency list as a string str in annotations (questbox = 1 -> open questdlg)
%   
%    str = ['$$TestFreq$$ cellfreqlist{1} '$$' cellfreqlist{2} '$$' ... '$$' cellfreqlist{n} '$$']

pathfilename = varargin{2};
cellfreqlist = varargin{3};
if nargin>3, questbox = varargin{4}; else, questbox = 0; end

status = []; errtype = '';
 
if ~isempty(pathfilename)&iscell(cellfreqlist)        
    if questbox==1
        button = questdlg('Frequency selection list was modified : do you want to save it ?',...
                          'Frequency list','Yes','No','Yes');
    else
        button = 'Yes';
    end
    if strcmp(button,'Yes')
        str = '$$TestFreq$$'; 
        [status,errtype] = HDFfileAScontrol('delannotations',pathfilename,'file_desc',str);
        for ii=1:length(cellfreqlist), str = [str cellfreqlist{ii} '$$']; end 
        [status,errtype] = HDFfileAScontrol('setannotations',pathfilename,'file_desc',str);
        %elseif strcmp(button,'Cancel')
        %status==[]; errtype = 'Cancel';
    elseif strcmp(button,'No')
        status = 0; errtype = 'no executed';        % jf corrected on 17/01/06:  status==0; errtype = 'no executed';
    end
end
if nargin>0, varargout{1} = status; end
if nargin>1, varargout{2} = errtype; end

%--------------------------------------------------------------------------------------------------------------------
case 'loadfreqlist'
%
%      cellfreqlist = HDFinterface('loadfreqlist',pathfilename)
%
%  load frequency list from annotations and put in cellfreqlist
%   

pathfilename = varargin{2};

cellfreqlist = []; 

if ~isempty(pathfilename)        
    TestFreq = HDFfileAScontrol('getannotations',pathfilename,'file_desc','$$TestFreq$$');
    if ~isempty(TestFreq)                  % put in the FreqList
        freqlist = TestFreq{1}(13:end); 
        indfr = findstr(freqlist,'$$');
        if ~isempty(indfr)
            %
            %  read str = ['$$TestFreq$$ cellfreqlist{1} '$$' cellfreqlist{2} '$$' ... '$$' cellfreqlist{n} '$$']
            %
            cellfreqlist = {''};
            Ncell = length(indfr);
            indfr = [-1 indfr];
            for ii=2:Ncell+1, cellfreqlist{ii-1} = freqlist(indfr(ii-1)+2:indfr(ii)-1); end
        end    
    end 
end

varargout{1} = cellfreqlist; 

%--------------------------------------------------------------------------------------------------------------------
case 'savedrawlist'
%
%  [status,errtype] = HDFinterface('savedrawlist',pathfilename,drawstruct) ...,questbox)
%
%  store drawing structure drawstruct as a string in annotations (desc)
%   
%    str = ['$$TestDraw$$linecolor$$' mat2str(drawstruct.linecolor) '$$linewidth$$' num2str(drawstruct.linewidth) ...
%           '$$lines$$' mat2str(drawstruct.lines) '$$... ']

pathfilename = varargin{2};
drawstruct = varargin{3};
if nargin>3, questbox = varargin{4}; else, questbox = 0; end

status = []; errtype = '';
 
if ~isempty(pathfilename)&isstruct(drawstruct)        
    if questbox==1
        button = questdlg('Drawing was modified : do you want to save it ?',...
                          'Drawing list','Yes','No','Cancel','Cancel');
    else
        button = 'Yes';
    end
    if strcmp(button,'Yes')
        str = '$$TestDraw$$'; 
        [status,errtype] = HDFfileAScontrol('delannotations',pathfilename,'file_desc',str);
        drawfieldnames = fieldnames(drawstruct);
        for ii=1:length(drawfieldnames)
            if isfield(drawstruct,drawfieldnames{ii})
                str = [str drawfieldnames{ii} '$$' mat2str(getfield(drawstruct,drawfieldnames{ii})) '$$'];
            end
        end
        [status,errtype] = HDFfileAScontrol('setannotations',pathfilename,'file_desc',str);
    elseif strcmp(button,'Cancel')
        status = []; errtype = 'Cancel';                  % jf corrected on 17/01/06 : status==[]; errtype = 'Cancel';
    else
        status = 0; errtype = 'no executed';              % jf corrected on 17/01/06 : status==0; errtype = 'no executed';
    end
end

if nargin>0, varargout{1} = status; end
if nargin>1, varargout{2} = errtype; end

%--------------------------------------------------------------------------------------------------------------------
case 'loaddrawlist'
%
%      drawstruct = HDFinterface('loaddrawlist',pathfilename)
%
%  load drawing string from annotations and put in structure drawstruct
%

pathfilename = varargin{2};

drawstruct = []; 

if ~isempty(pathfilename)        
    TestDraw = HDFfileAScontrol('getannotations',pathfilename,'file_desc','$$TestDraw$$');
    if ~isempty(TestDraw)                  
        drawlist = TestDraw{1}(13:end); 
        indfr = findstr(drawlist,'$$');
        if ~isempty(indfr)
            Ncell = length(indfr);
            indfr = [-1 indfr];
            for ii=1:2:Ncell 
                drawstruct = setfield(drawstruct,drawlist(indfr(ii)+2:indfr(ii+1)-1), ...
                                                         str2num(drawlist(indfr(ii+1)+2:indfr(ii+2)-1))); 
            end
        end    
    end 
end

varargout{1} = drawstruct; 
%--------------------------------------------------------------------------------------------------------------------
case 'micropnts'
%
%  [micropnts,coordsys,arraysys] = HDFinterface('micropnts',pathfilename)
%  return coordinates of probes and informations on theirs arangements
%  micropnts    -> array of probe coordinates [Nprb,3]
%  coordsys     -> system of ccordinates 'cartesian'/'cylindrical'/ etc.
%  arraysys     -> system of array 'rectangular' or 'cluster'
%  pathfilename -> file name (in optional with path)

pathfilename = varargin{2};
micropnts = [];
arraysys = [];
coordsys = [];

%-- read info
HDF_struct = HDFfileAScontrol('getinfo',pathfilename); 
if ~isempty(HDF_struct)
    Nchild = length(HDF_struct.children);
    % read coordsys attribute
    coordsys = [];
    for ii=1:Nchild
        if strcmp(HDF_struct.children(ii).properties.name,'Coordinates')
            coordsys = HDF_struct.children(ii).attributes.coordsys;
        end
    end
    % read arraysys attribute
    arraysys = [];
    for ii=1:Nchild
        if strcmp(HDF_struct.children(ii).properties.name,'IndexMap')
            arraysys = HDF_struct.children(ii).attributes.arraysys;  
        end
    end
    % determine micropnts array
    Coordinates = HDFfileAScontrol('getdata',pathfilename,'Coordinates',':',':',':');
    if ndims(Coordinates)>2, Nseq = size(Coordinates,3); else, Nseq = 1; end
    Nprb = size(Coordinates,2);
    micropnts = zeros(Nprb*Nseq,3);
    for ii=1:Nseq
        ndeb = Nprb*(ii-1) + 1;
        micropnts(ndeb:ndeb-1+Nprb,:) = Coordinates(:,:,ii)';
    end
            
else
    badfilewarn;    
end
varargout{1} = micropnts;
if nargout>1, varargout{2} = coordsys; end
if nargout>2, varargout{3} = arraysys; end

%--------------------------------------------------------------------------------------------------------------------
case 'averageddata'
%
% [spectrum/quadsignal,frequency/time,axename] = HDFinterface('averageddata',pathfilename)
%  return the averaged spectrum/quadsignal data on the array and the corresponding frequency/time vector
% HDFinterface('averageddata',pathfilename) with no output arguments : plot the averaged data
%
%  spectrum       -> spectrum averaged on the array in the case of 'MeasSpct'
%  quadsignal     -> quadratic signal averaged on the array in the case of 'TimeSpct'
%  frequency/time -> axis
%  axename        -> name of the axis : 'Frequency' or 'Time'
%  pathfilename   -> file name (in optional with path)

pathfilename = varargin{2};
averageddata = []; 
axevector = [];
axename = [];

%-- read info
HDF_struct = HDFfileAScontrol('getinfo',pathfilename); 
if ~isempty(HDF_struct)
    Nchild = length(HDF_struct.children);
    Nprb = HDF_struct.attributes.probe_num;
    Nseq = HDF_struct.attributes.seq_num;
    Nref = HDF_struct.attributes.ref_num;
    % read the frequency or time vector
    if strcmp(HDF_struct.attributes.AStype,'MeasSpct')
        axename = 'Frequency'; 
        if any(strcmp(HDF_struct.properties.datasetnames,'Aspct_prb'))
            [Spectrum,axevector,axeprb,axeseq] = HDFfileAScontrol('getdata',pathfilename,'Aspct_prb',':',':',':');
        elseif Nprb==Nref
            [Xspectrum,axevector,axeorder,axeseq] = HDFfileAScontrol('getdata',pathfilename,'Refarray',':',':',':');
            Spectrum = TransRefarray('diagmat',Xspectrum,'multiple');
        else
            error('there are not datasets named ''Aspct_prb'' or ''Refarray'''),
        end
        averageddata = sum(sum(Spectrum,3),2)/Nprb/Nseq;

    elseif strcmp(HDF_struct.attributes.AStype,'MeasTime')
        axename = 'Time';
        [sig,axevector,axeprb,axeseq] = HDFfileAScontrol('getdata',pathfilename,'Time_prb',':',':',':');
        averageddata = sum(sum(sig.^2,3),2)/Nprb/Nseq;
    else
        error('incorrect AS type of HDF file'),
    end
else
    badfilewarn;    
end
if nargout==0&~isempty(averageddata)    
    [pathstr,filename,ext] = fileparts(pathfilename);
    if strcmp(axename,'Frequency')
        semilogy(axevector,averageddata), 
        title([filename ': MeasSpct averaged autospectra']), xlabel([axename '  [Hz]']),
    else
        plot(axevector,averageddata), set(gca,'Xlim',[axevector(1) axevector(end)]),
        title([filename ': MeasTime averaged quadratic signals']), xlabel([axename '  [s]']),
    end
else    
    varargout{1} = averageddata;
    if nargout>1, varargout{2} = axevector; end
    if nargout>2, varargout{3} = axename; end
end

%--------------------------------------------------------------------------------------------------------------------
case 'specgram'
%
% [specto,frequency,time] = HDFinterface('specgram',pathfilename) ...,prbvect) ...,seqvect) ...,nfft) ...,window) ...,numoverlap)
%  return the specto matrix and the corresponding frequency and time vectors (only with 'MeasTime' AStype)
% HDFinterface('specgram',...) with no output arguments : plot the specgram
%
%  specto      -> spectrum averaged on the array in the case of 'MeasSpct'
%  frequency   -> system of array 'rectangular' or 'cluster'
%  pathfilename -> file name (in optional with path)

pathfilename = varargin{2};
specto = []; 
frequency = [];
time = [];

%-- read info
HDF_struct = HDFfileAScontrol('getinfo',pathfilename); 
if ~isempty(HDF_struct)
    Nchild = length(HDF_struct.children);
    Nprb = HDF_struct.attributes.probe_num;
    Nseq = HDF_struct.attributes.seq_num;
    Fs = HDF_struct.attributes.sampling_frequency;
    if nargin>2, prbvect = varargin{3}; else, prbvect = []; end
    if nargin>3, seqvect = varargin{4}; else, seqvect = []; end
    if nargin>4, nfft = varargin{5}; else, nfft = []; end
    if nargin>5, window = varargin{6}; else, window = []; end
    if nargin>6, numoverlap = varargin{7}; else, numoverlap = []; end
    
    % read 'Time_prb'
    if strcmp(HDF_struct.attributes.AStype,'MeasTime')
        [sig,time,axeprb,axeseq] = HDFfileAScontrol('getdata',pathfilename,'Time_prb',':',':',':');
        if isempty(Fs), Fs = (time(end)-time(1))/(length(time)-1); end
        if isempty(prbvect), prbvect = 1:Nprb; end        
        if isempty(seqvect), seqvect = 1:Nseq; end

        for ii=prbvect
            for jj=seqvect
                if ii==prbvect(1)&jj==seqvect(1) 
                    [B,frequency] = specgram(sig(:,ii,jj),nfft,Fs,window,numoverlap); 
                    specto = abs(B).^2;
                else
                    B = specgram(sig(:,ii,jj),nfft,Fs,window,numoverlap);
                end
                specto = specto + abs(B).^2;
            end
        end
        specto = sqrt(specto/(length(prbvect)*length(seqvect)));
    else
        error('incorrect AS type of HDF file'),
    end
else
    badfilewarn;    
end
if nargout==0&~isempty(specto)    
    [pathstr,filename,ext] = fileparts(pathfilename);
    newplot;
    if length(time)==1 , time = [0 1/frequency(2)]; end
    imagesc(time,frequency,20*log10(abs(specto)+eps)); axis xy; colormap(jet)
    xlabel('Time  [s]'), ylabel('Frequency  [Hz]'),
    title([filename ': MeasTime spectogram']), 
else    
    varargout{1} = specto;
    if nargout>1, varargout{2} = frequency; end
    if nargout>2, varargout{3} = time; end
end

%--------------------------------------------------------------------------------------------------------------------
end

%====================================================================================================================
function badfilewarn

        cmsgbox = cell(1,1);
        cmsgbox{1} = 'not HDF file or bad HDF format';
        CreateStruct.WindowStyle = 'modal';
        CreateStruct.Interpreter = 'tex';
        h = msgbox(cmsgbox,'Open HDF file','warn',CreateStruct);
        set(h,'Color',[0.85 0.85 0.85]);
        HDFfileAScontrol('closeall');
