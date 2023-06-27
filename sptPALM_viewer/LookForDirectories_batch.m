%*****************************
%
% LookForDirectories_spt.m

% ****************************
%
% JB Fiche
% June, 2023
% fiche@cbs.cnrs.fr
% -------------------------------------------------------------------------
% Purpose: Look for the all the folders containing at least one csv file
% (output by TrackMate analysis). 
% -------------------------------------------------------------------------
% Specific: 
% -------------------------------------------------------------------------
% To fix: 
% -------------------------------------------------------------------------
% Copyright Centre National de la Recherche Scientifique, 2020.


function FinalDirectories = LookForDirectories_batch(FileName)

dim = 0;
AllDirectories = {};
FinalDirectories = {};

% Select the root directory where the search will take place
% ----------------------------------------------------------
DirectoryName = uigetdir;
cd(DirectoryName)
dirinfo = dir();

% Look at the folders and subfolders inside the selected directory. This
% step will return a list of all the directories within the initial
% directly indicated by "DirectoryName". All the directories names are
% saved in the cell "AllDirectories".
% -----------------------------------

dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
NFolder = 0;

if length(dirinfo) > 2
    AllDirectories = cell(length(dirinfo)-2, 1);
    for k = 3 : length(dirinfo) % The two first are not directories '.' and '..'
        if isunix
            Path = strcat(DirectoryName, '/', dirinfo(k).name);
        else
            Path = strcat(DirectoryName, '\', dirinfo(k).name);
        end
        
        % For MatlabR2017 the function is isdir. For R2019a is dir does not
        % work anymore and needs to be replaced by isfolder.
        % --------------------------------------------------
        
        if isdir(Path)
            if isunix
                AllDirectories{k-2} = strcat(DirectoryName, '/', dirinfo(k).name);
            else
                AllDirectories{k-2} = strcat(DirectoryName, '\', dirinfo(k).name);
            end
            NFolder = NFolder + 1;
        end
    end
end

while NFolder > 0
    
    dim = dim + 1;
    NFolder = 0;
    
    for nFolder = 1 : size(AllDirectories,1)
        
        Path = AllDirectories{nFolder, dim};
        AllSubDirectories = {};
        
        for nSubFolder = 1 : size(Path, 1)
            if iscell(Path)
                dirinfo = dir(Path{nSubFolder});
            else
                dirinfo = dir(Path);
            end
            dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
            
            if length(dirinfo) > 2
                %                 AllSubDirectories = cell(length(dirinfo)-2, 1);
                for k = 3 : length(dirinfo) % The two first are not directories '.' and '..'
                    if iscell(Path)
                        if isunix
                            NewPath = strcat(Path{nSubFolder}, '/', dirinfo(k).name);
                        else
                            NewPath = strcat(Path{nSubFolder}, '\', dirinfo(k).name);
                        end
                    else
                        if isunix
                            NewPath = strcat(Path, '/', dirinfo(k).name);
                        else
                            NewPath = strcat(Path, '\', dirinfo(k).name);
                        end
                    end
                    if isfolder(NewPath)
                        %                         AllSubDirectories{k-2} = strcat(NewPath);
                        AllSubDirectories{end+1,1} = strcat(NewPath);
                        NFolder = NFolder + 1;
                    end
                end
            end
        end
        
        AllDirectories{nFolder, dim+1} = AllSubDirectories;
        
    end
end

% Concatenate all the names from the 2D cell "AllDirectories" into a single
% 1D cell "Directories".
% ----------------------

Directories = {};

for n = 1 : size(AllDirectories,1)
    for m = 1 : size(AllDirectories,2)
        
        Directories = cat(1, Directories, AllDirectories{n,m});
    end
end

% For each directory, check that they contain at least one file as 
% indicated by "FileName_format".
% ----------------------------------

for n = 1 : size(Directories,1)
    cd(Directories{n})
    
    dirinfo_ROI = dir(FileName);
    
    if ~isempty(dirinfo_ROI)
        FinalDirectories{end+1,1} = Directories{n};
    end
end

cd(DirectoryName)