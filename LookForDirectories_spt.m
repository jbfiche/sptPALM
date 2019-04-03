function FinalDirectories = LookForDirectories_spt(DirectoryName)

dim = 0;
AllDirectories = {};
FinalDirectories = {};
dirinfo = dir();

% Look for the folders corresponding to an ROI or the TL_Channel or
% containing a Kymo results file
% ------------------------------


dirinfo_ROI = dir('*tif.mat');

if ~isempty(dirinfo_ROI)
    if isunix
        FinalDirectories{end+1,1} = strcat(cd, '/', dirinfo_ROI.name);
    else
        FinalDirectories{end+1,1} = strcat(cd, '\', dirinfo_ROI.name);
    end
end

% Look at the folders inside the selected directory
% -------------------------------------------------

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
                Path{nSubFolder};
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
                    if isdir(NewPath)
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

Directories = {};

for n = 1 : size(AllDirectories,1)
    for m = 1 : size(AllDirectories,2)
        
        Directories = cat(1, Directories, AllDirectories{n,m});
    end
end

for n = 1 : size(Directories,1)
    cd(Directories{n})
    
    dirinfo_ROI = dir('*tif.mat');
    
    if ~isempty(dirinfo_ROI)
        if isunix
            FinalDirectories{end+1,1} = strcat(Directories{n}, '/', dirinfo_ROI.name);
        else
            FinalDirectories{end+1,1} = strcat(Directories{n}, '\', dirinfo_ROI.name);
        end
    end
end

cd(DirectoryName)