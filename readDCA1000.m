%%% This script is used to read the binary file produced by the DCA1000
%%% and Mmwave Studio
%%% Command to run in Matlab GUI - readDCA1000('<ADC capture bin file>') 
function [retVal] = readDCA1000(folder_locaion)
%% global variables
% change based on sensor config
numADCSamples = 128; % number of ADC samples per chirp
numRX = 4; % number of receivers
numTX = 3; % number of receivers
Loops = 96; % Number of transmissions per antenna
%% read file
% 获取目录下所有文件
files = dir(folder_locaion);

% 过滤出所有 .bin 文件
binFiles = files(endsWith({files.name}, '.bin'));

% 如果没有找到 .bin 文件，返回空
if isempty(binFiles)
    disp('No .bin files found in the specified directory.');
    return;
end

% 预分配 cell 数组存储数据，提高拼接效率
adcDataCell = cell(1, length(binFiles));

% 读取所有 .bin 文件
for index = 1:length(binFiles)
    file_name = binFiles(index).name;
    file_location = fullfile(folder_locaion, file_name);
    
    % 打开文件
    fid = fopen(file_location, 'r');
    if fid == -1
        warning(['Could not open file: ', file_name]);
        continue;
    end
    
    % 读取 .bin 数据
    adcDataCell{index} = fread(fid, 'uint16', 'l');
    
    % 关闭文件
    fclose(fid);
    
    disp(['Processed file: ', file_name]);
end

% 将 cell 转换为数组，提高拼接效率
adcData = vertcat(adcDataCell{:});

disp('All bin files have been read and accumulated.');

fileSize = size(adcData, 1);

numFrames = fileSize / (2*numADCSamples * numRX * numTX * Loops);

% data is captured over two LVDS lanes (parallel), 
% calculate how many parallel data, with shape (LVDStimes, 2).
LVDS = zeros(1, fileSize/2);
%combine real and imaginary part into complex data
%read in file: 2I is followed by 2Q
counter = 1;
for i=1:4:fileSize
    LVDS(1,counter)   = adcData(i) + 1i * adcData(i+2); 
    LVDS(1,counter+1) = adcData(i+1) + 1i * adcData(i+3); 
    counter = counter + 2;
end

Samples_Chirp = numADCSamples * numRX; % 128*4 = 512
numChirps  = numTX * Loops; % 3*96 = 288
Samples_Frame = Samples_Chirp * numChirps; % 512*288 = 147456

%organize data
adcData = zeros(numFrames,numChirps,numRX,numADCSamples);

for i = 1:numFrames
    frame = LVDS(1, (i-1)*Samples_Frame+1:i*Samples_Frame); % samples of one frame
    for j = 1:numChirps % 288
        chirp = frame(1, (j-1)*Samples_Chirp+1:j*Samples_Chirp); % samples of one chirp
        for k = 1:numRX % 4
            rx = chirp(1, (k-1)*numADCSamples+1:k*numADCSamples); % samples of one rx
            adcData(i,j,k,:) = rx; % 128
        end
    end
end

% return receiver data
retVal = adcData;