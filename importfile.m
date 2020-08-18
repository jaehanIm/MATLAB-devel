function gdLog = importfile(filename, startRow, endRow)
%IMPORTFILE1 텍스트 파일의 숫자형 데이터를 행렬로 가져옵니다.
%   GDLOG181031195314 = IMPORTFILE1(FILENAME) 디폴트 선택 항목의 텍스트 파일 FILENAME에서
%   데이터를 읽습니다.
%
%   GDLOG181031195314 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) 텍스트 파일
%   FILENAME의 STARTROW 행에서 ENDROW 행까지 데이터를 읽습니다.
%
% Example:
%   gdLog181031195314 = importfile1('gdLog_181031_195314.csv', 2, 8899);
%
%    TEXTSCAN도 참조하십시오.

% MATLAB에서 다음 날짜에 자동 생성됨: 2018/10/31 20:50:05

%% 변수를 초기화합니다.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% 각 텍스트 라인의 형식:
%   열1: double (%f)
%	열2: double (%f)
%   열3: double (%f)
%	열4: double (%f)
%   열5: double (%f)
%	열6: double (%f)
%   열7: double (%f)
%	열8: double (%f)
%   열9: double (%f)
%	열10: double (%f)
%   열11: double (%f)
%	열12: double (%f)
%   열13: double (%f)
%	열14: double (%f)
%   열15: double (%f)
%	열16: double (%f)
%   열17: double (%f)
%	열18: double (%f)
%   열19: double (%f)
%	열20: double (%f)
%   열21: double (%f)
%	열22: double (%f)
%   열23: double (%f)
%	열24: double (%f)
%   열25: double (%f)
%	열26: double (%f)
% 자세한 내용은 도움말 문서에서 TEXTSCAN을 참조하십시오.

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f %[^\n\r]'; % total 140 variables

%% 텍스트 파일을 엽니다.
fid = fopen(filename, 'r+');
fseek(fid, -1, 'eof');
fe = fread(fid, 1);
if fe ~= 13
    fseek(fid, -1, 'eof');
    fwrite(fid, 13);
end
fclose(fid);

fileID = fopen(filename,'r');

%% 형식에 따라 데이터 열을 읽습니다.
% 이 호출은 이 코드를 생성하는 데 사용되는 파일의 구조체를 기반으로 합니다. 다른 파일에 대한 오류가 발생하는 경우 가져오기 툴에서
% 코드를 다시 생성하십시오.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% 텍스트 파일을 닫습니다.
fclose(fileID);

%% 가져올 수 없는 데이터에 대한 사후 처리 중입니다.
% 가져오기 과정에서 가져올 수 없는 데이터에 규칙이 적용되지 않았으므로 사후 처리 코드가 포함되지 않았습니다. 가져올 수 없는
% 데이터에 사용할 코드를 생성하려면 파일에서 가져올 수 없는 셀을 선택하고 스크립트를 다시 생성하십시오.

%% 출력 변수 만들기
gdLog = table(dataArray{1:end-1}, 'VariableNames', {
    'rosTime','flightMode','ctrlDeviceStatus', ...
    'fcMcMode','nSat','gpsFix', 'jobSeq', ...
    'velNedGps_0','velNedGps_1','velNedGps_2', ...
    'posNed_0','posNed_1','posNed_2', ...
    'velNed_0','velNed_1','velNed_2', ...
    'rpy_0','rpy_1','rpy_2','ySpType', ... % 20
    'ctrlUser','ctrlStruct','ctrlSetpointType','ctrlOutputType', ...
    'ctrlSp_0','ctrlSp_1','ctrlSp_2', ...
    'ctrlOp_0','ctrlOp_1','ctrlOp_2','ySp', ...
    'rcRoll', 'rcPitch', 'rcYaw', 'rcThrottle', ... % 20+15
    'GpsNSV', 'RtkHealthFlag', 'GpsFusedNSV', 'GpHealth', ...
    'posGPS_0', 'posGPS_1', 'posGPS_2', ...
    'posRTK_0', 'posRTK_1', 'posRTK_2', ...
    'posGpsFused_0', 'posGpsFused_1', 'posGpsFused_2', ...
    'posGp_0', 'posGp_1', 'posGp_2', ...
    'errLatMix', 'errLatVis', 'errLatLid', 'cmdLatVelIgain', 'cmdLatVelMix', ...
    'errLatMixDt', 'errLatMixCov00', 'errLatMixCov11', ... % 20+15+24
    'vbx', 'vby', 'vbz', ...
    'AcWarnStat', 'AcHorWarnAc', 'AcVerWarnAc', ...
    'AcXRel', 'AcYRel', 'AcZRel', ...
    'AcHorWarnRange', 'AcHorWarnAngle', 'AcVerWarnRange', 'AcVerWarnAngle', ...
    'LidarDist', 'LidarAngle', ... // 20+15+24+15
    'LidarRaw_0', 'LidarRaw_1', 'LidarRaw_2', 'LidarRaw_3', ...
    'LidarRaw_4', 'LidarRaw_5', 'LidarRaw_6', 'LidarRaw_7', ...
    'LongVelCmd', 'LatVelCmd', 'HeaveVelCmd', ...
    'velCtrlI_u', 'velCtrlI_v', 'velCtrlI_d', ...
    'posCtrlI_N', 'posCtrlI_E', 'posCtrlI_D', ...
    'gimbalRPY_0', 'gimbalRPY_1', 'gimbalRPY_2', ... // 20+15+24+15+20
    'windStatus', 'windSpeed', 'windAngle', 'windQueryTime', 'windResponseTime', ...
    'acousticTemp', 'tempQueryTime', 'tempResponseTime', ...
    'accBody_0', 'accBody_1', 'accBody_2', ... // 20+15+24+15+20+11
    'trajUnitVectorT_0', 'trajUnitVectorT_1', 'trajUnitVectorT_2', ...
    'trajUnitVectorN_0', 'trajUnitVectorN_1', 'trajUnitVectorN_2', ...
    'trajUnitVectorB_0', 'trajUnitVectorB_1', 'trajUnitVectorB_2', ...
    'trajCtrlI_T', 'trajCtrlI_N', 'trajCtrlI_B', ...
    'StdJobLongPidErr', 'StdJobLongPidRate', 'StdJobLongPidIgain', ...
    'GuideModeLongPidErr', 'GuideModeLongPidRate', 'GuideModeLongPidIgain', ... % 20+15+24+15+20+11+18
    'pqr_0', 'pqr_1', 'pqr_2', ... 
    'rpdCmd_0', 'rpdCmd_1', 'rpdCmd_2', ... 
    'velCmdNav_0', 'velCmdNav_1', 'velCmdNav_2', ...
    'posCmdNed_0', 'posCmdNed_1', 'posCmdNed_2', ...
    'missionType', 'jobType', ...
    'bladeTravelDistance', 'trajTimeCur', 'trajTimeMax'}); % 20+15+24+15+20+11+18+17 = 140

