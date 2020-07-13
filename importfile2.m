function gdLog = importfile2(filename, startRow, endRow)
%IMPORTFILE1 �ؽ�Ʈ ������ ������ �����͸� ��ķ� �����ɴϴ�.
%   GDLOG181031195314 = IMPORTFILE1(FILENAME) ����Ʈ ���� �׸��� �ؽ�Ʈ ���� FILENAME����
%   �����͸� �н��ϴ�.
%
%   GDLOG181031195314 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) �ؽ�Ʈ ����
%   FILENAME�� STARTROW �࿡�� ENDROW ����� �����͸� �н��ϴ�.
%
% Example:
%   gdLog181031195314 = importfile1('gdLog_181031_195314.csv', 2, 8899);
%
%    TEXTSCAN�� �����Ͻʽÿ�.

% MATLAB���� ���� ��¥�� �ڵ� ������: 2018/10/31 20:50:05

%% ������ �ʱ�ȭ�մϴ�.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% �� �ؽ�Ʈ ������ ����:
%   ��1: double (%f)
%	��2: double (%f)
%   ��3: double (%f)
%	��4: double (%f)
%   ��5: double (%f)
%	��6: double (%f)
%   ��7: double (%f)
%	��8: double (%f)
%   ��9: double (%f)
%	��10: double (%f)
%   ��11: double (%f)
%	��12: double (%f)
%   ��13: double (%f)
%	��14: double (%f)
%   ��15: double (%f)
%	��16: double (%f)
%   ��17: double (%f)
%	��18: double (%f)
%   ��19: double (%f)
%	��20: double (%f)
%   ��21: double (%f)
%	��22: double (%f)
%   ��23: double (%f)
%	��24: double (%f)
%   ��25: double (%f)
%	��26: double (%f)
% �ڼ��� ������ ���� �������� TEXTSCAN�� �����Ͻʽÿ�.

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f %[^\n\r]'; % total 128 variables

%% �ؽ�Ʈ ������ ���ϴ�.
fid = fopen(filename, 'r+');
fseek(fid, -1, 'eof');
fe = fread(fid, 1);
if fe ~= 13
    fseek(fid, -1, 'eof');
    fwrite(fid, 13);
end
fclose(fid);

fileID = fopen(filename,'r');

%% ���Ŀ� ���� ������ ���� �н��ϴ�.
% �� ȣ���� �� �ڵ带 �����ϴ� �� ���Ǵ� ������ ����ü�� ������� �մϴ�. �ٸ� ���Ͽ� ���� ������ �߻��ϴ� ��� �������� ������
% �ڵ带 �ٽ� �����Ͻʽÿ�.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% �ؽ�Ʈ ������ �ݽ��ϴ�.
fclose(fileID);

%% ������ �� ���� �����Ϳ� ���� ���� ó�� ���Դϴ�.
% �������� �������� ������ �� ���� �����Ϳ� ��Ģ�� ������� �ʾ����Ƿ� ���� ó�� �ڵ尡 ���Ե��� �ʾҽ��ϴ�. ������ �� ����
% �����Ϳ� ����� �ڵ带 �����Ϸ��� ���Ͽ��� ������ �� ���� ���� �����ϰ� ��ũ��Ʈ�� �ٽ� �����Ͻʽÿ�.

%% ��� ���� �����
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
    'LidarRaw_0', 'LidarRaw_1', 'LidarRaw_2', 'LidarRaw_3', 'LidarRaw_4', 'LidarRaw_5', 'LidarRaw_6', 'LidarRaw_7',... // 20+15+24+15+8
    'LongVelCmd', 'LatVelCmd', 'HeaveVelCmd', ...
    'velCtrlI_u', 'velCtrlI_v', 'velCtrlI_d', ...
    'posCtrlI_N', 'posCtrlI_E', 'posCtrlI_D', ...
    'gimbalRPY_0', 'gimbalRPY_1', 'gimbalRPY_2', ... // 20+15+24+15+8+12
    'windStatus', 'windSpeed', 'windAngle', 'windQueryTime', 'windResponseTime', ...
    'acousticTemp', 'tempQueryTime', 'tempResponseTime', ...
    'accBody_0', 'accBody_1', 'accBody_2', ... // 20+15+24+15+8+12+11
    'trajUnitVectorT_0', 'trajUnitVectorT_1', 'trajUnitVectorT_2', ...
    'trajUnitVectorN_0', 'trajUnitVectorN_1', 'trajUnitVectorN_2', ...
    'trajUnitVectorB_0', 'trajUnitVectorB_1', 'trajUnitVectorB_2', ...
    'trajCtrlI_T', 'trajCtrlI_N', 'trajCtrlI_B', ...
    'StdJobLongPidErr', 'StdJobLongPidRate', 'StdJobLongPidIgain', ...
    'GuideModeLongPidErr', 'GuideModeLongPidRate', 'GuideModeLongPidIgain', ... % 20+15+24+15+8+12+11+18
    'pqr_0', 'pqr_1', 'pqr_2', ... 
    'rpdCmd_0', 'rpdCmd_1', 'rpdCmd_2', ... 
    'velCmdNav_0', 'velCmdNav_1', 'velCmdNav_2', ...
    'posCmdNed_0', 'posCmdNed_1', 'posCmdNed_2', ...
    'missionType', 'jobType', 'bladeTravelDistance', ...
    'trajTimeCur', 'trajTimeMax'}); % 20+15+24+15+8+12+11+18+17 = 140
