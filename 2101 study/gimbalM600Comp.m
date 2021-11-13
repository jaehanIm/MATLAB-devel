cmdTemp = zeros(1,7000);
cmdTemp(cmdIndex(1).DataIndex:cmdIndex(2).DataIndex) = 3;
cmdTemp(cmdIndex(3).DataIndex+10:cmdIndex(4).DataIndex+10) = 5;
cmdTemp(cmdIndex(5).DataIndex+10:cmdIndex(6).DataIndex+10) = 10;

data_1 = readtable('/home/jaehan/log/210830/gdLog_210831_114817.csv');
time = seconds(data_1.rosTime - data_1.rosTime(1));

figure(1)
clf
plot(time,data_1.gimbalRpy_deg_1)
hold on
plot(time,data_1.gimbalRpyCmd_deg_1)
grid on
plot(time,data_1.jobSeq)

data_2 = readtable('/home/jaehan/log/210902/gdLog_210902_193651.csv');
time = seconds(data_2.rosTime - data_2.rosTime(1));
figure(2)
clf
plot(time,data_2.gimbalRpy_deg_1)
hold on
grid on
plot(time, cmdTemp(1:length(time)))




cmd{1,1} = data_1.gimbalRpyCmd_deg_1(index1(1).DataIndex:index1(2).DataIndex);
cmd{1,2} = data_1.gimbalRpyCmd_deg_1(index1(3).DataIndex:index1(4).DataIndex);
cmd{1,3} = data_1.gimbalRpyCmd_deg_1(index1(5).DataIndex:index1(6).DataIndex);
rsp{1,1} = data_1.gimbalRpy_deg_1(index1(1).DataIndex:index1(2).DataIndex);
rsp{1,2} = data_1.gimbalRpy_deg_1(index1(3).DataIndex:index1(4).DataIndex);
rsp{1,3} = data_1.gimbalRpy_deg_1(index1(5).DataIndex:index1(6).DataIndex);
cmd{2,1} = cmdTemp(index2(1).DataIndex:index2(2).DataIndex)';
cmd{2,2} = cmdTemp(index2(3).DataIndex:index2(4).DataIndex)';
cmd{2,3} = cmdTemp(index2(5).DataIndex:index2(6).DataIndex)';
rsp{2,1} = data_2.gimbalRpy_deg_1(index2(1).DataIndex:index2(2).DataIndex);
rsp{2,2} = data_2.gimbalRpy_deg_1(index2(3).DataIndex:index2(4).DataIndex);
rsp{2,3} = data_2.gimbalRpy_deg_1(index2(5).DataIndex:index2(6).DataIndex);

%%

for i = 1:3
[Num, Den, delay, tfEq] = estimate_tf(cmd{1,i},rsp{1,i});
tfEq
figure()
step(tfEq)
stepinfo(tfEq)
end

for i = 1:3
[Num, Den, delay, tfEq] = estimate_tf(cmd{2,i},rsp{2,i});
tfEq
figure()
step(tfEq)
stepinfo(tfEq)
end

%% 
function [Num, Den, delay, tfEq] = estimate_tf(cmd,resp)
nD = 3;
nN = 2;

Cmd = cmd;
Response =  resp;

Cmd = detrend(Cmd,0);
Response = detrend(Response,0);

timeseriesSet = iddata(Response,Cmd,0.02);
sys = tfest(timeseriesSet,nD,nN,nan);
Num = sys.Numerator;
Den = sys.Denominator;
delay = sys.IODelay;

tfEq = tf(Num,Den,'Iodelay',delay);
end