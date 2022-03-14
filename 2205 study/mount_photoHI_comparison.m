%% yongdaeri
sample{1} = [2 4 12 16 24 28 30 32 40 67 92 105 109 110 111 112 116 121 122 131 138 140 169 174 178 182 192 194 204 206 210 230 232 234 244 248 261 263 265 268 270 272 274 281 283 286 288 305 307];
sample{2} = [2 6 12 16 32 121 122 169 174 175 178 180 182 184 200 203 230 232 258 259 260 261 265 266 268 269 270 272 274 281 283 286 290 295 310 311 314 316 318 320 324];
sample{3} = [2, 12, 16, 32, 121, 122, 169, 174, 178, 180, 182, 230, 281, 288, 295, 315, 328];

odd_single = 0; odd_double = 0; odd_triple = 0;
eve_single = 0; eve_double = 0; eve_triple = 0;
for i = 1:341
count = 0;
for j = 1:3
if ismember(i,sample{j})
count = count + 1;
end
end
if mod(i,2) == 0
if count == 1
eve_single = eve_single+1;
elseif count == 2
eve_double = eve_double+1;
elseif count == 3
eve_triple = eve_triple+1;
end
else
if count == 1
odd_single = odd_single+1;
elseif count == 2
odd_double = odd_double+1;
elseif count == 3
odd_triple = odd_triple+1;
end
end
end

result = [eve_single,eve_double,eve_triple;...
odd_single,odd_double,odd_triple];

figure(1)
clf
hold on
grid on
bar(result)
xticks = ([1 2]);
xticklabels({'stiff-soft','','soft-soft'})
legend('single','double','triple')
title('HI result')
ylabel('count')

%% oksang
s1{1} = [7 10 13 17 24 31 32 35 45 48 51 56 65 67 75 81 104 114 125 127 134 144 147 181 183 185 186 190 192 194 197 207];
s1{2} = [31, 45, 51, 56, 65, 67, 81, 104, 117, 134, 144, 147, 185, 186, 192, 194, 197, 207];
s1{3} = [7 31 51 56 61 65 67 75 81 104 114 147 181 183 185 186 192 194 195  197 207];
s2{1} = [4 11 20 34 35 38 58 64 68 89 96 102 106 113 132 134];
s2{2} = [4, 11, 20, 34, 35, 48, 50, 58, 64, 68, 102, 106, 113, 132, 134];
s2{3} = [3 11 20 34 35 38 47 48 50 58 64 68 96 102 106 113 132 134];

s1_single = []; s1_double = []; s1_triple = [];
s2_single = []; s2_double = []; s2_triple = [];

for i = 1:207
    count1 = 0;
    count2 = 0;
    for j = 1:3
        if ismember(i,s1{j})
            count1 = count1 + 1;
        end
        if ismember(i,s2{j})
            count2 = count2 + 1;
        end        
    end

    if count1 == 1
        s1_single = vertcat(s1_single,i);
    elseif count1 == 2
        s1_double = vertcat(s1_double,i);
    elseif count1 == 3
        s1_triple = vertcat(s1_triple,i);
    end

    if count2 == 1
        s2_single = vertcat(s2_single,i);
    elseif count2 == 2
        s2_double = vertcat(s2_double,i);
    elseif count2 == 3
        s2_triple = vertcat(s2_triple,i);
    end
end

s1_list = dir('/home/jaehan/Desktop/220303_image/maxGen/result_s1/');
s2_list = dir('/home/jaehan/Desktop/220303_image/maxGen/result_s2/');

s1_list = s1_list(3:end); s2_list = s2_list(3:end);

s1_data = zeros(size(s1_list,1),1); % FFTS score
s2_data = zeros(size(s2_list,1),1);
for i = 1:size(s1_list,1)
    fix = find(s1_list(i).name == '_');
    num = s1_list(i).name(1:fix-1);
    val = s1_list(i).name(fix+1:end-4);
    s1_data(str2num(num)) = str2num(val);
end
for i = 1:size(s2_list,1)
    fix = find(s2_list(i).name == '_');
    num = s2_list(i).name(1:fix-1);
    val = s2_list(i).name(fix+1:end-4);
    s2_data(str2num(num)) = str2num(val);
end

s1_score = zeros(size(s1_data,1),1); % HI
s2_score = zeros(size(s2_data,1),1);
for i = 1:size(s1_single,1)
    s1_score(s1_single(i)) = 1;
end
for i = 1:size(s1_double,1)
    s1_score(s1_double(i)) = 2;
end
for i = 1:size(s1_triple,1)
    s1_score(s1_triple(i)) = 3;
end
for i = 1:size(s2_single,1)
    s2_score(s2_single(i)) = 1;
end
for i = 1:size(s2_double,1)
    s2_score(s2_double(i)) = 2;
end
for i = 1:size(s2_triple,1)
    s2_score(s2_triple(i)) = 3;
end

s1_omitIdx = find(s1_data==0);
s2_omitIdx = find(s2_data==0);

s1_data(s1_omitIdx) = [];
s2_data(s2_omitIdx) = [];
s1_score(s1_omitIdx) = [];
s2_score(s2_omitIdx) = [];

figure(2)
clf
hold on
grid on
plot(s1_data,s1_score,'o')
plot(s2_data,s2_score,'o')

mission_FFTS = [13.2540, 14.1874, mean(s1_data), mean(s2_data)];
mission_score = [0.7806, 0.6122, mean(s1_score), mean(s2_score)];

figure(3)
clf
hold on
grid on
plot(mission_FFTS,mission_score,'r*','MarkerSize',10,'LineWidth',3);