function IMFresult = EMD(x,y)

% x = 1:0.1:500*pi;
% y = sin(x) + sin(x/2.7)*2;
% noise = randn(1,length(x))/100; 
% % noise = 0
% yOrig = y + noise;
% y = y + noise;

IMFresult = [];

for i = 1:100
    [IMF,residue] = EMDiter(x,y);
    if isempty(IMF)
        break
    end
    y = residue;
    IMFresult(i,:) = IMF;
end

% figure(1)
% clf
% hold on
% grid on
% L = size(IMFresult,1);
% subplot(L+1,1,1)
% plot(x,yOrig);
% for i = 1:size(IMFresult,1)
%     subplot(L+1,1,i+1)
%     plot(x,IMFresult(i,:))
% end
% legend()

end